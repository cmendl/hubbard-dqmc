#include "monte_carlo.h"
#include "greens_func.h"
#include "util.h"
#include "dupio.h"
#include "profiler.h"
#include "checkpoint.h"
#include "progress.h"
#include <mkl.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>
#include <stdint.h>


//________________________________________________________________________________________________________________________
///
/// \brief Perform a Determinant Quantum Monte Carlo (DQMC) iteration
///
/// \param kinetic              matrix exponential of kinetic energy operator
/// \param stratonovich_params  precomputed Hubbard-Stratonovich parameters
/// \param nwraps               number of "time slice wraps" before recomputing the Green's function; must be a multiple of 'prodBlen'
/// \param seed                 random number "seed" structure, will be updated during function call
/// \param s                    Hubbard-Stratonovich field; will be updated
/// \param tsm_u                time step matrices for the spin-up Green's function, will be updated
/// \param tsm_d                time step matrices for the spin-down Green's function, will be updated
/// \param Gu                   spin-up   Green's function, must have been computed on input and will be updated
/// \param Gd                   spin-down Green's function, must have been computed on input and will be updated
/// \param neqlt                perform an equal time measurement every 'neqlt' time slices
/// \param meas_data            measurement data
///
void DQMCIteration(const kinetic_t *restrict kinetic, const stratonovich_params_t *restrict stratonovich_params, const int nwraps,
	randseed_t *restrict seed, spin_field_t *restrict s, time_step_matrices_t *restrict tsm_u, time_step_matrices_t *restrict tsm_d,
	greens_func_t *restrict Gu, greens_func_t *restrict Gd, const int neqlt, measurement_data_t *restrict meas_data)
{
	Profile_Begin("DQMCIter");
	__assume_aligned(s, MEM_DATA_ALIGN);

	// dimension consistency checks
	assert(tsm_u->N == kinetic->Ncell * kinetic->Norb);
	assert(tsm_u->N == tsm_d->N);
	assert(tsm_u->L == tsm_d->L);
	const int N = tsm_u->N;
	const int Ncell = kinetic->Ncell;
	const int L = tsm_u->L;

	assert(tsm_u->prodBlen == tsm_d->prodBlen);
	assert(nwraps % tsm_u->prodBlen == 0);	// must be a multiple of 'prodBlen'

	// store Green's functions before recomputing them to estimate error
	#if defined(DEBUG) | defined(_DEBUG)
	greens_func_t Gu_old, Gd_old;
	AllocateGreensFunction(N, &Gu_old);
	AllocateGreensFunction(N, &Gd_old);
	__assume_aligned(Gu_old.mat, MEM_DATA_ALIGN);
	__assume_aligned(Gd_old.mat, MEM_DATA_ALIGN);
	#endif

	// random shuffle of lattice cells and orbitals
	int *orb_cell_order = MKL_malloc(N * sizeof(int), MEM_DATA_ALIGN);
	__assume_aligned(orb_cell_order, MEM_DATA_ALIGN);

	// iterate over time slices
	int l;
	for (l = 0; l < L; l++)
	{
		Profile_Begin("DQMCIter_Wraps");
		#pragma omp parallel sections
		{
			#pragma omp section
			GreenTimeSliceWrap(N, tsm_u->B[l], tsm_u->invB[l], Gu->mat);
			#pragma omp section
			GreenTimeSliceWrap(N, tsm_d->B[l], tsm_d->invB[l], Gd->mat);
		}
		Profile_End("DQMCIter_Wraps");

		// iterate over lattice sites in random order, updating the Hubbard-Stratonovich field
		Profile_Begin("DQMCIter_SiteUpdate");
		Random_Shuffle(seed, N, orb_cell_order);
		int j;
		for (j = 0; j < N; j++)
		{
			const int i = orb_cell_order[j];
			const int o = i / Ncell;	// orbital index
			assert(0 <= i && i < N);
			assert(0 <= o && o < kinetic->Norb);

			// Eq. (13)
			// suggest flipping s_{i,l}
			const double du = 1 + (1 - Gu->mat[i + i*N]) * stratonovich_params->delta[  s[i + l*N]][o];
			const double dd = 1 + (1 - Gd->mat[i + i*N]) * stratonovich_params->delta[1-s[i + l*N]][o];
			if (Random_GetUniform(seed) < fabs(du*dd))
			{
				// Eq. (15)
				#pragma omp parallel sections
				{
					#pragma omp section
					GreenShermanMorrisonUpdate(stratonovich_params->delta[  s[i + l*N]][o], N, i, Gu->mat);
					#pragma omp section
					GreenShermanMorrisonUpdate(stratonovich_params->delta[1-s[i + l*N]][o], N, i, Gd->mat);
				}
				// correspondingly update determinants
				Gu->logdet -= log(fabs(du));
				Gd->logdet -= log(fabs(dd));
				if (du < 0) { Gu->sgndet = -Gu->sgndet; }
				if (dd < 0) { Gd->sgndet = -Gd->sgndet; }

				// actually flip spin of Hubbard-Stratonovich field entry
				s[i + l*N] = 1 - s[i + l*N];
			}
		}
		Profile_End("DQMCIter_SiteUpdate");

		// re-compute corresponding B matrices
		Profile_Begin("DQMCIter_Brecomp");
		#pragma omp parallel sections
		{
			#pragma omp section
			UpdateTimeStepMatrices(kinetic, stratonovich_params->expVu, s, l, tsm_u);
			#pragma omp section
			UpdateTimeStepMatrices(kinetic, stratonovich_params->expVd, s, l, tsm_d);
		}
		Profile_End("DQMCIter_Brecomp");

		// recompute Green's function after several time slice "wraps"
		if ((l + 1) % nwraps == 0)
		{
			Profile_Begin("DQMCIter_Grecomp");
			// store current Green's function matrices to compare with newly constructed ones
			#if defined(DEBUG) | defined(_DEBUG)
			CopyGreensFunction(Gu, &Gu_old);
			CopyGreensFunction(Gd, &Gd_old);
			#endif

			#pragma omp parallel sections
			{
				#pragma omp section
				GreenConstruct(tsm_u, (l + 1) % L, Gu);
				#pragma omp section
				GreenConstruct(tsm_d, (l + 1) % L, Gd);
			}

			#if defined(DEBUG) | defined(_DEBUG)
			// deviation of matrix entries
			double err_u = UniformDistance(N*N, Gu_old.mat, Gu->mat);
			double err_d = UniformDistance(N*N, Gd_old.mat, Gd->mat);
			if (err_u > 1e-8*N || err_d > 1e-8*N) {
				duprintf("Warning: after calling 'GreenConstruct()', largest entrywise distance between previous and new Green's functions is %g (up) and %g (down).\n", err_u, err_d);
			}
			// deviation of matrix determinants
			err_u = fabs(Gu_old.logdet - Gu->logdet);
			err_d = fabs(Gd_old.logdet - Gd->logdet);
			if (err_u > 1e-8 || err_d > 1e-8) {
				duprintf("Warning: after calling 'GreenConstruct()', largest distance between logarithm of previous and new Green's function determinants is %g (up) and %g (down).\n", err_u, err_d);
			}
			if (Gu_old.sgndet != Gu->sgndet || Gd_old.sgndet != Gd->sgndet) {
				duprintf("Warning: after calling 'GreenConstruct()', determinant sign has changed.\n");
			}
			#endif
			Profile_End("DQMCIter_Grecomp");
		}

		if (neqlt > 0 && (l + 1) % neqlt == 0)
		{
			// accumulate equal time "measurement" data
			Profile_Begin("DQMCIter_AccumulateEqMeas");
			AccumulateMeasurement(Gu, Gd, meas_data);
			Profile_End("DQMCIter_AccumulateEqMeas");
		}
	}

	// clean up
	MKL_free(orb_cell_order);
	#if defined(DEBUG) | defined(_DEBUG)
	DeleteGreensFunction(&Gd_old);
	DeleteGreensFunction(&Gu_old);
	#endif
	Profile_End("DQMCIter");
}



//________________________________________________________________________________________________________________________
///
/// \brief Perform block updates of the phonon field
///
/// \param dt                   imaginary-time step
/// \param kinetic              matrix exponential of kinetic energy operator
/// \param stratonovich_params  precomputed Hubbard-Stratonovich parameters
/// \param phonon_params        phonon field parameters
/// \param seed                 random number "seed" structure, will be updated during function call
/// \param s                    Hubbard-Stratonovich field, remains constant during phonon block updates
/// \param X                    phonon field, will be updated
/// \param expX                 entrywise exponential of the phonon field: exp(-dt g X)
/// \param tsm_u                time step matrices for the spin-up Green's function, will be updated
/// \param tsm_d                time step matrices for the spin-down Green's function, will be updated
/// \param Gu                   spin-up Green's function, will be updated
/// \param Gd                   spin-down Green's function, will be updated
///
void PhononBlockUpdates(const double dt, const kinetic_t *restrict kinetic, const stratonovich_params_t *restrict stratonovich_params,
	const phonon_params_t *restrict phonon_params, randseed_t *restrict seed, const spin_field_t *restrict s, double *restrict X, double *restrict expX,
	time_step_matrices_t *restrict tsm_u, time_step_matrices_t *restrict tsm_d, greens_func_t *restrict Gu, greens_func_t *restrict Gd)
{
	__assume_aligned(   X, MEM_DATA_ALIGN);
	__assume_aligned(expX, MEM_DATA_ALIGN);

	// dimension consistency checks
	assert(tsm_u->N == kinetic->Ncell * kinetic->Norb);
	assert(tsm_u->N == tsm_d->N);
	assert(tsm_u->L == tsm_d->L);
	const int N = tsm_u->N;
	const int Ncell = kinetic->Ncell;
	const int L = tsm_u->L;

	// fast return for zero block updates
	if (phonon_params->nblock_updates <= 0) {
		return;
	}

	// store X_{i,l} and corresponding exponential for all 'l'
	double *X_i    = (double *)MKL_malloc(L * sizeof(double), MEM_DATA_ALIGN);
	double *expX_i = (double *)MKL_malloc(L * sizeof(double), MEM_DATA_ALIGN);
	__assume_aligned(X_i,    MEM_DATA_ALIGN);
	__assume_aligned(expX_i, MEM_DATA_ALIGN);

	// backup storage for Green's functions
	greens_func_t Gu_ref, Gd_ref;
	AllocateGreensFunction(N, &Gu_ref);
	AllocateGreensFunction(N, &Gd_ref);

	// new time step B matrices after the phonon field update
	time_step_matrices_t tsm_u_new;
	time_step_matrices_t tsm_d_new;
	AllocateTimeStepMatrices(N, L, tsm_u->prodBlen, &tsm_u_new);
	AllocateTimeStepMatrices(N, L, tsm_d->prodBlen, &tsm_d_new);

	int nblock_accept = 0;		// number of accepted block updates
	int nblock_reject = 0;		// number of rejected block updates

	int n;
	for (n = 0; n < phonon_params->nblock_updates; n++)
	{
		// randomly select a lattice site
		int i = (int)(Random_GetBoundedUint(seed, N));
		int o = i / Ncell;
		// ignore sites without phonon coupling
		if (phonon_params->g[o] == 0)
		{
			continue;
		}

		// backup X_{i,l} and corresponding exponential for all time slices
		int l;
		for (l = 0; l < L; l++)
		{
			   X_i[l] =    X[i + l*N];
			expX_i[l] = expX[i + l*N];
		}

		// backup Green's functions
		CopyGreensFunction(Gu, &Gu_ref);
		CopyGreensFunction(Gd, &Gd_ref);

		// suggest a simultaneous shift of X_{i,l} for all 'l'
		const double dx = (Random_GetUniform(seed) - 0.5) * phonon_params->box_width;

		// calculate change of the phonon (lattice) energy
		double dEph = 0;
		for (l = 0; l < L; l++)
		{
			dEph += (dx + 2*X[i + l*N]);
		}
		dEph *= 0.5*square(phonon_params->omega[o]) * dx;

		// actually shift phonon field entries
		for (l = 0; l < L; l++)
		{
			   X[i + l*N] += dx;
			expX[i + l*N] = exp(-dt*phonon_params->g[o] * X[i + l*N]);
		}

		// calculate new time step matrices
		#pragma omp parallel sections
		{
			#pragma omp section
			InitPhononTimeStepMatrices(kinetic, stratonovich_params->expVu, s, expX, &tsm_u_new);
			#pragma omp section
			InitPhononTimeStepMatrices(kinetic, stratonovich_params->expVd, s, expX, &tsm_d_new);
		}

		// calculate new Green's functions
		#pragma omp parallel sections
		{
			#pragma omp section
			GreenConstruct(&tsm_u_new, 0, Gu);
			#pragma omp section
			GreenConstruct(&tsm_d_new, 0, Gd);
		}

		// decide whether block update is accepted; note that det(G) = 1/det(M)
		if (Random_GetUniform(seed) < exp((Gu_ref.logdet + Gd_ref.logdet) - (Gu->logdet + Gd->logdet) - dt * dEph))
		{
			nblock_accept++;

			// copy new time step matrices
			CopyTimeStepMatrices(&tsm_u_new, tsm_u);
			CopyTimeStepMatrices(&tsm_d_new, tsm_d);
		}
		else
		{
			nblock_reject++;

			// undo changes
			int l;
			for (l = 0; l < L; l++)
			{
				   X[i + l*N] =    X_i[l];
				expX[i + l*N] = expX_i[l];
			}
			CopyGreensFunction(&Gu_ref, Gu);
			CopyGreensFunction(&Gd_ref, Gd);
		}
	}

	// clean up
	DeleteTimeStepMatrices(&tsm_d_new);
	DeleteTimeStepMatrices(&tsm_u_new);
	DeleteGreensFunction(&Gd_ref);
	DeleteGreensFunction(&Gu_ref);
	MKL_free(expX_i);
	MKL_free(X_i);
}


//________________________________________________________________________________________________________________________
///
/// \brief Perform a Determinant Quantum Monte Carlo (DQMC) iteration, taking phonons into account
///
/// \param dt                   imaginary-time step
/// \param kinetic              matrix exponential of kinetic energy operator
/// \param noHS                 set to true to skip updating the Hubbard-Stratonovich field
/// \param stratonovich_params  precomputed Hubbard-Stratonovich parameters
/// \param phonon_params        phonon field parameters
/// \param nwraps               number of "time slice wraps" before recomputing the Green's function
/// \param seed                 random number "seed" structure, will be updated during function call
/// \param s                    Hubbard-Stratonovich field; will be updated
/// \param X                    phonon field; will be updated
/// \param expX                 entrywise exponential of the phonon field: exp(-dt g X)
/// \param tsm_u                time step matrices for the spin-up Green's function, will be updated
/// \param tsm_d                time step matrices for the spin-down Green's function, will be updated
/// \param Gu                   spin-up   Green's function, must have been computed on input and will be updated
/// \param Gd                   spin-down Green's function, must have been computed on input and will be updated
/// \param neqlt                perform an equal time measurement every 'neqlt' time slices
/// \param meas_data            basic measurement data
/// \param meas_data_phonon     phonon measurement data
///
void DQMCPhononIteration(const double dt, const kinetic_t *restrict kinetic, const int noHS, const stratonovich_params_t *restrict stratonovich_params, const phonon_params_t *restrict phonon_params,
	const int nwraps, randseed_t *restrict seed, spin_field_t *restrict s, double *restrict X, double *restrict expX,
	time_step_matrices_t *restrict tsm_u, time_step_matrices_t *restrict tsm_d, greens_func_t *restrict Gu, greens_func_t *restrict Gd,
	const int neqlt, measurement_data_t *restrict meas_data, measurement_data_phonon_t *restrict meas_data_phonon)
{
	Profile_Begin("DQMCIter");
	__assume_aligned(s, MEM_DATA_ALIGN);

	// dimension consistency checks
	assert(tsm_u->N == kinetic->Ncell * kinetic->Norb);
	assert(tsm_u->N == tsm_d->N);
	assert(tsm_u->N == Gu->N);
	assert(tsm_u->N == Gd->N);
	assert(tsm_u->L == tsm_d->L);
	const int N = tsm_u->N;
	const int Ncell = kinetic->Ncell;
	const int L = tsm_u->L;

	assert(tsm_u->prodBlen == tsm_d->prodBlen);
	assert(nwraps % tsm_u->prodBlen == 0);	// must be a multiple of 'prodBlen'

	// store Green's functions before recomputing them to estimate error
	#if defined(DEBUG) | defined(_DEBUG)
	greens_func_t Gu_old, Gd_old;
	AllocateGreensFunction(N, &Gu_old);
	AllocateGreensFunction(N, &Gd_old);
	__assume_aligned(Gu_old.mat, MEM_DATA_ALIGN);
	__assume_aligned(Gd_old.mat, MEM_DATA_ALIGN);
	#endif

	// pre-compute 1/dt^2
	const double inv_dt_sq = 1.0 / square(dt);

	// random shuffle of lattice cells and orbitals
	int *orb_cell_order = MKL_malloc(N * sizeof(int), MEM_DATA_ALIGN);
	 __assume_aligned(orb_cell_order, MEM_DATA_ALIGN);

	// iterate over time slices
	int l;
	for (l = 0; l < L; l++)
	{
		Profile_Begin("DQMCIter_Wraps");
		#pragma omp parallel sections
		{
			#pragma omp section
			GreenTimeSliceWrap(N, tsm_u->B[l], tsm_u->invB[l], Gu->mat);
			#pragma omp section
			GreenTimeSliceWrap(N, tsm_d->B[l], tsm_d->invB[l], Gd->mat);
		}
		Profile_End("DQMCIter_Wraps");

		if (!noHS)
		{
			// iterate over lattice sites randomly, updating the Hubbard-Stratonovich field
			Profile_Begin("DQMCIter_HSUpdate");
			Random_Shuffle(seed, N, orb_cell_order);
			int j;
			for (j = 0; j < N; j++)
			{
				const int i = orb_cell_order[j];
				const int o = i / Ncell;	// orbital index
				assert(0 <= i && i < N);
				assert(0 <= o && o < kinetic->Norb);

				// Eq. (13)
				// suggest flipping s_{i,l}
				const double du = 1 + (1 - Gu->mat[i + i*N]) * stratonovich_params->delta[  s[i + l*N]][o];
				const double dd = 1 + (1 - Gd->mat[i + i*N]) * stratonovich_params->delta[1-s[i + l*N]][o];
				if (Random_GetUniform(seed) < fabs(du*dd))
				{
					// Eq. (15)
					#pragma omp parallel sections
					{
						#pragma omp section
						GreenShermanMorrisonUpdate(stratonovich_params->delta[  s[i + l*N]][o], N, i, Gu->mat);
						#pragma omp section
						GreenShermanMorrisonUpdate(stratonovich_params->delta[1-s[i + l*N]][o], N, i, Gd->mat);
					}
					// correspondingly update determinants
					Gu->logdet -= log(fabs(du));
					Gd->logdet -= log(fabs(dd));
					if (du < 0) { Gu->sgndet = -Gu->sgndet; }
					if (dd < 0) { Gd->sgndet = -Gd->sgndet; }

					// actually flip spin
					s[i + l*N] = 1 - s[i + l*N];
				}
			}
			Profile_End("DQMCIter_HSUpdate");
		}

		// next and previous time slices; required for the phonon field with periodic boundary conditions
		const int l_next = (l + 1    ) % L;
		const int l_prev = (l + L - 1) % L;

		// iterate over lattice sites, updating the phonon field
		Profile_Begin("DQMCIter_XUpdate");
		Random_Shuffle(seed, N, orb_cell_order);
		int j;
		for (j = 0; j < N; j++)
		{
			const int i = orb_cell_order[j];
			const int o = i / Ncell;	// orbital index
			assert(0 <= i && i < N);
			assert(0 <= o && o < kinetic->Norb);

			// skip orbitals without phonon coupling
			if (phonon_params->g[o] == 0)
			{
				continue;
			}

			// suggest a shift of X_{i,l}
			const double dx = (Random_GetUniform(seed) - 0.5) * phonon_params->box_width;

			// Eq. (19) in PRB 87, 235133 (2013)
			const double delta = expm1(-dt*phonon_params->g[o] * dx);
			const double du = 1 + (1 - Gu->mat[i + i*N]) * delta;
			const double dd = 1 + (1 - Gd->mat[i + i*N]) * delta;

			// change of the phonon (lattice) energy
			const double dEph = dx * (0.5*square(phonon_params->omega[o])*(dx + 2*X[i + l*N]) + inv_dt_sq*(dx - (X[i + l_next*N] - 2*X[i + l*N] + X[i + l_prev*N])));

			if (Random_GetUniform(seed) < fabs(du*dd) * exp(-dt * dEph))
			{
				// Eq. (15)
				#pragma omp parallel sections
				{
					#pragma omp section
					GreenShermanMorrisonUpdate(delta, N, i, Gu->mat);
					#pragma omp section
					GreenShermanMorrisonUpdate(delta, N, i, Gd->mat);
				}
				// correspondingly update determinants
				Gu->logdet -= log(fabs(du));
				Gd->logdet -= log(fabs(dd));
				if (du < 0) { Gu->sgndet = -Gu->sgndet; }
				if (dd < 0) { Gd->sgndet = -Gd->sgndet; }

				// actually update the phonon field
				   X[i + l*N] += dx;
				expX[i + l*N] = exp(-dt*phonon_params->g[o] * X[i + l*N]);
			}
		}
		Profile_End("DQMCIter_XUpdate");

		// re-compute corresponding B matrices
		Profile_Begin("DQMCIter_Brecomp");
		#pragma omp parallel sections
		{
			#pragma omp section
			UpdatePhononTimeStepMatrices(kinetic, stratonovich_params->expVu, s, expX, l, tsm_u);
			#pragma omp section
			UpdatePhononTimeStepMatrices(kinetic, stratonovich_params->expVd, s, expX, l, tsm_d);
		}
		Profile_End("DQMCIter_Brecomp");

		// recompute Green's function after several time slice "wraps"
		if ((l + 1) % nwraps == 0)
		{
			Profile_Begin("DQMCIter_Grecomp");
			// store current Green's function matrices to compare with newly constructed ones
			#if defined(DEBUG) | defined(_DEBUG)
			CopyGreensFunction(Gu, &Gu_old);
			CopyGreensFunction(Gd, &Gd_old);
			#endif

			#pragma omp parallel sections
			{
				#pragma omp section
				GreenConstruct(tsm_u, (l + 1) % L, Gu);
				#pragma omp section
				GreenConstruct(tsm_d, (l + 1) % L, Gd);
			}

			#if defined(DEBUG) | defined(_DEBUG)
			// deviation of matrix entries
			double err_u = UniformDistance(N*N, Gu_old.mat, Gu->mat);
			double err_d = UniformDistance(N*N, Gd_old.mat, Gd->mat);
			if (err_u > 1e-8*N || err_d > 1e-8*N) {
				duprintf("Warning: after calling 'GreenConstruct()', largest entrywise distance between previous and new Green's functions is %g (up) and %g (down).\n", err_u, err_d);
			}
			// deviation of matrix determinants
			err_u = fabs(Gu_old.logdet - Gu->logdet);
			err_d = fabs(Gd_old.logdet - Gd->logdet);
			if (err_u > 1e-8 || err_d > 1e-8) {
				duprintf("Warning: after calling 'GreenConstruct()', largest distance between logarithm of previous and new Green's function determinants is %g (up) and %g (down).\n", err_u, err_d);
			}
			if (Gu_old.sgndet != Gu->sgndet || Gd_old.sgndet != Gd->sgndet) {
				duprintf("Warning: after calling 'GreenConstruct()', determinant sign has changed.\n");
			}
			#endif
			Profile_End("DQMCIter_Grecomp");
		}

		if (neqlt > 0 && (l + 1) % neqlt == 0)
		{
			// accumulate equal time "measurement" data
			Profile_Begin("DQMCIter_AccumulateEqMeas");
			AccumulateMeasurement(Gu, Gd, meas_data);
			Profile_End("DQMCIter_AccumulateEqMeas");

			// accumulate phonon data
			Profile_Begin("DQMCIter_AccumulatePhonon");
			AccumulatePhononData(Gu, Gd, X, meas_data_phonon);
			Profile_End("DQMCIter_AccumulatePhonon");
		}
	}

	// perform block updates
	Profile_Begin("DQMCIter_PhononBlock");
	PhononBlockUpdates(dt, kinetic, stratonovich_params, phonon_params, seed, s, X, expX, tsm_u, tsm_d, Gu, Gd);
	Profile_End("DQMCIter_PhononBlock");
	// clean up
	MKL_free(orb_cell_order);
	#if defined(DEBUG) | defined(_DEBUG)
	DeleteGreensFunction(&Gd_old);
	DeleteGreensFunction(&Gu_old);
	#endif
	Profile_End("DQMCIter");
}


//________________________________________________________________________________________________________________________
///
/// \brief Perform a Determinant Quantum Monte Carlo (DQMC) simulation
///
/// \param params               simulation parameters
/// \param meas_data            measurement data structure for accumulating measurements
/// \param meas_data_uneqlt     unequal time measurement data structure
/// \param meas_data_phonon     phonon measurement data structure
/// \param iteration            pointer to iteration counter (i.e. number of iterations completed from previous runs)
/// \param seed                 random number generator seed
/// \param s                    Hubbard-Stratonovich field; will be updated
/// \param X                    phonon field; only accessed if params->use_phonons is true
/// \param expX                 entrywise exponential of the phonon field: exp(-dt g X)
///
void DQMCSimulation(const sim_params_t *restrict params,
	measurement_data_t *restrict meas_data, measurement_data_unequal_time_t *restrict meas_data_uneqlt, measurement_data_phonon_t *restrict meas_data_phonon,
	int *restrict iteration, randseed_t *restrict seed, spin_field_t *restrict s, double *restrict X, double *restrict expX)
{
	const int Norb  = params->Norb;
	const int Ncell = params->Nx * params->Ny;
	const int N     = Norb * Ncell;

	// get the time (in ticks) for automatic stopping and checkpointing
	const uint64_t t_end = GetTicks() + params->max_time * GetTickRes();

	// calculate matrix exponential of the kinetic nearest neighbor hopping matrix
	kinetic_t kinetic;
	RectangularKineticExponential(params, &kinetic);

	// pre-calculate some stuff related to the Hubbard-Stratonovich field, for every orbital
	int noHS = 1; // flag to disable H-S updates if all U == 0
	int o;
	for (o = 0; o < Norb; o++)
	{
		if (params->U[o] != 0)
		{
			noHS = 0;
			break;
		}
	}
	stratonovich_params_t stratonovich_params;
	FillStratonovichParameters(Norb, params->U, params->dt, &stratonovich_params);

	// time step matrices
	time_step_matrices_t tsm_u;
	time_step_matrices_t tsm_d;
	AllocateTimeStepMatrices(N, params->L, params->prodBlen, &tsm_u);
	AllocateTimeStepMatrices(N, params->L, params->prodBlen, &tsm_d);
	if (params->use_phonons)
	{
		#pragma omp parallel sections
		{
			#pragma omp section
			InitPhononTimeStepMatrices(&kinetic, stratonovich_params.expVu, s, expX, &tsm_u);
			#pragma omp section
			InitPhononTimeStepMatrices(&kinetic, stratonovich_params.expVd, s, expX, &tsm_d);
		}
	}
	else
	{
		#pragma omp parallel sections
		{
			#pragma omp section
			InitTimeStepMatrices(&kinetic, stratonovich_params.expVu, s, &tsm_u);
			#pragma omp section
			InitTimeStepMatrices(&kinetic, stratonovich_params.expVd, s, &tsm_d);
		}
	}

	// allocate Green's functions
	greens_func_t Gu, Gd;
	AllocateGreensFunction(N, &Gu);
	AllocateGreensFunction(N, &Gd);

	// construct initial Green's functions
	#pragma omp parallel sections
	{
		#pragma omp section
		GreenConstruct(&tsm_u, 0, &Gu);
		#pragma omp section
		GreenConstruct(&tsm_d, 0, &Gd);
	}

	// register a signal handler to stop simulation on SIGINT
	StopOnSIGINT();

	// first call to UpdateProgress to indicate completion of initialization phase
	UpdateProgress();

	// perform dqmc iterations
	duprintf("Starting DQMC iterations...\n");
	for (; *iteration < params->nequil + params->nsampl; (*iteration)++)
	{
		if (params->max_time > 0 && GetTicks() >= t_end)
		{
			duprintf("Reached time limit of %d seconds.\n", params->max_time);
			stopped = 1;
		}

		if (stopped == 1) // either the above happened or SIGINT was received
		{
			duprintf("Stopping DQMC iterations early.\n");
			break;
		}

		// set neqlt to 0 in equilibration stage so that no measurements are made
		const int neqlt = (*iteration >= params->nequil) ? params->neqlt : 0;

		if (params->use_phonons)
		{
			DQMCPhononIteration(params->dt, &kinetic, noHS, &stratonovich_params, &params->phonon_params, params->nwraps, seed, s, X, expX, &tsm_u, &tsm_d, &Gu, &Gd, neqlt, meas_data, meas_data_phonon);
		}
		else
		{
			DQMCIteration(&kinetic, &stratonovich_params, params->nwraps, seed, s, &tsm_u, &tsm_d, &Gu, &Gd, neqlt, meas_data);
		}

		// accumulate unequal time "measurement" data
		if (*iteration >= params->nequil && params->nuneqlt > 0 && (*iteration % params->nuneqlt) == 0)
		{
			Profile_Begin("DQMCSim_AccumulateUneqMeas");
			AccumulateUnequalTimeMeasurement((double)(Gu.sgndet * Gd.sgndet), &tsm_u, &tsm_d, meas_data_uneqlt);
			Profile_End("DQMCSim_AccumulateUneqMeas");
		}

		UpdateProgress();
	}

	// clean up
	DeleteGreensFunction(&Gd);
	DeleteGreensFunction(&Gu);
	DeleteTimeStepMatrices(&tsm_d);
	DeleteTimeStepMatrices(&tsm_u);
	DeleteStratonovichParameters(&stratonovich_params);
	DeleteKineticExponential(&kinetic);
}
