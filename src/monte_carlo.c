#include "monte_carlo.h"
#include "greens_func.h"
#include "util.h"
#include "dupio.h"
#include "profiler.h"
#include <mkl.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>


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
/// \param Gu                   spin-up Green's function, will be updated
/// \param Gd                   spin-down Green's function, will be updated
///
void DQMCIteration(const kinetic_t *restrict kinetic, const stratonovich_params_t *restrict stratonovich_params, const int nwraps,
	randseed_t *restrict seed, spin_field_t *restrict s, time_step_matrices_t *restrict tsm_u, time_step_matrices_t *restrict tsm_d, greens_func_t *restrict Gu, greens_func_t *restrict Gd)
{
	PROFILE_BEGIN(DQMCIter);
	__assume_aligned(s, MEM_DATA_ALIGN);

	// dimension consistency checks
	assert(tsm_u->N == kinetic->N);
	assert(tsm_u->N == tsm_d->N);
	assert(tsm_u->L == tsm_d->L);
	const int N = kinetic->N;

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

	// random shuffle of lattice sites
	int *site_order = MKL_malloc(N * sizeof(int), MEM_DATA_ALIGN);
	__assume_aligned(site_order, MEM_DATA_ALIGN);

	// iterate over time slices
	int l;
	for (l = 0; l < tsm_u->L; l++)
	{
		// recompute Green's function after several time slice "wraps"
		if ((l % nwraps) == 0)
		{
			PROFILE_BEGIN(DQMCIter_Grecomp);
			// store current Green's function matrices to compare with newly constructed ones
			#if defined(DEBUG) | defined(_DEBUG)
			CopyGreensFunction(Gu, &Gu_old);
			CopyGreensFunction(Gd, &Gd_old);
			#endif

			#pragma omp parallel sections
			{
				#pragma omp section
				GreenConstruct(tsm_u, l, Gu);
				#pragma omp section
				GreenConstruct(tsm_d, l, Gd);
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
			PROFILE_END(DQMCIter_Grecomp);
		}

		PROFILE_BEGIN(DQMCIter_Wraps);
		#pragma omp parallel sections
		{
			#pragma omp section
			GreenTimeSliceWrap(N, tsm_u->B[l], tsm_u->invB[l], Gu->mat);
			#pragma omp section
			GreenTimeSliceWrap(N, tsm_d->B[l], tsm_d->invB[l], Gd->mat);
		}
		PROFILE_END(DQMCIter_Wraps);

		// iterate over lattice sites in random order, updating the Hubbard-Stratonovich field
		Random_Shuffle(seed, N, site_order);
		int j;
		for (j = 0; j < N; j++)
		{
			PROFILE_BEGIN(DQMCIter_SiteUpdate);
			int i = site_order[j];
			assert(0 <= i && i < N);
			// Eq. (13)
			// suggest flipping s_{i,l}
			const double du = 1 + (1 - Gu->mat[i + i*N]) * stratonovich_params->delta[  s[i + l*N]];
			const double dd = 1 + (1 - Gd->mat[i + i*N]) * stratonovich_params->delta[1-s[i + l*N]];
			if (Random_GetUniform(seed) < fabs(du*dd))
			{
				// Eq. (15)
				#pragma omp parallel sections
				{
					#pragma omp section
					GreenShermanMorrisonUpdate(stratonovich_params->delta[  s[i + l*N]], N, i, Gu->mat);
					#pragma omp section
					GreenShermanMorrisonUpdate(stratonovich_params->delta[1-s[i + l*N]], N, i, Gd->mat);
				}
				// correspondingly update determinants
				Gu->logdet -= log(fabs(du));
				Gd->logdet -= log(fabs(dd));
				if (du < 0) { Gu->sgndet = -Gu->sgndet; }
				if (dd < 0) { Gd->sgndet = -Gd->sgndet; }

				// actually flip spin of Hubbard-Stratonovich field entry
				s[i + l*N] = 1 - s[i + l*N];
			}
			PROFILE_END(DQMCIter_SiteUpdate);
		}

		PROFILE_BEGIN(DQMCIter_Brecomp);
		// re-compute corresponding B matrices
		#pragma omp parallel sections
		{
			#pragma omp section
			UpdateTimeStepMatrices(kinetic, stratonovich_params->expVu, s, l, tsm_u);
			#pragma omp section
			UpdateTimeStepMatrices(kinetic, stratonovich_params->expVd, s, l, tsm_d);
		}
		PROFILE_END(DQMCIter_Brecomp);
	}

	// clean up
	MKL_free(site_order);
	#if defined(DEBUG) | defined(_DEBUG)
	DeleteGreensFunction(&Gd_old);
	DeleteGreensFunction(&Gu_old);
	#endif
	PROFILE_END(DQMCIter);
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
	assert(tsm_u->N == kinetic->N);
	assert(tsm_u->N == tsm_d->N);
	assert(tsm_u->L == tsm_d->L);
	const int N = kinetic->N;
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
		int i = (int)(Random_GetUint(seed) % (uint64_t)N);

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
		dEph *= 0.5*square(phonon_params->omega) * dx;

		// actually shift phonon field entries
		for (l = 0; l < L; l++)
		{
			   X[i + l*N] += dx;
			expX[i + l*N] = exp(-dt*phonon_params->g * X[i + l*N]);
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
/// \param stratonovich_params  precomputed Hubbard-Stratonovich parameters
/// \param phonon_params        phonon field parameters
/// \param nwraps               number of "time slice wraps" before recomputing the Green's function
/// \param seed                 random number "seed" structure, will be updated during function call
/// \param s                    Hubbard-Stratonovich field; will be updated
/// \param X                    phonon field; will be updated
/// \param expX                 entrywise exponential of the phonon field: exp(-dt g X)
/// \param tsm_u                time step matrices for the spin-up Green's function, will be updated
/// \param tsm_d                time step matrices for the spin-down Green's function, will be updated
/// \param Gu                   spin-up Green's function, will be updated
/// \param Gd                   spin-down Green's function, will be updated
///
void DQMCPhononIteration(const double dt, const kinetic_t *restrict kinetic, const stratonovich_params_t *restrict stratonovich_params, const phonon_params_t *restrict phonon_params,
	const int nwraps, randseed_t *restrict seed, spin_field_t *restrict s, double *restrict X, double *restrict expX,
	time_step_matrices_t *restrict tsm_u, time_step_matrices_t *restrict tsm_d, greens_func_t *restrict Gu, greens_func_t *restrict Gd)
{
	// dimension consistency checks
	assert(tsm_u->N == kinetic->N);
	assert(tsm_u->N == tsm_d->N);
	assert(tsm_u->N == Gu->N);
	assert(tsm_u->N == Gd->N);
	assert(tsm_u->L == tsm_d->L);
	const int N = kinetic->N;
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

	// pre-compute some constants
	const double omega_ph_sqhalf = 0.5*square(phonon_params->omega);
	const double inv_dt_sq = 1.0 / square(dt);

	// random shuffle of lattice sites
	int *site_order = MKL_malloc(N * sizeof(int), MEM_DATA_ALIGN);
	 __assume_aligned(site_order, MEM_DATA_ALIGN);

	// iterate over time slices
	int l;
	for (l = 0; l < L; l++)
	{
		// recompute Green's function after several time slice "wraps"
		if ((l % nwraps) == 0)
		{
			// store current Green's function matrices to compare with newly constructed ones
			#if defined(DEBUG) | defined(_DEBUG)
			CopyGreensFunction(Gu, &Gu_old);
			CopyGreensFunction(Gd, &Gd_old);
			#endif

			#pragma omp parallel sections
			{
				#pragma omp section
				GreenConstruct(tsm_u, l, Gu);
				#pragma omp section
				GreenConstruct(tsm_d, l, Gd);
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
		}

		#pragma omp parallel sections
		{
			#pragma omp section
			GreenTimeSliceWrap(N, tsm_u->B[l], tsm_u->invB[l], Gu->mat);
			#pragma omp section
			GreenTimeSliceWrap(N, tsm_d->B[l], tsm_d->invB[l], Gd->mat);
		}

		// iterate over lattice sites randomly, updating the Hubbard-Stratonovich field
		Random_Shuffle(seed, N, site_order);
		int j;
		for (j = 0; j < N; j++)
		{
			int i = site_order[j];
			assert(0 <= i && i < N);
			// Eq. (13)
			// suggest flipping s_{i,l}
			const double du = 1 + (1 - Gu->mat[i + i*N]) * stratonovich_params->delta[  s[i + l*N]];
			const double dd = 1 + (1 - Gd->mat[i + i*N]) * stratonovich_params->delta[1-s[i + l*N]];
			if (Random_GetUniform(seed) < fabs(du*dd))
			{
				// Eq. (15)
				#pragma omp parallel sections
				{
					#pragma omp section
					GreenShermanMorrisonUpdate(stratonovich_params->delta[  s[i + l*N]], N, i, Gu->mat);
					#pragma omp section
					GreenShermanMorrisonUpdate(stratonovich_params->delta[1-s[i + l*N]], N, i, Gd->mat);
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

		// next and previous time slices; required for the phonon field with periodic boundary conditions
		const int l_next = (l + 1    ) % L;
		const int l_prev = (l + L - 1) % L;

		// iterate over lattice sites, updating the phonon field
		Random_Shuffle(seed, N, site_order);
		for (j = 0; j < N; j++)
		{
			int i = site_order[j];
			assert(0 <= i && i < N);
			// suggest a shift of X_{i,l}
			const double dx = (Random_GetUniform(seed) - 0.5) * phonon_params->box_width;

			// Eq. (19) in PRB 87, 235133 (2013)
			const double delta = expm1(-dt*phonon_params->g * dx);
			const double du = 1 + (1 - Gu->mat[i + i*N]) * delta;
			const double dd = 1 + (1 - Gd->mat[i + i*N]) * delta;

			// change of the phonon (lattice) energy
			const double dEph = dx * (omega_ph_sqhalf*(dx + 2*X[i + l*N]) + inv_dt_sq*(dx - (X[i + l_next*N] - 2*X[i + l*N] + X[i + l_prev*N])));

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
				expX[i + l*N] = exp(-dt*phonon_params->g * X[i + l*N]);
			}
		}

		// re-compute corresponding B matrices
		#pragma omp parallel sections
		{
			#pragma omp section
			UpdatePhononTimeStepMatrices(kinetic, stratonovich_params->expVu, s, expX, l, tsm_u);
			#pragma omp section
			UpdatePhononTimeStepMatrices(kinetic, stratonovich_params->expVd, s, expX, l, tsm_d);
		}
	}

	// perform block updates
	PhononBlockUpdates(dt, kinetic, stratonovich_params, phonon_params, seed, s, X, expX, tsm_u, tsm_d, Gu, Gd);

	// clean up
	MKL_free(site_order);
	#if defined(DEBUG) | defined(_DEBUG)
	DeleteGreensFunction(&Gd_old);
	DeleteGreensFunction(&Gu_old);
	#endif
}


//________________________________________________________________________________________________________________________
///
/// \brief Perform a Determinant Quantum Monte Carlo (DQMC) simulation
///
/// \param U                Coulomb coupling constant in the Hubbard hamiltonian
/// \param dt               imaginary-time step
/// \param L                number of time steps
/// \param kinetic          matrix exponential of kinetic energy operator
/// \param prodBlen         largest number of B_l matrices multiplied together before performing a QR decomposition
/// \param nwraps           number of "time slice wraps" before recomputing the Green's function
/// \param nequil           number of equilibration iterations
/// \param nsampl           number of sample iterations
/// \param seed             random number "seed" structure, will be updated during function call
/// \param meas_data        measurement data structure for accumulating measurements
///
void DQMCSimulation(const double U, const double dt, const int L, const kinetic_t *restrict kinetic, const int prodBlen, const int nwraps,
	const int nequil, const int nsampl, randseed_t *restrict seed, measurement_data_t *restrict meas_data)
{
	int i;
	const int N = kinetic->N;

	// Hubbard-Stratonovich parameters
	stratonovich_params_t stratonovich_params;
	FillStratonovichParameters(U, dt, &stratonovich_params);

	// random initial Hubbard-Stratonovich field
	spin_field_t *s = (spin_field_t *)MKL_malloc(L*N * sizeof(spin_field_t), MEM_DATA_ALIGN);
	for (i = 0; i < L*N; i++)
	{
		s[i] = (Random_GetUniform(seed) < 0.5 ? 0 : 1);
	}

	// time step matrices
	time_step_matrices_t tsm_u;
	time_step_matrices_t tsm_d;
	AllocateTimeStepMatrices(N, L, prodBlen, &tsm_u);
	AllocateTimeStepMatrices(N, L, prodBlen, &tsm_d);
	#pragma omp parallel sections
	{
		#pragma omp section
		InitTimeStepMatrices(kinetic, stratonovich_params.expVu, s, &tsm_u);
		#pragma omp section
		InitTimeStepMatrices(kinetic, stratonovich_params.expVd, s, &tsm_d);
	}

	// allocate and construct initial Green's functions
	greens_func_t Gu, Gd;
	AllocateGreensFunction(N, &Gu);
	AllocateGreensFunction(N, &Gd);
	#pragma omp parallel sections
	{
		#pragma omp section
		GreenConstruct(&tsm_u, 0, &Gu);
		#pragma omp section
		GreenConstruct(&tsm_d, 0, &Gd);
	}

	// perform equilibration
	duprintf("Starting equilibration iterations...\n");
	for (i = 0; i < nequil; i++)
	{
		DQMCIteration(kinetic, &stratonovich_params, nwraps, seed, s, &tsm_u, &tsm_d, &Gu, &Gd);
	}

	// perform measurement iterations
	duprintf("Starting measurement iterations...\n");
	for (i = 0; i < nsampl; i++)
	{
		DQMCIteration(kinetic, &stratonovich_params, nwraps, seed, s, &tsm_u, &tsm_d, &Gu, &Gd);

		// accumulate "measurement" data
		PROFILE_BEGIN(DQMCSim_AccumulateMeasurements);
		AccumulateMeasurement(&Gu, &Gd, meas_data);
		PROFILE_END(DQMCSim_AccumulateMeasurements);
	}

	// clean up
	DeleteGreensFunction(&Gd);
	DeleteGreensFunction(&Gu);
	DeleteTimeStepMatrices(&tsm_d);
	DeleteTimeStepMatrices(&tsm_u);
	MKL_free(s);
}


//________________________________________________________________________________________________________________________
///
/// \brief Perform a Determinant Quantum Monte Carlo (DQMC) simulation, taking phonons into account
///
/// \param U                Coulomb coupling constant in the Hubbard hamiltonian
/// \param dt               imaginary-time step
/// \param L                number of time steps
/// \param kinetic          matrix exponential of kinetic energy operator
/// \param prodBlen         largest number of B_l matrices multiplied together before performing a QR decomposition
/// \param nwraps           number of "time slice wraps" before recomputing the Green's function
/// \param phonon_params    phonon field parameters
/// \param nequil           number of equilibration iterations
/// \param nsampl           number of sample iterations
/// \param seed             random number "seed" structure, will be updated during function call
/// \param meas_data        measurement data structure for accumulating measurements
///
void DQMCPhononSimulation(const double U, const double dt, const int L, const kinetic_t *restrict kinetic, const int prodBlen, const int nwraps,
	const phonon_params_t *restrict phonon_params, const int nequil, const int nsampl, randseed_t *restrict seed, measurement_data_t *restrict meas_data)
{
	int i;
	const int N = kinetic->N;

	// Hubbard-Stratonovich parameters
	stratonovich_params_t stratonovich_params;
	FillStratonovichParameters(U, dt, &stratonovich_params);

	// random initial Hubbard-Stratonovich field
	spin_field_t *s = (spin_field_t *)MKL_malloc(L*N * sizeof(spin_field_t), MEM_DATA_ALIGN);
	for (i = 0; i < L*N; i++)
	{
		s[i] = (Random_GetUniform(seed) < 0.5 ? 0 : 1);
	}

	// random initial phonon field
	double *X    = (double *)MKL_malloc(L*N * sizeof(double), MEM_DATA_ALIGN);
	double *expX = (double *)MKL_malloc(L*N * sizeof(double), MEM_DATA_ALIGN);
	for (i = 0; i < L*N; i++)
	{
		X[i] = (Random_GetUniform(seed) - 0.5) * phonon_params->box_width;
		expX[i] = exp(-dt*phonon_params->g * X[i]);
	}

	// time step matrices
	time_step_matrices_t tsm_u;
	time_step_matrices_t tsm_d;
	AllocateTimeStepMatrices(N, L, prodBlen, &tsm_u);
	AllocateTimeStepMatrices(N, L, prodBlen, &tsm_d);
	#pragma omp parallel sections
	{
		#pragma omp section
		InitPhononTimeStepMatrices(kinetic, stratonovich_params.expVu, s, expX, &tsm_u);
		#pragma omp section
		InitPhononTimeStepMatrices(kinetic, stratonovich_params.expVd, s, expX, &tsm_d);
	}

	// allocate and construct initial Green's functions
	greens_func_t Gu, Gd;
	AllocateGreensFunction(N, &Gu);
	AllocateGreensFunction(N, &Gd);
	#pragma omp parallel sections
	{
		#pragma omp section
		GreenConstruct(&tsm_u, 0, &Gu);
		#pragma omp section
		GreenConstruct(&tsm_d, 0, &Gd);
	}

	// perform equilibration
	duprintf("Starting equilibration iterations (including phonons)...\n");
	for (i = 0; i < nequil; i++)
	{
		DQMCPhononIteration(dt, kinetic, &stratonovich_params, phonon_params, nwraps, seed, s, X, expX, &tsm_u, &tsm_d, &Gu, &Gd);
	}

	// perform measurement iterations
	duprintf("Starting measurement iterations (including phonons)...\n");
	for (i = 0; i < nsampl; i++)
	{
		DQMCPhononIteration(dt, kinetic, &stratonovich_params, phonon_params, nwraps, seed, s, X, expX, &tsm_u, &tsm_d, &Gu, &Gd);

		// accumulate "measurement" data
		AccumulateMeasurement(&Gu, &Gd, meas_data);
	}

	// clean up
	DeleteGreensFunction(&Gd);
	DeleteGreensFunction(&Gu);
	DeleteTimeStepMatrices(&tsm_d);
	DeleteTimeStepMatrices(&tsm_u);
	MKL_free(expX);
	MKL_free(X);
	MKL_free(s);
}
