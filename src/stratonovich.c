#include "stratonovich.h"
#include "util.h"
#include <math.h>


//________________________________________________________________________________________________________________________
///
/// \brief Fill and allocate the Stratonovich parameter structure
///
void FillStratonovichParameters(const int Norb, const double *U, const double dt, stratonovich_params_t *params)
{
	params->Norb = Norb;

	params->expVu[0] = algn_malloc(Norb * sizeof(double));
	params->expVu[1] = algn_malloc(Norb * sizeof(double));
	params->expVd[0] = algn_malloc(Norb * sizeof(double));
	params->expVd[1] = algn_malloc(Norb * sizeof(double));
	params->delta[0] = algn_malloc(Norb * sizeof(double));
	params->delta[1] = algn_malloc(Norb * sizeof(double));

	int o;
	for (o = 0; o < Norb; o++)
	{
		const double x = 0.5*dt*U[o];
		const double exp_x = exp(x);
		const double y = expm1(x) * (exp_x + 1);
		const double sqrt_y = sqrt(y);

		// pre-compute exponentials of lambda = acosh(exp(dt*U/2))
		const double exp_pos_lambda = exp_x + sqrt_y;
		const double exp_neg_lambda = 1.0 / exp_pos_lambda;

		// spin-up Green's function
		params->expVu[0][o] = exp_neg_lambda;
		params->expVu[1][o] = exp_pos_lambda;

		// spin-down Green's function
		params->expVd[0][o] = exp_pos_lambda;
		params->expVd[1][o] = exp_neg_lambda;

		// delta parameter exp(2*lambda) - 1 and exp(-2*lambda) - 1
		params->delta[0][o] = 2 * (y + exp_x*sqrt_y);
		params->delta[1][o] = 1.0 / (params->delta[0][o] + 1) - 1;
	}
}


//________________________________________________________________________________________________________________________
///
/// \brief Free memory of the Stratonovich parameter structure
///
void DeleteStratonovichParameters(stratonovich_params_t *params)
{
	algn_free(params->delta[1]);
	algn_free(params->delta[0]);
	algn_free(params->expVd[1]);
	algn_free(params->expVd[0]);
	algn_free(params->expVu[1]);
	algn_free(params->expVu[0]);
}
