#include "stratonovich.h"
#include <math.h>


//________________________________________________________________________________________________________________________
///
/// \brief Fill the Stratonovich parameter structure
///
void FillStratonovichParameters(const double U, const double dt, stratonovich_params_t *params)
{
	const double x = 0.5*dt*U;
	const double exp_x = exp(x);
	const double y = expm1(x) * (exp_x + 1);
	const double sqrt_y = sqrt(y);

	// pre-compute exponentials of lambda = acosh(exp(dt*U/2))
	const double exp_pos_lambda = exp_x + sqrt_y;
	const double exp_neg_lambda = 1.0 / exp_pos_lambda;

	// spin-up Green's function
	params->expVu[0] = exp_neg_lambda;
	params->expVu[1] = exp_pos_lambda;

	// spin-down Green's function
	params->expVd[0] = exp_pos_lambda;
	params->expVd[1] = exp_neg_lambda;

	// delta parameter exp(2*lambda) - 1 and exp(-2*lambda) - 1
	params->delta[0] = 2*(y + exp_x*sqrt_y);
	params->delta[1] = 1.0 / (params->delta[0] + 1) - 1;
}
