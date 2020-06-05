/*  actuar: Actuarial Functions and Heavy Tailed Distributions
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the Gumbel distribution. See ../R/Gumbel.R for
 *  details.
 *
 *  We work with the density expressed as
 *
 *    e^(-u) * e^(-e^(-u)) / scale
 *
 *  with u = (x - alpha)/scale.
 *
 *  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dgumbel(double x, double alpha, double scale, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(alpha) || ISNAN(scale))
	return x + alpha + scale;
#endif
    if (!R_FINITE(scale))
	return ACT_D__0;
    if (!R_FINITE(x) && alpha == x)
	return R_NaN;		/* x - alpha is NaN */
    if (scale <= 0)
    {
	if (scale < 0) return R_NaN;
	/* scale == 0 */
	return (x == alpha) ? R_PosInf : ACT_D__0;
    }
    x = (x - alpha) / scale;

    if (!R_FINITE(x))
        return ACT_D__0;

    return ACT_D_exp(-(x + exp(-x) + log(scale)));
}

double pgumbel(double q, double alpha, double scale, int lower_tail,
                   int log_p)
{
#ifdef IEEE_754
    if (ISNAN(q) || ISNAN(alpha) || ISNAN(scale))
	return q + alpha + scale;
#endif
    if (!R_FINITE(q) && alpha == q)
	return R_NaN;		/* q - alpha is NaN */
    if (scale <= 0)
    {
	if (scale < 0) return R_NaN;
	/* scale == 0 : */
	return (q < alpha) ? ACT_DT_0 : ACT_DT_1;
    }
    double p = (q - alpha) / scale;
    if (!R_FINITE(p))
	return (q < alpha) ? ACT_DT_0 : ACT_DT_1;
    q = p;

    return ACT_DT_val(exp(-exp(-q)));
}

double qgumbel(double p, double alpha, double scale, int lower_tail,
                   int log_p)
{
#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(alpha) || ISNAN(scale))
	return p + alpha + scale;
#endif
    if (!R_FINITE(alpha) ||
        !R_FINITE(scale) ||
        scale <= 0.0)
        return R_NaN;;

    ACT_Q_P01_boundaries(p, R_NegInf, R_PosInf);
    p = ACT_DT_qIv(p);

    return alpha - scale * log(-log(p));
}

double rgumbel(double alpha, double scale)
{
    if (!R_FINITE(alpha) ||
        !R_FINITE(scale) ||
        scale <= 0.0)
        return R_NaN;;

    return alpha - scale * log(exp_rand());
}

#define EULER_CNST 0.577215664901532860606512090082

double mgumbel(double order, double alpha, double scale, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(order) || ISNAN(alpha) || ISNAN(scale))
	return order + alpha + scale;
#endif
    if (!R_FINITE(alpha) ||
        !R_FINITE(scale) ||
        !R_FINITE(order) ||
        scale <= 0.0 ||
        order <= 0.0 ||
	order > 2.0)
        return R_NaN;

    if (order == 1.0)
	return alpha + EULER_CNST * scale;
    if (order == 2.0)
	return R_pow_di(M_PI * scale, 2)/6 + R_pow_di(alpha + EULER_CNST * scale, 2);

    return R_NaN;		/* order != 1 or 2 */
}

double mgfgumbel(double t, double alpha, double scale, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(t) || ISNAN(alpha) || ISNAN(scale))
	return t + alpha + scale;
#endif
    if (!R_FINITE(alpha) ||
        !R_FINITE(scale) ||
        scale <= 0.0 ||
        scale * t < 1.0)
        return R_NaN;

    if (t == 0.0)
        return ACT_D__1;

    return ACT_D_exp(alpha * t + lgamma(1 - scale * t));
}
