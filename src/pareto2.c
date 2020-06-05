/*  actuar: Actuarial Functions and Heavy Tailed Distributions
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the Pareto (type) II distribution. See ../R/Pareto2.R for
 *  details.
 *
 *  We work with the density expressed as
 *
 *    shape * u^shape * (1 - u) / (x - min)
 *
 *  with u = 1/(1 + v), v = (x - min)/scale.
 *
 *  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"
#include "actuar.h"

double dpareto2(double x, double min, double shape, double scale,
		int give_log)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(min) || ISNAN(shape) || ISNAN(scale))
	return x + min + shape + scale;
#endif
    if (!R_FINITE(min) ||
        !R_FINITE(shape) ||
        shape <= 0.0 ||
        scale <= 0.0)
        return R_NaN;

    if (!R_FINITE(x) || x < min)
        return ACT_D__0;

    /* handle x == min separately */
    if (x == min) return ACT_D_val(shape / scale);

    double logv, logu, log1mu;

    logv = log(x - min) - log(scale);
    logu = - log1pexp(logv);
    log1mu = - log1pexp(-logv);

    return ACT_D_exp(log(shape) + shape * logu + log1mu - log(x - min));
}

double ppareto2(double q, double min, double shape, double scale,
		int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(q) || ISNAN(min) || ISNAN(shape) || ISNAN(scale))
	return q + min + shape + scale;
#endif
    if (!R_FINITE(min) ||
        !R_FINITE(shape) ||
        shape <= 0.0 ||
        scale <= 0.0)
        return R_NaN;

    if (q <= min)
        return ACT_DT_0;

    double u = exp(-log1pexp(log(q - min) - log(scale)));

    return ACT_DT_Cval(R_pow(u, shape));
}

double qpareto2(double p, double min, double shape, double scale,
		int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(min) || ISNAN(shape) || ISNAN(scale))
	return p + min + shape + scale;
#endif
    if (!R_FINITE(min) ||
        !R_FINITE(shape) ||
        !R_FINITE(scale) ||
        shape <= 0.0 ||
        scale <= 0.0)
        return R_NaN;

    ACT_Q_P01_boundaries(p, 0, R_PosInf);
    p = ACT_D_qIv(p);

    return min + scale * (R_pow(ACT_D_Cval(p), -1.0/shape) - 1.0);
}

double rpareto2(double min, double shape, double scale)
{
    if (!R_FINITE(min) ||
        !R_FINITE(shape) ||
        !R_FINITE(scale) ||
        shape <= 0.0 ||
        scale <= 0.0)
        return R_NaN;

    return min + scale * (R_pow(unif_rand(), -1.0/shape) - 1.0);
}

double mpareto2(double order, double min, double shape, double scale,
		int give_log)
{
#ifdef IEEE_754
    if (ISNAN(order) || ISNAN(min) || ISNAN(shape) || ISNAN(scale))
	return order + shape + scale;
#endif
    if (!R_FINITE(min) ||
        !R_FINITE(shape) ||
        !R_FINITE(scale) ||
        !R_FINITE(order) ||
        shape <= 0.0  ||
        scale <= 0.0)
        return R_NaN;

    /* The case min = 0 is a Pareto with a larger range of admissible
     * values for order: -1 < order < shape. */
    if (min == 0.0)
        return mpareto(order, shape, scale, give_log);

    /* From now on min != 0 and order must be a stricly non negative
     * integer < shape. */
    if (order < 0.0)
        return R_NaN;
    if (order >= shape)
        return R_PosInf;

    int i;
    double order0 = order;
    double sum, r = scale/min;
    double Ga = gammafn(shape);

    if (ACT_nonint(order))
    {
	order = ACT_forceint(order);
	warning(_("'order' (%.2f) must be integer, rounded to %.0f"),
		order0, order);
    }

    sum = Ga;			/* first term in the sum */
    for (i = 1; i <= order; i++)
    {
        sum += choose(order, i) * R_pow(r, i)
	    * gammafn(1.0 + i) * gammafn(shape - i);
    }

    return R_pow(min, order) * sum / Ga;
}

double levpareto2(double limit, double min, double shape, double scale, double order,
                 int give_log)
{
#ifdef IEEE_754
    if (ISNAN(limit) || ISNAN(min) || ISNAN(shape) || ISNAN(scale) || ISNAN(order))
	return limit + min + shape + scale + order;
#endif
    if (!R_FINITE(min) ||
        !R_FINITE(shape) ||
        !R_FINITE(scale) ||
        !R_FINITE(order) ||
        shape <= 0.0 ||
        scale <= 0.0)
        return R_NaN;

    if (limit <= min)
        return 0.0;

    /* The case min = 0 is a Pareto with a larger range of admissible
     * values for order: order > -1. */
    if (min == 0.0)
        return levpareto(limit, shape, scale, order, give_log);

    /* From now on min != 0 and order must be a stricly non negative
     * integer. */
    if (order < 0.0)
        return R_NaN;

    int i;
    double order0 = order;
    double logv, u, u1m;
    double sum, r = scale / min;

    logv = log(limit - min) - log(scale);
    u = exp(-log1pexp(logv));
    u1m = exp(-log1pexp(-logv));

    if (ACT_nonint(order))
    {
	order = ACT_forceint(order);
	warning(_("'order' (%.2f) must be integer, rounded to %.0f"),
		order0, order);
    }

    sum = betaint_raw(u1m, 1.0, shape, u); /* first term in the sum */
    for (i = 1; i <= order; i++)
    {
        sum += choose(order, i) * R_pow(r, i)
            * betaint_raw(u1m, 1.0 + i, shape - i, u);
    }

    return R_pow(min, order) * sum / gammafn(shape)
        + ACT_DLIM__0(limit, order) * R_pow(u, shape);
}
