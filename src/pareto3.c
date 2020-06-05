/*  actuar: Actuarial Functions and Heavy Tailed Distributions
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the Pareto (type) III distribution. See ../R/Pareto3.R for
 *  details.
 *
 *  We work with the density expressed as
 *
 *    shape * u * (1 - u) / (x - min)
 *
 *  with u = v/(1 + v), v = ((x - min)/scale)^shape.
 *
 *  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"
#include "actuar.h"

double dpareto3(double x, double min, double shape, double scale,
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
    if (x == min)
    {
	if (shape < 1) return R_PosInf;
	if (shape > 1) return ACT_D__0;
	/* else */
	return ACT_D_val(1.0/scale);
    }

    double logv, logu, log1mu;

    logv = shape * (log(x - min) - log(scale));
    logu = - log1pexp(-logv);
    log1mu = - log1pexp(logv);

    return ACT_D_exp(log(shape) + logu + log1mu - log(x - min));
}

double ppareto3(double q, double min, double shape, double scale,
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

    double u = exp(-log1pexp(shape * (log(scale) - log(q - min))));

    return ACT_DT_val(u);
}

double qpareto3(double p, double min, double shape, double scale,
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

    return min + scale * R_pow(1.0/ACT_D_Cval(p) - 1.0, 1.0/shape);
}

double rpareto3(double min, double shape, double scale)
{
    if (!R_FINITE(min) ||
        !R_FINITE(shape) ||
        !R_FINITE(scale) ||
        shape <= 0.0 ||
        scale <= 0.0)
        return R_NaN;

    return min + scale * R_pow(1.0/unif_rand() - 1.0, 1.0/shape);
}

double mpareto3(double order, double min, double shape, double scale,
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

    /* The case min = 0 is a loglogistic with a larger range of
     * admissible values for order: -shape < order < shape. */
    if (min == 0.0)
        return mllogis(order, shape, scale, give_log);

    /* From now on min != 0 and order must be a stricly non negative
     * integer < shape. */
    if (order < 0.0)
        return R_NaN;
    if (order >= shape)
        return R_PosInf;

    int i;
    double order0 = order;
    double tmp, sum, r = scale/min;

    if (ACT_nonint(order))
    {
	order = ACT_forceint(order);
	warning(_("'order' (%.2f) must be integer, rounded to %.0f"),
		order0, order);
    }

    sum = 1.0;			/* first term in the sum */
    for (i = 1; i <= order; i++)
    {
	tmp = i / shape;
        sum += choose(order, i) * R_pow(r, i)
	    * gammafn(1.0 + tmp) * gammafn(1.0 - tmp);
    }

    return R_pow(min, order) * sum;
}

double levpareto3(double limit, double min, double shape, double scale, double order,
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

    /* The case min = 0 is a loglogistic with a larger range of
     * admissible values for order: order > -shape. */
    if (min == 0.0)
        return levllogis(limit, shape, scale, order, give_log);

    /* From now on min != 0 and order must be a stricly non negative
     * integer. */
    if (order < 0.0)
        return R_NaN;

    int i;
    double order0 = order;
    double logv, u, u1m;
    double tmp, sum, r = scale / min;

    logv = shape * (log(limit - min) - log(scale));
    u = exp(-log1pexp(-logv));
    u1m = exp(-log1pexp(logv));

    if (ACT_nonint(order))
    {
	order = ACT_forceint(order);
	warning(_("'order' (%.2f) must be integer, rounded to %.0f"),
		order0, order);
    }

    sum = betaint_raw(u, 1.0, 1.0, u1m); /* first term in the sum */
    for (i = 1; i <= order; i++)
    {
	tmp = i / shape;
        sum += choose(order, i) * R_pow(r, i)
            * betaint_raw(u, 1.0 + tmp, 1.0 - tmp, u1m);
    }

    return R_pow(min, order) * sum
        + ACT_DLIM__0(limit, order) * (0.5 - u + 0.5);
}
