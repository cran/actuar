/*  actuar: Actuarial Functions and Heavy Tailed Distributions
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the inverse Burr distribution. See ../R/InverseBurr.R for details.
 *
 *  We work with the density expressed as
 *
 *    shape1 * shape2 * u^shape1 * (1 - u) / x
 *
 *  with u = v/(1 + v) = 1/(1 + 1/v), v = (x/scale)^shape2.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"
#include "actuar.h"

double dinvburr(double x, double shape1, double shape2, double scale,
                int give_log)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(shape1) || ISNAN(shape2) || ISNAN(scale))
	return x + shape1 + shape2 + scale;
#endif
    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        scale  <= 0.0)
        return R_NaN;

    if (!R_FINITE(x) || x < 0.0)
        return ACT_D__0;

    /* handle x == 0 separately */
    if (x == 0.0)
    {
	if (shape1 * shape2 < 1) return R_PosInf;
	if (shape1 * shape2 > 1) return ACT_D__0;
	/* else */
	return ACT_D_val(1.0/scale);
    }

    double logv, logu, log1mu;

    logv = shape2 * (log(x) - log(scale));
    logu = - log1pexp(-logv);
    log1mu = - log1pexp(logv);

    return ACT_D_exp(log(shape1) + log(shape2) + shape1 * logu
                   + log1mu - log(x));
}

double pinvburr(double q, double shape1, double shape2, double scale,
                int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(q) || ISNAN(shape1) || ISNAN(shape2) || ISNAN(scale))
	return q + shape1 + shape2 + scale;
#endif
    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        scale  <= 0.0)
        return R_NaN;

    if (q <= 0)
        return ACT_DT_0;

    double u = exp(-log1pexp(shape2 * (log(scale) - log(q))));

    return ACT_DT_val(R_pow(u, shape1));
}

double qinvburr(double p, double shape1, double shape2, double scale,
                int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(shape1) || ISNAN(shape2) || ISNAN(scale))
	return p + shape1 + shape2 + scale;
#endif
    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(scale)  ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        scale  <= 0.0)
        return R_NaN;

    ACT_Q_P01_boundaries(p, 0, R_PosInf);
    p = ACT_D_qIv(p);

    return scale * R_pow(R_pow(ACT_D_Lval(p), -1.0/shape1) - 1.0, -1.0/shape2);
}

double rinvburr(double shape1, double shape2, double scale)
{
    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(scale)  ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        scale  <= 0.0)
        return R_NaN;

    return scale * R_pow(R_pow(unif_rand(), -1.0/shape1) - 1.0, -1.0/shape2);
}

double minvburr(double order, double shape1, double shape2, double scale,
                int give_log)
{
#ifdef IEEE_754
    if (ISNAN(order) || ISNAN(shape1) || ISNAN(shape2) || ISNAN(scale))
	return order + shape1 + shape2 + scale;
#endif
    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(scale)  ||
        !R_FINITE(order)  ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        scale  <= 0.0)
        return R_NaN;

    if (order <= - shape1 * shape2 ||
        order >= shape2)
	return R_PosInf;

    double tmp = order / shape2;

    return R_pow(scale, order) * gammafn(shape1 + tmp) * gammafn(1.0 - tmp)
        / gammafn(shape1);
}

double levinvburr(double limit, double shape1, double shape2, double scale,
                  double order, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(limit) || ISNAN(shape1) || ISNAN(shape2) || ISNAN(scale) || ISNAN(order))
	return limit + shape1 + shape2 + scale + order;
#endif
    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(scale)  ||
        !R_FINITE(order) ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        scale  <= 0.0)
        return R_NaN;

    if (order  <= -shape1 * shape2)
	return R_PosInf;

    if (limit <= 0.0)
        return 0.0;

    double logv, u, u1m;
    double tmp = order / shape2;

    logv = shape2 * (log(limit) - log(scale));
    u = exp(-log1pexp(-logv));
    u1m = exp(-log1pexp(logv));

    return R_pow(scale, order)
	* betaint_raw(u, shape1 + tmp, 1.0 - tmp, u1m)
	/ gammafn(shape1)
	+ ACT_DLIM__0(limit, order) * (0.5 - R_pow(u, shape1) + 0.5);
}
