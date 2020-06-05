/*  actuar: Actuarial Functions and Heavy Tailed Distributions
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the transformed beta distribution. See ../R/TransformedBeta.R for
 *  details.
 *
 *  We work with the density expressed as
 *
 *    shape2 * u^shape3 * (1 - u)^shape1 / (x * beta(shape1, shape3))
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

double dtrbeta(double x, double shape1, double shape2, double shape3,
               double scale, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(shape1) || ISNAN(shape2) || ISNAN(shape3) || ISNAN(scale))
	return x + shape1 + shape2 + shape3 + scale;
#endif
    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(shape3) ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        shape3 <= 0.0 ||
        scale <= 0.0)
        return R_NaN;

    if (!R_FINITE(x) || x < 0.0)
        return ACT_D__0;

    /* handle x == 0 separately */
    if (x == 0.0)
    {
	if (shape2 * shape3 < 1) return R_PosInf;
	if (shape2 * shape3 > 1) return ACT_D__0;
	/* else */
	return give_log ?
	    log(shape2) - log(scale) - lbeta(shape3, shape1) :
	    shape2 / (scale * beta(shape3, shape1));
    }

    double logv, logu, log1mu;

    logv = shape2 * (log(x) - log(scale));
    logu = - log1pexp(-logv);
    log1mu = - log1pexp(logv);

    return ACT_D_exp(log(shape2) + shape3 * logu + shape1 * log1mu
		     - log(x) - lbeta(shape3, shape1));
}

double ptrbeta(double q, double shape1, double shape2, double shape3,
               double scale, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(q) || ISNAN(shape1) || ISNAN(shape2) || ISNAN(shape3) || ISNAN(scale))
	return q + shape1 + shape2 + shape3 + scale;
#endif
    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(shape3) ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        shape3 <= 0.0 ||
        scale  <= 0.0)
        return R_NaN;

    if (q <= 0)
        return ACT_DT_0;

    double logvm, u;

    logvm = shape2 * (log(scale) - log(q)); /* -log v */
    u = exp(-log1pexp(logvm));

    if (u > 0.5)
    {
        /* Compute (1 - x) accurately */
        double u1m = exp(-log1pexp(-logvm));
        return pbeta(u1m, shape1, shape3, 1 - lower_tail, log_p);
    }

    /* else u <= 0.5 */
    return pbeta(u, shape3, shape1, lower_tail, log_p);
}

double qtrbeta(double p, double shape1, double shape2, double shape3,
               double scale, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(shape1) || ISNAN(shape2) || ISNAN(shape3) || ISNAN(scale))
	return p + shape1 + shape2 + shape3 + scale;
#endif
    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(shape3) ||
        !R_FINITE(scale)  ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        shape3 <= 0.0 ||
        scale  <= 0.0)
        return R_NaN;

    ACT_Q_P01_boundaries(p, 0, R_PosInf);
    p = ACT_D_qIv(p);

    return scale * R_pow(1.0/qbeta(p, shape3, shape1, lower_tail, 0) - 1.0,
                         -1.0/shape2);
}

double rtrbeta(double shape1, double shape2, double shape3, double scale)
{
    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(shape3) ||
        !R_FINITE(scale) ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        shape3 <= 0.0 ||
        scale <= 0.0)
        return R_NaN;

    return scale * R_pow(1.0/rbeta(shape3, shape1) - 1.0, -1.0/shape2);
}

double mtrbeta(double order, double shape1, double shape2, double shape3,
               double scale, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(order) || ISNAN(shape1) || ISNAN(shape2) || ISNAN(shape3) || ISNAN(scale))
	return order + shape1 + shape2 + shape3 + scale;
#endif
    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(shape3) ||
        !R_FINITE(scale)  ||
        !R_FINITE(order)  ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        shape3 <= 0.0 ||
        scale  <= 0.0)
	return R_NaN;

    if (order <= - shape3 * shape2 ||
        order >= shape1 * shape2)
        return R_PosInf;

    double tmp = order / shape2;

    return R_pow(scale, order) * beta(shape3 + tmp, shape1 - tmp)
        / beta(shape1, shape3);
}

double levtrbeta(double limit, double shape1, double shape2, double shape3,
                 double scale, double order, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(limit) || ISNAN(shape1) || ISNAN(shape2) || ISNAN(shape3) ||
        ISNAN(scale) || ISNAN(order))
	return limit + shape1 + shape2 + shape3 + scale + order;
#endif
    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(shape3) ||
        !R_FINITE(scale) ||
        !R_FINITE(order) ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        shape3 <= 0.0 ||
        scale  <= 0.0)
        return R_NaN;

    if (order <= - shape3 * shape2)
        return R_PosInf;

    if (limit <= 0.0)
        return 0.0;

    double logv, u, u1m, Ix;
    double tmp = order / shape2;

    logv = shape2 * (log(limit) - log(scale));
    u = exp(-log1pexp(-logv));
    u1m = exp(-log1pexp(logv));

    Ix = (u > 0.5) ?
	pbeta(u1m, shape1, shape3, /*l._t.*/1, /*give_log*/0) :
        pbeta(u,   shape3, shape1, /*l._t.*/0, /*give_log*/0);

    return R_pow(scale, order)
	* betaint_raw(u, shape3 + tmp, shape1 - tmp, u1m)
	/ (gammafn(shape1) * gammafn(shape3))
	+ ACT_DLIM__0(limit, order) * Ix;
}
