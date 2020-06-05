/*  actuar: Actuarial Functions and Heavy Tailed Distributions
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the Feller-Pareto distribution. See ../R/FellerPareto.R for
 *  details.
 *
 *  We work with the density expressed as
 *
 *    shape2 * u^shape3 * (1 - u)^shape1 / ((x - min) * beta(shape1, shape3))
 *
 *  with u = v/(1 + v) = 1/(1 + 1/v), v = ((x - min)/scale)^shape2.
 *
 *  AUTHORS: Nicholas Langevin and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"
#include "actuar.h"

double dfpareto(double x, double min, double shape1, double shape2,
                double shape3, double scale, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(x)      || ISNAN(min) || ISNAN(shape1) || ISNAN(shape2) ||
        ISNAN(shape3) || ISNAN(scale))
	return x + min + shape1 + shape2 + shape3 + scale;
#endif
    if (!R_FINITE(min) ||
        !R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(shape3) ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        shape3 <= 0.0 ||
        scale <= 0.0)
        return R_NaN;

    if (!R_FINITE(x) || x < min)
        return ACT_D__0;

    /* handle (x - min) == 0 separately */
    if (x == min)
    {
	if (shape2 * shape3 < 1) return R_PosInf;
	if (shape2 * shape3 > 1) return ACT_D__0;
	/* else */
	return give_log ?
	    log(shape2) - log(scale) - lbeta(shape3, shape1) :
	    shape2 / (scale * beta(shape3, shape1));
    }

    double logv, logu, log1mu;

    logv = shape2 * (log(x - min) - log(scale));
    logu = - log1pexp(-logv);
    log1mu = - log1pexp(logv);

    return ACT_D_exp(log(shape2) + shape3 * logu + shape1 * log1mu
                   - log(x - min) - lbeta(shape3, shape1));
}

double pfpareto(double q, double min, double shape1, double shape2,
                double shape3, double scale, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(q)      || ISNAN(min)    || ISNAN(shape1) ||
        ISNAN(shape2) || ISNAN(shape3) || ISNAN(scale))
	return q + min + shape1 + shape2 + shape3 + scale;
#endif
    if (!R_FINITE(min) ||
        !R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(shape3) ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        shape3 <= 0.0 ||
        scale  <= 0.0)
        return R_NaN;

    if (q <= min)
        return ACT_DT_0;

    double logvm, u;

    logvm = shape2 * (log(scale) - log(q - min)); /* -log v */
    u = exp(-log1pexp(logvm));

    if (u > 0.5)
    {
        /* Compute (1 - u) accurately */
        double u1m = exp(-log1pexp(-logvm));
        return pbeta(u1m, shape1, shape3, 1 - lower_tail, log_p);
    }

    /* else u <= 0.5 */
    return pbeta(u, shape3, shape1, lower_tail, log_p);
}

double qfpareto(double p, double min, double shape1, double shape2,
                double shape3, double scale, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(p)      || ISNAN(min) ||  ISNAN(shape1) || ISNAN(shape2) ||
        ISNAN(shape3) || ISNAN(scale))
	return p + min + shape1 + shape2 + shape3 + scale;
#endif
    if (!R_FINITE(min) ||
        !R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(shape3) ||
        !R_FINITE(scale)  ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        shape3 <= 0.0 ||
        scale  <= 0.0)
        return R_NaN;

    ACT_Q_P01_boundaries(p, min, R_PosInf);
    p = ACT_D_qIv(p);

    return min + scale * R_pow(1.0
            / qbeta(p, shape3, shape1, lower_tail, 0) - 1.0, -1.0/shape2);
}

double rfpareto(double min, double shape1, double shape2, double shape3,
                double scale)
{
    if (!R_FINITE(min) ||
        !R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(shape3) ||
        !R_FINITE(scale) ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        shape3 <= 0.0 ||
        scale <= 0.0)
        return R_NaN;

    return min + scale * R_pow(1.0/rbeta(shape1, shape3) - 1.0, 1.0/shape2);
}

double mfpareto(double order, double min,  double shape1, double shape2,
                double shape3, double scale, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(order)  || ISNAN(min) ||  ISNAN(shape1) || ISNAN(shape2) ||
        ISNAN(shape3) || ISNAN(scale))
	return order + min + shape1 + shape2 + shape3 + scale;
#endif
    if (!R_FINITE(min) ||
        !R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(shape3) ||
        !R_FINITE(scale)  ||
        !R_FINITE(order)  ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        shape3 <= 0.0 ||
        scale  <= 0.0)
	return R_NaN;

    /* The case min = 0 is a Transformed Beta with a larger range of
     * admissible values for order: - shape3 * shape2 < order <
     * shape1 * shape2. */
    if (min == 0.0)
        return mtrbeta(order, shape1, shape2, shape3, scale, give_log);

    /* From now on min != 0 and order must be a stricly non negative
     * integer < shape1 * shape2. */
    if (order < 0.0)
        return R_NaN;
    if (order >= shape1 * shape2)
        return R_PosInf;

    int i;
    double order0 = order;
    double tmp, sum, r = scale/min;
    double Be = beta(shape1, shape3);

    if (ACT_nonint(order))
    {
	order = ACT_forceint(order);
	warning(_("'order' (%.2f) must be integer, rounded to %.0f"),
		order0, order);
    }

    sum = Be;			/* first term in the sum */
    for (i = 1; i <= order; i++)
    {
        tmp = i/shape2;
        sum += choose(order, i) * R_pow(r, i)
	    * beta(shape3 + tmp, shape1 - tmp);
    }

    return R_pow(min, order) * sum / Be;
}

double levfpareto(double limit, double min, double shape1, double shape2,
                  double shape3, double scale, double order, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(limit) || ISNAN(min) || ISNAN(shape1) || ISNAN(shape2) ||
        ISNAN(shape3) || ISNAN(scale) || ISNAN(order))
	return limit + min + shape1 + shape2 + shape3 + scale + order;
#endif
    if (!R_FINITE(min) ||
        !R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(shape3) ||
        !R_FINITE(scale) ||
        !R_FINITE(order) ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        shape3 <= 0.0 ||
        scale  <= 0.0)
        return R_NaN;

    if (limit <= min)
        return 0.0;

    /* The case min = 0 is a Transformed Beta with a larger range of
     * admissible values for order: order > - shape3 * shape2. */
    if (min == 0.0)
        return levtrbeta(limit, shape1, shape2, shape3, scale, order, give_log);

    /* From now on min != 0 and order must be a stricly non negative
     * integer. */
    if (order < 0.0)
        return R_NaN;

    int i;
    double order0 = order;
    double logv, u, u1m, Ix;
    double tmp, sum, r = scale / min;

    logv = shape2 * (log(limit - min) - log(scale));
    u = exp(-log1pexp(-logv));
    u1m = exp(-log1pexp(logv));

    if (ACT_nonint(order))
    {
	order = ACT_forceint(order);
	warning(_("'order' (%.2f) must be integer, rounded to %.0f"),
		order0, order);
    }

    sum = betaint_raw(u, shape3, shape1, u1m); /* first term in the sum */
    for (i = 1; i <= order; i++)
    {
        tmp = i / shape2;
        sum += choose(order, i) * R_pow(r, i)
            * betaint_raw(u, shape3 + tmp, shape1 - tmp, u1m);
    }

    Ix = (u > 0.5) ?
	pbeta(u1m, shape1, shape3, /*l._t.*/1, /*give_log*/0) :
        pbeta(u,   shape3, shape1, /*l._t.*/0, /*give_log*/0);

    return R_pow(min, order) * sum / (gammafn(shape1) * gammafn(shape3))
        + ACT_DLIM__0(limit, order) * Ix;
}
