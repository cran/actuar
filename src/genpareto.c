/*  actuar: Actuarial Functions and Heavy Tailed Distributions
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the Generalized Pareto distribution.. See ../R/GeneralizedPareto.R
 *  for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"
#include "actuar.h"

double dgenpareto(double x, double shape1, double shape2, double scale,
                  int give_log)
{
    /*  We work with the density expressed as
     *
     *  u^shape2 * (1 - u)^shape1 / (x * beta(shape1, shape2))
     *
     *  with u = v/(1 + v) = 1/(1 + 1/v), v = x/scale.
     */

#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(shape1) || ISNAN(shape2) || ISNAN(scale))
	return x + shape1 + shape2 + scale;
#endif
    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(scale)  ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        scale  <= 0.0)
        return R_NaN;

    if (!R_FINITE(x) || x < 0.0)
        return ACT_D__0;

    /* handle x == 0 separately */
    if (x == 0.0)
    {
	if (shape2 < 1) return R_PosInf;
	if (shape2 > 1) return ACT_D__0;
	/* else */
	return give_log ?
	    - log(scale) - lbeta(shape2, shape1) :
	    1 / (scale * beta(shape2, shape1));
    }

    double tmp, logu, log1mu;

    tmp = log(x) - log(scale);
    logu = - log1pexp(-tmp);
    log1mu = - log1pexp(tmp);

    return ACT_D_exp(shape2 * logu + shape1 * log1mu - log(x)
                   - lbeta(shape2, shape1));
}

double pgenpareto(double q, double shape1, double shape2, double scale,
                  int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(q) || ISNAN(shape1) || ISNAN(shape2) || ISNAN(scale))
	return q + shape1 + shape2 + scale;
#endif
    if (!R_FINITE(shape1) ||
        !R_FINITE(scale)  ||
        !R_FINITE(shape2) ||
        shape1 <= 0.0 ||
        scale <= 0.0 ||
        shape2 <= 0.0)
        return R_NaN;

    if (q <= 0)
        return ACT_DT_0;

    double u = exp(-log1pexp(log(scale) - log(q)));

    return pbeta(u, shape2, shape1, lower_tail, log_p);
}

double qgenpareto(double p, double shape1, double shape2, double scale,
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

    return scale / (1.0 / qbeta(p, shape2, shape1, lower_tail, 0) - 1.0);
}

double rgenpareto(double shape1, double shape2, double scale)
{
    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(scale)  ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        scale  <= 0.0)
        return R_NaN;

    return scale / (1.0 / rbeta(shape2, shape1) - 1.0);
}

double mgenpareto(double order, double shape1, double shape2, double scale,
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

    if (order <= -shape2 ||
        order >= shape1)
	return R_PosInf;

    return R_pow(scale, order) * beta(shape1 - order, shape2 + order)
        / beta(shape1, shape2);
}

double levgenpareto(double limit, double shape1, double shape2, double scale,
                    double order, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(limit) || ISNAN(shape1) || ISNAN(shape2) || ISNAN(scale) || ISNAN(order))
	return limit + shape1 + shape2 + scale + order;
#endif
    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(scale)  ||
        !R_FINITE(order)  ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        scale  <= 0.0)
        return R_NaN;

    if (order <= -shape2)
	return R_PosInf;

    if (limit <= 0.0)
        return 0.0;

    double u = exp(-log1pexp(log(scale) - log(limit)));

    return R_pow(scale, order)
	* betaint_raw(u, shape2 + order, shape1 - order)
	/ (gammafn(shape1) * gammafn(shape2))
        + ACT_DLIM__0(limit, order) * pbeta(u, shape2, shape1, 0, 0);
}
