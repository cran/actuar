/*  actuar: Actuarial Functions and Heavy Tailed Distributions
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the Pareto (type) IV distribution. See ../R/Pareto4.R for
 *  details.
 *
 *  We work with the density expressed as
 *
 *    shape1 * shape2 * u^shape1 * (1 - u) / (x - min)
 *
 *  with u = 1/(1 + v), v = ((x - min)/scale)^shape2.
 *
 *  AUTHORS: Nicholas Langevin and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"
#include "actuar.h"

double dpareto4(double x, double min, double shape1, double shape2,
                double scale, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(x)     || ISNAN(min) || ISNAN(shape1) || ISNAN(shape2) ||
	ISNAN(scale))
	return x + min + shape1 + shape2 + scale;
#endif
    if (!R_FINITE(min) ||
        !R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        scale <= 0.0)
        return R_NaN;

    if (!R_FINITE(x) || x < min)
        return ACT_D__0;

    /* handle (x - min) == 0 separately */
    if (x == min)
    {
	if (shape2 < 1) return R_PosInf;
	if (shape2 > 1) return ACT_D__0;
	/* else */
	return ACT_D_val(shape1 / scale);
    }

    double logv, logu, log1mu;

    logv = shape2 * (log(x - min) - log(scale));
    logu = - log1pexp(logv);
    log1mu = - log1pexp(-logv);

    return ACT_D_exp(log(shape1) + log(shape2) + shape1 * logu
		     + log1mu - log(x - min));
}

double ppareto4(double q, double min, double shape1, double shape2,
                double scale, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(q)     || ISNAN(min) || ISNAN(shape1) || ISNAN(shape2) ||
	ISNAN(scale))
	return q + min + shape1 + shape2 + scale;
#endif
    if (!R_FINITE(min) ||
        !R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        scale  <= 0.0)
        return R_NaN;

    if (q <= min)
        return ACT_DT_0;

    double u = exp(-log1pexp(shape2 * (log(q - min) - log(scale))));

    return ACT_DT_Cval(R_pow(u, shape1));
}

double qpareto4(double p, double min, double shape1, double shape2,
                double scale, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(p)     || ISNAN(min) ||  ISNAN(shape1) || ISNAN(shape2) ||
	ISNAN(scale))
	return p + min + shape1 + shape2 + scale;
#endif
    if (!R_FINITE(min) ||
        !R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(scale)  ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        scale  <= 0.0)
        return R_NaN;

    ACT_Q_P01_boundaries(p, min, R_PosInf);
    p = ACT_D_qIv(p);

    return min + scale * R_pow(R_pow(ACT_D_Cval(p), -1.0/shape1) - 1.0, 1.0/shape2);
}

double rpareto4(double min, double shape1, double shape2, double scale)
{
    if (!R_FINITE(min) ||
        !R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(scale) ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        scale <= 0.0)
        return R_NaN;

    return min + scale * R_pow(R_pow(unif_rand(), -1.0/shape1) - 1.0, 1.0/shape2);
}

double mpareto4(double order, double min,  double shape1, double shape2,
                double scale, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(order) || ISNAN(min) ||  ISNAN(shape1) || ISNAN(shape2) ||
	ISNAN(scale))
	return order + min + shape1 + shape2 + scale;
#endif
    if (!R_FINITE(min) ||
        !R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(scale)  ||
        !R_FINITE(order)  ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        scale  <= 0.0)
	return R_NaN;

    /* The case min = 0 is a Burr with a larger range of admissible
     * values for order: - shape2 < order < shape1 * shape2. */
    if (min == 0.0)
        return mburr(order, shape1, shape2, scale, give_log);

    /* From now on min != 0 and order must be a stricly non negative
     * integer < shape1 * shape2. */
    if (order < 0.0)
        return R_NaN;
    if (order >= shape1 * shape2)
        return R_PosInf;

    int i;
    double order0 = order;
    double tmp, sum, r = scale/min;
    double Ga = gammafn(shape1);

    if (ACT_nonint(order))
    {
	order = ACT_forceint(order);
	warning(_("'order' (%.2f) must be integer, rounded to %.0f"),
		order0, order);
    }

    sum = Ga;			/* first term in the sum */
    for (i = 1; i <= order; i++)
    {
        tmp = i/shape2;
        sum += choose(order, i) * R_pow(r, i)
	    * gammafn(1.0 + tmp) * gammafn(shape1 - tmp);
    }

    /* The first term of the sum is always min^order. */
    return R_pow(min, order) * sum / Ga;
}

double levpareto4(double limit, double min, double shape1, double shape2,
                  double scale, double order, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(limit) || ISNAN(min) || ISNAN(shape1) || ISNAN(shape2) ||
	ISNAN(scale) || ISNAN(order))
	return limit + min + shape1 + shape2 + scale + order;
#endif
    if (!R_FINITE(min) ||
        !R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(scale) ||
        !R_FINITE(order) ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        scale  <= 0.0)
        return R_NaN;

    if (limit <= min)
        return 0.0;

    /* The case min = 0 is a Burr with a larger range of admissible
     * values for order: order > - shape2. */
    if (min == 0.0)
        return levburr(limit, shape1, shape2, scale, order, give_log);

    /* From now on min != 0 and order must be a stricly non negative
     * integer. */
    if (order < 0.0)
        return R_NaN;

    int i;
    double order0 = order;
    double logv, u, u1m;
    double tmp, sum, r = scale / min;

    logv = shape2 * (log(limit - min) - log(scale));
    u = exp(-log1pexp(logv));
    u1m = exp(-log1pexp(-logv));

    if (ACT_nonint(order))
    {
	order = ACT_forceint(order);
	warning(_("'order' (%.2f) must be integer, rounded to %.0f"),
		order0, order);
    }

    sum = betaint_raw(u1m, 1.0, shape1, u); /* first term in the sum */
    for (i = 1; i <= order; i++)
    {
        tmp = i / shape2;
        sum += choose(order, i) * R_pow(r, i)
            * betaint_raw(u1m, 1.0 + tmp, shape1 - tmp, u);
    }

    return R_pow(min, order) * sum / gammafn(shape1)
        + ACT_DLIM__0(limit, order) * R_pow(u, shape1);
}
