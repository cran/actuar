/*  actuar: Actuarial Functions and Heavy Tailed Distributions
 *
 *  Functions to compute probability function, cumulative distribution
 *  and quantile functions, and to simulate random variates for the
 *  zero-modified geometric distribution. See
 *  ../R/ZeroModifiedGeometric.R for details.
 *
 *  Let X ~ Geometric(prob). The probability mass function of the
 *  zero-modified Geometric random variable Z is
 *
 *    Pr[Z = 0] = p0m
 *    Pr[Z = x] = (1 - p0m) * Pr[X = x]/(1 - prob), x = 1, 2, ...
 *
 *  The distribution function is, for all x = 0, 1, 2, ...,
 *
 *    Pr[Z <= x] = 1 - (1 - p0m) * (1 - Pr[X <= x])/(1 - prob).
 *
 *  Limiting case: prob == 1 has mass (1 - p0m) at x = 1.
 *
 *  AUTHOR: Jérémy Déraspe and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"
#include "actuar.h"

double dzmgeom(double x, double prob, double p0m, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(prob) || ISNAN(p0m))
	return x + prob + p0m;
#endif
    if (prob <= 0 || prob > 1 || p0m < 0 || p0m > 1) return R_NaN;

    if (x < 0 || !R_FINITE(x)) return ACT_D__0;
    if (x == 0) return ACT_D_val(p0m);
    /* NOTE: from now on x > 0 */

    /* limiting case as prob -> 1 is point mass (1 - p0m) at one */
    if (prob == 1) return (x == 1) ? ACT_D_Clog(p0m) : ACT_D__0;

    return ACT_D_val((1 - p0m) * dgeom(x - 1, prob, /*give_log*/0));
}

double pzmgeom(double x, double prob, double p0m, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(prob) || ISNAN(p0m))
	return x + prob + p0m;
#endif
    if (prob <= 0 || prob > 1 || p0m < 0 || p0m > 1) return R_NaN;

    if (x < 0) return ACT_DT_0;
    if (!R_FINITE(x)) return ACT_DT_1;
    if (x < 1) return ACT_DT_val(p0m);
    /* NOTE: from now on x >= 1 */

    /* limiting case as prob -> 1 is mass (1 - p0m) at one */
    if (prob == 1) return ACT_DT_1;

    /* working in log scale improves accuracy */
    return ACT_DT_CEval(log1p(-p0m)
			+ pgeom(x - 1, prob, /*l._t.*/0, /*log_p*/1));
}

double qzmgeom(double p, double prob, double p0m, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(prob) || ISNAN(p0m))
	return p + prob + p0m;
#endif
    if (prob <= 0 || prob > 1 || p0m < 0 || p0m > 1) return R_NaN;
    ACT_Q_P01_check(p);
    if (p0m == 1) return 0.0;
    /* limiting case as prob -> 1 is mass (1 - p0m) at one */
    if (prob == 1)
        return (ACT_DT_qIv(p) <= p0m) ? ACT_Q_p0lim(p0m) : 1.0;
    if (p == ACT_DT_0) return ACT_Q_p0lim(p0m);
    if (p == ACT_DT_1) return R_PosInf;

    p = ACT_DT_qIv(p);

    /* at this point 0 < p < 1, so p0m = 0 is not an issue */
    /* working in log scale improves accuracy */
    return qgeom(-expm1(log1p(-prob) - log1p(-p0m) + log1p(-p)),
		 prob, /*l._t.*/1, /*log_p*/0);

}

/* ALGORITHM FOR GENERATION OF RANDOM VARIATES
 *
 * 1. p0m >= p0: just simulate variates from the discrete mixture.
 *
 * 2. p0m < p0: fastest method depends p0m.
 *
 *    2.1 p0m < ACT_INVERSION: inversion method on a restricted range.
 *
 *    2.2 p0m >= ACT_INVERSION: simulate variates from discrete mixture
 *        with the corresponding zero truncated distribution.
 *
 * The threshold ACT_INVERSION is distribution specific.
 */

#define ACT_INVERSION 0.4

double rzmgeom(double prob, double p0m)
{
    if (!R_FINITE(prob) || prob <= 0 || prob > 1 || p0m < 0 || p0m > 1) return R_NaN;

    /* limiting case as p -> 1 is mass (1 - p0m) at one */
    if (prob == 1) return (unif_rand() <= p0m) ? 0.0 : 1.0;

    /* p0m >= prob: generate from mixture */
    if (p0m >= prob)
	return (unif_rand() * (1 - prob) < (1 - p0m)) ? rgeom(prob) : 0.0;

    /* inversion method */
    if (p0m < ACT_INVERSION)
    return qgeom(runif((prob - p0m)/(1 - p0m), 1), prob, 1, 0);

    /* generate from zero truncated mixture */
    return (unif_rand() <= p0m) ? 0.0 : 1 + rpois(exp_rand() * ((1 - prob) / prob));
}
