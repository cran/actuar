/*  actuar: Actuarial Functions and Heavy Tailed Distributions
 *
 *  Functions to compute probability function, cumulative distribution
 *  and quantile functions, and to simulate random variates for the
 *  zero-modified geometric distribution. See
 *  ../R/ZeroModifiedGeometric.R for details.
 *
 *  Zero-modified distributions are discrete mixtures between a
 *  degenerate distribution at zero and the corresponding,
 *  non-modified, distribution. As a mixture, they have density
 *
 *      Pr[Z = x] = [1 - (1 - p0m)/(1 - p0)] 1(x)
 *                  + [(1 - p0m)/(1 - p0)] Pr[X = 0],
 *
 *  where p0 = Pr[X = 0]. The density can also be expressed as
 *  Pr[Z = 0] = p0m and
 *
 *      Pr[Z = x] = (1 - p0m) * Pr[X = x]/(1 - Pr[X = 0]),
 *
 *  for x = 1, 2, ... The distribution function is, for all x,
 *
 *      Pr[Z <= x] = 1 - (1 - p0m) * (1 - Pr[X <= x])/(1 - p0).
 *
 *  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"
#include "actuar.h"

/* The geometric distribution has p0 = prob.
 *
 * Limiting case: prob == 1 is mass (1-p0m) at x = 1.
 */

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

    /* limiting case as prob approaches one is point mass (1-p0m) at one */
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

    /* limiting case as prob approaches one is mass (1-p0m) at one */
    if (prob == 1) return ACT_DT_1;

    /* working in log scale improves accuracy */
    return ACT_DT_CEval(log1p(-p0m)
			+ pgeom(x - 1, prob, /*l._t.*/0, /*log_p*/1));
}

double qzmgeom(double x, double prob, double p0m, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(prob) || ISNAN(p0m))
	return x + prob + p0m;
#endif
    if (prob <= 0 || prob > 1 || p0m < 0 || p0m > 1) return R_NaN;

    /* limiting case as prob approaches one is mass (1-p0m) at one */
    if (prob == 1)
    {
	/* simplified ACT_Q_P01_boundaries macro */
	if (log_p)
	{
	    if (x > 0)
		return R_NaN;
	    return (x <= log(p0m)) ? 0.0 : 1.0;
	}
	else /* !log_p */
	{
	    if (x < 0 || x > 1)
		return R_NaN;
	    return (x <= p0m) ? 0.0 : 1.0;
	}
    }

    ACT_Q_P01_boundaries(x, 1, R_PosInf);
    x = ACT_DT_qIv(x);

    /* working in log scale improves accuracy */
    return qgeom(-expm1(log1p(-prob) - log1p(-p0m) + log1p(-x)),
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

    /* limiting case as p approaches one is mass (1-p0m) at one */
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
