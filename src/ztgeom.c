/*  actuar: Actuarial Functions and Heavy Tailed Distributions
 *
 *  Functions to compute probability function, cumulative distribution
 *  and quantile functions, and to simulate random variates for the
 *  zero-truncated geometric distribution. See
 *  ../R/ZeroTruncatedGeometric.R for details.
 *
 *  Let X ~ Geometric(prob). The probability mass function of the
 *  zero-truncated Geometric random variable Z is
 *
 *    Pr[Z = 0] = 0
 *    Pr[Z = x] = Pr[X = x]/(1 - prob), x = 1, 2, ...
 *
 *  The distribution function is, for all x = 0, 1, 2, ...,
 *
 *    Pr[Z <= x] = (Pr[X <= x] - prob)/(1 - prob)
 *
 *  Limiting case: prob == 1 is point mass at x = 1.
 *
 *  AUTHOR: Jérémy Déraspe and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dztgeom(double x, double prob, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(prob))
	return x + prob;
#endif
    if (prob <= 0 || prob > 1) return R_NaN;

    if (x < 1 || !R_FINITE(x)) return ACT_D__0;

    /* limiting case as prob -> 1 is point mass at one */
    if (prob == 1) return (x == 1) ? ACT_D__1 : ACT_D__0;

    return ACT_D_val(dgeom(x - 1, prob, /*give_log*/0));
}

double pztgeom(double x, double prob, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(prob))
	return x + prob;
#endif
    if (prob <= 0 || prob > 1) return R_NaN;

    if (x < 1) return ACT_DT_0;
    if (!R_FINITE(x)) return ACT_DT_1;

    /* limiting case as prob -> 1 is point mass at one */
    if (prob == 1) return (x >= 1) ? ACT_DT_1 : ACT_DT_0;

    return ACT_DT_Cval(pgeom(x - 1, prob, /*l._t.*/0, /*log_p*/0));
}

double qztgeom(double p, double prob, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(prob))
	return p + prob;
#endif
    if (prob <= 0 || prob > 1) return R_NaN;
    ACT_Q_P01_check(p);
    /* limiting case as prob -> 1 is point mass at one */
    if (prob == 1) return 1.0;
    if (p == ACT_DT_0) return 1.0;
    if (p == ACT_DT_1) return R_PosInf;

    p = ACT_DT_qIv(p);

    return 1 + qgeom(p, prob, /*l._t.*/1, /*log_p*/0);
}

double rztgeom(double prob)
{
    if (!R_FINITE(prob) || prob <= 0 || prob > 1) return R_NaN;

    /* limiting case as p -> 1 is point mass at one */
    if (prob == 1) return 1.0;

    return 1 + rpois(exp_rand() * ((1 - prob) / prob));
}
