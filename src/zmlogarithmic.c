/*  actuar: Actuarial Functions and Heavy Tailed Distributions
 *
 *  Functions to compute probability function, cumulative distribution
 *  and quantile functions, and to simulate random variates for the
 *  zero-modified logarithmic distribution. See ../R/ZeroModifiedLogarithmic.R
 *  for details.
 *
 *  Let X ~ Logarithmic(prob). The probability mass function of the
 *  zero-modified Logarithmic random variable Z is
 *
 *    Pr[Z = 0] = p0m
 *    Pr[Z = x] = (1 - p0m) * Pr[X = x], x = 1, 2, ...
 *
 *  The distribution function is, for all x = 0, 1, 2, ...,
 *
 *    Pr[Z <= x] = 1 - (1 - p0m) * (1 - Pr[X <= x]).
 *
 *  Limiting case: prob == 0 has mass (1 - p0m) at x = 1.
 *
 *  AUTHOR: Jérémy Déraspe and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"
#include "actuar.h"

double dzmlogarithmic(double x, double prob, double p0m, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(prob) || ISNAN(p0m))
	return x + prob + p0m;
#endif
    if (prob < 0 || prob >= 1 || p0m < 0 || p0m > 1) return R_NaN;
    ACT_D_nonint_check(x);

    if (!R_FINITE(x) || x < 0) return ACT_D__0;
    if (x == 0) return ACT_D_val(p0m);
    /* NOTE: from now on x > 0 */

    /* limiting case as prob -> 0 is mass (1 - p0m) at one */
    if (prob == 0) return (x == 1) ? ACT_D_Clog(p0m) : ACT_D__0;

    x = ACT_forceint(x);

    double a = -1.0/log1p(-prob);

    return ACT_D_exp(log(a) + x * log(prob) + log1p(-p0m) - log(x));
}

double pzmlogarithmic(double x, double prob, double p0m, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(prob) || ISNAN(p0m))
	return x + prob + p0m;
#endif
    if (prob < 0 || prob >= 1 || p0m < 0 || p0m > 1) return R_NaN;

    if (x < 0) return ACT_DT_0;
    if (!R_FINITE(x)) return ACT_DT_1;
    if (x < 1) return ACT_DT_val(p0m);
    /* NOTE: from now on x >= 1 */

    /* simple case for all x >= 1 */
    if (p0m == 1) return ACT_DT_1;

    /* limiting case as prob -> 0 is mass (1 - p0m) at one. */
    if (prob == 0) return ACT_DT_1;

    return ACT_DT_Cval((1 - p0m) * plogarithmic(x, prob, /*l._t.*/0, /*log_p*/0));
}

double qzmlogarithmic(double p, double prob, double p0m, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(prob) || ISNAN(p0m))
	return p + prob + p0m;
#endif
    if (prob < 0 || prob >= 1 || p0m < 0 || p0m > 1) return R_NaN;
    ACT_Q_P01_check(p);
    if (p0m == 1) return 0.0;
    /* limiting case as prob -> 0 is mass (1 - p0m) at one. */
    if (prob == 0)
        return (ACT_DT_qIv(p) <= p0m) ? ACT_Q_p0lim(p0m) : 1.0;
    if (p == ACT_DT_0) return ACT_Q_p0lim(p0m);
    if (p == ACT_DT_1) return R_PosInf;

    p = ACT_DT_qIv(p);

    /* avoid rounding errors if p was given in log form */
    if (log_p)
	p0m = exp(log(p0m));

    /* avoid rounding errors if p was given as upper tail */
    if (!lower_tail)
	p0m = 0.5 - (0.5 - p0m + 0.5) + 0.5;

    /* at this point 0 < p < 1, so p0m = 0 is not an issue */
    return (p <= p0m) ? 0.0 : qlogarithmic((p - p0m)/(1 - p0m), prob, /*l._t.*/1, /*log_p*/0);
}

/* ALGORITHM FOR GENERATION OF RANDOM VARIATES
 *
 * Just simulate variates from the discrete mixture.
 *
 */
double rzmlogarithmic(double prob, double p0m)
{
    if (prob < 0 || prob >= 1 || p0m < 0 || p0m > 1) return R_NaN;

    return (unif_rand() < p0m) ? 0.0 : rlogarithmic(prob);
}
