/*  actuar: Actuarial Functions and Heavy Tailed Distributions
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, and to simulate random variates for the logarithmic
 *  discrete distribution. See ../R/Logarithmic.R for details.
 *
 *  We work with the probability mass function expressed as
 *
 *    a * p^x / x,  x = 1, 2, ...
 *
 *  with a = -1/log(1 - p).
 *
 *  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>	/* for R_CheckUserInterrupt() */
#include "locale.h"
#include "dpq.h"

double dlogarithmic(double x, double prob, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(prob))
	return x + prob;
#endif
    if (prob < 0 || prob >= 1) return R_NaN;
    ACT_D_nonint_check(x);

    if (!R_FINITE(x) || x < 1) return ACT_D__0;

    /* limiting case as prob approaches zero is point mass at one */
    if (prob == 0) return (x == 1) ? ACT_D__1 : ACT_D__0;

    x = ACT_forceint(x);

    double a = -1.0/log1p(-prob);

    return ACT_D_exp(log(a) + x * log(prob) - log(x));
}

/*  For plogarithmic(), there does not seem to be algorithms much more
 *  elaborate that successive computations of the probabilities using
 *  the recurrence relationship
 *
 *  P[X = x + 1] = p * x * Pr[X = x] / (x + 1), x = 2, 3, ...
 *
 *  with Pr[X = 1] = -p/log(1 - p). This is what is done here.
 */

double plogarithmic(double q, double prob, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(q) || ISNAN(prob))
	return q + prob;
#endif
    if (prob < 0 || prob >= 1) return R_NaN;

    if (q < 1) return ACT_DT_0;
    if (!R_FINITE(q)) return ACT_DT_1;

    /* limiting case as prob approaches zero is point mass at one. */
    if (prob == 0) return (q >= 1) ? ACT_DT_1 : ACT_DT_0;

    int k;
    double s, pk;

    pk = -prob/log1p(-prob);		      /* Pr[X = 1] */
    s = pk;

    if (q == 1) return ACT_DT_val(s); /* simple case */

    for (k = 1; k < q; k++)
    {
	pk *= prob * k/(k + 1.0);
	s += pk;
    }

    return ACT_DT_val(s);
}

/* For qlogarithmic() we mostly reuse the code from qnbinom() et al.
 * of R sources. From src/nmath/qnbinom.c:
 *
 *  METHOD
 *
 *	Uses the Cornish-Fisher Expansion to include a skewness
 *	correction to a normal approximation.  This gives an
 *	initial value which never seems to be off by more than
 *	1 or 2.	 A search is then conducted of values close to
 *	this initial start point.
 */

#define _thisDIST_ logarithmic
#define _dist_PARS_DECL_ double prob
#define _dist_PARS_      prob

#include "qDiscrete_search.h"	/* do_search() et al. */

double qlogarithmic(double p, double prob, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(prob))
	return p + prob;
#endif
    if (prob < 0 || prob >= 1) return R_NaN;

    /* limiting case as prob approaches zero is point mass at one */
    if (prob == 0)
    {
	/* simplified ACT_Q_P01_boundaries macro */
	if (log_p)
	{
	    if (p > 0)
		return R_NaN;
	    return 1.0;
	}
	else /* !log_p */
	{
	    if (p < 0 || p > 1)
		return R_NaN;
	    return 1.0;
	}
    }

    ACT_Q_P01_boundaries(p, 1.0, R_PosInf);

    double
	a = -1.0/log1p(-prob),
	P = a * prob,
	Q = 1.0/(0.5 - prob + 0.5),
	mu = P * Q,
	sigma = sqrt(mu * (Q - mu)),
	gamma = (P * (1 + prob - P*(3 + 2*P)) * R_pow_di(Q, 3))/R_pow_di(sigma, 3);

    /* q_DISCRETE_01_CHECKS(); */
    q_DISCRETE_DECL;
    q_DISCRETE_BODY();
}

/*  rlogarithmic() is an implementation with automatic selection of
 *  the LS and LK algorithms of:
 *
 *  Kemp, A. W. (1981), Efficient Generation of Logarithmically
 *  Distributed Pseudo-Random Variables, Journal of the Royal
 *  Statistical Society, Series C. Vol. 30, p. 249-253.
 *  URL http://www.jstor.org/stable/2346348
 *
 *  The algorithms are also discussed in chapter 10 of Devroye (1986).
 */

double rlogarithmic(double prob)
{
    if (prob < 0 || prob > 1) return R_NaN;

    /* limiting case as prob approaches zero is point mass at one. */
    if (prob == 0) return 1.0;

    /* Automatic selection between the LS and LK algorithms */
    if (prob < 0.95)
    {
	double s = -prob/log1p(-prob);
	double x = 1.0;
	double u = unif_rand();

	while (u > s)
	{
	    u -= s;
	    x += 1.0;
	    s *= prob * (x - 1.0)/x;
	}

	return x;
    }

    /* else (prob >= 0.95) */
    {
	double r = log1p(-prob);
	double v = unif_rand();

	if (v >= prob)       return 1.0;

	double u = unif_rand();
	double q = -expm1(r * u);

	if (v <= (q * q)) return floor(1.0 + log(v)/log(q));
	if (v <= q)       return 2.0; /* case q^2 < v <= q */
	return 1.0;		      /* case v > q */
    }
}
