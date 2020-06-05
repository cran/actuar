/*  actuar: Actuarial Functions and Heavy Tailed Distributions
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, and to simulate random variates for the Poisson-inverse
 *  gaussian distribution. See ../R/PoissonInverseGaussian.R for
 *  details.
 *
 *  We work with the density expressed as
 *
 *    p(x) = sqrt(1/phi) sqrt(1/(pi/2)) exp(1/(phi mu))/x!
 *         * [sqrt(2 phi (1 + (2 phi mu^2)^(-1)))]^(-(x - 0.5))
 *         * bessel_k(sqrt(2/phi (1 + (2 phi mu^2)^(-1))), x - 0.5)
 *
 *  or, is essence,
 *
 *    p(x) = A exp(1/(phi mu))/x! B^(-y) bessel_k(B/phi, y)
 *
 *  The limiting case mu = Inf is handled "automatically" with terms
 *  going to zero when mu is Inf. Specific code not worth it since the
 *  function should rarely be evaluated with mu = Inf in practice.
 *
 *  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"
#include "actuar.h"

double dpoisinvgauss_raw(double x, double mu, double phi, int give_log)
{
/* Here assume that x is integer, 0 < x < Inf, mu > 0, 0 < phi < Inf */

    int i;
    double p, pi1m, pi2m;
    double twophi = phi + phi;

    /* limiting case mu = Inf with simpler recursive formulas */
    if (!R_FINITE(mu))
    {
	p = -sqrt(2/phi);	      /* log p[0] */
	if (x == 0.0)
	    return ACT_D_exp(p);

	pi2m = exp(p);		      /* p[i - 2] = p[0]*/
	p = p - (M_LN2 + log(phi))/2; /* log p[1] */
	if (x == 1.0)
	    return ACT_D_exp(p);

	pi1m = exp(p);		/* p[i - 1] = p[1] */
	for (i = 2; i <= x; i++)
	{
	    p = (1 - 1.5/i) * pi1m + pi2m/twophi/(i * (i - 1));
	    pi2m = pi1m;
	    pi1m = p;
	}
	return ACT_D_val(p);
    }

    /* else: "standard" case with mu < Inf */
    double A, B;
    double mu2 = mu * mu;
    double twophimu2 = twophi * mu2;

    p = (1.0 - sqrt(1.0 + twophimu2))/phi/mu; /* log p[0] */
    if (x == 0.0)
	return ACT_D_exp(p);

    pi2m = exp(p);			      /* p[i - 2] = p[0]*/
    p = log(mu) + p - log1p(twophimu2)/2.0;   /* log p[1] */
    if (x == 1.0)
	return ACT_D_exp(p);

    pi1m = exp(p);		              /* p[i - 1] = p[1] */
    A = 1.0/(1.0 + 1.0/twophimu2);	      /* constant in first term */
    B = mu2/(1.0 + twophimu2);		      /* constant in second term */
    for (i = 2; i <= x; i++)
    {
	p = A * (1 - 1.5/i) * pi1m + (B * pi2m)/(i * (i - 1));
	pi2m = pi1m;
	pi1m = p;
    }

    return ACT_D_val(p);
}

double dpoisinvgauss(double x, double mu, double phi, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(mu) || ISNAN(phi))
	return x + mu + phi;
#endif
    if (mu <= 0.0 || phi <= 0.0)
        return R_NaN;
    ACT_D_nonint_check(x);

    if (!R_FINITE(x) || x < 0.0)
	return ACT_D__0;

    /* limiting case phi = Inf */
    if (!R_FINITE(phi))
	return (x == 0) ? ACT_D__1 : ACT_D__0;

    return dpoisinvgauss_raw(x, mu, phi, give_log);
}

/*  For ppoisinvgauss(), there does not seem to be algorithms much
 *  more elaborate than successive computations of the probabilities.
 */

double ppoisinvgauss(double q, double mu, double phi, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(q) || ISNAN(mu) || ISNAN(phi))
	return q + mu + phi;
#endif
    if (mu <= 0.0 || phi <= 0.0)
        return R_NaN;

    if (q < 0)
        return ACT_DT_0;

    /* limiting case phi = Inf */
    if (!R_FINITE(phi))
    	return ACT_DT_1;

    if (!R_FINITE(q))
	return ACT_DT_1;

    int x;
    double s = 0;

    for (x = 0; x <= q; x++)
	s += dpoisinvgauss_raw(x, mu, phi, /*give_log*/ 0);

    return ACT_DT_val(s);
}

/*  For qpoiinvgauss(), we mostly reuse the code for qnbinom() et al.
 *  in the R sources. From src/nmath/qnbinom.c:
 *
 *  METHOD
 *
 *	Uses the Cornish-Fisher Expansion to include a skewness
 *	correction to a normal approximation.  This gives an
 *	initial value which never seems to be off by more than
 *	1 or 2.	 A search is then conducted of values close to
 *	this initial start point.
 *
 *  For the limiting case mu = Inf (that has no finite moments), we
 *  use instead the quantile of an inverse chi-square distribution as
 *  starting point.
 */

static double
do_search(double y, double *z, double p, double mu, double phi, double incr)
{
    if(*z >= p) {	/* search to the left */
	for(;;) {
	    if(y == 0 ||
	       (*z = ppoisinvgauss(y - incr, mu, phi, /*l._t.*/1, /*log_p*/0)) < p)
		return y;
	    y = fmax2(0, y - incr);
	}
    }
    else {		/* search to the right */
	for(;;) {
	    y = y + incr;
	    if((*z = ppoisinvgauss(y, mu, phi, /*l._t.*/1, /*log_p*/0)) >= p)
		return y;
	}
    }
}

double qpoisinvgauss(double p, double mu, double phi, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(mu) || ISNAN(phi))
	return p + mu + phi;
#endif
    if (mu <= 0.0 || phi <= 0.0)
        return R_NaN;

    /* limiting case phi = Inf */
    if (!R_FINITE(phi))
    	return 0.0;

    ACT_Q_P01_boundaries(p, 0, R_PosInf);

    /* ## From R sources ##
     * Note : "same" code in qpois.c, qbinom.c, qnbinom.c --
     * FIXME: This is far from optimal [cancellation for p ~= 1, etc]: */
    if (!lower_tail || log_p)
    {
	p = ACT_DT_qIv(p); /* need check again (cancellation!): */
	if (p == ACT_DT_0) return 0;
	if (p == ACT_DT_1) return R_PosInf;
    }
    /* ## From R sources ##
     * temporary hack --- FIXME --- */
    if (p + 1.01 * DBL_EPSILON >= 1.0) return R_PosInf;

    double z, y;

    /* limiting case mu = Inf -> inverse chi-square as starting point*/
    if (!R_FINITE(mu))
	y = ACT_forceint(1/phi/qchisq(p, 1,  /*l._t.*/0, /*give_log*/0));
    /* other cases -> Corning-Fisher */
    else
    {
	double sigma, sigma2, gamma;
	double phim2 = phi * mu * mu;

	sigma2 = phim2 * mu + mu;
	sigma = sqrt(sigma2);
	gamma = (3 * phim2 * sigma2 + mu)/sigma2/sigma;

	z = qnorm(p, 0.0, 1.0, /*lower_tail*/1, /*log_p*/0);
	y = ACT_forceint(mu + sigma * (z + gamma * (z*z - 1)/6));
    }

    z = ppoisinvgauss(y, mu, phi, /*lower_tail*/1, /*log_p*/0);

    /* ## From R sources ##
     * fuzz to ensure left continuity: */
    p *= 1 - 64*DBL_EPSILON;

    /* ## From R sources ##
     * If the C-F value is not too large a simple search is OK */
    if (y < 1e5) return do_search(y, &z, p, mu, phi, 1);
    /* ## From R sources ##
     * Otherwise be a bit cleverer in the search */
    {
	double incr = floor(y * 0.001), oldincr;
	do {
	    oldincr = incr;
	    y = do_search(y, &z, p, mu, phi, incr);
	    incr = fmax2(1, floor(incr/100));
	} while(oldincr > 1 && incr > y*1e-15);
	return y;
    }
}

double rpoisinvgauss(double mu, double phi)
{
    if (mu <= 0.0 || phi <= 0.0)
        return R_NaN;

    return rpois(rinvgauss(mu, phi));
}
