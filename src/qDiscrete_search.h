/*  actuar: Actuarial Functions and Heavy Tailed Distributions
 *
 *  Find quantiles for discrete distributions using the Cornish-Fisher
 *  Expansion.
 *
 *  This file is a copy of src/nmath/qDiscrete_search.h of R sources,
 *  but with the following minor changes:
 *
 *  1. all debugging material is deleted;
 *  2. the macro function q_DISCRETE_01_CHECKS(), which does not serve
 *     any purpose without the debugging material is deleted;
 *  3. the declaration of variables in q_DISCRETE_BODY() is moved into
 *     a separate macro; see comments marked 'VG' below.
 *
 *  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 *          based on code from the R Core Team
 */

/* This is #included from ./logarithmic.c and ./poisinvgauss.c */

#define PST_0(a, b) a ## b
#define PASTE(a, b) PST_0(a, b)

#define _pDIST_  PASTE(p, _thisDIST_)
#define _qDIST_  PASTE(q, _thisDIST_)

#ifdef MATHLIB_STANDALONE
# define MAYBE_R_CheckUserInterrupt()
#else
# define MAYBE_R_CheckUserInterrupt() R_CheckUserInterrupt()
#endif

#define DO_SEARCH_FUN(...)					\
    do_search(double y, double *z, double p, __VA_ARGS__,	\
              double incr, int lower_tail, int log_p)

#define DO_SEARCH_(Y_, incr_, ...)				\
    do_search(Y_, &z, p, __VA_ARGS__, incr_, lower_tail, log_p)

#define P_DIST(Y_, ...) _pDIST_(Y_, __VA_ARGS__, lower_tail, log_p)

static double DO_SEARCH_FUN(_dist_PARS_DECL_)
{
    Rboolean left = (lower_tail ? (*z >= p) : (*z < p));

    if(left) {	// (lower_tail, *z >= p)  or  (upper tail, *z < p): search to the __left__
	for(int iter = 0; ; iter++) {
	    double newz = -1.; // -Wall
	    if(y > 0)
		newz = P_DIST(y - incr, _dist_PARS_);
	    else if(y < 0)
		y = 0;
	    // note that newz may be NaN because of remaining border line bugs in _pDIST_() {bug from pbeta()}
	    if(y == 0 || ISNAN(newz) || (lower_tail ? (newz < p) : (newz >= p))) {
		return y; // and previous *z
	    }
	    y = fmax2(0, y - incr);
	    *z = newz;
	}
    }
    else { // (lower_tail, *z < p)  or  (upper tail, *z >= p): search to the __right__
	for(int iter = 0; ; iter++) {
	    y += incr;
	    *z = P_DIST(y, _dist_PARS_);

	    if(ISNAN(*z) || (lower_tail ? (*z >= p) : (*z < p)))
	    {
		return y;
	    }
	}
    }
} // do_search()

#define q_DISCR_CHECK_BOUNDARY(Y_) if(Y_ < 0) Y_ = 0.

/* VG: the Poisson-inverse gaussian requires different declarations
 * for a limiting case. Therefore, the standard declaration is taken
 * out of the q_DISCRETE_BODY() macro found in R sources. */
#define q_DISCRETE_DECL							\
    double								\
    z = qnorm(p, 0., 1., lower_tail, log_p),				\
	y = ACT_forceint(mu + sigma * (z + gamma * (z*z - 1) / 6))

#define q_DISCRETE_BODY() do {						\
    q_DISCR_CHECK_BOUNDARY(y);						\
									\
    z = P_DIST(y, _dist_PARS_);						\
									\
    /* Algorithmic "tuning parameters", used to be hardwired; changed for speed &| precision */	\
    const double							\
	_pf_n_  = 8,     /* was hardwired to 64 */			\
	_pf_L_  = 2,     /* was hardwired to 64 */			\
	_yLarge_ = 4096, /* was hardwired to 1e5 */			\
	_incF_ = (1./64),/* was hardwired to 0.001 (= 1./1000 ) */	\
	_iShrink_ = 8,   /* was hardwired to  100 */			\
	_relTol_ = 1e-15,/* was hardwired to 1e-15 */			\
	_xf_ = 4; /* extra factor, *must* be >= 1 (new anyway) */	\
									\
    /* fuzz to ensure left continuity: do not loose too much (=> error in upper tail) */ \
    if(log_p) { /* <==> p \in [-Inf, 0]  different adjustment: "other sign" */ \
	double e = _pf_L_ * DBL_EPSILON;				\
	if(lower_tail && p > - DBL_MAX) /* prevent underflow to -Inf */	\
    	    p *= 1 + e;							\
    	else /* if(p < - DBL_MIN) // not too close to -0 */		\
    	    p *= 1 - e;							\
									\
    } else { /* not log scale */					\
	double e = _pf_n_ * DBL_EPSILON;				\
	if(lower_tail)							\
	    p *= 1 - e;							\
	else if(1 - p > _xf_*e) /* otherwise get p > 1 */		\
	    p *= 1 + e;							\
    }									\
									\
    /* If the C-F value  y is not too large a simple search is OK */	\
    if(y < _yLarge_) return DO_SEARCH_(y, 1, _dist_PARS_);		\
    /* Otherwise be a bit cleverer in the search: use larger increments, notably initially: */ \
    {  /* y >= _yLarge_ */						\
	double oldincr, incr = floor(y * _incF_);			\
	int qIt = 0;							\
									\
	do {								\
	    oldincr = incr;						\
	    y = DO_SEARCH_(y, incr, _dist_PARS_); /* also updating *z */ \
	    if(++qIt % 10000 == 0) MAYBE_R_CheckUserInterrupt();	\
	    incr = fmax2(1, floor(incr / _iShrink_));			\
	} while(oldincr > 1 && incr > y * _relTol_);			\
	return y;							\
    }									\
} while(0)
