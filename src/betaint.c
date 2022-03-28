/*  actuar: Actuarial Functions and Heavy Tailed Distributions
 *
 *  Function to compute the integral
 *
 *    B(a, b; x) = gammafn(a + b) int_0^x t^(a-1) (1-t)^(b-1) dt
 *
 *  for a > 0, b != -1, -2, ... and 0 < x < 1. When b > 0,
 *
 *    B(a, b; x) = gammafn(a) gammafn(b) pbeta(x, a, b).
 *
 *  When b < 0 and b != -1, -2, ... and a > 1 + floor(-b),
 *
 *    B(a, b; x)
 *    = -gammafn(a + b) {(x^(a-1) (1-x)^b)/b
 *        + [(a-1) x^(a-2) (1-x)^(b+1)]/[b(b+1)]
 *        + ...
 *        + [(a-1)...(a-r) x^(a-r-1) (1-x)^(b+r)]/[b(b+1)...(b+r)]}
 *      + [(a-1)...(a-r-1)]/[b(b+1)...(b+r)] gammafn(a-r-1)
 *        * gammafn(b+r+1) pbeta(x, a-r-1, b+r+1)
 *
 *  See Appendix A of Klugman, Panjer & Willmot, Loss Models,
 *  Fourth Edition, Wiley, 2012 for the formula.
 *
 *  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "actuar.h"
#include "dpq.h"
#include "locale.h"

double betaint_raw(double x, double a, double b, double x1m)
{
/* Here, assume that (x, a, b) are not NA, 0 < x < 1 and 0 < a < Inf. */

    if (b > 0)
    {
	/* I(x, a, b) = 1 - I(1 - x, b, a) */
	double Ix = (x > 0.5) ?
	    pbeta(x1m, b, a, /*l._t.*/0, /*give_log*/0) :
	    pbeta(x,   a, b, /*l._t.*/1, /*give_log*/0);
	return gammafn(a) * gammafn(b) * Ix;
    }

    double r = floor(-b);

    if (! (ACT_nonint(b) && a - r - 1 > 0))
	return R_NaN;

    /* There are two quantities to accumulate in order to compute the
     * final result: the alternating sum (to be stored in 'sum') and
     * the ratio [(a - 1) ... (a - r)]/[b(b + 1) ... (b + r)] (to be
     * stored in 'ratio'). Some calculations are done in the log
     * scale. */
    int i;
    double ap = a, bp = b;	/* copies of a and b */
    double lx = log(x);		/* log(x) */
    double lx1m = log(x1m);	/* log(1 - x) */
    double x1 = exp(lx1m - lx);	/* (1 - x)/x */
    double c, tmp, sum, ratio;

    /* Computation of the first term in the alternating sum. */
    ap--;			     /* a - 1 */
    c = exp(ap * lx + bp * lx1m)/bp; /* (x^(a - 1) (1 - x)^b) / b */
    sum = c;			     /* first term */
    ratio = 1/bp;		     /* 1 / b */
    bp++;			     /* b + 1 */

    /* Other terms in the alternating sum iff r > 0.
     * Relies on the fact that each new term in the sum is
     *
     *  [previous term] * (a - i - 1)(1 - x)/[(b + i + 1) x]
     *
     * for i = 0, ..., r - 1. We need to compute this value as
     *
     *  {[previous term] * [(1 - x)/x]} * [(a - i - 1)/(b + i + 1)]
     *
     * to preserve accuracy for very small values of x (near
     * DBL_MIN).
     */
    for (i = 0; i < r; i++)
    {
	tmp = ap/bp;		/* (a - i - 1)/(b + i + 1) */
	c = tmp * (c * x1);	/* new term in the sum  */
	sum += c;
	ratio *= tmp;
	ap--;
	bp++;
    }

    /* I(x, a, b) = 1 - I(1 - x, b, a) */
    double lIx = (x > 0.5) ?
	pbeta(x1m, bp, ap, /*l._t.*/0, /*give_log*/1) :
	pbeta(x,   ap, bp, /*l._t.*/1, /*give_log*/1);

    return(-gammafn(a + b) * sum
	   + (ratio * ap) * exp(lgammafn(ap) + lgammafn(bp) + lIx));
}


/* The frontend called by actuar_do_betaint() */
double betaint(double x, double a, double b)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(a) || ISNAN(b))
	return x + a + b;
#endif
    if (!R_FINITE(a))
	return(R_PosInf);
    if (a <= 0 || x <= 0 || x >= 1)
	return R_NaN;

    return betaint_raw(x, a, b, 0.5 - x + 0.5);
}


/*
 *  R TO C INTERFACE
 *
 *  This is a streamlined version of the scheme in dpq.c
 *
 */
#define mod_iterate2(n1, n2, n3, i1, i2, i3)    \
        for (i = i1 = i2 = i3 = 0; i < n;       \
             i1 = (++i1 == n1) ? 0 : i1,        \
             i2 = (++i2 == n2) ? 0 : i2,        \
             i3 = (++i3 == n3) ? 0 : i3,        \
             ++i)

/* Function called by .External() */
SEXP actuar_do_betaint(SEXP args)
{
    SEXP sx, sa, sb, sy;
    int i, ix, ia, ib, n, nx, na, nb;
    double xi, ai, bi, *x, *a, *b, *y;
    Rboolean naflag = FALSE;

    args = CDR(args);	       /* drop function name from arguments */

    if (!isNumeric(CAR(args))|| !isNumeric(CADR(args)) || !isNumeric(CADDR(args)))
        error(_("invalid arguments"));

    nx = LENGTH(CAR(args));
    na = LENGTH(CADR(args));
    nb = LENGTH(CADDR(args));
    if ((nx == 0) || (na == 0) || (nb == 0))
        return(allocVector(REALSXP, 0));

    n = nx;
    if (n < na) n = na;
    if (n < nb) n = nb;

    PROTECT(sx = coerceVector(CAR(args), REALSXP));
    PROTECT(sa = coerceVector(CADR(args), REALSXP));
    PROTECT(sb = coerceVector(CADDR(args), REALSXP));
    PROTECT(sy = allocVector(REALSXP, n));
    x = REAL(sx);
    a = REAL(sa);
    b = REAL(sb);
    y = REAL(sy);

    mod_iterate2(nx, na, nb, ix, ia, ib)
    {
        xi = x[ix];
        ai = a[ia];
        bi = b[ib];
        if (ISNA(xi) || ISNA(ai) || ISNA(bi))
	    y[i] = NA_REAL;
        else if (ISNAN(xi) || ISNAN(ai) || ISNAN(bi))
	    y[i] = R_NaN;
        else
        {
	    y[i] = betaint(xi, ai, bi);
	    if (ISNAN(y[i])) naflag = TRUE;
        }
    }

    if (naflag)
        warning(R_MSG_NA);

    if (n == nx)
    {
        SET_ATTRIB(sy, duplicate(ATTRIB(sx)));
        SET_OBJECT(sy, OBJECT(sx));
    }
    else if (n == na)
    {
        SET_ATTRIB(sy, duplicate(ATTRIB(sa)));
        SET_OBJECT(sy, OBJECT(sa));
    }
    else if (n == nb)
    {
        SET_ATTRIB(sy, duplicate(ATTRIB(sb)));
        SET_OBJECT(sy, OBJECT(sb));
    }

    UNPROTECT(4);

    return sy;
}
