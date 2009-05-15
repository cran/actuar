/*  ===== actuar: An R Package for Actuarial Science =====
 *
 *  Functions to calculate raw and limited moments for the Beta
 *  distribution. See ../R/BetaMoments.R for details.
 *
 *  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double mbeta(double order, double shape1, double shape2, int give_log)
{
    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(order)  ||
        shape1 <= 0.0 ||
        shape2 <= 0.0)
        return R_NaN;

    if (order <= -shape1)
	return R_PosInf;

    return beta(shape1 + order, shape2) / beta(shape1, shape2);
}

double levbeta(double limit, double shape1, double shape2, double order,
                int give_log)
{
    double tmp;

    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(order) ||
        shape1 <= 0.0 ||
        shape2 <= 0.0)
        return R_NaN;

    if (order <= -shape1)
	return R_PosInf;

    if (limit <= 0.0)
        return 0.0;

    tmp = order + shape1;

    return beta(tmp, shape2) / beta(shape1, shape2) *
        pbeta(limit, tmp, shape2, 1, 0) +
        R_VG__0(limit, order) * pbeta(limit, shape1, shape2, 0, 0);
}
