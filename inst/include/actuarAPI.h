/*  actuar: Actuarial Functions and Heavy Tailed Distributions
 *
 *  Support for exported functions at the C level.
 *
 *  This is derived from code in package zoo.
 *
 *  Copyright (C) 2020 Vincent Goulet
 *  Copyright (C) 2010 Jeffrey A. Ryan jeff.a.ryan @ gmail.com
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License as
 *  published by the Free Software Foundation; either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 *  02110-1301, USA.
 *
 *  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#ifdef  __cplusplus
extern "C" {
#endif

/* Special integrals */
double betaint(double x, double a, double b, int foo);
double betaint_raw(double x, double a, double b, double x1m);

/* One parameter distributions */
double mexp(double order, double scale, int give_log);
double levexp(double limit, double scale, double order, int give_log);
double mgfexp(double t, double scale, int give_log);

double dinvexp(double x, double scale, int give_log);
double pinvexp(double q, double scale, int lower_tail, int log_p);
double qinvexp(double p, double scale, int lower_tail, int log_p);
double rinvexp(double scale);
double minvexp(double order, double scale, int give_log);
double levinvexp(double limit, double scale, double order, int give_log);

double dlogarithmic(double x, double p, int give_log);
double plogarithmic(double x, double p, int lower_tail, int log_p);
double qlogarithmic(double x, double p, int lower_tail, int log_p);
double rlogarithmic(double p);

double dztpois(double x, double lambda, int give_log);
double pztpois(double q, double lambda, int lower_tail, int log_p);
double qztpois(double p, double lambda, int lower_tail, int log_p);
double rztpois(double lambda);

double dztgeom(double x, double prob, int give_log);
double pztgeom(double q, double prob, int lower_tail, int log_p);
double qztgeom(double p, double prob, int lower_tail, int log_p);
double rztgeom(double prob);

/* Two parameter distributions */
double munif(double order, double min, double max, int give_log);
double levunif(double limit, double min, double max, double order, int give_log);
double mgfunif(double t, double min, double max, int give_log);

double mnorm(double order, double mean, double sd, int give_log);
double mgfnorm(double t, double mean, double sd, int give_log);

double mbeta(double order, double shape1, double shape2, int give_log);
double levbeta(double limit, double shape1, double shape2, double order, int give_log);

double mgamma(double order, double shape, double scale, int give_log);
double levgamma(double limit, double shape, double scale, double order, int give_log);
double mgfgamma(double t, double shape, double scale, int give_log);

double mchisq(double order, double df, double ncp, int give_log);
double levchisq(double limit, double df, double ncp, double order, int give_log);
double mgfchisq(double t, double df, double ncp, int give_log);

double dinvgamma(double x, double scale, double shape, int give_log);
double pinvgamma(double q, double scale, double shape, int lower_tail, int log_p);
double qinvgamma(double p, double scale, double shape, int lower_tail, int log_p);
double rinvgamma(double scale, double shape);
double minvgamma(double order, double scale, double shape, int give_log);
double levinvgamma(double limit, double scale, double shape, double order, int give_log);
double mgfinvgamma(double t, double shape, double scale, int give_log);

double dinvparalogis(double x, double shape, double scale, int give_log);
double pinvparalogis(double q, double shape, double scale, int lower_tail, int log_p);
double qinvparalogis(double p, double shape, double scale, int lower_tail, int log_p);
double rinvparalogis(double shape, double scale);
double minvparalogis(double order, double shape, double scale, int give_log);
double levinvparalogis(double limit, double shape, double scale, double order, int give_log);

double dinvpareto(double x, double shape, double scale, int give_log);
double pinvpareto(double q, double shape, double scale, int lower_tail, int log_p);
double qinvpareto(double p, double shape, double scale, int lower_tail, int log_p);
double rinvpareto(double shape, double scale);
double minvpareto(double order, double shape, double scale, int give_log);
double levinvpareto(double limit, double shape, double scale, double order, int log_p);

double dinvweibull(double x, double scale, double shape, int give_log);
double pinvweibull(double q, double scale, double shape, int lower_tail, int log_p);
double qinvweibull(double p, double scale, double shape, int lower_tail, int log_p);
double rinvweibull(double scale, double shape);
double minvweibull(double order, double scale, double shape, int give_log);
double levinvweibull(double limit, double scale, double shape, double order, int give_log);

double dlgamma(double x, double shapelog, double ratelog, int give_log);
double plgamma(double q, double shapelog, double ratelog, int lower_tail, int log_p);
double qlgamma(double p, double shapelog, double ratelog, int lower_tail, int log_p);
double rlgamma(double ratelog, double shapelog);
double mlgamma(double order, double shapelog, double ratelog, int give_log);
double levlgamma(double limit, double shapelog, double ratelog, double order, int give_log);

double dllogis(double x, double shape, double scale, int give_log);
double pllogis(double q, double shape, double scale, int lower_tail, int log_p);
double qllogis(double p, double shape, double scale, int lower_tail, int log_p);
double rllogis(double shape, double scale);
double mllogis(double order, double shape, double scale, int give_log);
double levllogis(double limit, double shape, double scale, double order, int give_log);

double mlnorm(double order, double logmean, double logsd, int give_log);
double levlnorm(double limit, double logmean, double logsd, double order, int give_log);

double dparalogis(double x, double shape, double scale, int give_log);
double pparalogis(double q, double shape, double scale, int lower_tail, int log_p);
double qparalogis(double p, double shape, double scale, int lower_tail, int log_p);
double rparalogis(double shape, double scale);
double mparalogis(double order, double shape, double scale, int give_log);
double levparalogis(double limit, double shape, double scale, double order, int give_log);

double dpareto(double x, double shape, double scale, int give_log);
double ppareto(double q, double shape, double scale, int lower_tail, int log_p);
double qpareto(double p, double shape, double scale, int lower_tail, int log_p);
double rpareto(double shape, double scale);
double mpareto(double order, double shape, double scale, int give_log);
double levpareto(double limit, double shape, double scale, double order, int give_log);

double dpareto1(double x, double shape, double scale, int give_log);
double ppareto1(double q, double shape, double scale, int lower_tail, int log_p);
double qpareto1(double p, double shape, double scale, int lower_tail, int log_p);
double rpareto1(double shape, double scale);
double mpareto1(double order, double shape, double scale, int give_log);
double levpareto1(double limit, double shape, double scale, double order, int give_log);

double mweibull(double order, double scale, double shape, int give_log);
double levweibull(double limit, double scale, double shape, double order, int give_log);

double dgumbel(double x, double alpha, double beta, int give_log);
double pgumbel(double q, double alpha, double beta, int lower_tail, int log_p);
double qgumbel(double p, double alpha, double beta, int lower_tail, int log_p);
double rgumbel(double alpha, double beta);
double mgumbel(double order, double alpha, double beta, int give_log);
double mgfgumbel(double t, double alpha, double beta, int give_log);

double dinvgauss(double x, double mu, double phi, int give_log);
double pinvgauss(double q, double mu, double phi, int lower_tail, int log_p);
double qinvgauss(double q, double mu, double phi, int lower_tail, int log_p,
		 double tol, int maxit, int echo);
double rinvgauss(double mu, double phi);
double minvgauss(double order, double mean, double phi, int give_log);
double levinvgauss(double limit, double mean, double phi, double order, int give_log);
double mgfinvgauss(double t, double mean, double phi, int give_log);

double dztnbinom(double x, double size, double prob, int give_log);
double pztnbinom(double q, double size, double prob, int lower_tail, int log_p);
double qztnbinom(double p, double size, double prob, int lower_tail, int log_p);
double rztnbinom(double size, double prob);

double dztbinom(double x, double size, double prob, int give_log);
double pztbinom(double q, double size, double prob, int lower_tail, int log_p);
double qztbinom(double p, double size, double prob, int lower_tail, int log_p);
double rztbinom(double size, double prob);

double dzmlogarithmic(double x, double p, double p0m, int give_log);
double pzmlogarithmic(double x, double p, double p0m, int lower_tail, int log_p);
double qzmlogarithmic(double x, double p, double p0m, int lower_tail, int log_p);
double rzmlogarithmic(double p, double p0m);

double dzmpois(double x, double lambda, double p0m, int give_log);
double pzmpois(double q, double lambda, double p0m, int lower_tail, int log_p);
double qzmpois(double p, double lambda, double p0m, int lower_tail, int log_p);
double rzmpois(double lambda, double p0m);

double dzmgeom(double x, double prob, double p0m, int give_log);
double pzmgeom(double q, double prob, double p0m, int lower_tail, int log_p);
double qzmgeom(double p, double prob, double p0m, int lower_tail, int log_p);
double rzmgeom(double prob, double p0m);

double dpoisinvgauss(double x, double mu, double phi, int give_log);
double ppoisinvgauss(double q, double mu, double phi, int lower_tail, int log_p);
double qpoisinvgauss(double p, double mu, double phi, int lower_tail, int log_p);
double rpoisinvgauss(double mu, double phi);

/* Three parameter distributions */
double dburr(double x, double shape1, double shape2, double scale, int give_log);
double pburr(double q, double shape1, double shape2, double scale, int lower_tail, int log_p);
double qburr(double p, double shape1, double shape2, double scale, int lower_tail, int log_p);
double rburr(double shape1, double shape2, double scale);
double mburr(double order, double shape1, double shape2, double scale, int give_log);
double levburr(double limit, double shape1, double shape2, double scale, double order, int give_log);

double dgenpareto(double x, double shape1, double shape2, double scale, int give_log);
double pgenpareto(double q, double shape1, double shape2, double scale, int lower_tail, int log_p);
double qgenpareto(double p, double shape1, double shape2, double scale, int lower_tail, int log_p);
double rgenpareto(double shape1, double shape2, double scale);
double mgenpareto(double order, double shape1, double shape2, double scale, int give_log);
double levgenpareto(double limit, double shape1, double shape2, double scale, double order, int give_log);

double dinvburr(double x, double shape1, double shape2, double scale, int give_log);
double pinvburr(double q, double shape1, double shape2, double scale, int lower_tail, int log_p);
double qinvburr(double p, double shape1, double shape2, double scale, int lower_tail, int log_p);
double rinvburr(double shape1, double shape2, double scale);
double minvburr(double order, double shape1, double shape2, double scale, int give_log);
double levinvburr(double limit, double shape1, double shape2, double scale, double order, int give_log);

double dinvtrgamma(double x, double shape1, double shape2, double scale, int give_log);
double pinvtrgamma(double q, double shape1, double shape2, double scale, int lower_tail, int log_p);
double qinvtrgamma(double p, double shape1, double shape2, double scale, int lower_tail, int log_p);
double rinvtrgamma(double shape1, double shape2, double scale);
double minvtrgamma(double order, double shape1, double shape2, double scale, int give_log);
double levinvtrgamma(double limit, double shape1, double shape2, double scale, double order, int give_log);

double dtrgamma(double x, double shape1, double shape2, double scale, int give_log);
double ptrgamma(double q, double shape1, double shape2, double scale, int lower_tail, int log_p);
double qtrgamma(double p, double shape1, double shape2, double scale, int lower_tail, int log_p);
double rtrgamma(double shape1, double shape2, double scale);
double mtrgamma(double order, double shape1, double shape2, double scale, int give_log);
double levtrgamma(double limit, double shape1, double shape2, double scale, double order, int give_log);

double dpareto2(double x, double min, double shape, double scale, int give_log);
double ppareto2(double q, double min, double shape, double scale, int lower_tail, int log_p);
double qpareto2(double p, double min, double shape, double scale, int lower_tail, int log_p);
double rpareto2(double min, double shape, double scale);
double mpareto2(double order, double min, double shape, double scale, int give_log);
double levpareto2(double limit, double min, double shape, double scale, double order, int give_log);

double dpareto3(double x, double min, double shape, double scale, int give_log);
double ppareto3(double q, double min, double shape, double scale, int lower_tail, int log_p);
double qpareto3(double p, double min, double shape, double scale, int lower_tail, int log_p);
double rpareto3(double min, double shape, double scale);
double mpareto3(double order, double min, double shape, double scale, int give_log);
double levpareto3(double limit, double min, double shape, double scale, double order, int give_log);

double dzmnbinom(double x, double size, double prob, double p0m, int give_log);
double pzmnbinom(double q, double size, double prob, double p0m, int lower_tail, int log_p);
double qzmnbinom(double p, double size, double prob, double p0m, int lower_tail, int log_p);
double rzmnbinom(double size, double prob, double p0m);

double dzmbinom(double x, double size, double prob, double p0m, int give_log);
double pzmbinom(double q, double size, double prob, double p0m, int lower_tail, int log_p);
double qzmbinom(double p, double size, double prob, double p0m, int lower_tail, int log_p);
double rzmbinom(double size, double prob, double p0m);

/* Four parameter distributions */
double dgenbeta(double x, double shape1, double shape2, double shape3, double scale, int give_log);
double pgenbeta(double q, double shape1, double shape2, double shape3, double scale, int lower_tail, int log_p);
double qgenbeta(double p, double shape1, double shape2, double shape3, double scale, int lower_tail, int log_p);
double rgenbeta(double shape1, double shape2, double shape3, double scale);
double mgenbeta(double order, double shape1, double shape2, double shape3, double scale, int give_log);
double levgenbeta(double limit, double shape1, double shape2, double shape3, double scale, double order, int give_log);

double dtrbeta(double x, double shape1, double shape2, double shape3, double scale, int give_log);
double ptrbeta(double q, double shape1, double shape2, double shape3, double scale, int lower_tail, int log_p);
double qtrbeta(double p, double shape1, double shape2, double shape3, double scale, int lower_tail, int log_p);
double rtrbeta(double shape1, double shape2, double shape3, double scale);
double mtrbeta(double order, double shape1, double shape2, double shape3, double scale, int give_log);
double levtrbeta(double limit, double shape1, double shape2, double shape3, double scale, double order, int give_log);

double dpareto4(double x, double min, double shape1, double shape2, double scale, int give_log);
double ppareto4(double q, double min, double shape1, double shape2, double scale, int lower_tail, int log_p);
double qpareto4(double p, double min, double shape1, double shape2, double scale, int lower_tail, int log_p);
double rpareto4(double min, double shape1, double shape2, double scale);
double mpareto4(double order, double min, double shape1, double shape2, double scale, int give_log);
double levpareto4(double limit, double min, double shape1, double shape2, double scale, double order, int give_log);

/* Five parameter distributions */
double dfpareto(double x, double min, double shape1, double shape2, double shape3, double scale, int give_log);
double pfpareto(double q, double min, double shape1, double shape2, double shape3, double scale, int lower_tail, int log_p);
double qfpareto(double p, double min, double shape1, double shape2, double shape3, double scale, int lower_tail, int log_p);
double rfpareto(double min, double shape1, double shape2, double shape3, double scale);
double mfpareto(double order, double min,  double shape1, double shape2, double shape3, double scale, int give_log);
double levfpareto(double limit, double mu, double shape1, double shape2, double shape3, double scale, double order, int give_log);

/* Phase-type distributions */
double dphtype(double x, double *pi, double *T, int m, int give_log);
double pphtype(double x, double *pi, double *T, int m, int lower_tail, int log_p);
double rphtype(double *pi, double **Q, double *rates, int m);
double mphtype(double order, double *pi, double *T, int m, int give_log);
double mgfphtype(double x, double *pi, double *T, int m, int give_log);

#ifdef  __cplusplus
}
#endif
