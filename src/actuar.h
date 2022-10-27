#include <Rinternals.h>

/*Error messages */
#define R_MSG_NA        _("NaNs produced")

/* Interfaces to routines from package expint */
extern double(*actuar_gamma_inc)(double,double);

/* Functions accessed from .External() */
SEXP actuar_do_dpq(SEXP);
SEXP actuar_do_dpq0(int, SEXP);
SEXP actuar_do_dpq1(int, SEXP);
SEXP actuar_do_dpq2(int, SEXP);
SEXP actuar_do_dpq3(int, SEXP);
SEXP actuar_do_dpq4(int, SEXP);
SEXP actuar_do_dpq5(int, SEXP);
SEXP actuar_do_dpq6(int, SEXP);

SEXP actuar_do_random(SEXP);
SEXP actuar_do_random1(int, SEXP, SEXPTYPE);
SEXP actuar_do_random2(int, SEXP, SEXPTYPE);
SEXP actuar_do_random3(int, SEXP, SEXPTYPE);
SEXP actuar_do_random4(int, SEXP, SEXPTYPE);
SEXP actuar_do_random5(int, SEXP, SEXPTYPE);

SEXP actuar_do_dpqphtype(SEXP);
SEXP actuar_do_dpqphtype2(int, SEXP);

SEXP actuar_do_randomphtype(SEXP);
SEXP actuar_do_randomphtype2(int, SEXP, SEXPTYPE);

SEXP actuar_do_betaint(SEXP);

SEXP actuar_do_hierarc(SEXP);
SEXP actuar_do_panjer(SEXP);

/* Utility functions */
/*   Matrix algebra */
void actuar_expm(double *, int, double *);
double actuar_expmprod(double *, double *, double *, int);
void actuar_matpow(double *, int, int, double *);
void actuar_solve(double *, double *, int, int, double *);

/*   Special integrals */
double betaint(double, double, double);
double betaint_raw(double, double, double, double);

/*   Sampling */
int SampleSingleValue(int, double *);

/*   One parameter distributions */
double mexp(double, double, int);
double levexp(double, double, double, int);
double mgfexp(double, double, int);

double dinvexp(double, double, int);
double pinvexp(double, double, int, int);
double qinvexp(double, double, int, int);
double rinvexp(double);
double minvexp(double, double, int);
double levinvexp(double, double, double, int);

double dlogarithmic(double, double, int);
double plogarithmic(double, double, int, int);
double qlogarithmic(double, double, int, int);
double rlogarithmic(double);

double dztpois(double, double, int);
double pztpois(double, double, int, int);
double qztpois(double, double, int, int);
double rztpois(double);

double dztgeom(double, double, int);
double pztgeom(double, double, int, int);
double qztgeom(double, double, int, int);
double rztgeom(double);

/*   Two parameter distributions */
double munif(double, double, double, int);
double levunif(double, double, double, double, int);
double mgfunif(double, double, double, int);

double mnorm(double, double, double, int);
double mgfnorm(double, double, double, int);

double mbeta(double, double, double, int);
double levbeta(double, double, double, double, int);

double mgamma(double, double, double, int);
double levgamma(double, double, double, double, int);
double mgfgamma(double, double, double, int);

double mchisq(double, double, double, int);
double levchisq(double, double, double, double, int);
double mgfchisq(double, double, double, int);

double dinvgamma(double, double, double, int);
double pinvgamma(double, double, double, int, int);
double qinvgamma(double, double, double, int, int);
double rinvgamma(double, double);
double minvgamma(double, double, double, int);
double levinvgamma(double, double, double, double, int);
double mgfinvgamma(double, double, double, int);

double dinvparalogis(double, double, double, int);
double pinvparalogis(double, double, double, int, int);
double qinvparalogis(double, double, double, int, int);
double rinvparalogis(double, double);
double minvparalogis(double, double, double, int);
double levinvparalogis(double, double, double, double, int);

double dinvpareto(double, double, double, int);
double pinvpareto(double, double, double, int, int);
double qinvpareto(double, double, double, int, int);
double rinvpareto(double, double);
double minvpareto(double, double, double, int);
double levinvpareto(double, double, double, double, int);

double dinvweibull(double, double, double, int);
double pinvweibull(double, double, double, int, int);
double qinvweibull(double, double, double, int, int);
double rinvweibull(double, double);
double minvweibull(double, double, double, int);
double levinvweibull(double, double, double, double, int);

double dlgamma(double, double, double, int);
double plgamma(double, double, double, int, int);
double qlgamma(double, double, double, int, int);
double rlgamma(double, double);
double mlgamma(double, double, double, int);
double levlgamma(double, double, double, double, int);

double dllogis(double, double, double, int);
double pllogis(double, double, double, int, int);
double qllogis(double, double, double, int, int);
double rllogis(double, double);
double mllogis(double, double, double, int);
double levllogis(double, double, double, double, int);

double mlnorm(double, double, double, int);
double levlnorm(double, double, double, double, int);

double dparalogis(double, double, double, int);
double pparalogis(double, double, double, int, int);
double qparalogis(double, double, double, int, int);
double rparalogis(double, double);
double mparalogis(double, double, double, int);
double levparalogis(double, double, double, double, int);

double dpareto(double, double, double, int);
double ppareto(double, double, double, int, int);
double qpareto(double, double, double, int, int);
double rpareto(double, double);
double mpareto(double, double, double, int);
double levpareto(double, double, double, double, int);

double dpareto1(double, double, double, int);
double ppareto1(double, double, double, int, int);
double qpareto1(double, double, double, int, int);
double rpareto1(double, double);
double mpareto1(double, double, double, int);
double levpareto1(double, double, double, double, int);

double mweibull(double, double, double, int);
double levweibull(double, double, double, double, int);

double dgumbel(double, double, double, int);
double pgumbel(double, double, double, int, int);
double qgumbel(double, double, double, int, int);
double rgumbel(double, double);
double mgumbel(double, double, double, int);
double mgfgumbel(double, double, double, int);

double dinvgauss(double, double, double, int);
double pinvgauss(double, double, double, int, int);
double qinvgauss(double, double, double, int, int, double, int, int);
double rinvgauss(double, double);
double minvgauss(double, double, double, int);
double levinvgauss(double, double, double, double, int);
double mgfinvgauss(double, double, double, int);

double dztnbinom(double, double, double, int);
double pztnbinom(double, double, double, int, int);
double qztnbinom(double, double, double, int, int);
double rztnbinom(double, double);

double dztbinom(double, double, double, int);
double pztbinom(double, double, double, int, int);
double qztbinom(double, double, double, int, int);
double rztbinom(double, double);

double dzmlogarithmic(double, double, double, int);
double pzmlogarithmic(double, double, double, int, int);
double qzmlogarithmic(double, double, double, int, int);
double rzmlogarithmic(double, double);

double dzmpois(double, double, double, int);
double pzmpois(double, double, double, int, int);
double qzmpois(double, double, double, int, int);
double rzmpois(double, double);

double dzmgeom(double, double, double, int);
double pzmgeom(double, double, double, int, int);
double qzmgeom(double, double, double, int, int);
double rzmgeom(double, double);

double dpoisinvgauss(double, double, double, int);
double ppoisinvgauss(double, double, double, int, int);
double qpoisinvgauss(double, double, double, int, int);
double rpoisinvgauss(double, double);

/*   Three parameter distributions */
double dburr(double, double, double, double, int);
double pburr(double, double, double, double, int, int);
double qburr(double, double, double, double, int, int);
double rburr(double, double, double);
double mburr(double, double, double, double, int);
double levburr(double, double, double, double, double, int);

double dgenpareto(double, double, double, double, int);
double pgenpareto(double, double, double, double, int, int);
double qgenpareto(double, double, double, double, int, int);
double rgenpareto(double, double, double);
double mgenpareto(double, double, double, double, int);
double levgenpareto(double, double, double, double, double, int);

double dinvburr(double, double, double, double, int);
double pinvburr(double, double, double, double, int, int);
double qinvburr(double, double, double, double, int, int);
double rinvburr(double, double, double);
double minvburr(double, double, double, double, int);
double levinvburr(double, double, double, double, double, int);

double dinvtrgamma(double, double, double, double, int);
double pinvtrgamma(double, double, double, double, int, int);
double qinvtrgamma(double, double, double, double, int, int);
double rinvtrgamma(double, double, double);
double minvtrgamma(double, double, double, double, int);
double levinvtrgamma(double, double, double, double, double, int);

double dtrgamma(double, double, double, double, int);
double ptrgamma(double, double, double, double, int, int);
double qtrgamma(double, double, double, double, int, int);
double rtrgamma(double, double, double);
double mtrgamma(double, double, double, double, int);
double levtrgamma(double, double, double, double, double, int);

double dpareto2(double, double, double, double, int);
double ppareto2(double, double, double, double, int, int);
double qpareto2(double, double, double, double, int, int);
double rpareto2(double, double, double);
double mpareto2(double, double, double, double, int);
double levpareto2(double, double, double, double, double, int);

double dpareto3(double, double, double, double, int);
double ppareto3(double, double, double, double, int, int);
double qpareto3(double, double, double, double, int, int);
double rpareto3(double, double, double);
double mpareto3(double, double, double, double, int);
double levpareto3(double, double, double, double, double, int);

double dzmnbinom(double, double, double, double, int);
double pzmnbinom(double, double, double, double, int, int);
double qzmnbinom(double, double, double, double, int, int);
double rzmnbinom(double, double, double);

double dzmbinom(double, double, double, double, int);
double pzmbinom(double, double, double, double, int, int);
double qzmbinom(double, double, double, double, int, int);
double rzmbinom(double, double, double);

/*   Four parameter distributions */
double dgenbeta(double, double, double, double, double, int);
double pgenbeta(double, double, double, double, double, int, int);
double qgenbeta(double, double, double, double, double, int, int);
double rgenbeta(double, double, double, double);
double mgenbeta(double, double, double, double, double, int);
double levgenbeta(double, double, double, double, double, double, int);

double dtrbeta(double, double, double, double, double, int);
double ptrbeta(double, double, double, double, double, int, int);
double qtrbeta(double, double, double, double, double, int, int);
double rtrbeta(double, double, double, double);
double mtrbeta(double, double, double, double, double, int);
double levtrbeta(double, double, double, double, double, double, int);

double dpareto4(double, double, double, double, double, int);
double ppareto4(double, double, double, double, double, int, int);
double qpareto4(double, double, double, double, double, int, int);
double rpareto4(double, double, double, double);
double mpareto4(double, double, double, double, double, int);
double levpareto4(double, double, double, double, double, double, int);

/*   Five parameter distributions */
double dfpareto(double, double, double, double, double, double, int);
double pfpareto(double, double, double, double, double, double, int, int);
double qfpareto(double, double, double, double, double, double, int, int);
double rfpareto(double, double, double, double, double);
double mfpareto(double, double,  double, double, double, double, int);
double levfpareto(double, double, double, double, double, double, double, int);

/*   Phase-type distributions */
double dphtype(double, double *, double *, int, int);
double pphtype(double, double *, double *, int, int, int);
double rphtype(double *, double **, double *, int);
double mphtype(double, double *, double *, int, int);
double mgfphtype(double, double *, double *, int, int);

/* Definitions for the tables linking the first group of functions to
 * the second one. Tables found in names.c. One table for
 * {d,p,q,m,lev} functions and one for the {r} functions since we
 * need one more argument: the type of the result. */
typedef struct {
    char *name;
    SEXP (*cfun)(int, SEXP);
    int code;
} dpq_tab_struct;
extern dpq_tab_struct dpq_tab[];

typedef struct {
    char *name;
    SEXP (*cfun)(int, SEXP, SEXPTYPE);
    int code;
    SEXPTYPE type;
} random_tab_struct;
extern random_tab_struct random_tab[];
