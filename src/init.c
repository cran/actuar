/*
 *  Native routines registration, as per "Writing R extensions" and
 *  definition of native interface to one routine exported by package
 *  expint.
 *
 */

#include <Rinternals.h>
#include <R_ext/Visibility.h>
#include <R_ext/Rdynload.h>
#include <expintAPI.h>		/* optional; for expressiveness */
#include "actuar.h"

static const R_ExternalMethodDef ExternalEntries[] = {
    {"actuar_do_random", (DL_FUNC) &actuar_do_random, -1},
    {"actuar_do_randomphtype", (DL_FUNC) &actuar_do_randomphtype, -1},
    {"actuar_do_dpq", (DL_FUNC) &actuar_do_dpq, -1},
    {"actuar_do_dpqphtype", (DL_FUNC) &actuar_do_dpqphtype, -1},
    {"actuar_do_betaint", (DL_FUNC) &actuar_do_betaint, -1},
    {"actuar_do_hierarc", (DL_FUNC) &actuar_do_hierarc, -1},
    {"actuar_do_panjer", (DL_FUNC) &actuar_do_panjer, -1},
    {NULL, NULL, 0}
};

void attribute_visible R_init_actuar(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, NULL, ExternalEntries);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);

    /* Native interface to routine from package expint */
    actuar_gamma_inc = (double(*)(double,double)) R_GetCCallable("expint", "gamma_inc");

    /* Registration of actuar exported functions */
    /*   one parameter distributions */
    R_RegisterCCallable("actuar", "mexp", (DL_FUNC) mexp);
    R_RegisterCCallable("actuar", "levexp", (DL_FUNC) levexp);
    R_RegisterCCallable("actuar", "mgfexp", (DL_FUNC) mgfexp);

    R_RegisterCCallable("actuar", "dinvexp", (DL_FUNC) dinvexp);
    R_RegisterCCallable("actuar", "pinvexp", (DL_FUNC) pinvexp);
    R_RegisterCCallable("actuar", "qinvexp", (DL_FUNC) qinvexp);
    R_RegisterCCallable("actuar", "rinvexp", (DL_FUNC) rinvexp);
    R_RegisterCCallable("actuar", "minvexp", (DL_FUNC) minvexp);
    R_RegisterCCallable("actuar", "levinvexp", (DL_FUNC) levinvexp);

    R_RegisterCCallable("actuar", "dlogarithmic", (DL_FUNC) dlogarithmic);
    R_RegisterCCallable("actuar", "plogarithmic", (DL_FUNC) plogarithmic);
    R_RegisterCCallable("actuar", "qlogarithmic", (DL_FUNC) qlogarithmic);
    R_RegisterCCallable("actuar", "rlogarithmic", (DL_FUNC) rlogarithmic);

    R_RegisterCCallable("actuar", "dztpois", (DL_FUNC) dztpois);
    R_RegisterCCallable("actuar", "pztpois", (DL_FUNC) pztpois);
    R_RegisterCCallable("actuar", "qztpois", (DL_FUNC) qztpois);
    R_RegisterCCallable("actuar", "rztpois", (DL_FUNC) rztpois);

    R_RegisterCCallable("actuar", "dztgeom", (DL_FUNC) dztgeom);
    R_RegisterCCallable("actuar", "pztgeom", (DL_FUNC) pztgeom);
    R_RegisterCCallable("actuar", "qztgeom", (DL_FUNC) qztgeom);
    R_RegisterCCallable("actuar", "rztgeom", (DL_FUNC) rztgeom);

    /*   two parameter distributions */
    R_RegisterCCallable("actuar", "munif", (DL_FUNC) munif);
    R_RegisterCCallable("actuar", "levunif", (DL_FUNC) levunif);
    R_RegisterCCallable("actuar", "mgfunif", (DL_FUNC) mgfunif);

    R_RegisterCCallable("actuar", "mnorm", (DL_FUNC) mnorm);
    R_RegisterCCallable("actuar", "mgfnorm", (DL_FUNC) mgfnorm);

    R_RegisterCCallable("actuar", "mbeta", (DL_FUNC) mbeta);
    R_RegisterCCallable("actuar", "levbeta", (DL_FUNC) levbeta);

    R_RegisterCCallable("actuar", "mgamma", (DL_FUNC) mgamma);
    R_RegisterCCallable("actuar", "levgamma", (DL_FUNC) levgamma);
    R_RegisterCCallable("actuar", "mgfgamma", (DL_FUNC) mgfgamma);

    R_RegisterCCallable("actuar", "mchisq", (DL_FUNC) mchisq);
    R_RegisterCCallable("actuar", "levchisq", (DL_FUNC) levchisq);
    R_RegisterCCallable("actuar", "mgfchisq", (DL_FUNC) mgfchisq);

    R_RegisterCCallable("actuar", "dinvgamma", (DL_FUNC) dinvgamma);
    R_RegisterCCallable("actuar", "pinvgamma", (DL_FUNC) pinvgamma);
    R_RegisterCCallable("actuar", "qinvgamma", (DL_FUNC) qinvgamma);
    R_RegisterCCallable("actuar", "rinvgamma", (DL_FUNC) rinvgamma);
    R_RegisterCCallable("actuar", "minvgamma", (DL_FUNC) minvgamma);
    R_RegisterCCallable("actuar", "levinvgamma", (DL_FUNC) levinvgamma);
    R_RegisterCCallable("actuar", "mgfinvgamma", (DL_FUNC) mgfinvgamma);

    R_RegisterCCallable("actuar", "dinvparalogis", (DL_FUNC) dinvparalogis);
    R_RegisterCCallable("actuar", "pinvparalogis", (DL_FUNC) pinvparalogis);
    R_RegisterCCallable("actuar", "qinvparalogis", (DL_FUNC) qinvparalogis);
    R_RegisterCCallable("actuar", "rinvparalogis", (DL_FUNC) rinvparalogis);
    R_RegisterCCallable("actuar", "minvparalogis", (DL_FUNC) minvparalogis);
    R_RegisterCCallable("actuar", "levinvparalogis", (DL_FUNC) levinvparalogis);

    R_RegisterCCallable("actuar", "dinvpareto", (DL_FUNC) dinvpareto);
    R_RegisterCCallable("actuar", "pinvpareto", (DL_FUNC) pinvpareto);
    R_RegisterCCallable("actuar", "qinvpareto", (DL_FUNC) qinvpareto);
    R_RegisterCCallable("actuar", "rinvpareto", (DL_FUNC) rinvpareto);
    R_RegisterCCallable("actuar", "minvpareto", (DL_FUNC) minvpareto);
    R_RegisterCCallable("actuar", "levinvpareto", (DL_FUNC) levinvpareto);

    R_RegisterCCallable("actuar", "dinvweibull", (DL_FUNC) dinvweibull);
    R_RegisterCCallable("actuar", "pinvweibull", (DL_FUNC) pinvweibull);
    R_RegisterCCallable("actuar", "qinvweibull", (DL_FUNC) qinvweibull);
    R_RegisterCCallable("actuar", "rinvweibull", (DL_FUNC) rinvweibull);
    R_RegisterCCallable("actuar", "minvweibull", (DL_FUNC) minvweibull);
    R_RegisterCCallable("actuar", "levinvweibull", (DL_FUNC) levinvweibull);

    R_RegisterCCallable("actuar", "dlgamma", (DL_FUNC) dlgamma);
    R_RegisterCCallable("actuar", "plgamma", (DL_FUNC) plgamma);
    R_RegisterCCallable("actuar", "qlgamma", (DL_FUNC) qlgamma);
    R_RegisterCCallable("actuar", "rlgamma", (DL_FUNC) rlgamma);
    R_RegisterCCallable("actuar", "mlgamma", (DL_FUNC) mlgamma);
    R_RegisterCCallable("actuar", "levlgamma", (DL_FUNC) levlgamma);

    R_RegisterCCallable("actuar", "dllogis", (DL_FUNC) dllogis);
    R_RegisterCCallable("actuar", "pllogis", (DL_FUNC) pllogis);
    R_RegisterCCallable("actuar", "qllogis", (DL_FUNC) qllogis);
    R_RegisterCCallable("actuar", "rllogis", (DL_FUNC) rllogis);
    R_RegisterCCallable("actuar", "mllogis", (DL_FUNC) mllogis);
    R_RegisterCCallable("actuar", "levllogis", (DL_FUNC) levllogis);

    R_RegisterCCallable("actuar", "mlnorm", (DL_FUNC) mlnorm);
    R_RegisterCCallable("actuar", "levlnorm", (DL_FUNC) levlnorm);

    R_RegisterCCallable("actuar", "dparalogis", (DL_FUNC) dparalogis);
    R_RegisterCCallable("actuar", "pparalogis", (DL_FUNC) pparalogis);
    R_RegisterCCallable("actuar", "qparalogis", (DL_FUNC) qparalogis);
    R_RegisterCCallable("actuar", "rparalogis", (DL_FUNC) rparalogis);
    R_RegisterCCallable("actuar", "mparalogis", (DL_FUNC) mparalogis);
    R_RegisterCCallable("actuar", "levparalogis", (DL_FUNC) levparalogis);

    R_RegisterCCallable("actuar", "dpareto", (DL_FUNC) dpareto);
    R_RegisterCCallable("actuar", "ppareto", (DL_FUNC) ppareto);
    R_RegisterCCallable("actuar", "qpareto", (DL_FUNC) qpareto);
    R_RegisterCCallable("actuar", "rpareto", (DL_FUNC) rpareto);
    R_RegisterCCallable("actuar", "mpareto", (DL_FUNC) mpareto);
    R_RegisterCCallable("actuar", "levpareto", (DL_FUNC) levpareto);

    R_RegisterCCallable("actuar", "dpareto1", (DL_FUNC) dpareto1);
    R_RegisterCCallable("actuar", "ppareto1", (DL_FUNC) ppareto1);
    R_RegisterCCallable("actuar", "qpareto1", (DL_FUNC) qpareto1);
    R_RegisterCCallable("actuar", "rpareto1", (DL_FUNC) rpareto1);
    R_RegisterCCallable("actuar", "mpareto1", (DL_FUNC) mpareto1);
    R_RegisterCCallable("actuar", "levpareto1", (DL_FUNC) levpareto1);

    R_RegisterCCallable("actuar", "mweibull", (DL_FUNC) mweibull);
    R_RegisterCCallable("actuar", "levweibull", (DL_FUNC) levweibull);

    /* R_RegisterCCallable("actuar", "minvGauss", (DL_FUNC) minvGauss); [defunct v3.0-0] */
    /* R_RegisterCCallable("actuar", "levinvGauss", (DL_FUNC) levinvGauss); [idem] */
    /* R_RegisterCCallable("actuar", "mgfinvGauss", (DL_FUNC) mgfinvGauss); [idem] */

    R_RegisterCCallable("actuar", "dgumbel", (DL_FUNC) dgumbel);
    R_RegisterCCallable("actuar", "pgumbel", (DL_FUNC) pgumbel);
    R_RegisterCCallable("actuar", "qgumbel", (DL_FUNC) qgumbel);
    R_RegisterCCallable("actuar", "rgumbel", (DL_FUNC) rgumbel);
    R_RegisterCCallable("actuar", "mgumbel", (DL_FUNC) mgumbel);
    R_RegisterCCallable("actuar", "mgfgumbel", (DL_FUNC) mgfgumbel);

    R_RegisterCCallable("actuar", "dinvgauss", (DL_FUNC) dinvgauss);
    R_RegisterCCallable("actuar", "pinvgauss", (DL_FUNC) pinvgauss);
    R_RegisterCCallable("actuar", "qinvgauss", (DL_FUNC) qinvgauss);
    R_RegisterCCallable("actuar", "rinvgauss", (DL_FUNC) rinvgauss);
    R_RegisterCCallable("actuar", "minvgauss", (DL_FUNC) minvgauss);
    R_RegisterCCallable("actuar", "levinvgauss", (DL_FUNC) levinvgauss);
    R_RegisterCCallable("actuar", "mgfinvgauss", (DL_FUNC) mgfinvgauss);

    R_RegisterCCallable("actuar", "dztnbinom", (DL_FUNC) dztnbinom);
    R_RegisterCCallable("actuar", "pztnbinom", (DL_FUNC) pztnbinom);
    R_RegisterCCallable("actuar", "qztnbinom", (DL_FUNC) qztnbinom);
    R_RegisterCCallable("actuar", "rztnbinom", (DL_FUNC) rztnbinom);

    R_RegisterCCallable("actuar", "dztbinom", (DL_FUNC) dztbinom);
    R_RegisterCCallable("actuar", "pztbinom", (DL_FUNC) pztbinom);
    R_RegisterCCallable("actuar", "qztbinom", (DL_FUNC) qztbinom);
    R_RegisterCCallable("actuar", "rztbinom", (DL_FUNC) rztbinom);

    R_RegisterCCallable("actuar", "dzmlogarithmic", (DL_FUNC) dzmlogarithmic);
    R_RegisterCCallable("actuar", "pzmlogarithmic", (DL_FUNC) pzmlogarithmic);
    R_RegisterCCallable("actuar", "qzmlogarithmic", (DL_FUNC) qzmlogarithmic);
    R_RegisterCCallable("actuar", "rzmlogarithmic", (DL_FUNC) rzmlogarithmic);

    R_RegisterCCallable("actuar", "dzmpois", (DL_FUNC) dzmpois);
    R_RegisterCCallable("actuar", "pzmpois", (DL_FUNC) pzmpois);
    R_RegisterCCallable("actuar", "qzmpois", (DL_FUNC) qzmpois);
    R_RegisterCCallable("actuar", "rzmpois", (DL_FUNC) rzmpois);

    R_RegisterCCallable("actuar", "dzmgeom", (DL_FUNC) dzmgeom);
    R_RegisterCCallable("actuar", "pzmgeom", (DL_FUNC) pzmgeom);
    R_RegisterCCallable("actuar", "qzmgeom", (DL_FUNC) qzmgeom);
    R_RegisterCCallable("actuar", "rzmgeom", (DL_FUNC) rzmgeom);

    R_RegisterCCallable("actuar", "dpoisinvgauss", (DL_FUNC) dpoisinvgauss);
    R_RegisterCCallable("actuar", "ppoisinvgauss", (DL_FUNC) ppoisinvgauss);
    R_RegisterCCallable("actuar", "qpoisinvgauss", (DL_FUNC) qpoisinvgauss);
    R_RegisterCCallable("actuar", "rpoisinvgauss", (DL_FUNC) rpoisinvgauss);

    /*   three parameter distributions */
    R_RegisterCCallable("actuar", "dburr", (DL_FUNC) dburr);
    R_RegisterCCallable("actuar", "pburr", (DL_FUNC) pburr);
    R_RegisterCCallable("actuar", "qburr", (DL_FUNC) qburr);
    R_RegisterCCallable("actuar", "rburr", (DL_FUNC) rburr);
    R_RegisterCCallable("actuar", "mburr", (DL_FUNC) mburr);
    R_RegisterCCallable("actuar", "levburr", (DL_FUNC) levburr);

    R_RegisterCCallable("actuar", "dgenpareto", (DL_FUNC) dgenpareto);
    R_RegisterCCallable("actuar", "pgenpareto", (DL_FUNC) pgenpareto);
    R_RegisterCCallable("actuar", "qgenpareto", (DL_FUNC) qgenpareto);
    R_RegisterCCallable("actuar", "rgenpareto", (DL_FUNC) rgenpareto);
    R_RegisterCCallable("actuar", "mgenpareto", (DL_FUNC) mgenpareto);
    R_RegisterCCallable("actuar", "levgenpareto", (DL_FUNC) levgenpareto);

    R_RegisterCCallable("actuar", "dinvburr", (DL_FUNC) dinvburr);
    R_RegisterCCallable("actuar", "pinvburr", (DL_FUNC) pinvburr);
    R_RegisterCCallable("actuar", "qinvburr", (DL_FUNC) qinvburr);
    R_RegisterCCallable("actuar", "rinvburr", (DL_FUNC) rinvburr);
    R_RegisterCCallable("actuar", "minvburr", (DL_FUNC) minvburr);
    R_RegisterCCallable("actuar", "levinvburr", (DL_FUNC) levinvburr);

    R_RegisterCCallable("actuar", "dinvtrgamma", (DL_FUNC) dinvtrgamma);
    R_RegisterCCallable("actuar", "pinvtrgamma", (DL_FUNC) pinvtrgamma);
    R_RegisterCCallable("actuar", "qinvtrgamma", (DL_FUNC) qinvtrgamma);
    R_RegisterCCallable("actuar", "rinvtrgamma", (DL_FUNC) rinvtrgamma);
    R_RegisterCCallable("actuar", "minvtrgamma", (DL_FUNC) minvtrgamma);
    R_RegisterCCallable("actuar", "levinvtrgamma", (DL_FUNC) levinvtrgamma);

    R_RegisterCCallable("actuar", "dtrgamma", (DL_FUNC) dtrgamma);
    R_RegisterCCallable("actuar", "ptrgamma", (DL_FUNC) ptrgamma);
    R_RegisterCCallable("actuar", "qtrgamma", (DL_FUNC) qtrgamma);
    R_RegisterCCallable("actuar", "rtrgamma", (DL_FUNC) rtrgamma);
    R_RegisterCCallable("actuar", "mtrgamma", (DL_FUNC) mtrgamma);
    R_RegisterCCallable("actuar", "levtrgamma", (DL_FUNC) levtrgamma);

    R_RegisterCCallable("actuar", "dpareto2", (DL_FUNC) dpareto2);
    R_RegisterCCallable("actuar", "ppareto2", (DL_FUNC) ppareto2);
    R_RegisterCCallable("actuar", "qpareto2", (DL_FUNC) qpareto2);
    R_RegisterCCallable("actuar", "rpareto2", (DL_FUNC) rpareto2);
    R_RegisterCCallable("actuar", "mpareto2", (DL_FUNC) mpareto2);
    R_RegisterCCallable("actuar", "levpareto2", (DL_FUNC) levpareto2);

    R_RegisterCCallable("actuar", "dpareto3", (DL_FUNC) dpareto3);
    R_RegisterCCallable("actuar", "ppareto3", (DL_FUNC) ppareto3);
    R_RegisterCCallable("actuar", "qpareto3", (DL_FUNC) qpareto3);
    R_RegisterCCallable("actuar", "rpareto3", (DL_FUNC) rpareto3);
    R_RegisterCCallable("actuar", "mpareto3", (DL_FUNC) mpareto3);
    R_RegisterCCallable("actuar", "levpareto3", (DL_FUNC) levpareto3);

    R_RegisterCCallable("actuar", "dzmnbinom", (DL_FUNC) dzmnbinom);
    R_RegisterCCallable("actuar", "pzmnbinom", (DL_FUNC) pzmnbinom);
    R_RegisterCCallable("actuar", "qzmnbinom", (DL_FUNC) qzmnbinom);
    R_RegisterCCallable("actuar", "rzmnbinom", (DL_FUNC) rzmnbinom);

    R_RegisterCCallable("actuar", "dzmbinom", (DL_FUNC) dzmbinom);
    R_RegisterCCallable("actuar", "pzmbinom", (DL_FUNC) pzmbinom);
    R_RegisterCCallable("actuar", "qzmbinom", (DL_FUNC) qzmbinom);
    R_RegisterCCallable("actuar", "rzmbinom", (DL_FUNC) rzmbinom);

    /*   four parameter distributions */
    R_RegisterCCallable("actuar", "dtrbeta", (DL_FUNC) dtrbeta);
    R_RegisterCCallable("actuar", "ptrbeta", (DL_FUNC) ptrbeta);
    R_RegisterCCallable("actuar", "qtrbeta", (DL_FUNC) qtrbeta);
    R_RegisterCCallable("actuar", "rtrbeta", (DL_FUNC) rtrbeta);
    R_RegisterCCallable("actuar", "mtrbeta", (DL_FUNC) mtrbeta);
    R_RegisterCCallable("actuar", "levtrbeta", (DL_FUNC) levtrbeta);

    R_RegisterCCallable("actuar", "dgenbeta", (DL_FUNC) dgenbeta);
    R_RegisterCCallable("actuar", "pgenbeta", (DL_FUNC) pgenbeta);
    R_RegisterCCallable("actuar", "qgenbeta", (DL_FUNC) qgenbeta);
    R_RegisterCCallable("actuar", "rgenbeta", (DL_FUNC) rgenbeta);
    R_RegisterCCallable("actuar", "mgenbeta", (DL_FUNC) mgenbeta);
    R_RegisterCCallable("actuar", "levgenbeta", (DL_FUNC) levgenbeta);

    R_RegisterCCallable("actuar", "dpareto4", (DL_FUNC) dpareto4);
    R_RegisterCCallable("actuar", "ppareto4", (DL_FUNC) ppareto4);
    R_RegisterCCallable("actuar", "qpareto4", (DL_FUNC) qpareto4);
    R_RegisterCCallable("actuar", "rpareto4", (DL_FUNC) rpareto4);
    R_RegisterCCallable("actuar", "mpareto4", (DL_FUNC) mpareto4);
    R_RegisterCCallable("actuar", "levpareto4", (DL_FUNC) levpareto4);

    /*   five parameter distributions */
    R_RegisterCCallable("actuar", "dfpareto", (DL_FUNC) dfpareto);
    R_RegisterCCallable("actuar", "pfpareto", (DL_FUNC) pfpareto);
    R_RegisterCCallable("actuar", "qfpareto", (DL_FUNC) qfpareto);
    R_RegisterCCallable("actuar", "rfpareto", (DL_FUNC) rfpareto);
    R_RegisterCCallable("actuar", "mfpareto", (DL_FUNC) mfpareto);
    R_RegisterCCallable("actuar", "levfpareto", (DL_FUNC) levfpareto);

    /*   phase-type distributions */
    R_RegisterCCallable("actuar", "dphtype", (DL_FUNC) dphtype);
    R_RegisterCCallable("actuar", "pphtype", (DL_FUNC) pphtype);
    R_RegisterCCallable("actuar", "rphtype", (DL_FUNC) rphtype);
    R_RegisterCCallable("actuar", "mphtype", (DL_FUNC) mphtype);
    R_RegisterCCallable("actuar", "mgfphtype", (DL_FUNC) mgfphtype);

    /*   special integrals */
    R_RegisterCCallable("actuar", "betaint", (DL_FUNC) betaint);
}

/* Define imports from package expint */
double(*actuar_gamma_inc)(double,double);
