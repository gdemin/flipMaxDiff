#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP flipMaxDiff_densitiesP(SEXP);
extern SEXP flipMaxDiff_gradientBestWorst(SEXP, SEXP, SEXP, SEXP);
extern SEXP flipMaxDiff_logDensitiesBestWorst(SEXP, SEXP);
extern SEXP flipMaxDiff_logDensityBestWorst(SEXP, SEXP);
extern SEXP flipMaxDiff_logKernels(SEXP, SEXP, SEXP);
extern SEXP flipMaxDiff_sumWeightedOuterProducts(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"flipMaxDiff_densitiesP",               (DL_FUNC) &flipMaxDiff_densitiesP,               1},
    {"flipMaxDiff_gradientBestWorst",        (DL_FUNC) &flipMaxDiff_gradientBestWorst,        4},
    {"flipMaxDiff_logDensitiesBestWorst",    (DL_FUNC) &flipMaxDiff_logDensitiesBestWorst,    2},
    {"flipMaxDiff_logDensityBestWorst",      (DL_FUNC) &flipMaxDiff_logDensityBestWorst,      2},
    {"flipMaxDiff_logKernels",               (DL_FUNC) &flipMaxDiff_logKernels,               3},
    {"flipMaxDiff_sumWeightedOuterProducts", (DL_FUNC) &flipMaxDiff_sumWeightedOuterProducts, 2},
    {NULL, NULL, 0}
};

void R_init_flipMaxDiff(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
