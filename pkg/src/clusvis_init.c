#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP ClusVis_computeGradientCPP(SEXP, SEXP, SEXP, SEXP);
extern SEXP ClusVis_computeLikelihoodCPP(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"ClusVis_computeGradientCPP",   (DL_FUNC) &ClusVis_computeGradientCPP,   4},
    {"ClusVis_computeLikelihoodCPP", (DL_FUNC) &ClusVis_computeLikelihoodCPP, 4},
    {NULL, NULL, 0}
};

void R_init_ClusVis(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
