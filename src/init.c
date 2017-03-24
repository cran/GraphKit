#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP GraphKit_bootstrap(SEXP, SEXP, SEXP);
extern SEXP GraphKit_debias(SEXP, SEXP);
extern SEXP GraphKit_skipDownBipar(SEXP, SEXP, SEXP, SEXP);
extern SEXP GraphKit_skipDownChain(SEXP, SEXP, SEXP, SEXP);
extern SEXP GraphKit_skipDownChainCI(SEXP, SEXP, SEXP);
extern SEXP GraphKit_skipDownChrom(SEXP, SEXP, SEXP, SEXP);
extern SEXP GraphKit_skipDownChromCI(SEXP, SEXP, SEXP);
extern SEXP GraphKit_skipDownClique(SEXP, SEXP, SEXP, SEXP);
extern SEXP GraphKit_skipDownCliqueCI(SEXP, SEXP, SEXP);
extern SEXP GraphKit_skipDownConn(SEXP, SEXP, SEXP, SEXP);
extern SEXP GraphKit_skipDownConnCI(SEXP, SEXP, SEXP);
extern SEXP GraphKit_skipDownCycle(SEXP, SEXP, SEXP, SEXP);
extern SEXP GraphKit_skipDownDeg(SEXP, SEXP, SEXP, SEXP);
extern SEXP GraphKit_skipDownDegCI(SEXP, SEXP, SEXP);
extern SEXP GraphKit_skipDownGirth(SEXP, SEXP, SEXP, SEXP);
extern SEXP GraphKit_skipDownGirthCI(SEXP, SEXP, SEXP);
extern SEXP GraphKit_skipDownMatch(SEXP, SEXP, SEXP, SEXP);
extern SEXP GraphKit_skipDownPlan(SEXP, SEXP, SEXP, SEXP);
extern SEXP GraphKit_skipDownSingle(SEXP, SEXP, SEXP, SEXP);
extern SEXP GraphKit_skipDownSingleCI(SEXP, SEXP, SEXP);
extern SEXP GraphKit_stepDown(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"GraphKit_bootstrap",        (DL_FUNC) &GraphKit_bootstrap,        3},
    {"GraphKit_debias",           (DL_FUNC) &GraphKit_debias,           2},
    {"GraphKit_skipDownBipar",    (DL_FUNC) &GraphKit_skipDownBipar,    4},
    {"GraphKit_skipDownChain",    (DL_FUNC) &GraphKit_skipDownChain,    4},
    {"GraphKit_skipDownChainCI",  (DL_FUNC) &GraphKit_skipDownChainCI,  3},
    {"GraphKit_skipDownChrom",    (DL_FUNC) &GraphKit_skipDownChrom,    4},
    {"GraphKit_skipDownChromCI",  (DL_FUNC) &GraphKit_skipDownChromCI,  3},
    {"GraphKit_skipDownClique",   (DL_FUNC) &GraphKit_skipDownClique,   4},
    {"GraphKit_skipDownCliqueCI", (DL_FUNC) &GraphKit_skipDownCliqueCI, 3},
    {"GraphKit_skipDownConn",     (DL_FUNC) &GraphKit_skipDownConn,     4},
    {"GraphKit_skipDownConnCI",   (DL_FUNC) &GraphKit_skipDownConnCI,   3},
    {"GraphKit_skipDownCycle",    (DL_FUNC) &GraphKit_skipDownCycle,    3},
    {"GraphKit_skipDownDeg",      (DL_FUNC) &GraphKit_skipDownDeg,      4},
    {"GraphKit_skipDownDegCI",    (DL_FUNC) &GraphKit_skipDownDegCI,    3},
    {"GraphKit_skipDownGirth",    (DL_FUNC) &GraphKit_skipDownGirth,    4},
    {"GraphKit_skipDownGirthCI",  (DL_FUNC) &GraphKit_skipDownGirthCI,  3},
    {"GraphKit_skipDownMatch",    (DL_FUNC) &GraphKit_skipDownMatch,    4},
    {"GraphKit_skipDownPlan",     (DL_FUNC) &GraphKit_skipDownPlan,     4},
    {"GraphKit_skipDownSingle",   (DL_FUNC) &GraphKit_skipDownSingle,   4},
    {"GraphKit_skipDownSingleCI", (DL_FUNC) &GraphKit_skipDownSingleCI, 3},
    {"GraphKit_stepDown",         (DL_FUNC) &GraphKit_stepDown,         3},
    {NULL, NULL, 0}
};

void R_init_GraphKit(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
