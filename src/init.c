#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP diceR_connectivity_matrix(SEXP);
extern SEXP diceR_indicator_matrix(SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"diceR_connectivity_matrix", (DL_FUNC) &diceR_connectivity_matrix, 1},
  {"diceR_indicator_matrix",    (DL_FUNC) &diceR_indicator_matrix,    1},
  {NULL, NULL, 0}
};

void R_init_diceR(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
