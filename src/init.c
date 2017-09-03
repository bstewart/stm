#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP stm_gradcpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP stm_hpbcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP stm_lhoodcpp(SEXP, SEXP, SEXP, SEXP, SEXP);
