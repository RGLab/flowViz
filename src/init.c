/* Copyright Bioconductor Foundation of NA, 2007, all rights reserved */

#include "R.h"
#include "R_ext/Rdynload.h"

#include "flowViz.h"

static const R_CallMethodDef CallEntries[] = {
    {"binHex", (DL_FUNC) &binHex, 8},
    {NULL, NULL, 0}
};

void R_init_flowViz(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
}

