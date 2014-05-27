/*
	taken from teh util.c file from the xcms package
http://www.bioconductor.org/packages/release/bioc/src/contrib/xcms_1.38.0.tar.gz

	this file remains as long as the DescendMin function is not exported in xcms.

	Christian Panse <cp@fgcz.ethz.ch>
	Sat Mar 29 10:04:24 CET 2014

	
*/

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void DescendMin(double *yvals, int *numin, int *istart,
                int *ilower, int *iupper) {

    int i;

    for (i = *istart; i > 0; i--)
        if (yvals[i-1] >= yvals[i])
            break;
    *ilower = i;

    for (i = *istart; i < *numin-1; i++)
        if (yvals[i+1] >= yvals[i])
            break;
    *iupper = i;
}


//static const R_CallMethodDef cMethods[] = { 
static const R_CMethodDef CEntries[] = { 
    {"DescendMin", (DL_FUNC) &DescendMin, 5},
    {NULL, NULL, 0} 
};

void R_init_cosmiq(DllInfo *info)
{
       R_registerRoutines(info, CEntries, NULL, NULL, NULL);
}

