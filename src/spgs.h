/*
This file is part of the source code for
SPGS: an R package for identifying statistical patterns in genomic sequences.
Copyright (C) 2015  Universidad de Chile and INRIA-Chile

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of Version 2 of the GNU Public License is available in the 
share/licenses/gpl-2 file in the R installation directory or from 
http://www.R-project.org/Licenses/GPL-2.
*/



// spgs.h
// Header file for the C++ routines used by the SPGS package

#ifndef __spgs_h__
#define __spgs_h__

//#include "sharedef.h"
#include "uint64.h" // for uint64 type and UINT64MAX define
#include <R_ext/Rdynload.h>
#include <Rdefines.h>
#include <R_ext/Visibility.h>	

#ifdef __cplusplus
extern "C" {
#endif

	// Functions
void attribute_hidden GenerateMarkovSamplePath(double *dTransMat, double *dInitDist, 
int *iNumStates, double *u, int *iSamples, int *iSamplePath);

void attribute_hidden PairCounts(int *pSeq, int *pn, int *piUniqueElements, 
int *piCircular, int *piPairCounts);
void attribute_hidden TripleCounts(int *pSeq, int *pn, int *piUniqueElements, 
int *piCircular, int *piTripleCounts);
void attribute_hidden QuadrupleCounts(int *pSeq, int *pn, int *piUniqueElements, 
int *piCircular, int *piQuadrupleCounts);
void attribute_hidden CylinderCounts(int *pSeq, int *piN, int *piLags, int *piNLags,
int *piUniqueElements, int *piCircular, int *counts);
void attribute_hidden Cyl2lag2Counts(int *pSeq, int *piN, int *piLags, 
int *piCircular, int *piCounts);

void attribute_hidden ProbabilityNormalise(double *pdRexp, int *piN, int *piRows, 
	int *piCols, double *pdRes);
void attribute_hidden ComputeEta1Statistic(double *pdRexp, int *piN, double *pdEpsilon, double *pdRes);
void attribute_hidden ComputeEta2Statistic(double *pdRexp, int *piN, double *pdEpsilon, double *pdRes);
void attribute_hidden CountIncreasingPairs(double *pdSeries, int *piN, 
double *pdScratch, double *pdRes, int *piOverflow);
uint64 attribute_hidden CountIncreasingPairsHelper(double *pdSeries, size_t left, size_t right, 
double *pdScratch, int *piOverflow);


	// R registration
	
static R_NativePrimitiveArgType GMSP_T[] = {
    REALSXP, REALSXP, INTSXP, REALSXP, INTSXP, INTSXP
};

static R_NativePrimitiveArgType INTS5_T[] = {
    INTSXP, INTSXP, INTSXP, INTSXP, INTSXP
};

static R_NativePrimitiveArgType INTS7_T[] = {
    INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP
};

static R_NativePrimitiveArgType ProbabilityNormalise_T[] = {
    REALSXP, INTSXP, INTSXP, INTSXP, REALSXP
};

static R_NativePrimitiveArgType ComputeEta_T[] = {
    REALSXP, INTSXP, REALSXP, REALSXP
};

static R_NativePrimitiveArgType CountIncreasingPairs_T[] = {
    REALSXP, INTSXP, REALSXP, REALSXP, INTSXP
};

static const R_CMethodDef cMethods[] = {
  {"GenerateMarkovSamplePath", (DL_FUNC)GenerateMarkovSamplePath, 6, GMSP_T},
  {"PairCounts", (DL_FUNC)PairCounts, 5, INTS5_T},
  {"TripleCounts", (DL_FUNC)TripleCounts, 5, INTS5_T},
  {"QuadrupleCounts", (DL_FUNC)QuadrupleCounts, 5, INTS5_T},
  {"CylinderCounts", (DL_FUNC)CylinderCounts, 7, INTS7_T},
  {"Cyl2lag2Counts", (DL_FUNC)Cyl2lag2Counts, 5, INTS5_T},
  {"ProbabilityNormalise", (DL_FUNC)ProbabilityNormalise, 5, ProbabilityNormalise_T},
  {"ComputeEta1Statistic", (DL_FUNC)ComputeEta1Statistic, 4, ComputeEta_T},
  {"ComputeEta2Statistic", (DL_FUNC)ComputeEta2Statistic, 4, ComputeEta_T},
  {"CountIncreasingPairs", (DL_FUNC)CountIncreasingPairs, 5, CountIncreasingPairs_T},
  {NULL, NULL, 0, NULL}
};

void attribute_visible R_init_spgs(DllInfo *dll)
{
    R_registerRoutines(dll, cMethods, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
	}
  
void attribute_hidden R_unload_spgs(DllInfo *info)
{
// Nothing needs to be done to unload this dll
}

#ifdef __cplusplus
} // End of extern "C" block
#endif

#endif // __SPGS_h__
