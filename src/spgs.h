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


// Include shared object defines
#include "sharedef.h"
#include "uint64.h" // for uint64 type and UINT64MAX define

#ifndef __spgs_h__
#define __spgs_h__

#ifdef SPGS_DLL // defined if being compiled as a dll
  #ifdef SPGS_EXPORT // defined if building the dll
    #define SPGS_API DLL_EXPORT
  #else // otherwise, using the dll
    #define SPGS_API DLL_IMPORT
  #endif // SPGS_EXPORT
  #define SPGS_LOCAL DLL_LOCAL
#else // SPGS_DLL is not defined: this means it is a static lib.
  #define SPGS_API
  #define SPGS_LOCAL
#endif // SPGS_DLL


#ifdef __cplusplus
extern "C" {
#endif

	// Functions
SPGS_API void GenerateMarkovSamplePath(double *dTransMat, double *dInitDist, 
int *iNumStates, double *u, int *iSamples, int *iSamplePath);

SPGS_API void PairCounts(int *pSeq, int *pn, int *piUniqueElements, 
int *piCircular, int *piPairCounts);
SPGS_API void TripleCounts(int *pSeq, int *pn, int *piUniqueElements, 
int *piCircular, int *piTripleCounts);
SPGS_API void QuadrupleCounts(int *pSeq, int *pn, int *piUniqueElements, 
int *piCircular, int *piQuadrupleCounts);
SPGS_API void CylinderCounts(int *pSeq, int *piN, int *piLags, int *piNLags,
int *piUniqueElements, int *piCircular, int *counts);
SPGS_API void Cyl2lag2Counts(int *pSeq, int *piN, int *piLags, 
int *piCircular, int *piCounts);

SPGS_API void ProbabilityNormalise(double *pdRexp, int *piN, int *piRows, 
	int *piCols, double *pdRes);
SPGS_API void ComputeEta1Statistic(double *pdRexp, int *piN, double *pdEpsilon, double *pdRes);
SPGS_API void ComputeEta2Statistic(double *pdRexp, int *piN, double *pdEpsilon, double *pdRes);

SPGS_API void CountIncreasingPairs(double *pdSeries, int *piN, 
double *pdScratch, double *pdRes, int *piOverflow);
SPGS_LOCAL uint64 CountIncreasingPairsHelper(double *pdSeries, size_t left, size_t right, 
double *pdScratch, int *piOverflow);

#ifdef __cplusplus
} // End of extern "C" block
#endif

#endif // __SPGS_h__
