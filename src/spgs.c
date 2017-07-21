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



// spgs.c
// Source code to C routines used by the SPGS package


#include <stdio.h>
#include <stdlib.h> // for malloc and free
#include <math.h> // for fabs
#include <stdbool.h>
#include "spgs.h"


// Macro defines
#define INDEX(i,j) ((i)+(j)*s)
#define PROBS(i,j) (*(pdCurMat+(i-1)*4 + j-1))



void attribute_hidden GenerateMarkovSamplePath(double *dTransMat, double *dInitDist, 
int *iNumStates, double *dU, int *iSamples, int *iSamplePath)
{
	int i, j;
	int s = *iNumStates;
	double *dInitCdf = (double*)malloc(sizeof(double)*(size_t)s); // reserve space for initial cdf
	double *dTransCdf = (double*)malloc(sizeof(double)*(size_t)(s*s)); // reserve space for state-dependent transition cdfs

	// Prepare cdf of initial distribution and conditional state transitions
	dInitCdf[0] = 0.0;
	for (i=0; i<s; i++)
		dTransCdf[INDEX(i,0)] = 0.0;
	for (j=1; j<s; j++)
	{
		dInitCdf[j] = dInitCdf[j-1]+dInitDist[j-1];
		for (i=0; i<s; i++)
 dTransCdf[INDEX(i,j)] = dTransCdf[INDEX(i,j-1)]+dTransMat[INDEX(i,j-1)];
	} // for j

// Generate initial state
		j = 1;
		while(j<s && dInitCdf[j]<dU[0]) j++;
		iSamplePath[0] = j;

// Generate rest of sample path
	for (i=1; i<*iSamples; i++)
	{
		j = 1;
		while(j<s && dTransCdf[INDEX(iSamplePath[i-1]-1, j)]<dU[i]) j++;
		iSamplePath[i] = j;
	} // for i

// dispose of memory reserved for the various cdfs
	free(dTransCdf);
	free(dInitCdf);
} //function


void attribute_hidden PairCounts(int *pSeq, int *pn, int *piUniqueElements, int *piCircular, int *piPairCounts)
{
// Declare variables
	int *ps, *psn;
	int i, ind1, ind2;
	int n=*pn, iUniqueElements=*piUniqueElements, adjN;
	int iTotalElements = iUniqueElements*iUniqueElements;
	bool circular = *piCircular!=0;

// Initialise count array
	for (i=0; i<iTotalElements; i++) piPairCounts[i] = 0;
	psn = pSeq+n;
// Adjust the apparent length according to whether or not the sequence is being treated as scircular
	adjN = circular ? n : n-1;

// Count cylinders in the sequence
	ind2 = (*pSeq)-1; // get first element of sequence in preparation for iteration
	for (i=0, ps=pSeq+1; i<adjN; i++, ps++)
	{
		ind1 = ind2;
  		if (ps>=psn) ps -= n;
		ind2 = (*ps)-1;
		piPairCounts[ind1+ind2*iUniqueElements]++;
	} // for i
} // function


void attribute_hidden TripleCounts(int *pSeq, int *pn, int *piUniqueElements, 
int *piCircular, int*piTripleCounts)
{
// Declare variables
  int *ps, *psn;
  int i, ind1, ind2, ind3;
  int n=*pn, iUniqueElements=*piUniqueElements, adjN;
  int iUniqueElements2=iUniqueElements*iUniqueElements;
  int iTotalElements = iUniqueElements*iUniqueElements2;
  bool circular = *piCircular!=0;

// Initialise count array
  for (i=0; i<iTotalElements; i++) piTripleCounts[i] = 0;
  psn = pSeq+n;
// Adjust the apparent length according to whether or not the sequence is being treated as scircular
	adjN = circular ? n : n-2;

  // Count cylinders in the sequence
  ind2 = (*pSeq)-1;
  ind3 = (*(pSeq+1))-1;
	for (i=0, ps=pSeq+2; i<adjN; i++, ps++)
	{
		ind1 = ind2;
		ind2 = ind3;
		if (ps>=psn) ps -= n;
		ind3 = (*ps)-1;
		piTripleCounts[ind1 + ind2*iUniqueElements+ind3*iUniqueElements2]++;
	} // for i
} // function


void attribute_hidden QuadrupleCounts(int *pSeq, int *pn, int *piUniqueElements, 
int *piCircular, int *piQuadrupleCounts)
{
// Declare variables
  int *ps, *psn;
  int i, ind1, ind2, ind3, ind4;
  int n=*pn, iUniqueElements=*piUniqueElements, adjN;
  int iUniqueElements2=iUniqueElements*iUniqueElements;
  int iUniqueElements3 = iUniqueElements*iUniqueElements2;
  int iTotalElements = iUniqueElements*iUniqueElements3;
  bool circular = *piCircular!=0;

// Initialise count array
  for (i=0; i<iTotalElements; i++) piQuadrupleCounts[i] = 0;
  psn = pSeq+n;
// Adjust the apparent length according to whether or not the sequence is being treated as scircular
	adjN = circular ? n : n-3;

  // Count cylinders in the sequence
  ind2 = (*pSeq)-1;
  ind3 = (*(pSeq+1))-1;
  ind4 = (*(pSeq+2))-1;
	for (i=0, ps=pSeq+3; i<adjN; i++, ps++)
	{
		ind1 = ind2;
		ind2 = ind3;
		ind3 = ind4;
		if (ps>=psn) ps -= n;
		ind4 = (*ps)-1;
		piQuadrupleCounts[ind1 + ind2*iUniqueElements+ind3*iUniqueElements2+ind4*iUniqueElements3]++;
	} // for i
} // function


void attribute_hidden CylinderCounts(int *pSeq, int *piN, int *piLags, int *piNLags,
int *piUniqueElements, int *piCircular, int *counts)
{
// Declare variables
  int *psn;
  int i, j, iVal, index;
  int n=*piN, adjN, nLags=*piNLags;
  int *offset=(int*)malloc(sizeof(int)*(size_t)nLags);
  offset[0] = 1;
  for (j=1; j<nLags; j++) offset[j] = offset[j-1] * *piUniqueElements;
  int iTotalElements = offset[nLags-1] * *piUniqueElements;
  int **ps = (int**)malloc(sizeof(int*)*(size_t)nLags);
  bool circular = *piCircular!=0;

// Initialise count array
  for (j=0; j<iTotalElements; j++) counts[j] = 0;
  psn = pSeq+n;
// Adjust the apparent length according to whether or not the sequence is being treated as scircular
	adjN = circular ? n : n-piLags[nLags-1]+1;
  if (piLags[nLags-1]>adjN) return;

  // Count cylinders in the sequence
  	for (j=0; j<nLags; j++)
  ps[j] = pSeq+piLags[j]-1;
	for (i=0; i<adjN; i++)
	{
		index = 0;
		for (j=0; j<nLags; j++)
		{
			iVal = (*ps[j])-1;
			ps[j]++;
			if (ps[j]>=psn) ps[j] -= n;
				index += iVal*offset[j];
		} // for j
		counts[index]++;
	} // for i
	
// Clean up dynamically reserved temporary storage
	free(ps);
	free(offset);
} // function


void attribute_hidden Cyl2lag2Counts(int *pSeq, int *piN, int *piLags, 
int *piCircular, int *piCounts)
{
// Declare variables
  int *ps1, *ps2, *ps3, *ps4, *psn;
  int i, ind1, ind2, ind3, ind4, lag, lags;
  int n = *piN, adjN;
  int workSpace = 256*((*piLags)+1);
  bool circular = *piCircular!=0;

// Initialise count array
	for (i=0; i<workSpace; i++) piCounts[i] = 0;
	psn = pSeq+n;
// Adjust the apparent length according to whether or not the sequence is being treated as scircular
	adjN = circular ? n : n-3;

// Count cylinders in the sequence
	for (i=0, ps1=pSeq, ps2=pSeq+1; i<adjN; i++, ps1++, ps2++)
	{
		if (ps2>=psn) ps2 -= n;
		ind1 = (*ps1)-1;
		ind2 = (*ps2)-1;
		lags = *piLags;
		if (!circular && lags>adjN+1-i) lags = adjN+1-i;
		for (lag=0, ps3=ps1, ps4=ps1+1; lag<=lags; lag++, ps3++, ps4++)
		{
			if (ps3>=psn) ps3 -= n;
			if (ps4>=psn) ps4 -= n;
			ind3 = (*ps3)-1;
			ind4 = (*ps4)-1;
			piCounts[ind1+4*ind2+16*ind3+64*ind4+256*lag]++;
		} // for lag
	} // for i
} // function


void attribute_hidden ProbabilityNormalise(double *pdRexp, int *piN, int *piRows, 
	int *piCols, double *pdRes)
{
		int r, c, n;
		double sum;
		int N = *piN;
		int iRows = *piRows;
		int iCols = *piCols;
		int Nr = N*iRows;
		for (n=0; n<N; n++)
		for (r=0; r<iRows; r++)
		{
			sum = 0.0;
			for (c=0; c<iCols; c++) sum += pdRexp[n+N*r+Nr*c];
			for (c=0; c<iCols; c++) pdRes[n+N*r+Nr*c] = pdRexp[n+N*r+Nr*c]/sum;
		} // for r, n
} // function


void attribute_hidden ComputeEta1Statistic(double *pdRexp, int *piN, double *pdEpsilon, double *pdRes)
{
	int r, c, n; // counters
 double sum, stat; // for accumulating row sums and holding the test statistic
 double f1, f2; // temp vars for computing the test statistic
	double *pdRow, *pdCol; // pointers to current row and element within row for probability-normalising the rows of each 4 X 4 matrix
	double *pdCurMat = pdRexp; // pointer to current 4 X 4 matrix
	int iLessEpsilon = 0; // for counting the number of simulated test statistics that a smaller than epsilon
// Do some argument checking
	if (pdRes==NULL) return;
	if (*piN<=0) {*pdRes = 0.0; return;}
	if (*pdEpsilon<=0.0) {*pdRes = 0.0; return;}
	if (*pdEpsilon>=1.0) {*pdRes = 1.0; return;}
// Simulate n test statistics
	for (n=0; n<*piN; n++)
	{
	// Generate a stochastic matrix (i.e. one whose rows all sum to  unity)
		pdRow = pdCurMat; // get first row of matrix
		for (r=0; r<4; r++)
		{
			pdCol = pdRow; // get first element of row
			sum = 0.0; // set row sum to 0
			for (c=0; c<4; c++) {sum += *pdCol; pdCol++;}
			pdCol = pdRow; // reset to first element of row
			if (sum==0.0)
				for (c=0; c<4; c++) {*pdCol = 0.0; pdCol++;}
			else
				for (c=0; c<4; c++) {*pdCol /= sum; pdCol++;}
			pdRow += 4;
		} // for r
	// Compute the statistic
		if (PROBS(2,1)>0.0 || PROBS(3,1)>0.0)
		{
			f1 = (1-(PROBS(1,1)+PROBS(4,1)))*(1-(PROBS(2,2)+PROBS(3,2)))/(PROBS(2,1)+PROBS(3,1))-PROBS(1,2);
			f2 = (1-(PROBS(2,3)+PROBS(3,3)))/(1-(PROBS(2,2)+PROBS(3,2))) * (PROBS(1,2)+f1) - PROBS(1,3);
			if (PROBS(1,1)+PROBS(4,1)<1.0 && PROBS(2,2)+PROBS(3,2)<1.0 && PROBS(2,3)+PROBS(3,3)<1.0
				&& f1>0.0 && f2>0.0 && PROBS(4,1)+f1+f2<1.0)
			{
				f1 -= PROBS(4,2);
				f2 -= PROBS(4,3);
				stat = fabs(f1);
				stat = ((sum=fabs(f2))>stat) ? sum : stat;
			}
			else // 
				stat = 1.0;
		}
		else // PROBS(2,1)<=0.0 || PROBS(3,1)<=0.0
			stat = 1.0;
		pdCurMat += 16; // move to next matrix
	// if the test statistic <= epsilon, bump iLessEpsilon
		if (stat<=*pdEpsilon) iLessEpsilon++;
	} // for n
// Compute the estimated p-value
	*pdRes = (double)iLessEpsilon/(double)*piN;
} // function


void attribute_hidden ComputeEta2Statistic(double *pdRexp, int *piN, double *pdEpsilon, double *pdRes)
{
	int r, c, n; // counters
 double sum, stat; // for accumulating row sums and holding the test statistic
 double p24, f[5]; // temp vars for computing the test statistic
	double *pdRow, *pdCol; // pointers to current row and element within row for probability-normalising the rows of each 4 X 4 matrix
	double *pdCurMat = pdRexp; // pointer to current 4 X 4 matrix
	int iLessEpsilon = 0; // for counting the number of simulated test statistics that a smaller than epsilon
// Do some argument checking
	if (pdRes==NULL) return;
	if (*piN<=0) {*pdRes = 0.0; return;}
	if (*pdEpsilon<=0.0) {*pdRes = 0.0; return;}
	if (*pdEpsilon>=1.0) {*pdRes = 1.0; return;}
// Simulate n test statistics
	for (n=0; n<*piN; n++)
	{
	// Generate a stochastic matrix (i.e. one whose rows all sum to  unity)
		pdRow = pdCurMat; // get first row of matrix
		for (r=0; r<4; r++)
		{
			pdCol = pdRow; // get first element of row
			sum = 0.0; // set row sum to 0
			for (c=0; c<4; c++) {sum += *pdCol; pdCol++;}
			pdCol = pdRow; // reset to first element of row
			if (sum==0.0)
				for (c=0; c<4; c++) {*pdCol = 0.0; pdCol++;}
			else
				for (c=0; c<4; c++) {*pdCol /= sum; pdCol++;}
			pdRow += 4;
		} // for r
	// Compute the statistic
		p24 = 1.0 - (PROBS(2,1)+PROBS(2,2)+PROBS(2,3));
		if (p24>0.0 && PROBS(1,3)>0.0)
		{
			f[0] = PROBS(2,2);
			f[1] = 1.0 - PROBS(3,1) - PROBS(2,2) - p24*PROBS(1,2)/PROBS(1,3);
			f[2] = PROBS(2,1)*PROBS(1,3)/p24;
			f[3] = PROBS(3,1)*PROBS(1,3)/p24;
			f[4] = 1.0 - (PROBS(1,1) + f[2] + f[3]); //-PROBS(1,1) - (PROBS(2,1)+PROBS(3,1))*PROBS(1,3)/p24;
			if (f[1]>0.0 && f[4]>0.0)
			{
				f[0] -= PROBS(3,3);
				f[1] -= PROBS(3,2);
				f[2] -= PROBS(4,3);
				f[3] -= PROBS(4,2);
				f[4] -= PROBS(4,1);
				stat = fabs(f[0]);
				for (c=1; c<=4; c++)
					stat = ((sum=fabs(f[c]))>stat) ? sum : stat;
			}
			else // f[1]==0.0 || f[4]==0.0
				stat = 1.0;
		}
		else // p24==0.0 || PROBS(1,3)==0.0
			stat = 1.0;
		pdCurMat += 16; // move to next matrix
	// if the test statistic <= epsilon, bump iLessEpsilon
		if (stat<=*pdEpsilon) iLessEpsilon++;
 } // for n
// Compute the estimated p-value
	*pdRes = (double)iLessEpsilon/(double)*piN;
} // function


uint64 attribute_hidden CountIncreasingPairsHelper(double *pdSeries, size_t left, size_t right, 
double *pdScratch, int *piOverflow)
{
	*piOverflow = 0; // initialise overflow flag to 0
	if (left+1==right)
		return(0); // base case: no swaps can occur with one element
	else
	{
		uint64 inOrder, mergeInOrder;
		int i = 0;
		int length = right-left;
		int mid = length/2;
		int l = left, r = left + mid;

// Sort the left and right subarrays
		inOrder = CountIncreasingPairsHelper(pdSeries, left, left + mid, pdScratch, piOverflow);
		mergeInOrder = CountIncreasingPairsHelper(pdSeries, left+mid, right, pdScratch, piOverflow);
		if (inOrder>UINT64MAX-mergeInOrder)
		{
			*piOverflow = 1; // set overflow state
			return(UINT64MAX); // return max value
		}
		inOrder += mergeInOrder; // no overflow, just add

// Merge the left and right subarrays together using pdScratch as temporary storage
		mergeInOrder = 0;
		for (i=0; i<length; i++)
		{
			if (l<left+mid && (r==right || pdSeries[l]<=pdSeries[r]))
			{
				pdScratch[i] = pdSeries[l];
				l++;
				if (inOrder>UINT64MAX-mergeInOrder)
				{
				 *piOverflow = 1; // set overflow state
				 return(UINT64MAX); // return max value
				}
				mergeInOrder += right-r; // no overflow, just add
			}
			else
			{
				pdScratch[i] = pdSeries[r];
				r++;
			} // if
		} // for i
		if (inOrder>UINT64MAX-mergeInOrder)
		{
			*piOverflow = 1; // set overflow state
			return(UINT64MAX); // return max value
		}
		inOrder += mergeInOrder; // no overflow, just add

// Copy the result of the merge back to pdSeries
		for (i=left; i<right; i++) pdSeries[i] = pdScratch[i-left];
		return(inOrder);
	} // if
} // function

void attribute_hidden CountIncreasingPairs(double *pdSeries, int *piN, 
double *pdScratch, double *pdRes, int *piOverflow)
{
	*pdRes = (double)CountIncreasingPairsHelper(pdSeries, 0, *piN, pdScratch, piOverflow);
} // function

