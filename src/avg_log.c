/************************************************************************
 **
 ** avg_log.c
 **
 ** created by: B. M. Bolstad   <bolstad@stat.berkeley.edu>
 ** created on: Jan 7, 2002  (but based on earlier work from Nov 2002)
 **
 ** Copyright (C) 2002-2003 Ben Bolstad
 **
 ** last modified: Jan 7, 2003
 **
 ** License: GPL V2 or later (same as the rest of the Affy package)
 **
 ** General discussion
 **
 ** Implement avgerage log2 pm summarization, with or without normalization
 **
 ** Nov 2, 2002 - modify so that it will work efficently with affy2
 ** Jan 3, 2003 - Clean up commenting, prepare for integration in AffyExtensions
 ** Jan 7, 2003 - Make function standalone, to prepare for later combination into
 **               a more general framework.
 ** Jul 23, 2003 - add parameter for computing SE and SE implemented
 **
 ************************************************************************/

#include "qnorm.h"

#include <R.h> 
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/***************************************************************************
 **
 ** double AvgLog(double *x, int length)
 **
 ** double *x - a vector of PM intensities 
 ** int length - length of *x
 **
 ** take the average of log2 PM intensities.
 **
 ***************************************************************************/

double AvgLog(double *x, int length){
  int i;
  double sum = 0.0;
  double mean = 0.0;

  for (i=0; i < length; i++){
    sum = sum + x[i];
  }
  
  mean = sum/(double)length;

  return (mean);    
}

/***************************************************************************
 **
 ** double AvgLog(double *x, int length)
 **
 ** double *x - a vector of PM intensities 
 ** int length - length of *x
 **
 ** take the average of log2 PM intensities.
 **
 ***************************************************************************/

double AvgLogSE(double *x, double mean, int length){
  int i;
  double sum = 0.0;

  for (i=0; i < length; i++){
    sum = sum + (x[i]- mean)*(x[i] - mean);
  }
  
  sum = sqrt(sum/(double)length);

  return (sum);    
}


/***************************************************************************
 **
 ** double AverageLog(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes)
 **
 ** aim: given a data matrix of probe intensities, and a list of rows in the matrix 
 **      corresponding to a single probeset, compute average log2 expression measure. 
 **      Note that data is a probes by chips matrix.
 **
 ** double *data - Probe intensity matrix
 ** int rows - number of rows in matrix *data (probes)
 ** int cols - number of cols in matrix *data (chips)
 ** int *cur_rows - indicies of rows corresponding to current probeset
 ** double *results - already allocated location to store expression measures (cols length)
 ** int nprobes - number of probes in current probeset.
 **
 ***************************************************************************/

void AverageLog(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes, double *resultsSE){
  int i,j;
  double *z = Calloc(nprobes*cols,double);

  for (j = 0; j < cols; j++){
    for (i =0; i < nprobes; i++){
      z[j*nprobes + i] = log(data[j*rows + cur_rows[i]])/log(2.0);  
    }
  } 
  
  for (j=0; j < cols; j++){
    results[j] = AvgLog(&z[j*nprobes],nprobes);
    resultsSE[j] = AvgLogSE(&z[j*nprobes],results[j],nprobes);
  }




  Free(z);
}


/***************************************************************************
 **
 ** double AverageLog(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes)
 **
 ** aim: given a data matrix of probe intensities, and a list of rows in the matrix 
 **      corresponding to a single probeset, compute average log2 expression measure. 
 **      Note that data is a probes by chips matrix.
 **
 ** double *data - Probe intensity matrix
 ** int rows - number of rows in matrix *data (probes)
 ** int cols - number of cols in matrix *data (chips)
 ** int *cur_rows - indicies of rows corresponding to current probeset
 ** double *results - already allocated location to store expression measures (cols length)
 ** int nprobes - number of probes in current probeset.
 **
 ***************************************************************************/

void AverageLog_noSE(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes){
  int i,j;
  double *z = Calloc(nprobes*cols,double);

  for (j = 0; j < cols; j++){
    for (i =0; i < nprobes; i++){
      z[j*nprobes + i] = log(data[j*rows + cur_rows[i]])/log(2.0);  
    }
  } 
  
  for (j=0; j < cols; j++){
    results[j] = AvgLog(&z[j*nprobes],nprobes);
  }
  Free(z);
}

