/*********************************************************************
 **
 ** file: scaling.c
 **
 ** Aim: Implement the scaling normalization 
 **      
 ** Copyright (C) 2003 Ben Bolstad
 **
 ** created by: B. M. Bolstad <bolstad@stat.berkeley.edu>
 ** 
 ** created on: Jul 24, 2003
 **
 ** Last modified: Jul 24, 2003
 **
 ** This file implements scaling normalization.
 ** 
 **
 ** Modification history
 **
 ** Jul 24, 2003 - Initial version
 ** Jul 27, 2003 - can now calculate trimmed means
 ** Aug 22, 2003 - can call the function from R .Call using 
 **                R_normalize_scaling
 **                Remove a pesky debug printf
 **
 *********************************************************************/

#include "scaling.h"
#include "threestep_common.h"

#include "rma_common.h"

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>
 
/*****************************************************************
 **
 ** static double mean_trim(double *x, int length, double trim)
 **
 ** double *x
 ** int length length
 ** double trim - remove what fraction from each tail
 **
 ********************************************************************/

static double mean_trim(double *x, int length, double trim){
  
  int i;
  int low, high; /* where cutoffs are */
  double sum =0.0;
  double *buffer = malloc(length*sizeof(double));

  if (trim < 0.0 || trim >= 0.5){
    error("Trying to trim the mean to much or negative value");
  }
  
  if (trim == 0.0){
    for (i=0; i < length; i++){
      sum+=x[i];  
    }    
    return (sum/(double)length);
  } else {
    
    for (i = 0; i < length; i++)
      buffer[i] = x[i];
    qsort(buffer,length,sizeof(double), (int(*)(const void*, const void*))sort_double);

    low = (int)(trim*length);
    high = length - low -1;
    
    for (i= low; i < high; i++){
      sum+=buffer[i];
    }
    free(buffer);
    return (sum/(double)(high - low +1));
  }

}


/*********************************************************************
 **
 ** void scaling_norm(double *data, int rows, int cols, double trim, int baseline)
 **
 ** double *data - data matrix, columns of which we will normalize
 ** int rows, cols - dimensions of matrix
 ** double trim  - fraction of data to trim from each tail
 ** int baseline - index of array to be used as baseline. 
 **                this will be 0..cols-1, if  it is 
 **                -1 pick array with median overall intensity as baseline
 **                -2 pick array with median median as baseline
 **                -3 generate a probewise median array for baseline
 **                -4 generate a probewise mean array for baseline
 **
 ** this function implements the baseline method of normalizing arrays
 **
 ********************************************************************/


void scaling_norm(double *data, int rows, int cols, double trim, int baseline){

  int i,j;

  double beta = 0.0;
  double mean_baseline = 0.0;
  double mean_treatment = 0.0;
  double med_intensity = 0.0;

  double *buffer;
  double *row_buffer;

  if (baseline == -1){
    /* pick the array with median overall intensity as baseline */
    buffer = Calloc(cols,double);
    for (j=0; j < cols; j++){
      for(i=0; i < rows; i++){
	buffer[j]+=data[j*rows + i];
      }
    }
    
    med_intensity = median_low(buffer, cols);
    for (j = 0; j < cols; j++){
      
      if (buffer[j] == med_intensity){
	baseline = j;
	break;
      }
    }
    Free(buffer);
    mean_baseline = mean_trim(&data[baseline*rows],rows,trim);
    
  } else if (baseline == -2){
    /* pick array with with median median intensity as baseline */
    buffer = Calloc(cols,double);
    
    for (j =0; j < cols; j++){
      buffer[j] = median(&data[j*rows], rows);
    }
    
    med_intensity = median_low(buffer,cols);

    for (j = 0; j < cols; j++){
      if (buffer[j] == med_intensity){
	baseline = j;
	break;
      }
    }    
    Free(buffer);
    mean_baseline = mean_trim(&data[baseline*rows],rows,trim);

  } else if (baseline == -3){
    /* build a synthetic array using probewise medians */
     buffer = Calloc(rows, double);
     row_buffer = Calloc(cols, double);
     for(i=0; i < rows; i++){
       for (j=0; j < cols; j++){
  	 row_buffer[j]=data[j*rows + i];
       }
       buffer[i] = median(row_buffer,cols);
     }
     mean_baseline = mean_trim(buffer,rows,trim);
     Free(buffer);
     
  } else if (baseline == -4){
    /* build a synthetic array using probewise means */
    buffer = Calloc(rows, double);

    for (i = 0; i < rows; i++){
      for (j=0; j < cols; j++){
	buffer[i] += data[j*rows + i];
      }
      buffer[i]/=(double)cols;
    }
    mean_baseline = mean_trim(buffer,rows,trim);
    Free(buffer);

  } else {

    /* the conventional method */
    
    mean_baseline = mean_trim(&data[baseline*rows],rows,trim);
  }
  

  for (j =0; j < cols; j++){
    if (j !=baseline){
      mean_treatment = mean_trim(&data[j*rows],rows,trim);
      beta = mean_baseline/mean_treatment;
      for (i = 0; i < rows; i++){
	data[j*rows + i]*=beta;
      }
    }
  }
  
}


/*********************************************************************
 **
 ** SEXP R_normalize_scaling(SEXP X,SEXP trim,SEXP baseline)
 **
 ** SEXP X          - matrix to be scale normalized
 ** SEXP trim       - fraction to trim 
 ** SEXP baseline   - index of baseline array (0 .. n-1) or a number 
 **                -1 pick array with median overall intensity as baseline
 **                -2 pick array with median median as baseline
 **                -3 generate a probewise median array for baseline
 **                -4 generate a probewise mean array for baseline
 **
 ** This function should be called from R using the .Call() interface.
 **
 ** Note that we will copy the matrix before normalizing it.
 **
 *********************************************************************/

SEXP R_normalize_scaling(SEXP X,SEXP trim,SEXP baseline){

  int rows, cols;
  SEXP Xcopy,dim1;
  double *Xptr;
  double trimvalue;
  int baselinearray;


  PROTECT(dim1 = getAttrib(X,R_DimSymbol));
  rows = INTEGER(dim1)[0];
  cols = INTEGER(dim1)[1];
  PROTECT(Xcopy = allocMatrix(REALSXP,rows,cols));
  copyMatrix(Xcopy,X,0);
  
  Xptr = NUMERIC_POINTER(AS_NUMERIC(Xcopy));
  
  trimvalue = asReal(trim);
  baselinearray = asInteger(baseline);
    

  scaling_norm(Xptr, rows, cols, trimvalue, baselinearray);
  
    
  UNPROTECT(2);
  return Xcopy;

}
