/************************************************************************
 **
 ** file: medianpolish.c
 **
 ** Copyright (C) 2002-2003 Ben Bolstad
 **
 ** created by: B. M. Bolstad   <bolstad@stat.berkeley.edu>
 ** created on: Jan 7, 2003 (but based on code dating back as far as June 2002)
 **
 ** last modified: Jan 7, 2003
 **
 ** License: GPL V2 or later (same as the rest of the Affy package)
 **
 ** Median polish summary measure (used in the RMA expression measure)
 **
 **
 ** History
 **
 ** Jan 7, 2003 - Initial version to fit into the three-step framework.
 ** Jan 13, 2003 - move median() into threestep_common.c
 ** Feb 24, 2003 - make maxiter get used.
 ** Jul 23, 2003 - add ability to accept SE parameter
 ** Sept 13, 2003 - introduced medianpolishPLM which returns 
 **                 most of what is required by the fitting
 **                 algorithm
 ** Oct 05, 2003 - added in summary_param
 ** Apr 5, 2004 - change malloc/free to Calloc/Free
 ** Nov 13, 2006 - make median calls to median_nocopy
 **
 ************************************************************************/


#include "medianpolish.h"
#include "rma_common.h"
#include "threestep_common.h"

#include <R.h> 
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/*******************************************************************************
 **
 ** double sum_abs(double *z, int rows, int cols)
 **
 ** double *z - matrix of doubles
 ** int rows - dimension of matrix
 ** int cols - dimension of matrix
 **
 ** returns the sum of the absolute values of elements of the matrix *z
 **
 ******************************************************************************/

double sum_abs(double *z, int rows, int cols){
 
  int i, j;
  double sum = 0.0;

  for (i=0; i < rows; i++)
    for (j=0; j < cols; j++)
      sum+=fabs(z[j*rows+i]);

  return sum;
}

/********************************************************************************
 **
 ** void get_row_median(double *z, double *rdelta, int rows, int cols)
 **
 ** double *z - matrix of dimension  rows*cols
 ** double *rdelta - on output will contain row medians (vector of length rows)
 ** int rows, cols - dimesion of matrix
 **
 ** get the row medians of a matrix 
 **
 ********************************************************************************/

void get_row_median(double *z, double *rdelta, int rows, int cols){
  int i,j;
  double *buffer = (double *)Calloc(cols,double);

  for (i = 0; i < rows; i++){ 
    for (j = 0; j < cols; j++){
      buffer[j] = z[j*rows + i];
    }
    rdelta[i] = median_nocopy(buffer,cols);
  }
  
  Free(buffer);
}

/********************************************************************************
 **
 ** void get_col_median(double *z, double *cdelta, int rows, int cols)
 **
 ** double *z - matrix of dimension  rows*cols
 ** double *cdelta - on output will contain col medians (vector of length cols)
 ** int rows, cols - dimesion of matrix
 **
 ** get the col medians of a matrix 
 **
 ********************************************************************************/

void get_col_median(double *z, double *cdelta, int rows, int cols){
  
  int i, j;
  
  double *buffer = (double *)Calloc(rows,double);
  for (j = 0; j < cols; j++){
    for (i = 0; i < rows; i++){  
      buffer[i] = z[j*rows + i];
    }
    cdelta[j] = median_nocopy(buffer,rows);
  }
  
  Free(buffer);

}

/***********************************************************************************
 **
 ** void subtract_by_row(double *z, double *rdelta, int rows, int cols)
 ** 
 ** double *z - matrix of dimension rows by cols
 ** double *rdelta - vector of length rows
 ** int rows, cols dimensions of matrix
 **
 ** subtract the elements of *rdelta off each row of *z
 **
 ***********************************************************************************/

void subtract_by_row(double *z, double *rdelta, int rows, int cols){
  
  int i,j;

  for (i = 0; i < rows; i++){
    for (j = 0; j < cols; j++){
      z[j*rows +i]-= rdelta[i];
    }
  }
}


/***********************************************************************************
 **
 ** void subtract_by_col(double *z, double *cdelta, int rows, int cols)
 ** 
 ** double *z - matrix of dimension rows by cols
 ** double *cdelta - vector of length rows
 ** int rows, cols dimensions of matrix
 **
 ** subtract the elements of *cdelta off each col of *z
 **
 ***********************************************************************************/

void subtract_by_col(double *z, double *cdelta, int rows, int cols){
  
  int i,j;
  for (j = 0; j < cols; j++){
    for (i = 0; i < rows; i++){
      z[j*rows +i]-= cdelta[j];
    }
  }

}

/***********************************************************************************
 **
 ** void rmod(double *r, double *rdelta, int rows)
 ** 
 ** double *r - vector of length rows
 ** double *rdelta - vector of length rows
 ** int rows, cols dimensions of matrix
 **
 ** add elementwise *rdelta to *r
 **
 ***********************************************************************************/


void rmod(double *r, double *rdelta, int rows){
  int i;

  for (i = 0; i < rows; i++){
    r[i]= r[i] + rdelta[i];
  }
}

/***********************************************************************************
 **
 ** void cmod(double *c, double *cdelta, int cols)
 ** 
 ** double *c - vector of length rows
 ** double *cdelta - vector of length rows
 ** int cols length of vector
 **
 ** add elementwise *cdelta to *c
 **
 ***********************************************************************************/

void cmod(double *c, double *cdelta, int cols){
  int j;

  for (j = 0; j < cols; j++){
    c[j]= c[j] + cdelta[j];
  }
}


/*************************************************************************************
 **
 ** void median_polish(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes)
 **
 ** double *data - a data matrix of dimension rows by cols (the entire PM matrix)
 ** int rows, cols - rows and columns dimensions of matrix
 ** int cur_rows - vector of length nprobes containg row indicies of *data matrix which apply for a 
 **                particular probeset
 ** double *results - a vector of length cols already allocated. on output contains expression values
 ** int nprobes - number of probes in current probeset.
 **
 ** a function to do median polish expression summary.
 **
 *************************************************************************************/

void median_polish(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes, double *resultsSE, summary_plist *summary_param){

  int i,j,iter;
  int maxiter = 10;
  double eps=0.01;
  double oldsum = 0.0,newsum = 0.0;
  double t = 0.0;
  double delta;
  double *rdelta = Calloc(nprobes,double);
  double *cdelta = Calloc(cols,double);
  
  double *r = Calloc(nprobes,double);
  double *c = Calloc(cols,double);
  double *z = Calloc(nprobes*cols,double);

  for (j = 0; j < cols; j++){
    for (i =0; i < nprobes; i++){
      z[j*nprobes + i] = log(data[j*rows + cur_rows[i]])/log(2.0);  
    }
  } 
  
  
  for (iter = 1; iter <= maxiter; iter++){
    get_row_median(z,rdelta,nprobes,cols);
    subtract_by_row(z,rdelta,nprobes,cols);
    rmod(r,rdelta,nprobes);
    delta = median(c,cols);
    for (j = 0; j < cols; j++){
      c[j] = c[j] - delta;
    }
    t = t + delta;
    get_col_median(z,cdelta,nprobes,cols);
    subtract_by_col(z,cdelta,nprobes,cols);
    cmod(c,cdelta,cols);
    delta = median(r,nprobes);
    for (i =0; i < nprobes; i ++){
      r[i] = r[i] - delta;
    }
    t = t+delta;
    newsum = sum_abs(z,nprobes,cols);
    if (newsum == 0.0 || fabs(1.0 - oldsum/newsum) < eps)
      break;
    oldsum = newsum;
  }
  
  for (j=0; j < cols; j++){
    results[j] =  t + c[j]; 
    resultsSE[j] = R_NaReal;
  }
  
  Free(rdelta);
  Free(cdelta);
  Free(r);
  Free(c);
  Free(z); 
}

void median_polishPLM(double *data, int rows, int cols, int *cur_rows, double *probe_param, double *chip_param, double *intercept_param, int nprobes, double *residuals){

  int i,j,iter;
  int maxiter = 10;
  double eps=0.01;
  double oldsum = 0.0,newsum = 0.0;
  double t = 0.0;
  double delta;
  double *rdelta = Calloc(nprobes,double);
  double *cdelta = Calloc(cols,double);
  
  double *r = Calloc(nprobes,double);
  double *c = Calloc(cols,double);
  double *z = Calloc(nprobes*cols,double);

  for (j = 0; j < cols; j++){
    for (i =0; i < nprobes; i++){
      z[j*nprobes + i] = log(data[j*rows + cur_rows[i]])/log(2.0);  
    }
  } 
  
  
  for (iter = 1; iter <= maxiter; iter++){
    get_row_median(z,rdelta,nprobes,cols);
    subtract_by_row(z,rdelta,nprobes,cols);
    rmod(r,rdelta,nprobes);
    delta = median(c,cols);
    for (j = 0; j < cols; j++){
      c[j] = c[j] - delta;
    }
    t = t + delta;
    get_col_median(z,cdelta,nprobes,cols);
    subtract_by_col(z,cdelta,nprobes,cols);
    cmod(c,cdelta,cols);
    delta = median(r,nprobes);
    for (i =0; i < nprobes; i ++){
      r[i] = r[i] - delta;
    }
    t = t+delta;
    newsum = sum_abs(z,nprobes,cols);
    if (newsum == 0.0 || fabs(1.0 - oldsum/newsum) < eps)
      break;
    oldsum = newsum;
  }

  for (i=0; i < nprobes; i++){
    probe_param[i] = r[i];
  }

  
  for (j=0; j < cols; j++){
    chip_param[j] =  c[j];
  }

  intercept_param[0] = t;

  for (j =0; j < cols; j++){
    for (i=0; i < nprobes; i++){
      residuals[j*nprobes +i] = z[j*nprobes +i];
    }
  }


  
  Free(rdelta);
  Free(cdelta);
  Free(r);
  Free(c);
  Free(z); 
}





void median_polish_threestep_PLM(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes, double *resultsSE, double *residuals, summary_plist *summary_param){

  int i,j,iter;
  int maxiter = 10;
  double eps=0.01;
  double oldsum = 0.0,newsum = 0.0;
  double t = 0.0;
  double delta;
  double *rdelta = Calloc(nprobes,double);
  double *cdelta = Calloc(cols,double);
  
  double *r = Calloc(nprobes,double);
  double *c = Calloc(cols,double);
  double *z = Calloc(nprobes*cols,double);

  for (j = 0; j < cols; j++){
    for (i =0; i < nprobes; i++){
      z[j*nprobes + i] = log(data[j*rows + cur_rows[i]])/log(2.0);  
    }
  } 
  
  
  for (iter = 1; iter <= maxiter; iter++){
    get_row_median(z,rdelta,nprobes,cols);
    subtract_by_row(z,rdelta,nprobes,cols);
    rmod(r,rdelta,nprobes);
    delta = median(c,cols);
    for (j = 0; j < cols; j++){
      c[j] = c[j] - delta;
    }
    t = t + delta;
    get_col_median(z,cdelta,nprobes,cols);
    subtract_by_col(z,cdelta,nprobes,cols);
    cmod(c,cdelta,cols);
    delta = median(r,nprobes);
    for (i =0; i < nprobes; i ++){
      r[i] = r[i] - delta;
    }
    t = t+delta;
    newsum = sum_abs(z,nprobes,cols);
    if (newsum == 0.0 || fabs(1.0 - oldsum/newsum) < eps)
      break;
    oldsum = newsum;
  }
  
  for (j=0; j < cols; j++){
    results[j] =  t + c[j]; 
    resultsSE[j] = R_NaReal;
  }

  for (j = 0; j < cols; j++){
    for (i =0; i < nprobes; i++){
      residuals[j*nprobes+i] = z[j*nprobes + i];
    }
  } 


  
  Free(rdelta);
  Free(cdelta);
  Free(r);
  Free(c);
  Free(z); 
}



