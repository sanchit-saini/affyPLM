/*********************************************************************
 **
 ** file: do_PLMrlm.c
 **
 ** Aim: fit robust linear models for the PLMset object.
 **
 ** Copyright (C) 2003 Ben Bolstad
 **
 ** created by: B. M. Bolstad <bolstad@stat.berkeley.edu>
 ** 
 ** created on: Jan 21, 2003
 **
 ** Last modified: Feb 15, 2003
 **
 ** the aim here will be to fit specified robust linear models
 ** to affy data. this file actually contains the code that 
 ** takes preprocessed data and actually carries out the fitting of the models
 ** 
 ** in particular it breaks the data into "blocks" where a block is
 ** a probes by chips matrix of preprocessed PM intensites. to each block it fits 
 ** the specified model.
 **
 ** the probe effects part of the design matrix is handled here.
 **
 ** Note that the function calls here suffer a little bit from "fortran-itis", ie
 ** lots of parameters passed in and out, should probably to rewritten in a 
 ** cleaner form at somepoint.
 **
 ** Modification history
 **
 ** Jan 21, 2003 - Initial version.
 ** Jan 24, 2003 - Better handling of more general models.
 ** Jan 27, 2003 - Computation of standard errors
 ** Jan 28, 2003 - ability to select between se methods
 ** Jan 31, 2003 - add in a missing #include "rlm.h"
 **                note that we should be able to optimize routines, by testing 
 **                whether we really need to reallocate and initialize the X matrix 
 **                each iteration, typically almost all the probesets on a
 **                chip have same number of probesets, and X is not changed by the RLM
 **                fitting routine so we could actually keep it between blocks only changing it when
 **                the number of probes had changed.  <- TO DO WHEN WE WORK ON OPTIMIZATIONS
 **                Documentation Updates.
 **                Changed return type of rlm_PLM_block to void (At some point perhaps this should be an integer error flag)
 ** Feb 1, 2003 -  rework the design matrix code so that it is only reallocated/reinitialised if required.
 **                this results in a new function that actually allocates/assigns data values for the X matrix
 **                and changes in the parameter list passed to rlm_PLM_block
 ** Feb 11, 2003 - remove a loop that did nothing. Clean up some unused code.
 **                Optimizations in rlm_PLM_block to remove unnecessary duplication/copying of data.
 **                Clean-ups in the commenting.
 ** Feb 15, 2003 - generalize the class of models that can be fit by altering the design matrix
 **                realloc function. We also break the results copying parts of do_PLMrlm()
 **                out into a new function copy_PLM_results (this is for the model parameters,
 **                weights are still handled in the main function
 **                n,p definitions need cleaning up in this case.
 ** Feb 17, 2003 - Continue clean-up of n,p. add in the ability to fit models with no probe parameters
 ** Feb 18. 2003 - squash seg-faulting bug in rlm_designamtrix_realloc..... for method = 21
 ** Feb 23, 2003 - add implementations for method=10,11 into rlm_design_matrix_realloc(), copy_PLM_results()
 ** Feb 24, 2003 - comment out some unused variables. Trying to reduce compiler warnings.
 ** Apr 04, 2003  - make the number of rows in the PM matrix be dynamic.
 ** Jun 04, 2003  - add mechanism for more general psi functions.
 ** Jun 11, 2003  - make sure that the standard error call allows other psi functions.
 ** Jul 23, 2003 - remove one last compiler warning.
 **
 *********************************************************************/


#include "do_PLMrlm.h"
#include "rlm_se.h"
#include "rlm.h"
#include "psi_fns.h"
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>



/**********************************************************************
 **
 ** void rlm_design_matrix_realloc(double *X, int nprobes, int cols, int nchipparams, double *chipcovariates, int method)
 **
 ** double *X - on output contains new design matrix, should already allocated on entry
 ** int nprobes - number of probes 
 ** int cols - number of chips 
 ** int nchipparams - number of chip level factors/covariates
 ** double *chipcovariates - the chip covariates
 ** int method - an integer indicating what should be done to with the probe parameters
 **             we will code it so that 
 **              0 is no intercept with sum to zero constraint on probe effects.
 **              1 is intercept with sum to zero contraint on probe effects.
 **              10 is no intercept with the endpoint constraint on probe effects.
 **              11 is intercept with the endpoint constraint on probe effects.
 **              21 will fit model with no probe effects but an intercept
 **              20 or anything else will yield no probe effect, no intercept                   
 ** 
 ** this function sets the correct design matrix, it whould be called everytime
 ** the number of probes in the current probeset changes from the previous probeset.
 **
 **
 **********************************************************************/

void rlm_design_matrix_realloc(double *X, int nprobes, int cols, int p, double *chipcovariates, int method){

  int i, j, row,curcol,currow;
  int n = (nprobes*cols);
  //int p = (nchipparams + (nprobes-1));
  
  //printf("Reallocating %d %d matrix\n",n,p);

  
  /* Since we are reallocating the designmatrix we must also reintialize the entire matrix */
  for (i=0; i < n; i++){
    for (j=0; j < p; j++){
      X[j*n + i] = 0.0;
    }
  }

  
  /*
    Now make an X matrix. First columns for probe effect then chip effect columns
    Calloc puts everything to zero, so need change only non zero elements.
    
    If none of the method codes are matched then a no intercept, no probecoef model is fit.

  */

  if (method == 0){
    for (row =0; row < n; row++){
      curcol = row%nprobes;
      if (curcol == nprobes -1){
	for (j=0; j < (nprobes-1); j++){
	  X[j*n + row] = -1.0;
	}
      } else {
	X[curcol*n + row] = 1.0;
      }
    }
     /* DO chip effects by copying across from chip covariate matrix */
    for (row = 0; row < n; row++){
      currow = row/nprobes;
      for (curcol= (nprobes-1); curcol < p; curcol++){
	X[curcol*n +row] = chipcovariates[(curcol-(nprobes-1))*cols + currow];
      }
    }
  } else if (method == 1){
    for (row =0; row < n; row++){
      X[row] = 1.0;  /* Intercept term */
      curcol = row%nprobes;
      if (curcol == nprobes -1){
	for (j=1; j < nprobes; j++){
	  X[j*n + row] = -1.0;
	}
      } else {
	X[(curcol+1)*n + row] = 1.0;
      }
    } 

    /* DO chip effects by copying across from chip covariate matrix */
  
    for (row = 0; row < n; row++){
      currow = row/nprobes;
      for (curcol= (nprobes); curcol < p; curcol++){
	X[curcol*n +row] = chipcovariates[(curcol-(nprobes))*cols + currow];
      }
    }
    /* for (row =0; row < n; row++){
      for (curcol=0; curcol <p; curcol++){
	printf("%2.1f ",X[curcol*n + row]);
      }
      printf("\n");
      } */
  } else if (method == 10){
    

     for (row =0; row < n; row++){
      curcol = row%nprobes;
      if (curcol == 0){
	for (j=0; j < (nprobes-1); j++){
	  X[j*n + row] = 0.0;
	}
      } else {
	X[(curcol-1)*n + row] = 1.0;
      }
    }
     /* DO chip effects by copying across from chip covariate matrix */
    for (row = 0; row < n; row++){
      currow = row/nprobes;
      for (curcol= (nprobes-1); curcol < p; curcol++){
	X[curcol*n +row] = chipcovariates[(curcol-(nprobes-1))*cols + currow];
      }
    }
    /*for (row =0; row < n; row++){
      for (curcol=0; curcol <p; curcol++){
	printf("%2.1f ",X[curcol*n + row]);
      }
      printf("\n");
      } */

  } else if (method == 11){
    for (row =0; row < n; row++){
      X[row] = 1.0;  /* Intercept term */
      curcol = row%nprobes;
      if (curcol == 0){
	for (j=1; j < nprobes; j++){
	  X[j*n + row] = 0.0;
	}
      } else {
	X[(curcol)*n + row] = 1.0;
      }
    } 

    /* DO chip effects by copying across from chip covariate matrix */
  
    for (row = 0; row < n; row++){
      currow = row/nprobes;
      for (curcol= (nprobes); curcol < p; curcol++){
	X[curcol*n +row] = chipcovariates[(curcol-(nprobes))*cols + currow];
      }
    }
    /* for (row =0; row < n; row++){
      for (curcol=0; curcol <p; curcol++){
	printf("%2.1f ",X[curcol*n + row]);
      }
      printf("\n");
      } */
    
  } else if (method == 21){
     for (row =0; row < n; row++){
      X[row] = 1.0;  /* Intercept term */
     } 
    for (row = 0; row < n; row++){
      currow = row/nprobes;
      for (curcol= 0; curcol < (p-1); curcol++){
	X[(curcol+1)*n +row] = chipcovariates[(curcol)*cols + currow];
      }
    }

    /* for (row =0; row < n; row++){
      for (curcol=0; curcol <p; curcol++){
	printf("%2.1f ",X[curcol*n + row]);
      }
      printf("\n");
      } */ 
  } else if (method == 20){
    for (row = 0; row < n; row++){
      currow = row/nprobes;
      for (curcol= 0; curcol < p; curcol++){
	X[curcol*n +row] = chipcovariates[(curcol)*cols + currow];
      }
    }
  }
  


}



/**********************************************************************
 **
 ** void rlm_PLM_block(double *data, int rows, int cols, int *cur_rows, int nchipparams, double *chipcovariates, 
 **                      double *weights, double *params, double *se_estimates, int nprobes, int method,int se_method)
 **
 ** double *data - a matrix of preprocessed PM intensities: dimension rows * cols
 ** int rows - dimension matrix of *data (this will typically be the number of PM probes)
 ** int cols - dimension matrix of *data (this will typically be the number of chips)
 ** int *currows - a vector of row indices indicating which rows belong to the 
 **                current probeset: length nprobes
 ** int nchipparams - number of chip level parameters that will be fit by the model
 ** double *chipcovariates - a matrix of chip level covariates/factors (basically part 
 **                          of the designmatrix): dimension is cols*nchipparams
 ** double *weights - a matrix where the rlm weights based on the current fit will 
 **                   be stored: dimension nprobes*cols
 ** double *params - a vector where the parameter estimates are stored for the current 
 **                  fit are stored: length nchipparams + nprobes -1
 ** double *se_estimates - a vector where standard error estimates are to be 
 **                        stored:  length nchipparams + nprobes -1
 ** int nprobes - the number of probes in the current block
 ** int method - an integer code which will indicate what type of model the probe effects should follow
 ** int se_method - an integer code indicating with standard error method to use.
 **
 **
 ** this function fits a specified linear model to the probes from all chips for a particular probeset
 **
 **
 **
 **
 **
 *********************************************************************/

void rlm_PLM_block(double *data, int rows, int cols, int *cur_rows, int p, double *X, double *weights, double *params, double *se_estimates, int nprobes, int method,int se_method, int psitype,double psi_k){

  int i, j; /* row,curcol,currow; */

  int n = (nprobes*cols);

  double *Y = Calloc(n,double);
  double *out_resids=Calloc(n,double);
  
  double lg2 = log(2.0); /* Cache hopefully for speed :) */

  /* log2 transform and create Y vector */

  for (j = 0; j < cols; j++){
    for (i =0; i < nprobes; i++){
      Y[j*nprobes + i] = log(data[j*rows + cur_rows[i]])/lg2;
    }
  }


  rlm_fit(X,Y, n, p, params, out_resids, weights, PsiFunc(psitype),psi_k);
  rlm_compute_se(X,Y, n, p, params, out_resids, weights, se_estimates,se_method, PsiFunc(psitype),psi_k);


  Free(out_resids);
  Free(Y);
}


void copy_PLM_results(double *cur_params, double *cur_se_estimates, double *out_probeparams, double *out_chipparams, double *out_constparams,double *out_probeSE, double *out_chipSE, double *out_constSE, int rows, int cols, int nchipparams, int nprobes,int nps, int j, int i, int method){

  int k,l;

  if (method == 0){
    if (j == (rows -1)){
      out_probeparams[j] = 0.0;
      for (l=0; l < nprobes -1; l++){
	out_probeparams[j +1 - (nprobes) + l] = cur_params[l];
	out_probeparams[j]-=cur_params[l];
	out_probeSE[j +1 - (nprobes) + l] = cur_se_estimates[l];
      }
      out_probeSE[j] = 0.0/0.0;
    } else {
      out_probeparams[j -1] = 0.0;
      for (l=0; l < nprobes -1; l++){
	out_probeparams[j - (nprobes) + l] = cur_params[l];
	out_probeparams[j-1]-=cur_params[l];
	out_probeSE[j - (nprobes) + l] = cur_se_estimates[l]; 
      }
      out_probeSE[j-1] = 0.0/0.0;
    }
  
    for (k= 0; k < nchipparams ; k++){
      out_chipparams[k*nps +i] = cur_params[k+(nprobes-1)];
      out_chipSE[k*nps +i] = cur_se_estimates[k+(nprobes-1)];
    }
  } else if (method ==1){
    /* first element is  constcoef */
    
    if (j == (rows -1)){
      out_probeparams[j] = 0.0;
      for (l=0; l < nprobes -1; l++){
	out_probeparams[j +1 - (nprobes) + l] = cur_params[l+1];
	out_probeparams[j]-=cur_params[l+1];
	out_probeSE[j +1 - (nprobes) + l] = cur_se_estimates[l+1];
      }
      out_probeSE[j] = 0.0/0.0;
    } else {
      out_probeparams[j -1] = 0.0;
      for (l=0; l < nprobes -1; l++){
	out_probeparams[j - (nprobes) + l] = cur_params[l+1];
	out_probeparams[j-1]-=cur_params[l+1];
	out_probeSE[j - (nprobes) + l] = cur_se_estimates[l+1]; 
      }
      out_probeSE[j-1] = 0.0/0.0;
    }
    for (k= 0; k < nchipparams ; k++){
      out_chipparams[k*nps +i] = cur_params[k+(nprobes)];
      out_chipSE[k*nps +i] = cur_se_estimates[k+(nprobes)];
    }
    out_constparams[i] = cur_params[0];
    out_constSE[i] = cur_se_estimates[0];
  } else if (method == 20){
    for (k= 0; k < nchipparams ; k++){
      out_chipparams[k*nps +i] = cur_params[k];
      out_chipSE[k*nps +i] = cur_se_estimates[k];
    }
  } else if (method == 21){
    out_constparams[i] = cur_params[0];
    out_constSE[i] = cur_se_estimates[0];
    for (k= 0; k < nchipparams ; k++){
      out_chipparams[k*nps +i] = cur_params[k+1];
      out_chipSE[k*nps +i] = cur_se_estimates[k+1];
    }
  } else if (method == 10){
    
    if (j == (rows -1)){
      out_probeparams[j+1-nprobes] = 0.0;
      for (l=1; l < nprobes; l++){
	out_probeparams[j +1 - (nprobes) + l] = cur_params[l-1];
	//out_probeparams[j]-=cur_params[l];
	out_probeSE[j +1 - (nprobes) + l] = cur_se_estimates[l-1];
      }
      out_probeSE[j+1-nprobes] = 0.0;
    } else {
      out_probeparams[j-nprobes] = 0.0;
      for (l=1; l < nprobes; l++){
	out_probeparams[j - (nprobes) + l] = cur_params[l-1];
	//out_probeparams[j-1]-=cur_params[l];
	out_probeSE[j - (nprobes) + l] = cur_se_estimates[l-1]; 
      }
      out_probeSE[j-nprobes] = 0.0;
    }
    
    for (k= 0; k < nchipparams ; k++){
      out_chipparams[k*nps +i] = cur_params[k+(nprobes-1)];
      out_chipSE[k*nps +i] = cur_se_estimates[k+(nprobes-1)];
    }
    
    
  } else if (method == 11){
    /* first element is  constcoef */
    
    if (j == (rows -1)){
      out_probeparams[j+1-nprobes] = 0.0;
      for (l=1; l < nprobes; l++){
	out_probeparams[j +1 - (nprobes) + l] = cur_params[l];
	//out_probeparams[j]-=cur_params[l+1];
	out_probeSE[j +1 - (nprobes) + l] = cur_se_estimates[l];
      }
      out_probeSE[j+1-nprobes] = 0.0;
    } else {
      out_probeparams[j -nprobes] = 0.0;
      for (l=1; l < nprobes; l++){
	out_probeparams[j - (nprobes) + l] = cur_params[l];
	//out_probeparams[j-1]-=cur_params[l+1];
	out_probeSE[j - (nprobes) + l] = cur_se_estimates[l]; 
      }
      out_probeSE[j-nprobes] = 0.0;
    }
    for (k= 0; k < nchipparams ; k++){
      out_chipparams[k*nps +i] = cur_params[k+(nprobes)];
      out_chipSE[k*nps +i] = cur_se_estimates[k+(nprobes)];
    }
    out_constparams[i] = cur_params[0];
    out_constSE[i] = cur_se_estimates[0];

  }
  
}











/********************************************************************************************
 **
 ** void do_PLMrlm(double *PM, char **ProbeNames, int *rows, int *cols, int nps, int method, 
 **                double *chipcovariates, char **outNames, double *out_weights, double *out_probeparams, 
 **                double *out_chipparams, double *out_constparams)
 **
 ** double *PM - a matrix of Preprocessed (N and B) PM probe intensities: dimension *rows by *cols
 ** char **ProbeNames - a vector of character strings, each element is the name of
 **                     the corresponding probeset that this element should belong to
 **                     ie if Probenames[i] = "ABCDEFG" then the PM values in row i of 
 **                     the PM matrix belong to probeset "ABCDEFG"
 ** int *rows - dimension of *PM
 ** int *cols - dimension of *PM
 ** int nps - number of probesets on chip
 ** int method - integer value indicating what method should be used to form the design matrix in particular in 
 **              relation to the probe effects
 ** double *chipcovariates - A matrix of chip level  factor/covariates, basically one part of the designmatrix
 ** char **outNames - a place to put the names of each probeset processed
 ** double *out_weights - a matrix into which one dumps fitted model rlm weights
 ** double *out_probeparams - a matrix to store probe effect coefficents
 ** double *out_chipparams - a matrix to store chip parameter estimates
 ** double *out_constparams - a matrix to store constant parameters (if fitted)
 ** double *out_probeSE, double *out_chipSE, double *out_constSE
 **
 **
 **
 **
 *******************************************************************************************/


void do_PLMrlm(double *PM, char **ProbeNames, int *rows, int *cols, int nps, int nchipparams, int method,int se_method, double *chipcovariates, char **outNames, double *out_weights, double *out_probeparams, double *out_chipparams, double *out_constparams,double *out_probeSE, double *out_chipSE, double *out_constSE, int psitype, double psi_k){
 int j = 0;
  int i = 0;
  int k = 0,l=0;
  int n, p=0;
  int size;
  char *first;
  int first_ind;
  int max_nrows = 1000;

  /* buffers of size 200 should be enough. */

  int *cur_rows=Calloc(max_nrows,int);
  int nprobes=0;
  int old_nprobes =0;

  double *cur_weights = Calloc(*cols,double);
  double *cur_params = Calloc(*cols+100,double);
  double *cur_se_estimates = Calloc(*cols+100,double);
  /* double *OLDPM = NULL; */

  double *X = Calloc(10,double);

  first = ProbeNames[0];
  first_ind = 0;
  i =0;
  nprobes = 1;
  for (j = 1; j < *rows; j++){
    if ((strcmp(first,ProbeNames[j]) != 0) | (j == (*rows -1))){
      if (j == (*rows -1)){
        nprobes++;
        for (k = 0; k < nprobes; k++){
	  if (k >= max_nrows){
	    max_nrows = 2*max_nrows;
	    cur_rows = Realloc(cur_rows, max_nrows, int);
	  }
          cur_rows[k] = (j+1 - nprobes)+k;
        }
      } else {
        for (k = 0; k < nprobes; k++){
	  if (k >= max_nrows){
	    max_nrows = 2*max_nrows;
	    cur_rows = Realloc(cur_rows, max_nrows, int);
	  }
          cur_rows[k] = (j - nprobes)+k;
	}
      }

      /* Check last number of probes and only Realloc when needed */
      if (old_nprobes != nprobes){
	n = nprobes*(*cols);
	if (method % 10 == 1){
	  if (method == 21){
	    p = nchipparams+1;
	  } else {
	    p = nprobes + nchipparams;
	  }
	} else {
	  if (method == 20){
	    p = nchipparams;
	  } else {
	    p = nprobes + nchipparams-1;
	  }
	} 
	cur_weights = Realloc(cur_weights,n,double);
	cur_params  = Realloc(cur_params,p,double);
	cur_se_estimates  = Realloc(cur_se_estimates,p,double);
	X = Realloc(X,n*p,double);
	rlm_design_matrix_realloc(X, nprobes, *cols, p, chipcovariates, method);
	old_nprobes = nprobes;
      }


      rlm_PLM_block(PM, *rows, *cols, cur_rows, p, X , cur_weights, cur_params, cur_se_estimates, nprobes, method,se_method,psitype,psi_k);
      
      copy_PLM_results(cur_params, cur_se_estimates, out_probeparams, out_chipparams, out_constparams, out_probeSE, out_chipSE, out_constSE, *rows, *cols, nchipparams, nprobes, nps,j,i,method);
      
      /*      if (j == (*rows -1)){
	out_probeparams[j] = 0.0;
	for (l=0; l < nprobes -1; l++){
	  out_probeparams[j +1 - (nprobes) + l] = cur_params[l];
	  out_probeparams[j]-=cur_params[l];
	  out_probeSE[j +1 - (nprobes) + l] = cur_se_estimates[l];
	}
	out_probeSE[j] = 0.0/0.0;
      } else {
	out_probeparams[j -1] = 0.0;
	for (l=0; l < nprobes -1; l++){
	  out_probeparams[j - (nprobes) + l] = cur_params[l];
	  out_probeparams[j-1]-=cur_params[l];
	  out_probeSE[j - (nprobes) + l] = cur_se_estimates[l]; 
	}
	out_probeSE[j-1] = 0.0/0.0;
      }
	
      for (k= 0; k < nchipparams ; k++){
	  out_chipparams[k*nps +i] = cur_params[k+(nprobes-1)];
	  out_chipSE[k*nps +i] = cur_se_estimates[k+(nprobes-1)];
	  } */
      
      if (j == (*rows -1)){
	for(k=0; k < *cols; k++){
	  for (l=0; l < nprobes; l++){
	    out_weights[k*(*rows) + (j+1 - (nprobes) + l)] = cur_weights[k*(nprobes) + l];
	  }
	  //printf("\n");
	}
      } else {
	for(k=0; k < *cols; k++){
	  for (l=0; l < nprobes; l++){
	    out_weights[k*(*rows) + (j - (nprobes) + l)] = cur_weights[k*(nprobes) + l];
	    //printf("%d ",(j - (nprobes) + l));
	    //  printf("% f",cur_weights[k*(nprobes) + l]);
	  }
	  //printf("\n");
	}
      }
      
      size = strlen(first);
      outNames[i] = Calloc(size+1,char);
      strcpy(outNames[i],first);
      i++;
      first = ProbeNames[j];
      first_ind = j;
      nprobes = 0;
    }
    nprobes++;
  }

  Free(X);
  Free(cur_se_estimates);
  Free(cur_params);
  Free(cur_weights);
  Free(cur_rows);
}
