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
 ** cleaner form at somepoint. (THIS HAS BEEN DONE)
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
 ** Sep 02, 2003 - residuals are now stored
 ** Sep 05, 2003 - clean up in how parameters are passed to  do_PLMrlm
 ** Sep 06, 2003 - introduced the struct modelfit. It is used to group together
 **                items related to the current model (for each probeset) being
 **                fitted. More clean up in how parameters are passed between
 **                functions. residSE now outputted.
 ** Sep 07, 2003 - chip level part of variance covariance matrix returned.
 ** Sep 08, 2003 - more work on returning variance covariance
 ** Sep 13, 2003 - Modify rlm_PLM_block so it handles number of iterations
 **                and initialization method
 ** Oct 12, 2003 - fixed declaration order error                     
 **
 *********************************************************************/


#include "do_PLMrlm.h"
#include "rlm_se.h"
#include "rlm.h"
#include "psi_fns.h"
#include "common_types.h"
#include "medianpolish.h"

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/**********************************************************************
 **
 ** the modelfit struct is used for storing information about the 
 ** current model (the one being fitted individually to each probeset.
 **
 **********************************************************************/

typedef struct{
  double *cur_params;            /* storage */
  double *cur_se_estimates;
  double *cur_weights;
  double *cur_resids;
  double *cur_varcov;
  double *cur_residSE;
  int *cur_rows;  /* indices in the data matrix to use for current model */
  double *X;      /* design matrix */
  int n;          /* number of observations */
  int p;          /* number of parameters */
  int nprobes;    /* number of probes in current probeset */
  
} modelfit;








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
 **  rlm_PLM_block(const Datagroup *data, const PLMmodelparam *model, modelfit *current)
 **
 **
 **********************************************************************
 **
 ** HERE AND BELOW IS OUTDATED. LEFT HERE FOR HISTORICAL INTEREST
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

void rlm_PLM_block(const Datagroup *data, const PLMmodelparam *model, modelfit *current){
  
  int i, j;
  //  int n = (current->nprobes*data->cols);

  double *Y = Calloc(current->n,double);
  double lg2 = log(2.0); /* Cache hopefully for speed :) */

  double *probeparam;
  double *chipparam;
  double constparam;


  /* log2 transform and create Y vector */

  for (j = 0; j < data->cols; j++){
    for (i =0; i < current->nprobes; i++){
      Y[j*current->nprobes + i] = log(data->PM[j*data->rows + current->cur_rows[i]])/lg2;
    }
  }

  if (model->init_method == 1){
    /* median polish */
    probeparam = Calloc(current->nprobes,double);
    chipparam = Calloc(data->cols,double);
    median_polishPLM(data->PM,data->rows, data->cols, current->cur_rows, probeparam, chipparam, &constparam, current->nprobes, current->cur_resids); 
    //   if (model->method == 0){
    //  for (i =0; i < (current->nprobes-1); i++){
    //	current->cur_params[i] = probeparam[i];
    //  }
    //  for (i = 0; i < data->cols; i++){
    //	current->cur_params[i+(current->nprobes-1)] = constparam + chipparam[i];
    //  }
    // }
    Free(probeparam);
    Free(chipparam);
  } else if (model->init_method ==2){
    /* fully iterated huber regression */
    
    rlm_fit(current->X,Y, current->n, current->p, current->cur_params, current->cur_resids, current->cur_weights, PsiFunc(0),1.345,20,0);

  }
  /* Least Squares */
  rlm_fit(current->X,Y, current->n, current->p, current->cur_params, current->cur_resids, current->cur_weights, PsiFunc(model->psi_code),model->psi_k,model->n_rlm_iterations,model->init_method);
  rlm_compute_se(current->X,Y, current->n, current->p, current->cur_params, current->cur_resids, current->cur_weights, current->cur_se_estimates,current->cur_varcov, current->cur_residSE, model->se_method, PsiFunc(model->psi_code),model->psi_k);
  
  Free(Y);
}

/*********************************************************************
 **
 ** void copy_PLM_results(modelfit *current, PLMoutput *output, Datagroup *data,const PLMmodelparam *model, int j, int i)
 **
 ** This function should copy the results of the current fit into 
 ** the appropiate output areas.
 **
 ********************************************************************/

void copy_PLM_results(modelfit *current, PLMoutput *output, Datagroup *data,const PLMmodelparam *model, const outputsettings *store, int j, int i){

  int k,l;
  int offset;

  /* depending on what model was fit copy the parameter estimates and standard errors to the appropriate places */
 
  
  if (model->method == 0){
    if (j == (data->rows -1)){
      output->out_probeparams[j] = 0.0;
      for (l=0; l < current->nprobes -1; l++){
	output->out_probeparams[j +1 - (current->nprobes) + l] = current->cur_params[l];
	output->out_probeparams[j]-=current->cur_params[l];
	output->out_probe_SE[j +1 - (current->nprobes) + l] = current->cur_se_estimates[l];
      }
      output->out_probe_SE[j] = 0.0/0.0;
    } else {
      output->out_probeparams[j -1] = 0.0;
      for (l=0; l < current->nprobes -1; l++){
	output->out_probeparams[j - (current->nprobes) + l] = current->cur_params[l];
	output->out_probeparams[j-1]-=current->cur_params[l];
	output->out_probe_SE[j - (current->nprobes) + l] = current->cur_se_estimates[l]; 
      }
      output->out_probe_SE[j-1] = 0.0/0.0;
    }
  
    for (k= 0; k < model->nchipparams ; k++){
      output->out_chipparams[k*data->nprobesets +i] = current->cur_params[k+(current->nprobes-1)];
      output->out_chip_SE[k*data->nprobesets +i] = current->cur_se_estimates[k+(current->nprobes-1)];
    }
  } else if (model->method ==1){
    /* first element is  constcoef */
    
    if (j == (data->rows -1)){
      output->out_probeparams[j] = 0.0;
      for (l=0; l < current->nprobes -1; l++){
	output->out_probeparams[j +1 - (current->nprobes) + l] = current->cur_params[l+1];
	output->out_probeparams[j]-=current->cur_params[l+1];
	output->out_probe_SE[j +1 - (current->nprobes) + l] = current->cur_se_estimates[l+1];
      }
      output->out_probe_SE[j] = 0.0/0.0;
    } else {
      output->out_probeparams[j -1] = 0.0;
      for (l=0; l < current->nprobes -1; l++){
	output->out_probeparams[j - (current->nprobes) + l] = current->cur_params[l+1];
	output->out_probeparams[j-1]-=current->cur_params[l+1];
	output->out_probe_SE[j - (current->nprobes) + l] = current->cur_se_estimates[l+1]; 
      }
      output->out_probe_SE[j-1] = 0.0/0.0;
    }
    for (k= 0; k < model->nchipparams ; k++){
      output->out_chipparams[k*data->nprobesets +i] = current->cur_params[k+(current->nprobes)];
      output->out_chip_SE[k*data->nprobesets +i] = current->cur_se_estimates[k+(current->nprobes)];
    }
    output->out_constparams[i] = current->cur_params[0];
    output->out_const_SE[i] = current->cur_se_estimates[0];
  } else if (model->method == 20){
    for (k= 0; k < model->nchipparams ; k++){
      output->out_chipparams[k*data->nprobesets +i] = current->cur_params[k];
      output->out_chip_SE[k*data->nprobesets +i] = current->cur_se_estimates[k];
    }
  } else if (model->method == 21){
    output->out_constparams[i] = current->cur_params[0];
    output->out_const_SE[i] = current->cur_se_estimates[0];
    for (k= 0; k < model->nchipparams ; k++){
      output->out_chipparams[k*data->nprobesets +i] = current->cur_params[k+1];
      output->out_chip_SE[k*data->nprobesets +i] = current->cur_se_estimates[k+1];
    }
  } else if (model->method == 10){
    
    if (j == (data->rows -1)){
      output->out_probeparams[j+1-current->nprobes] = 0.0;
      for (l=1; l < current->nprobes; l++){
	output->out_probeparams[j +1 - (current->nprobes) + l] = current->cur_params[l-1];
	//out_probeparams[j]-=current->cur_params[l];
	output->out_probe_SE[j +1 - (current->nprobes) + l] = current->cur_se_estimates[l-1];
      }
      output->out_probe_SE[j+1-current->nprobes] = 0.0;
    } else {
      output->out_probeparams[j-current->nprobes] = 0.0;
      for (l=1; l < current->nprobes; l++){
	output->out_probeparams[j - (current->nprobes) + l] = current->cur_params[l-1];
	//out_probeparams[j-1]-=current->cur_params[l];
	output->out_probe_SE[j - (current->nprobes) + l] = current->cur_se_estimates[l-1]; 
      }
      output->out_probe_SE[j-current->nprobes] = 0.0;
    }
    
    for (k= 0; k < model->nchipparams ; k++){
      output->out_chipparams[k*data->nprobesets +i] = current->cur_params[k+(current->nprobes-1)];
      output->out_chip_SE[k*data->nprobesets +i] = current->cur_se_estimates[k+(current->nprobes-1)];
    }
    
    
  } else if (model->method == 11){
    /* first element is  constcoef */
    
    if (j == (data->rows -1)){
      output->out_probeparams[j+1-current->nprobes] = 0.0;
      for (l=1; l < current->nprobes; l++){
	output->out_probeparams[j +1 - (current->nprobes) + l] = current->cur_params[l];
	//out_probeparams[j]-=current->cur_params[l+1];
	output->out_probe_SE[j +1 - (current->nprobes) + l] = current->cur_se_estimates[l];
      }
      output->out_probe_SE[j+1-current->nprobes] = 0.0;
    } else {
      output->out_probeparams[j -current->nprobes] = 0.0;
      for (l=1; l < current->nprobes; l++){
	output->out_probeparams[j - (current->nprobes) + l] = current->cur_params[l];
	//out_probeparams[j-1]-=current->cur_params[l+1];
	output->out_probe_SE[j - (current->nprobes) + l] = current->cur_se_estimates[l]; 
      }
      output->out_probe_SE[j-current->nprobes] = 0.0;
    }
    for (k= 0; k < model->nchipparams ; k++){
      output->out_chipparams[k*data->nprobesets +i] = current->cur_params[k+(current->nprobes)];
      output->out_chip_SE[k*data->nprobesets +i] = current->cur_se_estimates[k+(current->nprobes)];
    }
    output->out_constparams[i] = current->cur_params[0];
    output->out_const_SE[i] = current->cur_se_estimates[0];

  }



  /* copy out the variance/covariance matrix (or parts thereof) if required */


  if (store->varcov){
    offset = current->p-model->nchipparams;
    if (store->varcov == 1){
      /**              0 is no intercept with sum to zero constraint on probe effects.
       **              1 is intercept with sum to zero contraint on probe effects.
       **              10 is no intercept with the endpoint constraint on probe effects.
       **              11 is intercept with the endpoint constraint on probe effects.
       **              21 will fit model with no probe effects but an intercept
       **              20 or anything else will yield no probe effect, no intercept  
       **/
      for (k = 0; k < model->nchipparams; k++){
	for (l = 0; l <= k ; l++){
	  output->out_varcov[i][k*model->nchipparams + l] = current->cur_varcov[(k+offset)*current->p + (l+offset)];
	  output->out_varcov[i][l*model->nchipparams + k] = output->out_varcov[i][k*model->nchipparams + l];
	}
      }
    } else {
      error("varcov option all not currently supported");
      
    }
  }
  
  /* copy the weights and residuals into output */
  /* note that we use the values in "store"
     to determine whether to save what has been returned
     for everything that follows                       */


  
  
  if (store->weights){
    if (j == (data->rows -1)){
      for(k=0; k < data->cols; k++){
	for (l=0; l < current->nprobes; l++){
	  output->out_weights[k*(data->rows) + (j+1 - (current->nprobes) + l)] = current->cur_weights[k*(current->nprobes) + l];
	}
	//printf("\n");
      }
    } else {
      for(k=0; k < data->cols; k++){
	for (l=0; l < current->nprobes; l++){
	  output->out_weights[k*(data->rows) + (j - (current->nprobes) + l)] = current->cur_weights[k*(current->nprobes) + l];
	  //printf("%d ",(j - (current->nprobes) + l));
	  //  printf("% f",current->cur_weights[k*(current->nprobes) + l]);
	}
	//printf("\n");
      }
    }
    
  }

  if (store->residuals){
    if (j == (data->rows -1)){
      for(k=0; k < data->cols; k++){
	for (l=0; l < current->nprobes; l++){
	  output->out_resids[k*(data->rows) + (j+1 - (current->nprobes) + l)] = current->cur_resids[k*(current->nprobes) + l];
	}
      }
    } else {
      for(k=0; k < data->cols; k++){
	for (l=0; l < current->nprobes; l++){
	  output->out_resids[k*(data->rows) + (j - (current->nprobes) + l)] = current->cur_resids[k*(current->nprobes) + l];
	}
      }
    }
    
  }

  if (store->residSE){
    output->out_residSE[i] = current->cur_residSE[0];
    output->out_residSE[data->nprobesets+i] = current->n - current->p;
  }
}










/********************************************************************************************
 ** 
 ** void do_PLMrlm(Datagroup *data,  PLMmodelparam * model, PLMoutput *output, 
 **                outputsettings *store) 
 **
 ** Datagroup *data - the data to which we will be fitting models
 ** PLMmodelparam *model - information about the model to be fitted
 ** PLMoutput *output - places to store various model outputs
 ** outputsettings *store - certain items are optional and won't
 **                         be stored unless user desired.
 **
 **
 ********************************************************************************************
 **  ANYTHING BELOW THIS LINE IS OUTDATED AND SHOULD BE IGNORED. IT IS LEFT HERE FOR
 **  HISTORICAL INTEREST ONLY.
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

void do_PLMrlm(Datagroup *data,  PLMmodelparam *model, PLMoutput *output, outputsettings *store){
  int i = 0,j=0,k=0;
  int size;
  char *first;
  int first_ind;
  int max_nrows = 1000;
  int old_nprobes =0;
 
  /* buffers of size 200 should be enough. */

  modelfit *current = malloc(sizeof(modelfit));

  current->cur_rows=Calloc(max_nrows,int);
  current->cur_weights = Calloc(data->cols,double);
  current->cur_params = Calloc(data->cols+100,double);
  current->cur_se_estimates = Calloc(data->cols+100,double);
  current->cur_resids = Calloc(data->cols,double);
  current->X = Calloc(10,double);
  current->p = 0;
  current->nprobes = 0;
  current->n = 0;
  current->cur_residSE = Calloc(2,double);
  current->cur_varcov = Calloc(4,double);

  first = data->ProbeNames[0];
  first_ind = 0;
  i =0;
  current->nprobes = 1;
  for (j = 1; j < data->rows; j++){
    if ((strcmp(first,data->ProbeNames[j]) != 0) | (j == (data->rows -1))){
      if (j == (data->rows -1)){
        current->nprobes++;
        for (k = 0; k < current->nprobes; k++){
	  if (k >= max_nrows){
	    max_nrows = 2*max_nrows;
	    current->cur_rows = Realloc(current->cur_rows, max_nrows, int);
	  }
          current->cur_rows[k] = (j+1 - current->nprobes)+k;
        }
      } else {
        for (k = 0; k < current->nprobes; k++){
	  if (k >= max_nrows){
	    max_nrows = 2*max_nrows;
	    current->cur_rows = Realloc(current->cur_rows, max_nrows, int);
	  }
          current->cur_rows[k] = (j - current->nprobes)+k;
	}
      }

      /* Check last number of probes and only Realloc when needed */
      if (old_nprobes != current->nprobes){
	current->n = current->nprobes*(data->cols);
	if (model->method % 10 == 1){
	  if (model->method == 21){
	    current->p = model->nchipparams+1;
	  } else {
	    current->p = current->nprobes + model->nchipparams;
	  }
	} else {
	  if (model->method == 20){
	    current->p = model->nchipparams;
	  } else {
	    current->p = current->nprobes + model->nchipparams-1;
	  }
	} 
	current->cur_weights = Realloc(current->cur_weights,current->n,double);
	current->cur_resids = Realloc(current->cur_resids,current->n,double);
	current->cur_params  = Realloc(current->cur_params,current->p,double);
	current->cur_se_estimates  = Realloc(current->cur_se_estimates,current->p,double);
	current->cur_varcov = Realloc(current->cur_varcov,current->p*current->p, double);
	current->X = Realloc(current->X,current->n*current->p,double);
	rlm_design_matrix_realloc(current->X, current->nprobes, data->cols, current->p, model->input_chipcovariates, model->method);
	old_nprobes = current->nprobes;
      }


      rlm_PLM_block(data, model, current);      
      copy_PLM_results(current, output, data, model, store, j,i);
      
    
      size = strlen(first);
      output->outnames[i] = Calloc(size+1,char);
      strcpy(output->outnames[i],first);  
      i++;
      first = data->ProbeNames[j];
      first_ind = j;
      current->nprobes = 0;
    }
    current->nprobes++;
  }


  Free(current->X);
  Free(current->cur_varcov);
  Free(current->cur_resids);
  Free(current->cur_se_estimates);
  Free(current->cur_params);
  Free(current->cur_weights);
  Free(current->cur_rows);
  free(current);
}
