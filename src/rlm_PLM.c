/*********************************************************************
 **
 ** file: rlm_PLM.c
 **
 ** Aim: fit robust linear models for the PLMset object.
 **
 ** Copyright (C) 2003 Ben Bolstad
 **
 ** created by: B. M. Bolstad <bolstad@stat.berkeley.edu>
 ** 
 ** created on: Jan 17, 2003
 **
 ** Last modified: Feb 15, 2003
 **
 ** the aim here will be to fit specified robust linear models
 ** to affy data. 
 **
 ** The code in this particular file takes the objects passed from R
 ** sends what is required to the preprocessing steps (Background and Normalization)
 ** and then takes that data, forms "C" data structures and calls routines to do the 
 ** actual robust model fitting
 **
 ** Modification history
 **
 ** Jan 17, 2003 - Initial version.
 ** Jan 18, 2003 - Better setup rlm procedure to take covariates and
 **                specify different models.
 ** Jan 19, 2003 - continued implementation
 ** Jan 20, 2003 - more implementation, clean up parameter passing methodology.
 ** Jan 24, 2003 - expand the range of models that can be fit by actually
 **                making use of the chipcovariates parameter.
 ** Jan 27, 2003 - Standard error calculation
 ** Jan 28, 2003 - Ability to select different types of standard error estimation
 ** Feb 1, 2003 - remove the row naming aspect on weights and probes, this will be handled by R 
 **               more documentation
 ** Feb 6, 2003 - Change printf("Fitting models ....") to an Rprintf
 ** Feb 15, 2003 - Add a mechanism for returning intercept parameters.
 ** Feb 17, 2003 - add in a free(outnames); free(ProbeNames) to rlmPLMset;
 ** Feb 24, 2003 - get rid of unused, but declared variables.
 ** Mar 21, 2003 - modify background for LESN methods
 ** Jun 4,  2003 - Add mechanism for different psi functions
 ** Jul 26, 2003 - add normalization options parameter
 **                add background options parameter
 ** Sep 02, 2003 - we now store residuals
 ** Sep 04, 2003 - a parameter which specifies what should be outputted is now
 **                passed. This item is an R list similar to normalization
 **                and background parameter lists.
 **                Considerable clean up of how parameters are passed to 
 **                do_PLMrlm
 ** Sep 06, 2003 - Make Storage allocation routine separate.
 ** Sep 07, 2003 - output varcov
 ** Sep 12, 2003 - remove psi, psi_k etc from arguments of functions
 **                they are now in model_params
 ** Sept 14, 2003 - can intialize M estimatation starting with a fully iterated
 **                 Huber regression
 ** Apr 5, 2004   - All malloc/free are now Calloc/Free
 ** May 3, 2004   - Fixed a subtle and small memory leak.
 ** May 27, 2004  - add a way to detect that default model is being fitted
 **
 *********************************************************************/

#include "preprocess.h"
#include "do_PLMrlm.h"
#include "common_types.h"

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>





/*********************************************************************
 **
 ** static void rlmPLM_alloc_space(PLMRoutput *Routput, PLMoutput *output,
 **                                outputsettings *store,Datagroup *data, 
 **                                PLMmodelparam *model)
 **
 ** 
 ** This function allocates all the space that is needed for storing
 ** user desired output from the PLM
 ** 
 **
 ********************************************************************/

static void rlmPLM_alloc_space(PLMRoutput *Routput, PLMoutput *output,outputsettings *store,Datagroup *data, PLMmodelparam *model){
  SEXP tmp;

  int i;
   
  Routput->nprotected = 0;

  
  output->outnames = (char **)Calloc(data->nprobesets,char *);

  if (store->weights){
    PROTECT(Routput->weights = allocMatrix(REALSXP, data->rows, data->cols));
  } else {
    PROTECT(Routput->weights = allocMatrix(REALSXP, 0, 0));
  }
  Routput->nprotected++;
  output->out_weights = NUMERIC_POINTER(Routput->weights);


  PROTECT(Routput->probe_coef = allocMatrix(REALSXP,data->rows,1));
  Routput->nprotected++;
  output->out_probeparams = NUMERIC_POINTER(Routput->probe_coef);


  PROTECT(Routput->chip_coef = allocMatrix(REALSXP, data->nprobesets, model->nchipparams));
  Routput->nprotected++;
  output->out_chipparams = NUMERIC_POINTER(Routput->chip_coef);
  
  PROTECT(Routput->const_coef = allocMatrix(REALSXP, data->nprobesets, 1));
  Routput->nprotected++;
  output->out_constparams = NUMERIC_POINTER(Routput->const_coef);

  PROTECT(Routput->chip_SE = allocMatrix(REALSXP, data->nprobesets, model->nchipparams));
  Routput->nprotected++;
  output->out_chip_SE = NUMERIC_POINTER(Routput->chip_SE);


  PROTECT(Routput->probe_SE = allocMatrix(REALSXP,data->rows,1));
  Routput->nprotected++;
  output->out_probe_SE = NUMERIC_POINTER(Routput->probe_SE);


  PROTECT(Routput->const_SE = allocMatrix(REALSXP, data->nprobesets, 1));
  Routput->nprotected++;
  output->out_const_SE = NUMERIC_POINTER(Routput->const_SE);


  if (store->residuals){
    PROTECT(Routput->residuals = allocMatrix(REALSXP, data->rows, data->cols));
  } else {
    PROTECT(Routput->residuals = allocMatrix(REALSXP, 0, 0));
  }
  Routput->nprotected++;
  output->out_resids = NUMERIC_POINTER(Routput->residuals); 
  

  if (store->residSE){
    PROTECT(Routput->residSE = allocMatrix(REALSXP,data->nprobesets, 2));
  } else {
    PROTECT(Routput->residSE = allocMatrix(REALSXP,0,0));
  }
  Routput->nprotected++;
  output->out_residSE = NUMERIC_POINTER(Routput->residSE);

  
  if (store->varcov == 0){
    PROTECT(Routput->varcov = allocVector(VECSXP,0));
    output->out_varcov= NULL;
  } else if (store->varcov == 1){
    PROTECT(Routput->varcov = allocVector(VECSXP,data->nprobesets));
    output->out_varcov = Calloc(data->nprobesets,double*);
    for (i =0; i < data->nprobesets; i++){
      PROTECT(tmp = allocMatrix(REALSXP,model->nchipparams,model->nchipparams));
      SET_VECTOR_ELT(Routput->varcov,i,tmp);
      UNPROTECT(1);
      output->out_varcov[i] = NUMERIC_POINTER(VECTOR_ELT(Routput->varcov,i));
    }
  }
  Routput->nprotected++;
  
  
  
  


}



/*********************************************************************
 **
 ** SEXP rlmPLMset(SEXP PMmat, SEXP MMmat, SEXP ProbeNamesVec,SEXP N_probes, SEXP rlm_model_type, SEXP chipcovariates)
 **
 **  
 ** SEXP PMmat - matrix with PM intensities - background corrected and normalized.
 ** SEXP MMmat - matrix with MM intensities - legacy, possibily for future use.
 ** SEXP ProbeNamesVec - vector containing names of Probeset related to each probe
 ** SEXP N_probes - number of probes
 ** SEXP rlm_model_type - codes indicating what model to be fit.
 ** SEXP chipcovariates - matrix with covariates for each chip, by convention the first column will
 **                       always be of type factor, specifying treatment group (or chip number or whatever)
 **
 ** returns an R list in order such that first element is chip parameters matrix,
 **         second element is probe parameters matrix
 **         third element is weight parameters matrix
 **         fourth, fifth are standard errors for chip and probe parameters respectively.
 **
 ** in particular  rlm_model_type will be an integer, the value
 ** will be a code describing which type of model is being fit as regards the 
 ** probe coefficients and intercept.
 **
 ** in particular for the probes (and intercept parameter) there are several options
 **
 ** 0 -  no slope parameter, probe parameters sum to zero
 ** 1 -  slope parameter, probeparameters sum to zero (sum to zero constraint)
 ** 2 -  slope parameter, probeparameters for first probe is zero (endpoint constraint)
 ** 
 ** if there is a slope parameter fitted, then there must be a constraint on the chip factor parameters.
 ** we will default to endpoint constraint, unless otherwise specified. This should be handled by the R code
 **
 ** this function will take the R objects and cast them in such away that they can be dealt with as c objects, allocate the memory
 ** needed for output from the model fitting procedure. Then a call is made to the actual routine which does all 
 ** the model fitting. After the model has been fit there is some naming applied to the output matrices 
 **
 ** as of Feb 1, 2003 probe.coef, se and weights naming will now occur in the R code, this is because we want distinct names for each probe
 ** and the probeNames accessor used in the calling of functionality here does not do this.
 **  
 **
 **
 **
 *********************************************************************/

SEXP rlmPLMset(SEXP PMmat, SEXP MMmat, SEXP ProbeNamesVec,SEXP N_probes, SEXP rlm_model_type,  SEXP chipcovariates, SEXP outputparam, SEXP modelparam){
  
  int i;

  outputsettings *store = (outputsettings *)Calloc(1,outputsettings);
  Datagroup *data = (Datagroup *)Calloc(1,Datagroup);
  PLMoutput *output = (PLMoutput *)Calloc(1,PLMoutput);
  PLMmodelparam *model = (PLMmodelparam *)Calloc(1,PLMmodelparam);
  PLMRoutput *Routput = (PLMRoutput *)Calloc(1,PLMRoutput);
  
  SEXP dim1,dim2;
  SEXP dimnames,names;
  SEXP output_list;
  SEXP param;




  /* organise data to be passed to model fitting routines */
  
  PROTECT(dim1 = getAttrib(PMmat,R_DimSymbol));
  
  data->rows = INTEGER(dim1)[0];
  data->cols = INTEGER(dim1)[1];
  
  PROTECT(dim2 = getAttrib(chipcovariates,R_DimSymbol));

  data->PM = NUMERIC_POINTER(AS_NUMERIC(PMmat));
  data->MM = NUMERIC_POINTER(AS_NUMERIC(MMmat));
  data->nprobesets = INTEGER(N_probes)[0];
  
  /* Get the names corresponding to each row */
    
  data->ProbeNames = (char **)Calloc(data->rows,char *);
  for (i =0; i < data->rows; i++){
    data->ProbeNames[i] = CHAR(VECTOR_ELT(ProbeNamesVec,i));
  }
  

  /* figure out what the covariate matrix should be based on the rlm_model_type information and chipcovariates */

  param = GetParameter(modelparam,"psi.type");
  model->psi_code = asInteger(param);
  model->method = asInteger(rlm_model_type);

  param = GetParameter(modelparam,"se.type");
  model->se_method = asInteger(param);
  
  param = GetParameter(modelparam,"psi.k");
  model->psi_k = asReal(param);
  model->input_chipcovariates = NUMERIC_POINTER(chipcovariates);
  model->nchipparams = INTEGER(dim2)[1];
  param = GetParameter(modelparam,"max.its");
  model->n_rlm_iterations = asInteger(param);
  param = GetParameter(modelparam,"init.method");
  if (strcmp(CHAR(VECTOR_ELT(param,0)),"ls") == 0){
    model->init_method = 0;
  } else if (strcmp(CHAR(VECTOR_ELT(param,0)),"median.polish") == 0){
    model->init_method = 1;
  } else if (strcmp(CHAR(VECTOR_ELT(param,0)),"Huber") == 0){
    model->init_method = 2;
  }
  param = GetParameter(modelparam,"isdefaultmodel");
  model->default_model = asInteger(param);

  param = GetParameter(modelparam,"MMorPM.covariate");
  model->mmorpm_covariate = asInteger(param);


  /* figure out what optional features we want to output */
  param = GetParameter(outputparam,"weights");
  store->weights = asInteger(param);
  param = GetParameter(outputparam,"residuals");
  store->residuals = asInteger(param);
  param = GetParameter(outputparam,"resid.SE");
  store->residSE = asInteger(param);
  param = GetParameter(outputparam,"varcov");
  
  if (strcmp(CHAR(VECTOR_ELT(param,0)),"none") == 0){
    store->varcov = 0;
  } else if (strcmp(CHAR(VECTOR_ELT(param,0)),"chiplevel") == 0){
    store->varcov =1;
  } else if (strcmp(CHAR(VECTOR_ELT(param,0)),"all") == 0){
    store->varcov =2;
  }


  /* Make space for output */
  rlmPLM_alloc_space(Routput, output, store, data, model);
  
  

  /* now go actually fit the model */

  Rprintf("Fitting models\n");

  do_PLMrlm(data, model, output,store);
  
  /* now lets put names on the output matrices */
  /* First the chip coef matrix */
  PROTECT(dimnames = allocVector(VECSXP,2));
  PROTECT(names = allocVector(STRSXP,data->nprobesets));
  for ( i =0; i < data->nprobesets; i++)
    SET_VECTOR_ELT(names,i,mkChar(output->outnames[i]));
  SET_VECTOR_ELT(dimnames,0,names);
  setAttrib(Routput->chip_coef, R_DimNamesSymbol, dimnames); 
  setAttrib(Routput->chip_SE,R_DimNamesSymbol, dimnames);

  
  /* Now the probe coef matrix */
  // setAttrib(probe_coef, R_DimNamesSymbol, dimnames);
  // setAttrib(probe_SE, R_DimNamesSymbol, dimnames);
   
  /*Now lets create the output_list */

  PROTECT(output_list = allocVector(VECSXP,10));

  SET_VECTOR_ELT(output_list,0,Routput->chip_coef);
  SET_VECTOR_ELT(output_list,1,Routput->probe_coef);
  SET_VECTOR_ELT(output_list,2,Routput->weights);
  SET_VECTOR_ELT(output_list,3,Routput->chip_SE);
  SET_VECTOR_ELT(output_list,4,Routput->probe_SE);
  SET_VECTOR_ELT(output_list,5,Routput->const_coef);
  SET_VECTOR_ELT(output_list,6,Routput->const_SE);
  SET_VECTOR_ELT(output_list,7,Routput->residuals);
  SET_VECTOR_ELT(output_list,8,Routput->residSE);
  SET_VECTOR_ELT(output_list,9,Routput->varcov);
  UNPROTECT(Routput->nprotected + 5);

  for ( i =0; i < data->nprobesets; i++)
    Free(output->outnames[i]);
  
  Free(output->outnames);
  Free(data->ProbeNames);
  Free(data);
  Free(output);
  Free(Routput);
  Free(store);
  Free(model);
  
  return output_list;
  
}



/*********************************************************************
 **
 ** SEXP R_rlmPLMset_c(SEXP PMmat, SEXP MMmat, SEXP ProbeNamesVec,SEXP N_probes,SEXP densfunc, 
 **                    SEXP rho,SEXP norm_flag, SEXP bg_flag, SEXP bg_type,SEXP norm_type, SEXP summary_type)
 **
 **
 **
 ** Interface to R should be called through .Call() interface. Fits a specified 
 ** robust linear model after doing specific background correction and normalization.
 **
 ** SEXP PMmat - matrix with PM intensities
 ** SEXP MMmat - matrix with MM intensities
 ** SEXP ProbeNamesVec - vector containing names of Probeset related to each probe
 ** SEXP N_probes - number of probes
 ** SEXP densfunc - R function for nonparametric density estimatation, used for 
 **                 RMA background procedure.
 ** SEXP rho - an R environment to work within when using densfunc
 ** SEXP norm_flag - TRUE if normalization to be carried out, false otherwise
 ** SEXP bg_flag - TRUE if background correction to be carried out, false otherwise
 ** SEXP bg_type - integer indicating which background method is to be used
 ** SEXP norm_type - integer indicating which normalization method is to be used
 ** SEXP rlm_model_type - code indicating which type of rlm to fit.
 ** SEXP chipcovariates - matrix with covariates for each chip, by convention the first column will
 **                       always be of type factor, specifying treatment group (or chip number or whatever)
 **
 ** returns an R list in order such that first element is chip parameters matrix,
 **         second element is probe parameters matrix
 **         third element is weight parameters matrix
 **         fourth, fifth are standard errors for chip and probe parameters respectively.
 **
 **
 ** SEE ALSO documentation for rlmPLMset
 **
 ********************************************************************/

SEXP R_rlmPLMset_c(SEXP PMmat, SEXP MMmat, SEXP ProbeNamesVec,SEXP N_probes,SEXP norm_flag, SEXP bg_flag, SEXP bg_type,SEXP norm_type, SEXP rlm_model_type, SEXP chipcovariates, SEXP background_parameters,SEXP norm_parameters, SEXP output_parameters, SEXP model_parameters){


  SEXP dim1,PMcopy,rlmPLMresults;
  int rows,cols;

  /*Create a copy matrix to work on. Allows us to modify data in background and normalization steps without affecting original data */
  PROTECT(dim1 = getAttrib(PMmat,R_DimSymbol));
  rows = INTEGER(dim1)[0];
  cols = INTEGER(dim1)[1];
  PROTECT(PMcopy = allocMatrix(REALSXP,rows,cols));
  copyMatrix(PMcopy,PMmat,0);

  /* If Background correction do it */
  if (INTEGER(bg_flag)[0]){
    PMcopy = pp_background(PMcopy, MMmat, ProbeNamesVec,N_probes,bg_type,background_parameters);
  }

  /* If Normalization do it */
  if (INTEGER(norm_flag)[0]){
    PMcopy = pp_normalize(PMcopy, MMmat, ProbeNamesVec,N_probes,norm_type, norm_parameters);
  }
  
  /* Now do RLM fit */

  rlmPLMresults = rlmPLMset(PMcopy, MMmat, ProbeNamesVec, N_probes, rlm_model_type, chipcovariates,output_parameters,model_parameters);
  
  UNPROTECT(2);

  return rlmPLMresults;
}
