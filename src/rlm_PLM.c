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
 **
 *********************************************************************/

#include "preprocess.h"
#include "do_PLMrlm.h"

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

SEXP rlmPLMset(SEXP PMmat, SEXP MMmat, SEXP ProbeNamesVec,SEXP N_probes, SEXP rlm_model_type, SEXP rlm_se_type, SEXP chipcovariates,SEXP psitype, SEXP psik){
  
  int rows, cols;
  double *out_weights, *out_probeparams, *out_chipparams, *out_constparams,*input_chipcovariates; /* *chipcovmatrix;*/
  double *out_chip_SE,*out_probe_SE, *out_const_SE;
  double *PM,*MM;
  char **outnames;
  char **ProbeNames;
  int i,nprobesets,nchipparams;
  int method,se_method,psi_code;
  double psi_k;


  SEXP dim1,dim2;
  SEXP weights, probe_coef,chip_coef, const_coef,chip_SE,probe_SE,const_SE;
  /* SEXP outnamesvec; */
  SEXP dimnames,names;
  SEXP output_list;

  PROTECT(dim1 = getAttrib(PMmat,R_DimSymbol));
  rows = INTEGER(dim1)[0];
  cols = INTEGER(dim1)[1];

  //  nchipparams = cols;

  PROTECT(dim2 = getAttrib(chipcovariates,R_DimSymbol));

  nchipparams = INTEGER(dim2)[1];


  PM = NUMERIC_POINTER(AS_NUMERIC(PMmat));
  MM = NUMERIC_POINTER(AS_NUMERIC(MMmat));

  nprobesets=INTEGER(N_probes)[0];

  /* figure out what the covariate matrix should be based on the rlm_model_type information and chipcovariates */
    
  method = asInteger(rlm_model_type);
  se_method = asInteger(rlm_se_type);
  psi_code = asInteger(psitype);
  psi_k = asReal(psik);

  /* Get the names corresponding to each row */
  
  ProbeNames =malloc(rows*(sizeof(char *)));
  for (i =0; i < rows; i++)
    ProbeNames[i] = CHAR(VECTOR_ELT(ProbeNamesVec,i));
  outnames = malloc(nprobesets*sizeof(char *));
  

  /* Make space for output */

  PROTECT(weights = allocMatrix(REALSXP, rows, cols));
  out_weights = NUMERIC_POINTER(weights);
  PROTECT(probe_coef = allocMatrix(REALSXP,rows,1));
  out_probeparams = NUMERIC_POINTER(probe_coef);
  PROTECT(chip_coef = allocMatrix(REALSXP, nprobesets, nchipparams));
  out_chipparams = NUMERIC_POINTER(chip_coef);
  PROTECT(const_coef = allocMatrix(REALSXP, nprobesets, 1));
  out_constparams = NUMERIC_POINTER(const_coef);
  PROTECT(chip_SE = allocMatrix(REALSXP, nprobesets, nchipparams));
  out_chip_SE = NUMERIC_POINTER(chip_SE);

  PROTECT(probe_SE = allocMatrix(REALSXP,rows,1));
  out_probe_SE = NUMERIC_POINTER(probe_SE);
  PROTECT(const_SE = allocMatrix(REALSXP, nprobesets, 1));
  out_const_SE = NUMERIC_POINTER(const_SE);

  /* now go fit the models */

  input_chipcovariates = NUMERIC_POINTER(chipcovariates);
  /* printf("Fitting models\n");*/
  Rprintf("Fitting models\n");

  do_PLMrlm(PM, ProbeNames, &rows, &cols, nprobesets, nchipparams, method, se_method,input_chipcovariates, outnames, out_weights, out_probeparams, out_chipparams, out_constparams,out_probe_SE, out_chip_SE, out_const_SE,psi_code,psi_k);
  //  UNPROTECT(5);

  /* now lets put names on the output matrices */

  /* First the chip coef matrix */
  PROTECT(dimnames = allocVector(VECSXP,2));
  PROTECT(names = allocVector(STRSXP,nprobesets));
  for ( i =0; i < nprobesets; i++)
      SET_VECTOR_ELT(names,i,mkChar(outnames[i]));
  SET_VECTOR_ELT(dimnames,0,names);
  setAttrib(chip_coef, R_DimNamesSymbol, dimnames); 
  setAttrib(chip_SE,R_DimNamesSymbol, dimnames);

  /* Now the weights matrix */

  /* the naming here will be handled in R as of Feb 1, 2003 */
  //SET_VECTOR_ELT(dimnames,0,ProbeNamesVec);
  // setAttrib(weights, R_DimNamesSymbol, dimnames);
  
  /* Now the probe coef matrix */
  // setAttrib(probe_coef, R_DimNamesSymbol, dimnames);
  // setAttrib(probe_SE, R_DimNamesSymbol, dimnames);
  UNPROTECT(11);

  /*Now lets create the output_list */

  PROTECT(output_list = allocVector(VECSXP,7));

  SET_VECTOR_ELT(output_list,0,chip_coef);
  SET_VECTOR_ELT(output_list,1,probe_coef);
  SET_VECTOR_ELT(output_list,2,weights);
  SET_VECTOR_ELT(output_list,3,chip_SE);
  SET_VECTOR_ELT(output_list,4,probe_SE);
  SET_VECTOR_ELT(output_list,5,const_coef);
  SET_VECTOR_ELT(output_list,6,const_SE);
  UNPROTECT(1);
  
  free(outnames);
  free(ProbeNames);
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

SEXP R_rlmPLMset_c(SEXP PMmat, SEXP MMmat, SEXP ProbeNamesVec,SEXP N_probes,SEXP norm_flag, SEXP bg_flag, SEXP bg_type,SEXP norm_type, SEXP rlm_model_type, SEXP rlm_se_type, SEXP chipcovariates,  SEXP psitype, SEXP psik, SEXP background_parameters,SEXP norm_parameters){


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

  rlmPLMresults = rlmPLMset(PMcopy, MMmat, ProbeNamesVec, N_probes, rlm_model_type, rlm_se_type, chipcovariates,psitype,psik);
  
  UNPROTECT(2);

  return rlmPLMresults;
}
