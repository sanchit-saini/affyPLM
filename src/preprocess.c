/*********************************************************************
 **
 ** file: preprocess.c
 **
 ** Aim: Have implementations of the background and normalization 
 **      steps 
 **
 ** Copyright (C) 2003 Ben Bolstad
 **
 ** created by: B. M. Bolstad <bolstad@stat.berkeley.edu>
 ** 
 ** created on: Jan 17, 2003
 **
 ** Last modified: Feb 24, 2003
 **
 ** This file contains the normalization and background preprocessing
 ** steps.
 **
 ** Modification history
 **
 ** Jan 17, 2003 - Initial version. Moved background and normalization
 **                steps from threestep.c to this file. renamed from
 **                threestep_* to pp_*
 **  
 ** Feb 6, 2003 - change printf to Rprintf
 ** Feb 12, 2003 - add a missing #include "preprocess.h"
 **                clean up documentation
 **                Put in some code, that should not yet be executed to allow
 **                a call to the MAS 5.0 style background.
 ** Feb 24, 2003 - cleanup declaration of variable that is not used
 ** Mar 21, 2003 - Add ability to call LESN background methods, this required
 **                adding additional parameters to pp_background function for 
 **                LESN parameters.
 ** Mar 23, 2003 - comment out printf in MAS background, this will be handled by R code.
 ** Jun 23, 2003 - Add quantiles_probeset in normalization options.
 ** Jul 24, 2003 - scaling normalization added as normalization option
 ** Jul 26, 2003 - more general framework for providing parameters to 
 **                normalization. The function GetParameter was introduced.
 **                A similar change was made to the background code
 **                so that SEXP densfunc, SEXP rho, SEXP LESN_param
 **                were all subsumed into bg_parameters
 ** Apr 5, 2004 - all malloc/free are now Calloc/Free
 **
 **
 *********************************************************************/

#include "rma_common.h" 
#include "rma_background2.h" 
#include "qnorm.h" 
#include "idealmismatch.h"
#include "LESN.h"
#include "preprocess.h"
#include "qnorm_probeset.h"
#include "scaling.h"

#include <R.h> 
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/********************************************************************************************
 **
 ** SEXP GetParameter(SEXP alist, char *param_name)
 **
 **
 **
 ** given a list find the parameter with the given name or halt because parameter not
 ** found.
 **
 ********************************************************************************************/

SEXP GetParameter(SEXP alist, char *param_name){

  int length,i;
  SEXP Names;

  if (isVectorList(alist) ==FALSE){
    error("Parameter list was not list.");
  }
  
  length = length(alist);
  Names = GET_NAMES(alist);

  if (length(Names) != length){
    error("Need names for all items in parameter list.");
  }

  for( i=0; i < length; i++){
    if (strcmp(CHAR(VECTOR_ELT(Names,i)),param_name) == 0){
      break;
    }
  }

  if (i == length){
    error("Did not find %s in parameter list.", param_name);
  }
  
  return VECTOR_ELT(alist,i);

  
}



/********************************************************************************************
 **
 ** void pp_background(SEXP PMmat, SEXP MMmat, SEXP ProbeNamesVec,SEXP N_probes,SEXP norm_type)
 **
 ** SEXP PMmat - matrix of Perfect-match values
 ** SEXP MMmat - matrix of Mismatch values
 ** SEXP ProbeNamesVec - vector containing names of probeset for each probe
 ** SEXP N_probes - number of PM/MM probes on an array
 ** SEXP bg_type  - an integer indicating the background method to be used.
 ** SEXP densfunc - an R function for computing non parameteric density in RMA background
 ** SEXP rho - an R environment used in computation of RMA background
 ** SEXP LESN_param - a vector of two elements containing parameters for LESN method
 **                   first parameter is p0, the second is theta  - ADDED Mar 21, 2003
 **
 ** Carry out the background/PM adjustment preprocessing step.
 **
 ** At some point this will be generalized so that it is even easier to add new
 ** and/or different background methods.
 **
 ** bg_type definitions
 **
 ** 1 - RMA version 1 background (this is what is used in the normalization paper)
 ** 2 - RMA version 2 background (used in the NAR paper)
 ** 3 - use Ideal Mismatch Idea of MAS 5.0 to subtract IM from PM
 ** 4 - use an MAS 5.0 style within grid background correction method
 ** 5 - use both MAS 5 style followed by Ideal Mismatch correction
 ** 6 - use LESN proposal 2, half gaussian correction
 ** 7 - use LESN proposal 1, exponential correction
 ** 8 - use LESN proposal 0, just shift intensities lower so that lowest value is small
 **
 ********************************************************************************************/

SEXP pp_background(SEXP PMmat, SEXP MMmat, SEXP ProbeNamesVec,SEXP N_probes,SEXP bg_type,SEXP background_param){
  int i;
  double *PM,*MM;
  
  double theta, baseline;
  int rows, cols;
  char **ProbeNames;
  
  SEXP dim1;
  SEXP densfunc;
  SEXP rho;
  SEXP LESN_param;




  /*printf("Background correcting\n");*/
  
  if (asInteger(bg_type) == 1){
    /* Version 1 RMA Background */
    densfunc = GetParameter(background_param,"densfun");
    rho = GetParameter(background_param,"rho");
    PMmat = bg_correct_c(PMmat, MMmat, densfunc,rho,bg_type);
  } else if (asInteger(bg_type) == 2){
    /* Version 2 RMA BAckground */
    densfunc = GetParameter(background_param,"densfun");
    rho = GetParameter(background_param,"rho");
    PMmat = bg_correct_c(PMmat, MMmat, densfunc,rho,bg_type);
  } else if (asInteger(bg_type) == 3){
    /* Correction using Ideal Mismatch */ 
    /* printf("Background correcting\n");*/
    Rprintf("Background correcting\n");

    PROTECT(dim1 = getAttrib(PMmat,R_DimSymbol)); 
    rows = INTEGER(dim1)[0];
    cols = INTEGER(dim1)[1]; 
    
    PM = NUMERIC_POINTER(AS_NUMERIC(PMmat));
    MM = NUMERIC_POINTER(AS_NUMERIC(MMmat));

    ProbeNames =(char **)Calloc(rows,char *);
    for (i =0; i < rows; i++)
      ProbeNames[i] = CHAR(VECTOR_ELT(ProbeNamesVec,i));

    IdealMM_correct(PM,MM, &rows, &cols,ProbeNames);
    Free(ProbeNames);
    UNPROTECT(1);

  } else if (asInteger(bg_type) == 4){
    /* Mas 5 location dependent background */
    /* MAS5 will be handled within R code */

    /* Rprintf("Background correcting - MAS Method currently not implemented\n"); */
    /* affy_background_adjust_R(probeintensity,x, y, *nprobes, *nchips, *rows, *cols, *grid_dim); */
  } else if (asInteger(bg_type) == 5){
    /* Mas 5 location dependent background 
       followed by IMM correction */
    /* MAS5 will be handled within R code */
    /* Rprintf("Background correcting - MAS Method currently not implemented, using only ideal mismatch\n"); */
    /* affy_background_adjust_R(probeintensity,x, y, *nprobes, *nchips, *rows, *cols, *grid_dim); */
    PROTECT(dim1 = getAttrib(PMmat,R_DimSymbol)); 
    rows = INTEGER(dim1)[0];
    cols = INTEGER(dim1)[1]; 
    
    PM = NUMERIC_POINTER(AS_NUMERIC(PMmat));
    MM = NUMERIC_POINTER(AS_NUMERIC(MMmat));

    ProbeNames = (char **)Calloc(rows,char *);
    for (i =0; i < rows; i++)
      ProbeNames[i] = CHAR(VECTOR_ELT(ProbeNamesVec,i));

    IdealMM_correct(PM,MM, &rows, &cols,ProbeNames);
    Free(ProbeNames);
    UNPROTECT(1);
    
  } else if (asInteger(bg_type) == 6){
    /* LESN proposal 2 - half gaussian weighting */
    LESN_param = GetParameter(background_param,"lesnparam");


    Rprintf("Background correcting\n");
    PROTECT(dim1 = getAttrib(PMmat,R_DimSymbol)); 
    rows = INTEGER(dim1)[0];
    cols = INTEGER(dim1)[1]; 
    PM = NUMERIC_POINTER(AS_NUMERIC(PMmat));
    baseline = NUMERIC_POINTER(LESN_param)[0];
    theta = NUMERIC_POINTER(LESN_param)[1];
    LESN_correct(PM, rows, cols, 2, baseline,theta);
    UNPROTECT(1);
  } else if (asInteger(bg_type) == 7){
    /* LESN proposal 1 - half gaussian weighting */
    Rprintf("Background correcting\n");
    LESN_param = GetParameter(background_param,"lesnparam");
     
    PROTECT(dim1 = getAttrib(PMmat,R_DimSymbol)); 
    rows = INTEGER(dim1)[0];
    cols = INTEGER(dim1)[1]; 
    PM = NUMERIC_POINTER(AS_NUMERIC(PMmat));
    baseline = NUMERIC_POINTER(LESN_param)[0];
    theta = NUMERIC_POINTER(LESN_param)[1];
    LESN_correct(PM, rows, cols, 1, baseline, theta);
    UNPROTECT(1);
  } else if (asInteger(bg_type) == 8){
    /*      LESN proposal 0 - shifting */
    Rprintf("Background correcting\n");
    LESN_param = GetParameter(background_param,"lesnparam");
    
    PROTECT(dim1 = getAttrib(PMmat,R_DimSymbol)); 
    rows = INTEGER(dim1)[0];
    cols = INTEGER(dim1)[1]; 
    PM = NUMERIC_POINTER(AS_NUMERIC(PMmat));
    baseline = NUMERIC_POINTER(LESN_param)[0];
    theta = NUMERIC_POINTER(LESN_param)[1];
    LESN_correct(PM, rows, cols, 0, baseline, theta);
    UNPROTECT(1);

  } else if (asInteger(bg_type) == 9){
    /* Specific biweight correction */
    Rprintf("Background correcting\n");

    PROTECT(dim1 = getAttrib(PMmat,R_DimSymbol)); 
    rows = INTEGER(dim1)[0];
    cols = INTEGER(dim1)[1]; 
    
    PM = NUMERIC_POINTER(AS_NUMERIC(PMmat));
    MM = NUMERIC_POINTER(AS_NUMERIC(MMmat));

    ProbeNames = (char **)Calloc(rows,char *);
    for (i =0; i < rows; i++)
      ProbeNames[i] = CHAR(VECTOR_ELT(ProbeNamesVec,i));

    SpecificBiweightCorrect(PM,MM, &rows, &cols,ProbeNames);
    Free(ProbeNames);
    UNPROTECT(1);


  }

  return PMmat;
}


/********************************************************************************************
 **
 ** void pp_normalize(SEXP PMmat, SEXP MMmat, SEXP ProbeNamesVec,SEXP N_probes,SEXP norm_type)
 **
 ** SEXP PMmat - matrix of Perfect-match values
 ** SEXP MMmat - matrix of Mismatch values
 ** SEXP ProbeNamesVec - vector containing names of probeset for each probe
 ** SEXP N_probes - number of PM/MM probes on an array
 ** SEXP norm_type  - an integer indicating the normalization method to be used.
 **
 ** a function to manipulate the data so that the R objects are converted
 ** into C data structures. The Normalization is then applied to the C data structures. 
 ** The summary statistic is applied only to PM data.
 **
 ** this function assumes any sort of background correction has been applied and only the
 ** PM probes need be normalized.
 **
 **
 ** norm_type  1 - Quantile Normalization
 ** norm_type  2 - Probeset Quantile Normalization
 ** norm_type  3 - Scaling normalization
 ** 
 *******************************************************************************************/

SEXP pp_normalize(SEXP PMmat, SEXP MMmat, SEXP ProbeNamesVec,SEXP N_probes,SEXP norm_type, SEXP norm_parameters){


  SEXP param;
  
  int rows, cols;
  /* double *outexpr; */
  double *PM,*MM;
  /*  char **outnames; */
  char **ProbeNames; 
  int i;
  int nprobesets;
  
  double trim;
  int baseline;

  int usemedian;
  int uselog2;


  SEXP dim1;
  /* SEXP outvec,outnamesvec;
     SEXP dimnames,names;*/
  
  PROTECT(dim1 = getAttrib(PMmat,R_DimSymbol)); 
  rows = INTEGER(dim1)[0];
  cols = INTEGER(dim1)[1]; 

  PM = NUMERIC_POINTER(AS_NUMERIC(PMmat));
  MM = NUMERIC_POINTER(AS_NUMERIC(MMmat));
  
  nprobesets=INTEGER(N_probes)[0];
  
  /*  printf("%i\n",nprobesets); */
  /* printf("%d ",INTEGER(norm_flag)[0]); */
  /* normalize PM using quantile normalization */
  /* printf("Normalizing\n"); */
  Rprintf("Normalizing\n");



  if (asInteger(norm_type) == 1){
    qnorm_c(PM,&rows,&cols);
  } else if (asInteger(norm_type) == 2) { 
    ProbeNames = (char **)Calloc(rows,char *);
    for (i =0; i < rows; i++)
      ProbeNames[i] = CHAR(VECTOR_ELT(ProbeNamesVec,i));
    
    param = GetParameter(norm_parameters,"use.median");
    usemedian = asInteger(param);
    param = GetParameter(norm_parameters,"use.log2");
    uselog2 = asInteger(param);
    qnorm_probeset_c(PM, rows, cols, nprobesets, ProbeNames, usemedian, uselog2);
    Free(ProbeNames);
  } else if (asInteger(norm_type) == 3) { 
    param = GetParameter(norm_parameters,"scaling.trim");
    trim = asReal(param);
    param = GetParameter(norm_parameters, "scaling.baseline");
    baseline = asInteger(param);
    scaling_norm(PM, rows, cols,trim, baseline);
  }
  UNPROTECT(1);
  return PMmat;
}

