/*********************************************************************
 **
 ** file: common_types.h
 **
 ** Aim: Define some structures that will be used for passing
 **      related data items to model fitting routines.
 **
 ** Copyright (C) 2003 Ben Bolstad
 **
 ** created by: B. M. Bolstad <bolstad@stat.berkeley.edu>
 ** 
 ** created on: Sept 04, 2003
 **
 ** History
 ** Sept 04, 2003 - Initial version (outputsettings, DataGroup, 
 **                 PLMoutput, RPLMoutput, PLMmodelparam)      
 ** Sept 13, 2003 - Added n_rlm_iterations and init_method to
 **                 PLMmodelparam (these will control how many
 **                 iterations of iteratively reweighted least
 **                 squares and how we intialize the IRLS: least 
 **                 squares or median polish.
 **
 *********************************************************************/

#ifndef COMMON_TYPES_H
#define COMMON_TYPES_H 1

#include "threestep_summary_methods.h"




#include <R.h>
#include <Rinternals.h>

/*******************************************************************
 **
 ** a structure for holding flags and settings on what should
 ** be stored as output from the PLM routines
 **
 **
 **
 ******************************************************************/

typedef struct {
  int weights;   /* Store Weights */
  int residuals;     /* Store Residuals */
  int residSE;   /* Store ResidSE */
  int pseudoSE;  /* Compute Pseudo SE in the case of an rmaPLM */
  int varcov;   /* Store varcov matrices what type. 0 means none, 1 means only chip level, 2 means all,*/
} outputsettings;



/*******************************************************************
 **
 ** a structure for holding probe intensity data as it gets
 ** used by the PLM routines
 **
 **
 **
 ******************************************************************/

typedef struct {
  double *PM;
  double *MM;
  int rows;
  int cols;
  int nprobesets;
  char **ProbeNames;
} Datagroup;

/*******************************************************************
 **
 ** a structure for holding computed quantities that result from
 ** using the PLM routines
 **
 **
 **
 ******************************************************************/

typedef struct {
  char **outnames;
  double *out_weights;
  double *out_probeparams;
  double *out_chipparams;
  double *out_constparams;
  double *out_probe_SE;
  double *out_chip_SE;
  double *out_const_SE;
  double *out_resids;
  double *out_residSE;
  double **out_varcov;

} PLMoutput;


/*******************************************************************
 **
 ** a structure for holding computed quantities as they will be 
 ** returned to R (ie in R data structures
 **
 **
 **
 ******************************************************************/

typedef struct {
  SEXP weights;
  SEXP probe_coef; 
  SEXP chip_coef;
  SEXP const_coef;
  SEXP chip_SE;
  SEXP probe_SE;
  SEXP const_SE;
  SEXP residuals;
  SEXP residSE;
  SEXP varcov;
  int nprotected;
} PLMRoutput;

/*******************************************************************
 **
 ** a structure for holding parameters and flags that describe the
 ** PLM model being fitted.
 ** 
 **
 **
 **
 ******************************************************************/

typedef struct{
  int nchipparams;
  int method;  
  int se_method;
  int psi_code;
  double psi_k;
  double *input_chipcovariates;
  int n_rlm_iterations;
  int init_method;
  pt2PLMSummary PLM3stepSummary;
} PLMmodelparam;





#endif
