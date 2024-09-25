#include "common_types.h"

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

PLM_modelfit *new_PLMmodelfit(){
  
  PLM_modelfit *x = R_Calloc(1,PLM_modelfit);
  x->cur_params = R_Calloc(1,double);
  x->cur_weights= R_Calloc(1,double);
  x->cur_resids= R_Calloc(1,double);
  x->cur_varcov= R_Calloc(1,double);
  x->cur_residSE= R_Calloc(1,double);
  x->X= R_Calloc(1,double);   
  x->n =0;
  x->p=0;
  x->nprobes=0;
    


  return x;

}

void free_PLMmodelfit(PLM_modelfit *x){
  
  R_Free(x->cur_params);
  R_Free(x->cur_weights);
  R_Free(x->cur_resids);
  R_Free(x->cur_varcov);
  R_Free(x->cur_residSE);
  R_Free(x->cur_se_estimates);
  R_Free(x->X);
  R_Free(x);
}



