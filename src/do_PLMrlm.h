#ifndef DO_PLMRLM_H
#define DO_PLMRLM_H 1

void do_PLMrlm(double *PM, char **ProbeNames, int *rows, int *cols, int nps, int nchipparams, int method,int se_method, double *chipcovariates, char **outNames, double *out_weights, double *out_probeparams, double *out_chipparams, double *out_constparams,double *out_probeSE, double *out_chipSE, double *out_constSE,int psitype,double psi_k);

#endif
