#ifndef RLM_SE_H
#define RLM_SE_H 1


void rlm_compute_se(double *X,double *Y, int n, int p, double *beta, double *resids,double *weights,double *se_estimates, int method,double (* PsiFn)(double, double, int), double psi_k);


#endif
