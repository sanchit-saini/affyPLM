#ifndef BIWEIGHT_H
#define BIWEIGHT_H 1

#include "threestep_summary_methods_param.h"

void tukeybiweight(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes, double *resultsSE, summary_plist *summary_param);
double Tukey_Biweight(double *x, int length);
void tukeybiweight_PLM(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes, double *resultsSE, double *residuals, summary_plist *summary_param);


#endif
