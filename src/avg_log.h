#ifndef AVG_LOG_H
#define AVG_LOG_H

#include "threestep_summary_methods_param.h"

void AverageLog(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes, double *resultsSE, summary_plist *summary_param);
void AverageLog_noSE(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes);

void AverageLog_PLM(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes, double *resultsSE, double *residuals, summary_plist *summary_param);


#endif
