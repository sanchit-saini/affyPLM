#ifndef MEDIAN_LOGPM_H
#define MEDIAN_LOGPM_H 1

void MedianLogPM(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes, double *resultsSE);
void MedianLogPM_noSE(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes);

#endif
