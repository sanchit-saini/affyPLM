#ifndef THREESTEP_SUMMARY_H 
#define THREESTEP_SUMMARY_H 1


void do_3summary(double *PM, char **ProbeNames, int *rows, int *cols, double *results, char **outNames, int nps,void (* SummaryMeth)(double*, int, int, int *,double *, int, double *),double *resultsSE);

#endif
