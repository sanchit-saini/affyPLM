/********************************************************************
 **
 ** file: threestep_summary_methods.c
 **
 ** written by: B. M. Bolstad <bolstad@stat.berkeley.edu>
 ** Created on: Jan 11, 2003
 **
 ** Copyright (C) 2003 Ben Bolstad
 **
 ** aim: general interface to threestep summary methods
 **      to add another summary method only the header file should need 
 **      be modified.
 **
 ** last modified: Feb 6, 2003
 **
 ** History
 ** 
 ** Jan 11, 2003 - Initial version
 ** Jan 13, 2003 - added rlm_threestep method.
 ** Feb 6, 2003 - added four new methods: LogAverage, LogMedianPM, MedianLogPM, LogNthLargestPM
 ** Jul 23, 2003 - a three step method should return a SE estimate
 **
 ********************************************************************/

#include "threestep_summary_methods.h"
#include "stdlib.h"

int number_summary_methods=9;

pt2Summary funcArr[9];

pt2Summary SummaryMethod(int code){
  
  funcArr[0] = &median_polish;
  funcArr[1] = &tukeybiweight;
  funcArr[2] = &AverageLog;
  funcArr[3] = &rlm_threestep;
  funcArr[4] = &LogAverage;
  funcArr[5] = &LogMedianPM;
  funcArr[6] = &MedianLogPM;
  funcArr[7] = &LogNthLargestPM;
  funcArr[8] = &lm_threestep;

  return funcArr[code];
}
