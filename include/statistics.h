#ifndef __RANDOM_AVERAGE_H__
#define __RANDOM_AVERAGE_H__

#include <math.h>

double average(double (*f)(), int n) ;
double mean(double *a, int n) ;
double meansq(double *a, int n) ;
double stddev(double *a, int n) ;

#endif
