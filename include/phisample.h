#ifndef __PHI_SAMPLE_H__
#define __PHI_SAMPLE_H__

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

// A data structure representing simulating requirements
typedef struct {
    double kbt;
    double kphi;
    double bondlength;
    double bondangle;
    size_t length;
    size_t capacity;
    int *samples;
    int *repeats;
} phisample;

double energy_function(double kphi, double phi);

// read requirements from file
phisample *read_phisample(FILE *f);

// print data of the sample
void print_data(phisample *samp);

// print header of the sample
void print_info_phisample();

#endif

