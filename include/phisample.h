#ifndef __PHI_SAMPLE_H__
#define __PHI_SAMPLE_H__

#include <stdio.h>
#include <stdlib.h>

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

phisample *read_phisample(FILE *f);

void print_data(phisample *samp);

void print_info_phisample();

#endif

