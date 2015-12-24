#ifndef __POLYMER_SAMPLE_H__
#define __POLYMER_SAMPLE_H__
#include <stdio.h>
#include <time.h>
typedef struct {
    int chainlength;
    double bondlength;
    double bondangle;
    double kbt;
    double kphi;
    double epsilon;
    double sigma;
    double alpha;
    size_t length;
    size_t capacity;
    int *samples;
    int *repeats;
} polymer_sample;

polymer_sample *read_polymer_sample(FILE *f);
void print_data_polymer(polymer_sample *samp);
void print_info_polymer();
#endif
