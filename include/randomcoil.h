#ifndef __RANDOMCOIL_H__
#define __RANDOMCOIL_H__

#include "zm.h"
#include <gsl/gsl_rng.h>
#include <stdio.h>

chain *randomcoil(double blength, double bangle, size_t length, gsl_rng *rng);
chain *randomcoil_from_file(FILE *f, gsl_rng *rng);

#endif
