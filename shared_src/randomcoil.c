#include "randomcoil.h"

#define NAME_SIZE 8
chain *randomcoil(double blength, double bangle, size_t length, gsl_rng *rng) {
    // create a coil, torsion angle is random
    // a GSL random number generator (gsl_rng) is required to generate random numbers
    // the result returned should be freed by free_chain().
    chain *ch = new_chain(length);
    for (int i = 0; i < length; i++) {
        char name[NAME_SIZE];
        snprintf(name, NAME_SIZE, "A%d", i + 1);
        chain_add(ch, name, blength, bangle, gsl_rng_uniform(rng) * 360. - 180.);
    }
    return ch;
}

#define BUF_SIZE 4096
chain *randomcoil_from_file(FILE *f, gsl_rng *rng) {
    char buf[BUF_SIZE];
    while (fgets(buf, BUF_SIZE, f) != NULL) {
        if (buf[0] == '#') {
            continue;
        }
        double blength, bangle;
        size_t length;
        if (sscanf(buf, "%lf %lf %zu", &blength, &bangle, &length) == 3) {
            return randomcoil(blength, bangle, length, rng);
        }
    }
    return NULL;
}
