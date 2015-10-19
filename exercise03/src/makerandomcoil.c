#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_rng.h>

#include "randomcoil.h"
#include "zm.h"
#include "xyz.h"
#include "zmxyz.h"

int main(void) {
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, time(NULL));
    chain *rc = randomcoil_from_file(stdin, rng);
    print_chain(rc, stdout);
    chain_xyz *rcx = new_chain_xyz(1);
    chain_xyz_from_zm(rcx, rc);
    print_chain_xyz(rcx, stderr);
    free_chain(rc);
    free_chain_xyz(rcx);
    gsl_rng_free(rng);
    return 0;
}
