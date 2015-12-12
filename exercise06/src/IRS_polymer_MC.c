#include "zm.h"
#include "xyz.h"
#include "zmxyz.h"
#include "randomcoil.h"
#include <gsl/gsl_rng.h>

double energy_vdw(chain *ch, double epsilon, double alpha, double sigma) {
    chain_xyz *chx = new_chain_xyz(ch->length);
    chain_xyz_from_zm(chx, ch);
    double evdw = 0;
    for (int i = 1; i <= ch->length; i++) {
        for (int j = i + 1; j <= ch->length; j++) {
            double r = vect_length(vect_minus(chx->coordinates[i], chx->coordinates[j]));
            if (r <= alpha * sigma) {
                evdw += 4 * epsilon * (pow(alpha, -12) - pow(alpha, -6));
            } else {
                evdw += 4 * epsilon * (pow(sigma / r, 12) - pow(sigma / r, 6));
            }
        }
    }
    free_chain_xyz(chx);
    return evdw;
}

double energy_tors(chain *ch, double kphi) {
    double et = 0;
    for (int i = 4; i <= ch->length; i++) {
        et += 0.5 * kphi * (1 + cos(3 * ch->torsionangles[i] * M_PI / 180));
    }
    return et;
}

void sample_irs(chain **chs, size_t n, double b, double th, size_t len, double kphi, double kbt, gsl_rng *rng1, gsl_rng *rng2) {
    for (int i = 0; i < n; i++) {
        for (;;) {
            chain *ch = randomcoil(b, th, len, rng1);
            double t = exp(-energy_tors(ch, kphi) / kbt);
            double x = gsl_rng_uniform(rng2);
            if (x < t) {
                chs[i] = ch;
            } else {
                free_chain(ch);
            }
        }
    }
}
