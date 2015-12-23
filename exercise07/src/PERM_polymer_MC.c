#include "zm.h"
#include "xyz.h"
#include "xyzbase.h"
#include "zmxyz.h"
#include "randomcoil.h"
#include "statistics.h"
#include "thermo.h"
#include <gsl/gsl_rng.h>
#include <time.h>

double energy_vdw_single(chain_xyz *chx, int i, double epsilon, double sigma) {
    double evdw = 0;
    for (int j = 0; j <= i - 4; j++) {
        double r = vect_length(vect_minus(chx->coordinates[i], chx->coordinates[j]));
        evdw += 4 * epsilon * (pow(sigma / r, 12) - pow(sigma / r, 6));
    }
    return evdw;
}

double energy_vdw(chain *ch, double epsilon, double sigma) {
    chain_xyz *chx;
    chx = new_chain_xyz(ch->length);
    chain_xyz_from_zm(chx, ch);
    double evdw = 0;
    for (int i = 4; i < ch->length; i++) {
        evdw += energy_vdw_single(chx, i, epsilon, sigma);
    }
    free_chain_xyz(chx);
    return evdw;
}

double energy_tors_single(double kphi, double phi) {
    return 0.5 * kphi * (1 + cos(3 * phi));
}

double energy_tors(chain *ch, double kphi) {
    double et = 0;
    for (int i = 3; i < ch->length; i++) {
        et += energy_tors_single(kphi, ch->torsionangles[i] * M_PI / 180);
    }
    return et;
}

typedef struct {
    chain *chain;
    double weight;
    int len;
} perm_step;

void polymer_sample_perm_step(chain *ch, double weight,
        double w_lower, double w_upper,
        int len, double b, double th,
        double kphi,
        double epsilon, double sigma,
        double kbt,
        int k, // length of current chain
        chain **samp, int *cur, int n,
        gsl_rng *rng) {
    if (*cur >= n) {
        free_chain(ch);
        return;
    }
    double phi = gsl_rng_uniform(rng) * 2 * M_PI - M_PI;
    chain_add(ch, "", b, th, phi * 180 / M_PI);
    chain_xyz *chx = new_chain_xyz(k + 1);
    chain_xyz_from_zm(chx, ch);
    double dw = exp(- (energy_tors_single(kphi, phi) + energy_vdw_single(chx, k, epsilon, sigma)));
    free_chain_xyz(chx);
    weight *= dw;
    //printf("step %d, %.12f\n", k, weight);
    if (k + 1 == len) {
        samp[*cur] = ch;
        printf("got %d polymers, %.17f\n", *cur + 1, weight);
        ++*cur;
        return;
    } else {
        if (weight < w_lower) {
            if (gsl_rng_uniform(rng) < 0.5) {
                printf("X\n");
                free_chain(ch);
                return;
            } else {
                printf(">");
                polymer_sample_perm_step(ch, 2 * weight,
                        w_lower, w_upper,
                        len, b, th,
                        kphi,
                        epsilon, sigma,
                        kbt,
                        k + 1,
                        samp, cur, n,
                        rng);
            }
        } else if (weight > w_upper) {
            chain *ch2 = chain_dup(ch);
            printf("O");
            polymer_sample_perm_step(ch, weight / 2,
                    w_lower, w_upper,
                    len, b, th,
                    kphi,
                    epsilon, sigma,
                    kbt,
                    k + 1,
                    samp, cur, n,
                    rng);
            printf("Q");
            polymer_sample_perm_step(ch2, weight / 2,
                    w_lower, w_upper,
                    len, b, th,
                    kphi,
                    epsilon, sigma,
                    kbt,
                    k + 1,
                    samp, cur, n,
                    rng);
        } else {
            printf("-");
            polymer_sample_perm_step(ch, weight,
                    w_lower, w_upper,
                    len, b, th,
                    kphi,
                    epsilon, sigma,
                    kbt,
                    k + 1,
                    samp, cur, n,
                    rng);
        }
    }
}

void polymer_sample_perm(chain **samp, int n, // samples
        int len, double b, double th,
        double kphi,
        double epsilon, double sigma,
        double kbt,
        double w_lower, double w_upper,
        gsl_rng *rng) {
    int cur = 0;
    while (cur < n) {
        printf("!!! %d \n#", cur + 1);
        chain *ch = new_chain(len);
        ch->begin.x = 0;
        ch->begin.y = 0;
        ch->begin.z = 0;
        chain_add(ch, "", 0, 0, 0); // first
        chain_add(ch, "", b, 0, 0); // second
        chain_add(ch, "", b, th, 0); // third
        polymer_sample_perm_step(ch, 1,
                w_lower, w_upper,
                len, b, th,
                kphi,
                epsilon, sigma,
                kbt,
                3,
                samp, &cur, n,
                rng);
    };
}

int main() {
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, 5000);
    chain *samp[5000];
    polymer_sample_perm(samp, 5000,
            20, 1.5, 109.5,
            2,
            0.1, 3,
            1.,
            1e-10, 1e-7,
            rng);
    return 0;
}
