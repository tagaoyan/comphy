#include "zm.h"
#include "xyz.h"
#include "xyzbase.h"
#include "zmxyz.h"
#include "randomcoil.h"
#include "statistics.h"
#include <gsl/gsl_rng.h>

#include "phisample.h"

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

// calculate them together to avoid redundant calculation
typedef struct {
    double partation;
    double energy;
} partation_and_energy;

partation_and_energy partation_and_energy_irs(chain **chs, size_t n, double epsilon, double alpha, double sigma, double kphi, double kbt) {
    double t[n];
    double vdw[n];
    double tors[n];
    double z, u;
    for (int i = 0; i < n; i++) {
        vdw[i] = energy_vdw(chs[i], epsilon, alpha, sigma);
        tors[i] = energy_tors(chs[i], kphi);
    }
    for (int i = 0; i < n; i++) {
        t[i] = exp(tors[i] / kbt);
    }
    z = pow(2 * M_PI, chs[0]->length - 3) / mean(t, n);
    for (int i = 0; i < n; i++) {
        t[i] = exp(- vdw[i] / kbt);
    }
    u = 1 / mean(t, n);
    z /= u;
    for (int i = 0; i < n; i++) {
        t[i] *= (tors[i] + vdw[i]);
    }
    u *= mean(t, n);
    partation_and_energy r;
    r.partation = z;
    r.energy = u;
    return r;
}

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

#define BUF_SIZE 4096
polymer_sample *read_polymer_sample(FILE *f) {
    char buf[BUF_SIZE];
    polymer_sample *samp = malloc(sizeof(polymer_sample));
    samp->length = 0;
    samp->capacity = 1;
    samp->samples = malloc(sizeof(int) * samp->capacity);
    samp->repeats = malloc(sizeof(int) * samp->capacity);
    samp->epsilon = -1; // mark for vdw parameters
    while (fgets(buf, BUF_SIZE, f) != NULL) {
        if (buf[0] == '#') {
            continue;
        }
        int cl;
        double b, th, kbt, kphi;
        double epsilon, sigma, alpha;
        int smp, rep;
        if (sscanf(buf, "%lf %lf %lf %lf", &b, &th, &kbt, &kphi) == 4) {
            samp->bondlength = b;
            samp->bondangle = th;
            samp->kbt = kbt;
            samp->kphi = kphi;
        } else if (sscanf(buf, "%lf %lf %lf", &epsilon, &sigma, &alpha) == 3) {
            // soft vdw
            samp->epsilon = epsilon;
            samp->sigma = sigma;
            samp->alpha = alpha;
        } else if (sscanf(buf, "%lf %lf", &epsilon, &sigma) == 2) {
            // vdw or sample data
            if (samp->epsilon < 0) {
                samp->epsilon = epsilon;
                samp->sigma = sigma;
            } else {
                sscanf(buf, "%d %d", &smp, &rep);
                if (samp->length >= samp->capacity) {
                    samp->capacity *= 2;
                    samp->samples = realloc(samp->samples, sizeof(int) * samp->capacity);
                    samp->repeats = realloc(samp->repeats, sizeof(int) * samp->capacity);
                }
                samp->samples[samp->length] = smp;
                samp->repeats[samp->length] = rep;
                samp->length++;
            }
        } else if (sscanf(buf, "%d", &cl) == 1) {
            samp->chainlength = cl;
        }
    }
    return samp;
}

void run_polymer_sample_irs(polymer_sample *samp, gsl_rng rng1, gsl_rng rng2) {
    for (int i = 0; i < samp->length; i++) {
        int n = samp->samples[i];
        int r = samp->repeats[i];
        double en[r], ent[r], fe[r]; // internal energy, entrophy, free energy
        for (int j = 0; j < r; j++) {
            double s[n]; // samples
            sample_bz(s, n, kphi, rng, rng2);
            double z = partation_bz(s, n, kphi, kbt);
            fe[j] = free_energy(z, kbt);
            en[j] = energy_bz(s, n, kphi);
            ent[j] = entrophy(en[j], fe[j], kbt);
        }
        printf("%d %d %.6f %.6f %.6f %.6f %.6f %.6f\n", r, n,
                mean(en, r), meansq(en, r),
                mean(fe, r), stddev(fe, r),
                mean(ent, r), stddev(ent, r));
    }
}

int main() {
    polymer_sample *sp = read_polymer_sample(stdin);
    run_polymer_sample_irs(sp);
    return 0;
}
