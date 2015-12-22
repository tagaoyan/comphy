#include "zm.h"
#include "xyz.h"
#include "xyzbase.h"
#include "zmxyz.h"
#include "randomcoil.h"
#include "statistics.h"
#include "thermo.h"
#include <gsl/gsl_rng.h>

#include "phisample.h"

double energy_vdw(chain *ch, double epsilon, double alpha, double sigma) {
    chain_xyz *chx;
    chx = new_chain_xyz(ch->length);
    chain_xyz_from_zm(chx, ch);
    double evdw = 0;
    for (int i = 4; i < ch->length; i++) {
        for (int j = 0; j <= i - 4; j++) {
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

void sample_irs(chain **chs, size_t n, double b, double th, size_t len, double kphi, double kbt, gsl_rng *rng1, gsl_rng *rng2) {
    for (int i = 0; i < n; i++) {
        chain *ch = new_chain(len);
        for (int j = 0; j < len; j++) {
            for (;;) {
                double phi = gsl_rng_uniform(rng1) * 2 * M_PI - M_PI;
                double t = exp(-energy_tors_single(kphi, phi) / kbt);
                double x = gsl_rng_uniform(rng2);
                if (x < t) {
                    chain_add(ch, "", b, th, phi * 180 / M_PI);
                    break;
                }
            }
        }
        chs[i] = ch;
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

void print_data_polymer(polymer_sample *samp) {
    puts("# chain parameters:");
    printf("# %s: %d\n", "chain length", samp->chainlength);
    printf("# %s: %g\n", "bond length", samp->bondlength);
    printf("# %s: %g\n", "bond angle", samp->bondangle);
    puts("# temperature:");
    printf("# %s: %g\n", "kbt", samp->kbt);
    puts("# torsional energy parameters:");
    printf("# %s: %g\n", "kphi", samp->kphi);
    puts("# vdw energy parameters:");
    printf("# %s: %g\n", "epsilon", samp->epsilon);
    printf("# %s: %g\n", "sigma", samp->sigma);
    printf("# %s: %g\n", "alpha", samp->alpha);
}

void print_info_polymer() {
    puts("# Result of polymer sample");
    char dstr[32];
    time_t t;
    struct tm tmp;
    t = time(NULL);
    localtime_r(&t, &tmp);
    strftime(dstr, 32, "%F %T", &tmp);
    printf("# %s\n", dstr);
    puts("# Format: There are 8 columns in one line. Comment lines begin with #");
    puts("# 1.<number of repeats> 2.<number of samples>");
    puts("# 3.<internal energy> 4.<error of internal energy_function> (kb * 1K)");
    puts("# 5.<free energy> 6.<error of free energy> (kb * 1K)");
    puts("# 7.<entrophy> 8.<error of entrophy> (kb)");
}

void run_polymer_sample_irs(polymer_sample *samp, gsl_rng *rng1, gsl_rng *rng2) {
    print_info_polymer();
    print_data_polymer(samp);
    for (int i = 0; i < samp->length; i++) {
        int n = samp->samples[i];
        int r = samp->repeats[i];
        double en[r], ent[r], fe[r]; // internal energy, entrophy, free energy
        for (int j = 0; j < r; j++) {
            chain *s[n]; // samples
            sample_irs(s, n,
                    samp->bondlength,
                    samp->bondangle,
                    samp->chainlength,
                    samp->kphi,
                    samp->kbt,
                    rng1,
                    rng2
                    );
            partation_and_energy zu = partation_and_energy_irs(s, n,
                    samp->epsilon,
                    samp->alpha,
                    samp->sigma,
                    samp->kphi,
                    samp->kbt
                    );
            double z = zu.partation;
            en[j] = zu.energy;
            fe[j] = free_energy(z, samp->kbt);
            ent[j] = entrophy(en[j], fe[j], samp->kbt);
            for (int k = 0; k < n; k++) {
                free_chain(s[k]);
            }
        }
        printf("%d %d %.6f %.6f %.6f %.6f %.6f %.6f\n", r, n,
                mean(en, r), stddev(en, r),
                mean(fe, r), stddev(fe, r),
                mean(ent, r), stddev(ent, r));
    }
}

int main() {
    polymer_sample *samp = read_polymer_sample(stdin);
    gsl_rng *rng1 = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng *rng2 = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng1, time(NULL));
    gsl_rng_set(rng2, time(NULL) / 2);
    run_polymer_sample_irs(samp, rng1, rng2);
    gsl_rng_free(rng1);
    gsl_rng_free(rng2);
    return 0;
}
