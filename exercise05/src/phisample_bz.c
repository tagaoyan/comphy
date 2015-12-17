#include <gsl/gsl_rng.h>
#include <math.h>
#include <time.h>

#include "statistics.h"
#include "phisample.h"


void sample_bz(double *a, int n, double kphi, double kbt, gsl_rng *rng, gsl_rng *rng2) {
    // get n samples that follows Boltzmann distribution
    for (int i = 0; i < n; i++) {
        for (; ;) {
            double x = gsl_rng_uniform(rng) * 2 * M_PI - M_PI;
            double t = gsl_rng_uniform(rng2);
            if (t < exp(- energy_function(kphi, x) / kbt)) {
                a[i] = x;
                break;
            }
        }
    }
}

double partation_bz(double *a, int n, double kphi, double kbt) {
    // calculate partation function of n samples
    double z[n];
    for (int i = 0; i < n; i++) {
        z[i] = exp(energy_function(kphi, a[i]) / kbt);
    }
    return 2 * M_PI / mean(z, n);
}

double energy_bz(double *a, int n, double kphi) {
    // calculate internal energy from n samples
    double u[n];
    for (int i = 0; i < n; i++) {
        u[i] = energy_function(kphi, a[i]);
    }
    return mean(u, n);
}

void run_phisample_bz(phisample *samp, gsl_rng *rng, gsl_rng *rng2) {
    double kphi = samp->kphi;
    double kbt = samp->kbt;
    print_data(samp);
    for (int i = 0; i < samp->length; i++) {
        int n = samp->samples[i];
        int r = samp->repeats[i];
        double en[r], ent[r], fe[r]; // internal energy, entrophy, free energy
        for (int j = 0; j < r; j++) {
            double s[n]; // samples
            sample_bz(s, n, kphi, kbt, rng, rng2);
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
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng *rng2 = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, time(NULL));
    gsl_rng_set(rng2, time(NULL) / 2);
    phisample *samp;
    print_info_phisample();
    samp = read_phisample(stdin);
    run_phisample_bz(samp, rng, rng2);
    free(samp);
    gsl_rng_free(rng);
    gsl_rng_free(rng2);
    return 0;
}
