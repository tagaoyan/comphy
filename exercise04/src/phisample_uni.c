#include <gsl/gsl_rng.h>
#include <time.h>
#include <math.h>

#include "phisample.h"
#include "statistics.h"


void sample(double *a, int n, gsl_rng *rng) {
    // uniformly sample n phi angle to a
    for (int i = 0; i < n; i++) {
        a[i] = gsl_rng_uniform(rng) * 2 * M_PI - M_PI;
    }
}

double partation(double *a, int n, double kphi, double kbt) {
    // calculate partation function of n samples
    double z[n];
    for (int i = 0; i < n; i++) {
        z[i] = exp(- energy_function(kphi, a[i]) / kbt);
    }
    return mean(z, n) * 2 * M_PI;
}


double energy(double *a, int n, double kphi, double kbt, double z) {
    // calculate internal energy from n samples
    double u[n];
    for (int i = 0; i < n; i++) {
        double e = energy_function(kphi, a[i]);
        u[i] = e * exp(-e/kbt);
    }
    return mean(u, n) / z * 2 * M_PI;
}

void run_phisample(phisample *samp, gsl_rng *rng) {
    double kphi = samp->kphi;
    double kbt = samp->kbt;
    printf("# bond length: %g\n", samp->bondlength);
    printf("# bond angle: %g\n", samp->bondangle);
    printf("# kphi: %g\n", samp->kphi);
    printf("# kbt: %g\n", samp->kbt);
    for (int i = 0; i < samp->length; i++) {
        int n = samp->samples[i];
        int r = samp->repeats[i];
        double en[r], ent[r], fe[r]; // internal energy, entrophy, free energy
        for (int j = 0; j < r; j++) {
            double s[n]; // samples
            sample(s, n, rng);
            double z = partation(s, n, kphi, kbt);
            fe[j] = free_energy(z, kbt);
            en[j] = energy(s, n, kphi, kbt, z);
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
    gsl_rng_set(rng, time(NULL));
    phisample *samp;
    print_info_phisample();
    samp = read_phisample(stdin);
    run_phisample(samp, rng);
    free(samp);
    gsl_rng_free(rng);
    return 0;
}
