#include <gsl/gsl_rng.h>
#include <time.h>
#include <math.h>

#include "statistics.h"

double energy_function(double kphi, double phi) {
    return 0.5 * kphi * (1 + cos(3 * phi));
}

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

double free_energy(double z, double kbt) {
    // calculate free energy from partation function
    return -kbt * log(z);
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

double entrophy(double u, double fe, double kbt) {
    // calculate entrophy from internal energy and free energy
    return (u - fe) / kbt;
}

// A data structure representing simulating requirements
typedef struct {
    double kbt;
    double kphi;
    double bondlength;
    double bondangle;
    size_t length;
    size_t capacity;
    int *samples;
    int *repeats;
} phisample;

#define BUF_SIZE 4096
phisample *read_phisample(FILE *f) {
    // read simulation requirements from file
    char buf[BUF_SIZE];
    phisample *samp = malloc(sizeof(phisample));
    samp->length = 0;
    samp->capacity = 1;
    samp->samples = malloc(sizeof(int) * samp->capacity);
    samp->repeats = malloc(sizeof(int) * samp->capacity);
    while (fgets(buf, BUF_SIZE, f) != NULL) {
        if (buf[0] == '#') {
            continue;
        }
        double b, th, kbt, kphi;
        int smp, rep;
        if (sscanf(buf, "%lf %lf %lf %lf", &b, &th, &kbt, &kphi) == 4) {
            samp->bondlength = b;
            samp->bondangle = th;
            samp->kbt = kbt;
            samp->kphi = kphi;
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
    }
    return samp;
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
void print_info() {
    puts("# Result of phisample");
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

int main() {
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, time(NULL));
    print_info();
    run_phisample(read_phisample(stdin), rng);
    return 0;
}
