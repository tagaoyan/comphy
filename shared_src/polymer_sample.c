#include "polymer_sample.h"
#include <stdio.h>
#include <stdlib.h>

#define BUF_SIZE 4096
polymer_sample *read_polymer_sample(FILE *f) {
    char buf[BUF_SIZE];
    polymer_sample *samp = malloc(sizeof(polymer_sample));
    samp->length = 0;
    samp->capacity = 1;
    samp->samples = malloc(sizeof(int) * samp->capacity);
    samp->repeats = malloc(sizeof(int) * samp->capacity);
    samp->epsilon = -1; // mark for vdw parameters
    samp->alpha = -1; // for non-soft vdw
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
    puts("# 9.<radius of gyration> 10.<error of Rg> (kb)");
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
    if (samp->alpha > 0) {
        printf("# %s: %g\n", "alpha", samp->alpha);
    }
}

double radius_of_gyration(chain *ch) {
    chain_xyz *chx;
    chx = new_chain_xyz(ch->length);
    chain_xyz_from_zm(chx, ch);
    double x[ch->length], y[ch->length], z[ch->length];
    for (int i = 0; i < ch->length; i++) {
        x[i] = chx->coordinates[i].x;
        y[i] = chx->coordinates[i].y;
        z[i] = chx->coordinates[i].z;
    }
    coordinate_xyz rm;
    rm.x = mean(x, ch->length);
    rm.y = mean(y, ch->length);
    rm.z = mean(z, ch->length);
    double rg = 0;
    for (int i = 0; i < ch->length; i++) {
        rg += vect_length_sq(vect_minus(rm, chx->coordinates[i]));
    }
    rg /= ch->length;
    free_chain_xyz(chx);
    return sqrt(rg);
}
