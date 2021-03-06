#include "phisample.h"

#define BUF_SIZE 4096

double energy_function(double kphi, double phi) {
    return 0.5 * kphi * (1 + cos(3 * phi));
}

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

void print_data(phisample *samp) {
    printf("# bond length: %g\n", samp->bondlength);
    printf("# bond angle: %g\n", samp->bondangle);
    printf("# kphi: %g\n", samp->kphi);
    printf("# kbt: %g\n", samp->kbt);
}

void print_info_phisample() {
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
    puts("# 3.<free energy> 4.<error of free energy> (kb * 1K)");
    puts("# 5.<internal energy> 6.<error of internal energy> (kb * 1K)");
    puts("# 7.<entrophy> 8.<error of entrophy> (kb)");
}
