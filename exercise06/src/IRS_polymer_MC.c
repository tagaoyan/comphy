#include "zm.h"
#include "xyz.h"
#include "xyzbase.h"
#include "zmxyz.h"
#include "randomcoil.h"
#include "statistics.h"
#include "thermo.h"
#include "polymer_sample.h"
#include <gsl/gsl_rng.h>
#include <time.h>

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
    return sqrt(rg);
}


typedef struct {
    double partation;
    double energy;
    double radius_of_gyration;
} sample_result;

sample_result sample_irs(size_t n, double b, double th, size_t len, double kphi, double epsilon, double alpha, double sigma, double kbt, gsl_rng *rng1, gsl_rng *rng2) {
    double t[n]; // temp
    double w[n]; // temp
    double vdw[n];
    double tors[n];
    double rg[n];
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
        vdw[i] = energy_vdw(ch, epsilon, alpha, sigma);
        tors[i] = energy_tors(ch, kphi);
        rg[i] = radius_of_gyration(ch);
        t[i] = exp(tors[i] / kbt);
        free_chain(ch);
    }
    double z = pow(2 * M_PI, len - 3) / mean(t, n);
    for (int i = 0; i < n; i++) {
        w[i] = t[i] = exp(- vdw[i] / kbt);
    }
    double u = 1 / mean(t, n);
    double rog = u;
    z /= u;
    for (int i = 0; i < n; i++) {
        t[i] *= tors[i] + vdw[i];
        w[i] *= rg[i];
    }
    u *= mean(t, n);
    rog *= mean(w, n);
    sample_result r;
    r.partation = z;
    r.energy = u;
    r.radius_of_gyration = rog;
    return r;
}

void run_polymer_sample_irs(polymer_sample *samp, gsl_rng *rng1, gsl_rng *rng2) {
    print_info_polymer();
    print_data_polymer(samp);
    for (int i = 0; i < samp->length; i++) {
        int n = samp->samples[i];
        int r = samp->repeats[i];
        double en[r], ent[r], fe[r], rg[r];
        for (int j = 0; j < r; j++) {
            sample_result rt = sample_irs(n,
                    samp->bondlength,
                    samp->bondangle,
                    samp->chainlength,
                    samp->kphi,
                    samp->epsilon,
                    samp->alpha,
                    samp->sigma,
                    samp->kbt,
                    rng1,
                    rng2
                    );
            double z = rt.partation;
            en[j] = rt.energy;
            fe[j] = free_energy(z, samp->kbt);
            ent[j] = entrophy(en[j], fe[j], samp->kbt);
            rg[j] = rt.radius_of_gyration;
        }
        printf("%d %d %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n", r, n,
                mean(en, r), stddev(en, r),
                mean(fe, r), stddev(fe, r),
                mean(ent, r), stddev(ent, r),
                mean(rg, r), stddev(rg, r));
    }
}

int main() {
    polymer_sample *samp = read_polymer_sample(stdin);
    gsl_rng *rng1 = gsl_rng_alloc(gsl_rng_taus2);
    gsl_rng *rng2 = gsl_rng_alloc(gsl_rng_taus2);
    gsl_rng_set(rng1, time(NULL));
    gsl_rng_set(rng2, time(NULL) / 2);
    run_polymer_sample_irs(samp, rng1, rng2);
    gsl_rng_free(rng1);
    gsl_rng_free(rng2);
    return 0;
}
