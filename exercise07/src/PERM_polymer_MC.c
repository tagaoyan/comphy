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


void polymer_sample_perm_step(chain *ch, double weight, double factor,
        double w_lower, double w_upper,
        int len, double b, double th,
        double kphi,
        double epsilon, double sigma,
        double kbt,
        int k, // length of current chain
        double *zs, double *us, double *rg, int *cur, int n, double *efn,
        gsl_rng *rng) {
    if (*cur >= n) {
        free_chain(ch);
        return;
    }
    double phi = gsl_rng_uniform(rng) * 2 * M_PI - M_PI;
    chain_add(ch, "", b, th, phi * 180 / M_PI);
    chain_xyz *chx = new_chain_xyz(k + 1);
    chain_xyz_from_zm(chx, ch);
    double dw = exp(- (energy_tors_single(kphi, phi) + energy_vdw_single(chx, k, epsilon, sigma)) / kbt);
    free_chain_xyz(chx);
    weight *= dw;
    if (k + 1 == len) {
        zs[*cur] = weight;
        us[*cur] = weight * (energy_vdw(ch, epsilon, sigma) + energy_tors(ch, kphi));
        rg[*cur] = weight * radius_of_gyration(ch);
        free_chain(ch);
        *efn += factor;
        ++*cur;
        return;
    } else {
        if (weight < w_lower) {
            if (gsl_rng_uniform(rng) < 0.5) {
                free_chain(ch);
                return;
            } else {
                polymer_sample_perm_step(ch, 2 * weight, 2 * factor,
                        w_lower, w_upper,
                        len, b, th,
                        kphi,
                        epsilon, sigma,
                        kbt,
                        k + 1,
                        zs, us, rg, cur, n, efn,
                        rng);
            }
        } else if (weight > w_upper) {
            chain *ch2 = chain_dup(ch);
            polymer_sample_perm_step(ch, weight / 2, factor / 2,
                    w_lower, w_upper,
                    len, b, th,
                    kphi,
                    epsilon, sigma,
                    kbt,
                    k + 1,
                    zs, us, rg, cur, n, efn,
                    rng);
            polymer_sample_perm_step(ch2, weight / 2, factor / 2,
                    w_lower, w_upper,
                    len, b, th,
                    kphi,
                    epsilon, sigma,
                    kbt,
                    k + 1,
                    zs, us, rg, cur, n, efn,
                    rng);
        } else {
            polymer_sample_perm_step(ch, weight, factor,
                    w_lower, w_upper,
                    len, b, th,
                    kphi,
                    epsilon, sigma,
                    kbt,
                    k + 1,
                    zs, us, rg, cur, n, efn,
                    rng);
        }
    }
}

sample_result polymer_sample_perm(int n, // samples
        int len, double b, double th,
        double kphi,
        double epsilon, double sigma,
        double kbt,
        double w_lower, double w_upper,
        gsl_rng *rng) {
    int cur = 0;
    double efn = 0;
    double zs[n], us[n], rg[n];
    while (cur < n) {
        chain *ch = new_chain(len);
        ch->begin.x = 0;
        ch->begin.y = 0;
        ch->begin.z = 0;
        chain_add(ch, "", 0, 0, 0); // first
        chain_add(ch, "", b, 0, 0); // second
        chain_add(ch, "", b, th, 0); // third
        polymer_sample_perm_step(ch, 1, 1,
                w_lower, w_upper,
                len, b, th,
                kphi,
                epsilon, sigma,
                kbt,
                3,
                zs, us, rg, &cur, n, &efn,
                rng);
    };
    sample_result rt;
    rt.partation = mean(zs, n) * n / efn * pow(2 * M_PI, len - 3);
    rt.energy = mean(us, n) / mean(zs, n);
    rt.radius_of_gyration = mean(rg, n) / mean(zs, n);
    return rt;
}

void run_polymer_sample_perm(polymer_sample *samp, gsl_rng *rng) {
    print_info_polymer();
    print_data_polymer(samp);
    for (int i = 0; i < samp->length; i++) {
        int n = samp->samples[i];
        int r = samp->repeats[i];
        double en[r], fe[r], ent[r], rg[r];
        for (int j = 0; j < r; j++) {
            sample_result rt = polymer_sample_perm(n,
                    samp->chainlength, samp->bondlength, samp->bondangle,
                    samp->kphi,
                    samp->epsilon, samp->sigma,
                    samp->kbt,
                    exp(-samp->chainlength * 0.75) / 10, exp(-samp->chainlength * 0.75) * 10,
                    rng
                    );
            fe[j] = free_energy(rt.partation, samp->kbt);
            en[j] = rt.energy;
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
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus2);
    gsl_rng_set(rng, time(NULL));
    polymer_sample *ss = read_polymer_sample(stdin);
    run_polymer_sample_perm(ss, rng);
    gsl_rng_free(rng);
    return 0;
}
