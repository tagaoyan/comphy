#include <gsl/gsl_rng.h>
#include <time.h>
#include <math.h>

double sample_uni(double (*f)(double *r), int dim, int n, gsl_rng *rng) {
    double sum = 0;
    for (int i = 0; i < n; i++) {
        double r[dim];
        for (int j = 0; j < dim; j++) {
            r[j] = gsl_rng_uniform(rng);
        }
        sum += f(r);
    }
    return sum / n;
}

double energy(double *r) {
    // a 1-dimensional random number is required
    double phi = r[0] * M_PI * 2 - M_PI;
    double res = 0.5 * (1 + cos(3 * phi));
    res *= exp(-res / 10.);
    return res;
}

int main() {
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, time(NULL));
    int n = 10000000;
    for (int i = 0; i < 10; i++) {
        printf("%g\n", sample_uni(energy, 1, n, rng));
    }
    gsl_rng_free(rng);
    return 0;
}
