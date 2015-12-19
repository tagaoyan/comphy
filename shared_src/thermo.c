#include "thermo.h"
#include "math.h"

double energy_function(double kphi, double phi) {
    return 0.5 * kphi * (1 + cos(3 * phi));
}

double free_energy(double z, double kbt) {
    // calculate free energy from partation function
    return -kbt * log(z);
}

double entrophy(double u, double fe, double kbt) {
    // calculate entrophy from internal energy and free energy
    return (u - fe) / kbt;
}
