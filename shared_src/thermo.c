#include "thermo.h"
#include "math.h"


double free_energy(double z, double kbt) {
    // calculate free energy from partation function
    return -kbt * log(z);
}

double entrophy(double u, double fe, double kbt) {
    // calculate entrophy from internal energy and free energy
    return (u - fe) / kbt;
}
