#include "xyzzm.h"

static double bondlength(coordinate_xyz a, coordinate_xyz b) {
    return vect_length(vect_minus(b, a));
}

static inline double rad2deg(double rad) {
    return rad * 180. / M_PI;
}

static double bondangle(coordinate_xyz a, coordinate_xyz m, coordinate_xyz b) {
    return rad2deg(acos(angle_cos(vect_minus(a, m), vect_minus(b, m))));
}

static double torsionangle(
        coordinate_xyz a,
        coordinate_xyz m,
        coordinate_xyz n,
        coordinate_xyz b) {
    coordinate_xyz ma = vect_minus(a, m);
    coordinate_xyz nb = vect_minus(b, n);
    coordinate_xyz mn = vect_minus(n, m);
    coordinate_xyz e1 = cross_prod(ma, mn);
    coordinate_xyz e2 = cross_prod(nb, mn);
    double tau = rad2deg(acos(angle_cos(e1, e2)));
    if (isnan(tau)) {
        // due to machine error in float calculation
        // cos value may be greater than 1
        tau = 180.;
    }
    if (inner_prod(cross_prod(e1, e2), mn) < 0) {
        tau *= -1;
    }
    return tau;
}

void chain_zm_from_xyz(chain *ch, chain_xyz *chx) {
    for (int i = 1; i <= chx->length; i++) {
        double blength = 0, bangle = 0, tangle = 0;
        if (i == 1) {
            ch->begin = chx->coordinates[i];
        }
        if (i >= 2) {
            blength = bondlength(
                    chx->coordinates[i - 1],
                    chx->coordinates[i]
                    );
        }
        if (i >= 3) {
            bangle = bondangle(
                    chx->coordinates[i - 2],
                    chx->coordinates[i - 1],
                    chx->coordinates[i]
                    );
        }
        if (i >= 4) {
            tangle = torsionangle(
                    chx->coordinates[i - 3],
                    chx->coordinates[i - 2],
                    chx->coordinates[i - 1],
                    chx->coordinates[i]
                    );
        }
        chain_add(ch, chx->atomnames[i], blength, bangle, tangle);
    }
}
