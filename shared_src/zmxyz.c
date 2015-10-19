#include "zmxyz.h"

// first atom, internal use
static coordinate_xyz zm2xyz_1(coordinate_xyz first) {
    return first;
}

// second atom, internal use
static coordinate_xyz zm2xyz_2(coordinate_xyz first, double l) {
    // move along z axis
    first.z += l;
    return first;
}

// internal use
static inline double deg2rad(double deg) {
    return deg * M_PI / 180.;
}

// third atom, internal use
static coordinate_xyz zm2xyz_3(coordinate_xyz second, double l, double alpha) {
    // bend to x axis
    alpha = deg2rad(alpha);
    second.x += l * sin(alpha); // cos(M_PI - alpha)
    second.z += l * -cos(alpha);
    return second;
}

// forth and after, internal use
static coordinate_xyz zm2xyz_n(
        coordinate_xyz a, // three atoms before
        coordinate_xyz b,
        coordinate_xyz c,
        double l, // bond length
        double alpha, // bond angle
        double tau // torsion angle
        ) {
    coordinate_xyz vk = unit_vect(vect_minus(c, b)); // k ~ c - b
    coordinate_xyz ba = vect_minus(a, b);
    coordinate_xyz ve = unit_vect(vect_minus(ba, vect_scale(vk, inner_prod(ba, vk))));
    coordinate_xyz vt = unit_vect(cross_prod(vk, ve));
    alpha = deg2rad(alpha);
    tau = deg2rad(tau);
    coordinate_xyz u = unit_vect(vect_plus(vect_scale(vk, -cos(alpha)), vect_plus(vect_scale(ve, sin(alpha) * cos(tau)), vect_scale(vt, sin(alpha) * sin(tau)))));
    return vect_plus(c, vect_scale(u, l));
}

void chain_xyz_from_zm(chain_xyz *chx, chain *ch) {
    for (int i = 1; i <= ch->length; i++) {
        coordinate_xyz p;
        switch (i) {
            case 1:
                p = zm2xyz_1(ch->begin);
                break;
            case 2:
                p = zm2xyz_2(chx->coordinates[1], ch->bondlengths[2]);
                break;
            case 3:
                p = zm2xyz_3(chx->coordinates[2], ch->bondlengths[3], ch->bondangles[3]);
                break;
            default:
                p = zm2xyz_n(
                        chx->coordinates[i - 3],
                        chx->coordinates[i - 2],
                        chx->coordinates[i - 1],
                        ch->bondlengths[i],
                        ch->bondangles[i],
                        ch->torsionangles[i]
                        );
                break;
        }
        chain_xyz_add(chx, ch->atomnames[i], p);
    }
}
