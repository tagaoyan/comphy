#include "xyzbase.h"

double inner_prod(coordinate_xyz a, coordinate_xyz b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

double vect_length_sq(coordinate_xyz a) {
    return inner_prod(a, a);
}

double vect_length(coordinate_xyz a) {
    return sqrt(inner_prod(a, a));
}

double angle_cos(coordinate_xyz a, coordinate_xyz b) {
    return inner_prod(a, b) / (vect_length(a) * vect_length(b));
}

coordinate_xyz vect_scale(coordinate_xyz a, double t) {
    coordinate_xyz p;
    p.x = a.x * t;
    p.y = a.y * t;
    p.z = a.z * t;
    return p;
}

coordinate_xyz unit_vect(coordinate_xyz a) {
    return vect_scale(a, 1 / vect_length(a));
}

coordinate_xyz vect_plus(coordinate_xyz a, coordinate_xyz b) {
    coordinate_xyz p;
    p.x = a.x + b.x;
    p.y = a.y + b.y;
    p.z = a.z + b.z;
    return p;
}

coordinate_xyz vect_minus(coordinate_xyz a, coordinate_xyz b) {
    coordinate_xyz p;
    p.x = a.x - b.x;
    p.y = a.y - b.y;
    p.z = a.z - b.z;
    return p;
}

coordinate_xyz cross_prod(coordinate_xyz a, coordinate_xyz b) {
    coordinate_xyz p;
    p.x = a.y * b.z - a.z * b.y;
    p.y = a.z * b.x - a.x * b.z;
    p.z = a.x * b.y - a.y * b.x;
    return p;
}
