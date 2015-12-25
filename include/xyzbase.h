#ifndef __XYZBASE_H__
#define __XYZBASE_H__

/*
 *  xyzbase.h
 *
 *  provides:
 *  coordinate_xyz: storing 3D point or vector
 *  3D vector calculations
 */


#include <math.h>

typedef struct {
    double x;
    double y;
    double z;
} coordinate_xyz;

double inner_prod(coordinate_xyz a, coordinate_xyz b);
double vect_length(coordinate_xyz a);
double vect_length_sq(coordinate_xyz a);
double angle_cos(coordinate_xyz a, coordinate_xyz b);
coordinate_xyz vect_scale(coordinate_xyz a, double t);
coordinate_xyz unit_vect(coordinate_xyz a);
coordinate_xyz vect_plus(coordinate_xyz a, coordinate_xyz b);
coordinate_xyz vect_minus(coordinate_xyz a, coordinate_xyz b);
coordinate_xyz cross_prod(coordinate_xyz a, coordinate_xyz b);

#endif
