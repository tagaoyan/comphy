#ifndef __XYZ_H__
#define __XYZ_H__
/*
 *  xyz.h provides
 *
 *  chain_xyz: storing molecular chain with xyz coordinate
 *  functions for creating, adding atoms to and destroying chain_xyz
 */


#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "xyzbase.h"

typedef struct {
    size_t length;
    size_t capacity;
    char **atomnames;
    coordinate_xyz *coordinates;
} chain_xyz;

chain_xyz* new_chain_xyz(size_t capacity);
void chain_xyz_add(chain_xyz *chx, char *aname, coordinate_xyz crd);
void free_chain_xyz(chain_xyz *chx);
void read_chain_xyz(chain_xyz *chx, FILE *f);
void print_chain_xyz(chain_xyz *chx, FILE *f);

#endif
