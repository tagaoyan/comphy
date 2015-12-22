#ifndef __ZM_H__
#define __ZM_H__

#include <stdio.h>
#include <stdlib.h>
#include "xyz.h"

typedef struct {
    // chain length
    size_t length;
    // max length of the chain
    size_t capacity;
    // begining of the chain in xyz coordinate
    coordinate_xyz begin;
    // array of atom names, the 0th element is skipped, making the first atom has the index 1
    char **atomnames;
    // array of bond lengths, index n means bond between n-th and (n-1)-th atom
    double *bondlengths;
    // array of bond angles, index n means bond angle of n-th and (n-2)-th atom
    double *bondangles;
    // array of torsion angles, index n means torsion angle of n-th and (n-3)-th atom
    double *torsionangles;
} chain;

typedef chain chain_zm;

chain *new_chain(size_t capacity);
void free_chain(chain *ch);
void chain_add(chain *ch, char *aname, double blength, double bangle, double tangle);
void chain_del(chain *ch);
void read_chain(chain *ch, FILE *f);
void print_chain(chain *ch, FILE *f);

#endif
