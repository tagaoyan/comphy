#ifndef __ZMXYZ_H__
#define __ZMXYZ_H__
#include "zm.h"
#include "xyz.h"
#include "xyzbase.h"
#include <math.h>
/*
 *  zmxyz.h provides
 *
 *  convert from zmatrix to xyz coordinate
 */

void chain_xyz_from_zm(chain_xyz *chx, chain *ch);
#endif
