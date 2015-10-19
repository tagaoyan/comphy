#include <stdio.h>
#include "zm.h"
#include "xyz.h"
#include "zmxyz.h"
#include "xyzbase.h"

int main(void) {
    chain *ch = new_chain(19);
    chain_xyz *chx = new_chain_xyz(19);
    read_chain(ch, stdin);
    chain_xyz_from_zm(chx, ch);
    print_chain_xyz(chx, stdout);
    free_chain(ch);
    free_chain_xyz(chx);
    return 0;
}
