#include "xyz.h"
#include "zm.h"
#include "xyzzm.h"
#include "zmxyz.h"
#include "xyzbase.h"

int main(void) {
    chain *ch = new_chain(1);
    chain_xyz *chx = new_chain_xyz(1);
    read_chain_xyz(chx, stdin);
    chain_zm_from_xyz(ch, chx);
    print_chain(ch, stdout);
    free_chain_xyz(chx);
    free_chain(ch);
    return 0;
}
