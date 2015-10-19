#include "zm.h"

chain *new_chain(size_t capacity) {
    chain *ch = malloc(sizeof(chain));
    ch->atomnames = malloc(sizeof(char *) * (capacity + 1));
    ch->bondlengths = malloc(sizeof(double) * (capacity + 1));
    ch->bondangles = malloc(sizeof(double) * (capacity + 1));
    ch->torsionangles = malloc(sizeof(double) * (capacity + 1));
    ch->capacity = capacity;
    ch->length = 0;
    return ch;
}

void free_chain(chain *ch) {
    for (int i = 0; i <= ch->capacity; i++) {
        free(ch->atomnames[i]);
    }
    free(ch->atomnames);
    free(ch->bondlengths);
    free(ch->bondangles);
    free(ch->torsionangles);
    free(ch);
    // manually set ch to NULL if necessary
}

void chain_add(chain *ch, char *aname, double blength, double bangle, double tangle) {
    // add an atom to the end of the chain,
    // some parameters will be ignored if the chain is too short
    int index = ch->length + 1;
    if (index > ch->capacity) {
        // double the capacity
        ch->capacity *= 2;
        ch->atomnames = realloc(ch->atomnames, sizeof(char *) * (ch->capacity + 1));
        ch->bondlengths = realloc(ch->bondlengths, sizeof(double) * (ch->capacity + 1));
        ch->bondangles = realloc(ch->bondangles, sizeof(double) * (ch->capacity + 1));
        ch->torsionangles = realloc(ch->torsionangles, sizeof(double) * (ch->capacity + 1));
    }
    ch->atomnames[index] = strdup(aname);
    ch->bondlengths[index] = blength;
    ch->bondangles[index] = bangle;
    ch->torsionangles[index] = tangle;
    ch->length++;
}

#define BUF_SIZE 4096
#define ANAME_SIZE 8
void read_chain(chain *ch, FILE *f) {
    char buf[BUF_SIZE];
    while (fgets(buf, BUF_SIZE, f) != NULL) {
        if (buf[0] == '#') {
            continue;
        }
        int id;
        sscanf(buf, "%d", &id);
        if (id > 1) {
            char aname[ANAME_SIZE];
            double blength, bangle, tangle;
            sscanf(buf, "%*d %s %*d %lf %*d %lf %*d %lf", aname, &blength, &bangle, &tangle);
            chain_add(ch, aname, blength, bangle, tangle);
        } else {
            char *aname = malloc(ANAME_SIZE * sizeof(char *));
            double x, y, z;
            sscanf(buf, "%*d %s %lf %lf %lf", aname, &x, &y, &z);
            chain_add(ch, aname, 0, 0, 0);
            ch->begin.x = x;
            ch->begin.y = y;
            ch->begin.z = z;
        }
    }
}

void print_chain(chain *ch, FILE *f) {
    for (int i = 1; i <= ch->length; i++) {
        switch(i) {
            case 1:
                fprintf(f, "%d %s\t%f %f %f\n",
                        1, ch->atomnames[1],
                        ch->begin.x, ch->begin.y, ch->begin.z
                        );
                break;
            case 2:
                fprintf(f, "%d %s\t%d %f\n",
                        i, ch->atomnames[i],
                        i - 1, ch->bondlengths[i]
                        );
                break;
            case 3:
                fprintf(f, "%d %s\t%d %f\t%d %f\n",
                        i, ch->atomnames[i],
                        i - 1, ch->bondlengths[i],
                        i - 2, ch->bondangles[i]
                        );
                break;
            default:
                fprintf(f, "%d %s\t%d %f\t%d %f\t%d %f\n",
                        i, ch->atomnames[i],
                        i - 1, ch->bondlengths[i],
                        i - 2, ch->bondangles[i],
                        i - 3, ch->torsionangles[i]
                        );
                break;
        }
    }
}
