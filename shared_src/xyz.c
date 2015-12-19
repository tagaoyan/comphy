#include "xyz.h"

chain_xyz* new_chain_xyz(size_t capacity) {
    // create a new chain in xyz coordinates
    chain_xyz *chx = malloc(sizeof(chain_xyz));
    chx->capacity = capacity;
    chx->length = 0;
    chx->atomnames = malloc(sizeof(char *) * capacity);
    chx->coordinates = malloc(sizeof(coordinate_xyz) * capacity);
    return chx;
}

void chain_xyz_add(chain_xyz *chx, char *aname, coordinate_xyz crd) {
    int index = chx->length;
    if (index >= chx->capacity) {
        // double the capacity
        chx->capacity *= 2;
        chx->atomnames = realloc(chx->atomnames, sizeof(char *) * chx->capacity);
        chx->coordinates = realloc(chx->coordinates, sizeof(coordinate_xyz) * chx->capacity);
    }
    chx->atomnames[index] = strdup(aname);
    chx->coordinates[index] = crd;
    chx->length++;
}

void free_chain_xyz(chain_xyz *chx) {
    for (int i = 0; i < chx->capacity; i++) {
        free(chx->atomnames[i]);
    }
    free(chx->coordinates);
    free(chx->atomnames);
    free(chx);
}

#define BUF_SIZE 4096
#define ANAME_SIZE 8
void read_chain_xyz(chain_xyz *chx, FILE *f) {
    char buf[BUF_SIZE];
    while (fgets(buf, BUF_SIZE, f) != NULL) {
        if (buf[0] == '#') {
            continue;
        }
        int id;
        char aname[ANAME_SIZE];
        coordinate_xyz p;
        sscanf(buf, "%d %s %lf %lf %lf", &id, aname, &p.x, &p.y, &p.z);
        chain_xyz_add(chx, aname, p);
    }
}

void print_chain_xyz(chain_xyz *chx, FILE *f) {
    for (int i = 0; i < chx->length; i++) {
        coordinate_xyz p = chx->coordinates[i];
        fprintf(f, "%d %s\t%.12f %.12f %.12f\n", i + 1, chx->atomnames[i], p.x, p.y, p.z);
    }
}
