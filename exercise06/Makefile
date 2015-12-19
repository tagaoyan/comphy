all: bin/IRS_polymer_MC

include ../Makefile.pre

LDLIBS = -lm -lgsl -lgslcblas

bin/IRS_polymer_MC: obj/zm.o obj/xyzbase.o  obj/xyz.o obj/zmxyz.o obj/randomcoil.o obj/statistics.o obj/phisample.o obj/thermo.o

