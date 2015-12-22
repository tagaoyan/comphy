all: bin/phisample_bz

include ../Makefile.pre

LDLIBS = -lm -lgsl -lgslcblas

bin/phisample_bz: obj/phisample_bz.o obj/statistics.o obj/phisample.o obj/thermo.o

run: bin/phisample_bz
	./bin/phisample_bz < phisample_in.dat > result/bz_sample_results.dat

