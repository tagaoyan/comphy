all: bin/phisample_uni

include ../Makefile.pre

LDLIBS = -lm -lgsl -lgslcblas

bin/phisample_uni: obj/phisample_uni.o obj/statistics.o obj/phisample.o obj/thermo.o

run: bin/phisample_uni
	./bin/phisample_uni < phisample_in.dat > result/uni_sample_results.dat

