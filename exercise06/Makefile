all: bin/IRS_polymer_MC

include ../Makefile.pre

LDLIBS = -lm -lgsl -lgslcblas

bin/IRS_polymer_MC: obj/zm.o obj/xyzbase.o  obj/xyz.o obj/zmxyz.o obj/randomcoil.o obj/statistics.o obj/polymer_sample.o obj/thermo.o

run: bin/IRS_polymer_MC
	@echo "This command may take 1 min to run..."
	./bin/IRS_polymer_MC < IRS_polymer_MC_in.dat > result/polymer_MC_results.dat

