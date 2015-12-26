all: bin/PERM_polymer_MC

include ../Makefile.pre

LDLIBS = -lm -lgsl -lgslcblas

bin/PERM_polymer_MC: obj/zm.o obj/xyzbase.o  obj/xyz.o obj/zmxyz.o obj/randomcoil.o obj/statistics.o obj/phisample.o obj/thermo.o obj/polymer_sample.o

run: bin/PERM_polymer_MC
	@echo "This command may take over 3 min to run..."
	./bin/PERM_polymer_MC < polymer_MC_in.dat > result/PERM_polymer_MC_results.dat
