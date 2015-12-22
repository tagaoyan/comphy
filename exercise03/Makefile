all: bin/makerandomcoil

include ../Makefile.pre

LDLIBS = -lm -lgsl -lgslcblas

bin/makerandomcoil: obj/makerandomcoil.o obj/xyzbase.o obj/zm.o obj/xyz.o obj/zmxyz.o obj/randomcoil.o

run: bin/makerandomcoil
	./bin/makerandomcoil < randomcoil_in.dat > result/rc_zmatrix.dat 2> result/rc_xyz.dat
