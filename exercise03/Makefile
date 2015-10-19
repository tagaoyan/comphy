all: bin/makerandomcoil

include ../Makefile.pre

LDLIBS = -lm -lgsl -lgslcblas

bin/makerandomcoil: obj/makerandomcoil.o obj/xyzbase.o obj/zm.o obj/xyz.o obj/zmxyz.o obj/randomcoil.o

