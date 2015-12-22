
all: bin/xyz2zm

include ../Makefile.pre

LDLIBS = -lm

bin/xyz2zm: obj/xyz2zm.o obj/zm.o obj/xyz.o obj/xyzbase.o obj/xyzzm.o obj/zmxyz.o

run: bin/xyz2zm
	./bin/xyz2zm < ../exercise01/result/xyz.dat > result/zmatrix.dat
