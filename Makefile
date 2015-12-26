.PHONY: all ex01 ex02 ex03 ex04 ex05 ex06 ex07 run

all: ex01 ex02 ex03 ex04 ex05 ex06 ex07

run: ex01 ex02 ex03 ex04 ex05 ex06 ex07
	cd exercise01 && make run
	cd exercise02 && make run
	cd exercise03 && make run
	cd exercise04 && make run
	cd exercise05 && make run
	cd exercise06 && make run
	cd exercise07 && make run

ex01:
	cd exercise01 && make

ex02:
	cd exercise02 && make

ex03:
	cd exercise03 && make

ex04:
	cd exercise04 && make

ex05:
	cd exercise05 && make

ex06:
	cd exercise06 && make

ex07:
	cd exercise07 && make
