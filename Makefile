COMPILER = g++
CCFLAGS = -Wno-deprecated -lglut -O3 -ffast-math  -march=pentium4

myprogram: main.cpp
	${COMPILER} ${CCFLAGS}  main.cpp MTRand.h -o Ising

