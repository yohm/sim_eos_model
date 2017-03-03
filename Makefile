CPP=g++
OPT=-O3 -Wall -std=c++11
INCLUDE=

all: eos.out

eos.out: eos.cpp
	$(CPP) $(OPT) -o eos.out eos.cpp

clean:
	rm -f *.out *~ *.bak *.o

