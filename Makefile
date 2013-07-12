
CC=gcc
CXX=g++
AR=ar
CFLAGS=-O2 -Wall --std=c99 -g 
CXXFLAGS=-O2 -Wall -g 
LIBS=

main: main.o hmm.o classifier.o
	$(CXX) $^ -o $@ $(CXXFLAGS) $(LIBS)

hmm.o: hmm.cpp hmm.h
classifier.o: classifier.cpp classifier.h
main.o: main.cpp

clean:
	rm -f *.o ./main



