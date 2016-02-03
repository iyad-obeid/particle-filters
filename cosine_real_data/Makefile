myFuncs: myFuncs.h myFuncs.cpp
	g++ -c myFuncs.cpp

particle: particle.h particle.cpp
	g++ -c particle.cpp
	
main: main.cpp
	g++ -c main.cpp
	
all: particle.o main.o myFuncs.o
	g++ main.o myFuncs.o particle.o -o main
	