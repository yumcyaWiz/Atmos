all:
	g++ -std=c++17 -fopenmp -g main.cpp

build:
	g++ -std=c++17 -fopenmp -O3 -mtune=native -march=native main.cpp
