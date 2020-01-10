all: test4.o
	g++ --std=c++11 -O0 -g test4.o -o test4 -lfst

test4.o: test4.cpp
	g++ --std=c++11 -O0 -g test4.cpp -c -lfst
