CC = g++
CFLAGS= -Wall -pipe -O2 -std=c++11 -ffast-math -funsafe-math-optimizations \
        -march=native -pedantic \
        -I/opt/local/include -I.

LDFLAGS= -L/opt/local/lib -larmadillo

OBJECTS= readMNIST.o main.o 

.PHONY : all clean

all: rateMNIST

rateMNIST : $(OBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OBJECTS) -o rateMNIST
main.o : main.cpp readMNIST.hpp
	$(CC) -c $(CFLAGS) main.cpp
readMNIST.o : readMNIST.cpp readMNIST.hpp
	$(CC) -c $(CFLAGS) readMNIST.cpp
clean: 
	rm -f *.o rateMNIST
