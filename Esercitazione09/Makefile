CC = g++
CFLAGS = -Wall -O3 --std=c++11

Genetic_TSP.exe : Genetic_TSP.o random.o
	$(CC) random.o Genetic_TSP.o -o Genetic_TSP.exe
Genetic_TSP.o : Genetic_TSP.cpp
	$(CC) -c Genetic_TSP.cpp -o Genetic_TSP.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Genetic_TSP.exe
