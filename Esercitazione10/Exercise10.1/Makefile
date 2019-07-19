CC = g++
CFLAGS = -Wall -O3 --std=c++11

SimulatedAnnealing_TSP.exe : SimulatedAnnealing_TSP.o random.o
	$(CC) random.o SimulatedAnnealing_TSP.o -o SimulatedAnnealing_TSP.exe
SimulatedAnnealing_TSP.o : SimulatedAnnealing_TSP.cpp
	$(CC) -c SimulatedAnnealing_TSP.cpp -o SimulatedAnnealing_TSP.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o SimulatedAnnealing_TSP.exe
