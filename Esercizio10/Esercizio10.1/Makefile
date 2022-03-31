CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o data_structures.o simulated_annealing.o
	$(CC) random.o main.o data_structures.o simulated_annealing.o -o main.exe

main.o : main.cpp random.h data_structures.h simulated_annealing.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
simulated_annealing.o: simulated_annealing.h simulated_annealing.cpp random.h data_structures.h
	$(CC) -c simulated_annealing.cpp -o simulated_annealing.o $(CFLAGS)
data_structures.o: data_structures.h data_structures.cpp random.h
	$(CC) -c data_structures.cpp -o data_structures.o $(CFLAGS)
random.o: random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean:
	rm *.o main.exe

esegui:
	./main.exe
