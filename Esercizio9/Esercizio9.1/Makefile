CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o data_structures.o genetic_algorithm.o
	$(CC) random.o main.o data_structures.o genetic_algorithm.o -o main.exe

main.o : main.cpp random.h data_structures.h genetic_algorithm.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
genetic_algorithm.o: genetic_algorithm.h genetic_algorithm.cpp random.h data_structures.h
	$(CC) -c genetic_algorithm.cpp -o genetic_algorithm.o $(CFLAGS)
data_structures.o: data_structures.h data_structures.cpp random.h
	$(CC) -c data_structures.cpp -o data_structures.o $(CFLAGS)
random.o: random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean:
	rm *.o main.exe

esegui:
	./main.exe
