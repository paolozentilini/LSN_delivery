CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o genetic_algorithm.o data_structures.o random.o
	$(CC) random.o genetic_algorithm.o data_structures.o main.o -o main.exe
main.o : main.cpp genetic_algorithm.h data_structures.h random.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
genetic_algorithm.o: genetic_algorithm.cpp genetic_algorithm.h random.h
	$(CC) -c genetic_algorithm.cpp -o genetic_algorithm.o $(CFLAGS)
data_structures.o: data_structures.cpp data_structures.h
	$(CC) -c data_structures.cpp -o data_structures.o $(CFLAGS)

clean:
	rm *.o main.exe *.dat

exe:
	mpiexec -np 4 ./main.exe
