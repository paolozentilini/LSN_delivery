CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11
LIBS =

main.exe: main.o random.o
	$(CC) $^ -o $@ $(LIBS)

%.o: %.cpp %.h
	$(CC) -c $< $(CFLAGS)

main.o : main.cpp random.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)

esegui: main.exe
	mpiexec -np 2 ./main.exe file.txt

clean :
	rm *.o main.exe
