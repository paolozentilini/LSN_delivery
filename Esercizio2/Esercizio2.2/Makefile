CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o statistic.o continous_random_walk.o discrete_random_walk.o
	$(CC) random.o main.o statistic.o continous_random_walk.o discrete_random_walk.o -o main.exe

main.o : main.cpp random.h statistic.h discrete_random_walk.h continous_random_walk.h random_walk.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
continous_random_walk.o: continous_random_walk.h continous_random_walk.cpp random_walk.h random.h
	$(CC) -c continous_random_walk.cpp -o continous_random_walk.o
discrete_random_walk.o: discrete_random_walk.h discrete_random_walk.cpp random_walk.h random.h
	$(CC) -c discrete_random_walk.cpp -o discrete_random_walk.o
statistic.o : statistic.cpp statistic.h random_walk.h
	$(CC) -c statistic.cpp -o statistic.o $(CFLAGS)
random.o: random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean:
	rm *.o main.exe

esegui:
	./main.exe
