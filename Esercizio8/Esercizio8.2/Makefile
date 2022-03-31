CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o uniform.o distribution_trial.o statistic.o metropolis.o integral.o
	$(CC) uniform.o main.o distribution_trial.o statistic.o metropolis.o integral.o -o main.exe

main.o : main.cpp random.h uniform.h distribution.h distribution_trial.h statistic.h metropolis.h integral.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
distribution_trial.o: distribution_trial.cpp distribution_trial.h distribution.h
	$(CC) -c distribution_trial.cpp -o distribution_trial.o $(CFLAGS)
statistic.o : statistic.cpp statistic.h integral.h
	$(CC) -c statistic.cpp -o statistic.o $(CFLAGS)
integral.o: integral.cpp integral.h metropolis.h
	$(CC) -c integral.cpp -o integral.o $(CFLAGS)
metropolis.o: metropolis.h metropolis.cpp random.h
	$(CC) -c metropolis.cpp -o metropolis.o $(CFLAGS)
uniform.o: uniform.cpp uniform.h random.h
	$(CC) -c uniform.cpp -o uniform.o $(CFLAGS)

clean:
	rm *.o main.exe

esegui:
	./main.exe
