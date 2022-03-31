CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o gaussian.o uniform.o distribution_1s.o distribution_2p.o statistic.o metropolis.o integral.o
	$(CC) uniform.o main.o gaussian.o distribution_1s.o distribution_2p.o statistic.o metropolis.o integral.o -o main.exe

main.o : main.cpp random.h gaussian.h uniform.h distribution.h distribution_1s.h distribution_2p.h statistic.h metropolis.h integral.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
distribution_1s.o: distribution_1s.cpp distribution_1s.h distribution.h
	$(CC) -c distribution_1s.cpp -o distribution_1s.o $(CFLAGS)
distribution_2p.o: distribution_2p.cpp distribution_2p.h distribution.h
	$(CC) -c distribution_2p.cpp -o distribution_2p.o $(CFLAGS)
statistic.o : statistic.cpp statistic.h integral.h
	$(CC) -c statistic.cpp -o statistic.o $(CFLAGS)
integral.o: integral.cpp integral.h metropolis.h
	$(CC) -c integral.cpp -o integral.o $(CFLAGS)
metropolis.o: metropolis.h metropolis.cpp random.h
	$(CC) -c metropolis.cpp -o metropolis.o $(CFLAGS)
gaussian.o: gaussian.cpp gaussian.h random.h
	$(CC) -c gaussian.cpp -o gaussian.o $(CFLAGS)
uniform.o: uniform.cpp uniform.h random.h
	$(CC) -c uniform.cpp -o uniform.o $(CFLAGS)

clean:
	rm *.o main.exe

esegui:
	./main.exe
