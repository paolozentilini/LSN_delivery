CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o statistic.o buffon_experiment.o
	$(CC) random.o main.o statistic.o buffon_experiment.o -o main.exe

main.o : main.cpp random.h statistic.h experiment.h buffon_experiment.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
buffon_experiment.o : buffon_experiment.cpp buffon_experiment.h random.h experiment.h
		$(CC) -c buffon_experiment.cpp -o buffon_experiment.o $(CFLAGS)
statistic.o : statistic.cpp statistic.h
	$(CC) -c statistic.cpp -o statistic.o $(CFLAGS)
random.o: random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean:
	rm *.o main.exe

esegui:
	./main.exe
