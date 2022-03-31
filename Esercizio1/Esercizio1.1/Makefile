CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o statistic.o integral_1.o integral.o
	$(CC) random.o main.o statistic.o integral_1.o integral.o -o main.exe

main.o : main.cpp random.h statistic.h integral.h experiment.h integral_1.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
integral.o : integral.cpp integral.h experiment.h
		$(CC) -c integral.cpp -o integral.o $(CFLAGS)
integral_1.o : integral_1.cpp integral_1.h experiment.h
				$(CC) -c integral_1.cpp -o integral_1.o $(CFLAGS)
statistic.o : statistic.cpp statistic.h
	$(CC) -c statistic.cpp -o statistic.o $(CFLAGS)
random.o: random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean:
	rm *.o main.exe

esegui:
	./main.exe
