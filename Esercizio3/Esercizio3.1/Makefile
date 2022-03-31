CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o statistic.o continous_wiener_process.o discrete_wiener_process.o plain_vanilla.o
	$(CC) random.o main.o statistic.o continous_wiener_process.o discrete_wiener_process.o plain_vanilla.o -o main.exe

main.o : main.cpp random.h statistic.h discrete_wiener_process.h continous_wiener_process.h wiener_process.h plain_vanilla.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
plain_vanilla.o: plain_vanilla.cpp plain_vanilla.h wiener_process.h
	$(CC) -c plain_vanilla.cpp -o plain_vanilla.o $(CFLAGS)
continous_wiener_process.o: continous_wiener_process.h continous_wiener_process.cpp wiener_process.h random.h
	$(CC) -c continous_wiener_process.cpp -o continous_wiener_process.o $(CFLAGS)
discrete_wiener_process.o: discrete_wiener_process.h discrete_wiener_process.cpp wiener_process.h random.h
	$(CC) -c discrete_wiener_process.cpp -o discrete_wiener_process.o $(CFLAGS)
statistic.o : statistic.cpp statistic.h plain_vanilla.h
	$(CC) -c statistic.cpp -o statistic.o $(CFLAGS)
random.o: random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean:
	rm *.o main.exe

esegui:
	./main.exe
