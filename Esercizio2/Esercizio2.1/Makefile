CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o statistic.o coseno.o uniform_sampling.o importance_sampling.o g_function.o f_function.o
	$(CC) random.o main.o statistic.o coseno.o uniform_sampling.o importance_sampling.o g_function.o f_function.o -o main.exe

main.o : main.cpp random.h statistic.h montecarlo.h funzionebase.h uniform_sampling.h importance_sampling.h g_function.h f_function.h 
	$(CC) -c main.cpp -o main.o $(CFLAGS)
uniform_sampling.o: uniform_sampling.h uniform_sampling.cpp montecarlo.h random.h
	$(CC) -c uniform_sampling.cpp -o uniform_sampling.o
importance_sampling.o: importance_sampling.h importance_sampling.cpp montecarlo.h random.h
	$(CC) -c importance_sampling.cpp -o importance_sampling.o
coseno.o : coseno.cpp coseno.h funzionebase.h
	$(CC) -c coseno.cpp -o coseno.o $(CFLAGS)
g_function.o: g_function.cpp g_function.h funzionebase.h
	$(CC) -c g_function.cpp -o g_function.o $(CFLAGS)
f_function.o: f_function.cpp f_function.h funzionebase.h
	$(CC) -c f_function.cpp -o f_function.o $(CFLAGS)
statistic.o : statistic.cpp statistic.h
	$(CC) -c statistic.cpp -o statistic.o $(CFLAGS)
random.o: random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean:
	rm *.o main.exe

esegui:
	./main.exe
