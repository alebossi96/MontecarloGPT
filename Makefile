CC = g++
CFLAGS = -Wall -Wextra -pedantic -Ofast#-g

main: main.o
	$(CC) main.o montecarlo.o -o main
checkdistribiution: main.o
	$(CC) checkdistribiution.o montecarlo.o -o checkdistribiution
checkdirection: main.o
	$(CC) checkdirection.o montecarlo.o -o checkdirection
python: montecarlo.cpp 
	$(CC) -c -I/usr/include/python3.8/ -I/usr/include/numpy/ -Ofast -Wall -fpic -pedantic -lm montecarlo.cpp montecarlomodule.cpp
	$(CC) -shared -lm -o montecarlomodule.so montecarlo.o montecarlomodule.o

main.o: main.cpp montecarlo.cpp
	$(CC) -c $(CFLAGS) main.cpp montecarlo.cpp

clean:
	rm *.o main *.so
