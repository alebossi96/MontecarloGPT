CC = g++
CFLAGS = -Wall -Wextra -pedantic -Ofast#-g

main: main.o
	$(CC) main.o montecarlo.o -o main
checkdistribiution: main.o
	$(CC) checkdistribiution.o montecarlo.o -o checkdistribiution
checkdirection: main.o
	$(CC) checkdirection.o montecarlo.o -o checkdirection
	
main.o: main.cpp montecarlo.cpp
	$(CC) -c $(CFLAGS) main.cpp montecarlo.cpp
clean:
	rm *.o main
