CC = g++
CFLAGS = -Wall -Wextra -pedantic -g

main: main.o
	$(CC) main.o -o main
main.o: main.cpp
	$(CC) -c $(CFLAGS) main.cpp 
clean:
	rm *.o main
