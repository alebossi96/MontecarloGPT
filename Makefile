CC = g++
CFLAGS = -Wall -pedantic -g

main: main.o
	$(CC) main.o -o main
main.o: main.cpp
	$(CC) main.cpp -c $(CFLAGS)
clean:
	rm *.o main
