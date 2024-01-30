CC = gcc
CFLAGS = -Wextra -Wall -std=c99 -Ofast 
EXECUTABLE = main

all: $(EXECUTABLE)

$(EXECUTABLE): main.o functions.o
	$(CC) $(CFLAGS) -o $(EXECUTABLE) main.o functions.o 

main.o: main.c functions.h
	$(CC) $(CFLAGS) -c main.c 

functions.o: functions.c functions.h
	$(CC) $(CFLAGS) -c functions.c

clean:
	rm -f $(EXECUTABLE)
