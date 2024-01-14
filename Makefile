CC = gcc
CFLAGS = -Wall -fopenmp

SRC = convexhull.c
OBJ = convexhull.o

all: convexhull

convexhull: $(OBJ)
	$(CC) $(CFLAGS) -o convexhull $(OBJ)

$(OBJ): $(SRC)
	$(CC) $(CFLAGS) -c $(SRC)

clean:
	rm -f convexhull $(OBJ)

