SRC = src
INC = inc

CC = gcc
CFLAGS = -I$(INC) -ggdb -Wall -DNDEBUG

SOURCES   = *.c
OBJ = *.o

all: libtree.a

libtree.a: $(SRC)/$(SOURCES)
	$(CC) $(CFLAGS) -c $(SRC)/$(SOURCES)
	ar rc libtree.a $(OBJ)
	ranlib libtree.a
	rm -f $(OBJ)

clean:
	rm -f libtree.a
