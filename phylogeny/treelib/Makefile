SRC = src
INC = inc

CC = gcc
CFLAGS = -I$(INC)

STRINGOBJ = strings.o lists.o global.o
NJTREEOBJ = align.o buildtree.o cluster.o distancemat.o sequence.o tree.o util.o
TREEOBJ = treelib.o

all: $(STRINGOBJ) $(NJTREEOBJ) $(TREEOBJ)
	ar rc libtree.a $(STRINGOBJ) $(NJTREEOBJ) $(TREEOBJ)
	ranlib libtree.a
	rm $(STRINGOBJ) $(NJTREEOBJ) $(TREEOBJ)

global.o : $(SRC)/global.c $(INC)/global.h
	$(CC) $(CFLAGS) -c $(SRC)/global.c

lists.o : $(SRC)/lists.c $(INC)/lists.h
	$(CC) $(CFLAGS) -c $(SRC)/lists.c

strings.o : $(SRC)/strings.c $(INC)/strings.h
	$(CC) $(CFLAGS) -c $(SRC)/strings.c 

align.o : $(SRC)/align.c $(INC)/align.h
	$(CC) $(CFLAGS) -c $(SRC)/align.c

buildtree.o : $(SRC)/buildtree.c $(INC)/buildtree.h
	$(CC) $(CFLAGS) -c $(SRC)/buildtree.c

cluster.o : $(SRC)/cluster.c $(INC)/cluster.h
	$(CC) $(CFLAGS) -c $(SRC)/cluster.c

distancemat.o : $(SRC)/distancemat.c $(INC)/distancemat.h
	$(CC) $(CFLAGS) -c $(SRC)/distancemat.c

sequence.o : $(SRC)/sequence.c $(INC)/sequence.h
	$(CC) $(CFLAGS) -c $(SRC)/sequence.c

tree.o : $(SRC)/tree.c $(INC)/tree.h
	$(CC) $(CFLAGS) -c $(SRC)/tree.c

util.o : $(SRC)/util.c $(INC)/util.h
	$(CC) $(CFLAGS) -c $(SRC)/util.c

treelib.o : $(SRC)/treelib.c $(INC)/treelib.h
	$(CC) $(CFLAGS) -c $(SRC)/treelib.c

clean:
	rm libtree.a
