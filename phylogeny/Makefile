include ../../../include.mk
binPath = ../../../bin
libPath = ../../../lib

cflags += ${tokyoCabinetIncl}

treeIncPath = treelib/inc

all : ${binPath}/cactus_phylogeny treelib/libtree.a

${binPath}/cactus_phylogeny : cactus_phylogeny.c reconcilliation.c phylogeny.h ${libPath}/sonLib.a ${libPath}/cactusLib.a treelib/libtree.a ${treeIncPath}/treelib.h
	${cxx} ${cflags} -I${libPath} -I${treeIncPath} -o ${binPath}/cactus_phylogeny cactus_phylogeny.c reconcilliation.c ${libPath}/sonLib.a ${libPath}/cactusLib.a treelib/libtree.a ${tokyoCabinetLib}

treelib/libtree.a:
	cd treelib && make

clean : 
	rm -f *.o
	rm -f ${binPath}/cactus_phylogeny
	rm treelib/libtree.a
