rootPath = ../
include ../include.mk

treeIncPath = treelib/inc
treeSrc = treelib/src/*.c

all : ${binPath}/cactus_phylogeny treelib/libtree.a

${binPath}/cactus_phylogeny : cactus_phylogeny.c reconcilliation.c phylogeny.h ${libPath}/cactusLib.a treelib/libtree.a ${treeIncPath}/treelib.h ${basicLibsDependencies}
	${cxx} ${cflags} -I${libPath} -I${treeIncPath} -o ${binPath}/cactus_phylogeny cactus_phylogeny.c reconcilliation.c treelib/libtree.a ${libPath}/cactusLib.a  ${basicLibs}

treelib/libtree.a: ${treeSrc} ${treeIncPath}/*.h
	cd treelib && ${MAKE}

clean : 
	rm -f *.o
	rm -f ${binPath}/cactus_phylogeny
	rm -f treelib/libtree.a
