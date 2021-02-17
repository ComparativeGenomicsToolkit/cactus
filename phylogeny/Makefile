rootPath = ..
include ${rootPath}/include.mk

treeIncPath = treelib/inc
treeSrc = treelib/src/*.c
CPPFLAGS += -I${treeIncPath}
LDLIBS := treelib/libtree.a ${cactusLibs} ${LDLIBS}

all: all_libs all_progs
all_libs: treelib/libtree.a
all_progs: all_libs


#${MAKE} ${BINDIR}/cactus_phylogeny

${BINDIR}/cactus_phylogeny : cactus_phylogeny.c reconcilliation.c phylogeny.h ${LIBDIR}/cactusLib.a treelib/libtree.a ${treeIncPath}/treelib.h ${LIBDEPENDS}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/cactus_phylogeny cactus_phylogeny.c reconcilliation.c ${LDLIBS}

treelib/libtree.a: ${treeSrc} ${treeIncPath}/*.h
	cd treelib && ${MAKE}

clean : 
	rm -f *.o
	rm -f ${BINDIR}/cactus_phylogeny
	rm -f treelib/libtree.a
