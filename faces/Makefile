rootPath = ..
include ${rootPath}/include.mk

all: all_libs all_progs
all_libs: 
all_progs: all_libs
	${MAKE} ${BINDIR}/cactus_fillAdjacencies

${BINDIR}/cactus_fillAdjacencies : cactus_fillAdjacencies.c cactus_buildFaces.c ${LIBDIR}/cactusLib.a ${LIBDEPENDS}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/cactus_fillAdjacencies cactus_fillAdjacencies.c cactus_buildFaces.c ${LIBDIR}/cactusLib.a ${LDLIBS}

clean : 
	rm -f *.o
	rm -f  ${BINDIR}/cactus_fillAdjacencies

