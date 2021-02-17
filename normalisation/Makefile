rootPath = ..
include ${rootPath}/include.mk

libSources = main.c normal.c
libHeaders = normal.h

all: all_libs all_progs
all_libs: 
all_progs: all_libs


#${MAKE} ${BINDIR}/cactus_normalisation

${BINDIR}/cactus_normalisation : ${libSources} ${libHeaders} ${LIBDIR}/cactusLib.a ${LIBDEPENDS}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/cactus_normalisation main.c normal.c ${LIBDIR}/cactusLib.a ${LDLIBS}

clean : 
	rm -f ${BINDIR}/cactus_normalisation
