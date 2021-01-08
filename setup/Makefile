rootPath = ..
include ${rootPath}/include.mk

libSources = impl/*.c
libHeaders = inc/*.h

CFLAGS += ${tokyoCabinetIncl} ${hiredisIncl}

all: all_libs all_progs
all_libs: ${LIBDIR}/stCactusSetup.a
all_progs: all_libs
	${MAKE} ${BINDIR}/cactus_setup

${LIBDIR}/stCactusSetup.a : ${libSources} ${libHeaders}
	${CC} ${CPPFLAGS} ${CFLAGS} -I inc -I ${LIBDIR}/ -c ${libSources}
	${AR} rc stCactusSetup.a *.o
	${RANLIB} stCactusSetup.a
	mv stCactusSetup.a ${LIBDIR}/

${BINDIR}/cactus_setup : cactus_setup.c ${LIBDEPENDS} ${LIBDIR}/cactusLib.a ${LIBDIR}/stCactusSetup.a
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/cactus_setup cactus_setup.c ${LIBDIR}/stCactusSetup.a ${LIBDIR}/cactusLib.a ${LDLIBS}

clean : 
	rm -f *.o
	rm -f ${LIBDIR}/stCactusSetup.a ${BINDIR}/cactus_setup
