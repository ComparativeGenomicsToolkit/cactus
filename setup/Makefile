rootPath = ..
include ${rootPath}/include.mk

libSources = impl/*.c
libHeaders = inc/*.h

CFLAGS += ${hiredisIncl}

all: all_libs
all_progs : all_libs
all_libs: ${LIBDIR}/stCactusSetup.a

${LIBDIR}/stCactusSetup.a : ${libSources} ${libHeaders}
	${CC} ${CPPFLAGS} ${CFLAGS} -I inc -I ${LIBDIR}/ -c ${libSources}
	${AR} rc stCactusSetup.a *.o
	${RANLIB} stCactusSetup.a
	mv stCactusSetup.a ${LIBDIR}/

clean : 
	rm -f *.o
	rm -f ${LIBDIR}/stCactusSetup.a
