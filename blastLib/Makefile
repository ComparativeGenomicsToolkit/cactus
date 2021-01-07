rootPath = ..
include ${rootPath}/include.mk

CFLAGS += ${tokyoCabinetIncl} ${hiredisIncl}

libSources = blastAlignmentLib.c
libHeaders = blastAlignmentLib.h

all: all_libs all_progs
all_libs: ${LIBDIR}/cactusBlastAlignment.a
all_progs: 

${LIBDIR}/cactusBlastAlignment.a : ${libSources} ${libHeaders}
	${CC} ${CPPFLAGS} ${CFLAGS} -I inc -I ${LIBDIR}/ -c ${libSources}
	${AR} rc cactusBlastAlignment.a *.o
	${RANLIB} cactusBlastAlignment.a 
	mv cactusBlastAlignment.a ${LIBDIR}/

clean : 
	rm -f *.o
	rm -f ${LIBDIR}/cactusBlastAlignment.a
