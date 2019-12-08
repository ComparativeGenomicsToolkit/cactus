rootPath = ..
include ${rootPath}/include.mk

cflags += ${tokyoCabinetIncl}

libSources = blastAlignmentLib.c
libHeaders = blastAlignmentLib.h

all: all_libs all_progs
all_libs: ${libPath}/cactusBlastAlignment.a
all_progs: 

${libPath}/cactusBlastAlignment.a : ${libSources} ${libHeaders}
	${cxx} ${cflags} -I inc -I ${libPath}/ -c ${libSources}
	ar rc cactusBlastAlignment.a *.o
	ranlib cactusBlastAlignment.a 
	rm *.o
	mv cactusBlastAlignment.a ${libPath}/

clean : 
	rm -f *.o
	rm -f ${libPath}/cactusBlastAlignment.a
