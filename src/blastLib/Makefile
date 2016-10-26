rootPath = ../
include ../include.mk

cflags += ${tokyoCabinetIncl}

libSources = blastAlignmentLib.c
libHeaders = blastAlignmentLib.h

all : ${libPath}/cactusBlastAlignment.a

${libPath}/cactusBlastAlignment.a : ${libSources} ${libHeaders} ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I inc -I ${libPath}/ -c ${libSources}
	ar rc cactusBlastAlignment.a *.o
	ranlib cactusBlastAlignment.a 
	rm *.o
	mv cactusBlastAlignment.a ${libPath}/
	cp ${libHeaders} ${libPath}/

clean : 
	rm -f *.o
	rm -f ${libPath}/cactusBlastAlignment.a