include ../../../include.mk
binPath = ../../../bin
libPath = ../../../lib

libSources = cactusGraph.c 
libHeaders = cactusGraph.h
cflags += ${tokyoCabinetIncl}

all : ${libPath}/cactusGraph.a 

${libPath}/cactusGraph.a : ${libSources} ${libHeaders} ${libPath}/pinchGraph.a ${libPath}/3EdgeConnected.a ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I ${libPath}/ -c ${libSources}
	ar rc cactusGraph.a *.o
	ranlib cactusGraph.a 
	rm *.o
	mv cactusGraph.a ${libPath}/
	cp ${libHeaders} ${libPath}/

clean : 
	rm -f ${libPath}/cactusGraph.a
