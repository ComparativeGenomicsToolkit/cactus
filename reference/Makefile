include ../../../include.mk
binPath = ../../../bin
libPath = ../../../lib

cflags += ${tokyoCabinetIncl}

libSources = cactus_reference.c reference.c pseudoAdjacencyBuilder.c balanceTangles.c
libHeaders = reference.h balanceTangles.h

all : ${binPath}/cactus_reference 

${binPath}/cactus_reference : ${libSources} ${libHeaders} ${libPath}/sonLib.a ${libPath}/cactusLib.a 
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_reference ${libSources} ${libPath}/cactusLib.a  ${libPath}/sonLib.a ${tokyoCabinetLib}

clean : 
	rm -f *.o
	rm -f ${binPath}/cactus_reference
