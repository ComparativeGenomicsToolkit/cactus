include ../../../include.mk
binPath = ../../../bin
libPath = ../../../lib

cflags += ${tokyoCabinetIncl}

all : ${binPath}/cactus_reference 

${binPath}/cactus_reference : cactus_reference.c reference.c pseudoAdjacencyBuilder.c reference.h ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_reference cactus_reference.c reference.c pseudoAdjacencyBuilder.c ${libPath}/cactusLib.a  ${libPath}/sonLib.a ${tokyoCabinetLib}

clean : 
	rm -f *.o
	rm -f ${binPath}/cactus_reference
