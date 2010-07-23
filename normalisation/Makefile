include ../../../include.mk
binPath = ../../../bin
libPath = ../../../lib

cflags += ${tokyoCabinetIncl}

libSources = main.c normal.c
libHeaders = normal.h

all : ${binPath}/cactus_normalisation

${binPath}/cactus_normalisation : ${libSources} ${libHeaders} ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_normalisation main.c normal.c ${libPath}/cactusLib.a ${libPath}/sonLib.a ${tokyoCabinetLib}

clean : 
	rm -f ${binPath}/cactus_normalisation
