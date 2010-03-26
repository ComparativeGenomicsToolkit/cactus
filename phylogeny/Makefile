include ../../../include.mk
binPath = ../../../bin
libPath = ../../../lib

cflags += ${tokyoCabinetIncl}

all : ${binPath}/cactus_phylogeny 

${binPath}/cactus_phylogeny : cactus_phylogeny.c ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_phylogeny cactus_phylogeny.c ${libPath}/sonLib.a ${libPath}/cactusLib.a ${tokyoCabinetLib}

clean : 
	rm -f *.o
	rm -f ${binPath}/cactus_phylogeny

