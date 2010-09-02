include ../../../include.mk
binPath = ../../../bin
libPath = ../../../lib

cflags += ${tokyoCabinetIncl}

all : ${binPath}/cactus_fillAdjacencies

${binPath}/cactus_fillAdjacencies : cactus_fillAdjacencies.c cactus_buildFaces.c ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_fillAdjacencies cactus_fillAdjacencies.c cactus_buildFaces.c ${libPath}/cactusLib.a ${libPath}/sonLib.a ${dblibs}

clean : 
	rm -f *.o
	rm -f  ${binPath}/cactus_fillAdjacencies

