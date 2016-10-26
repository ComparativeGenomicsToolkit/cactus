rootPath = ../
include ../include.mk

all : ${binPath}/cactus_fillAdjacencies

${binPath}/cactus_fillAdjacencies : cactus_fillAdjacencies.c cactus_buildFaces.c ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_fillAdjacencies cactus_fillAdjacencies.c cactus_buildFaces.c ${libPath}/cactusLib.a ${basicLibs}

clean : 
	rm -f *.o
	rm -f  ${binPath}/cactus_fillAdjacencies

