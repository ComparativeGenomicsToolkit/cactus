include ../../../include.mk
binPath = ../../../bin
libPath = ../../../lib

cflags += ${tokyoCabinetIncl}

libSources = cactus_core.c cactusFlowerFunctions.c 
libHeaders = cactus_core.h cactusFlowerFunctions.h

all : ${libPath}/cactusCore.a ${binPath}/cactus_core

${libPath}/cactusCore.a : ${libSources} ${libHeaders} ${libPath}/cactusGraph.a ${libPath}/pinchGraph.a ${libPath}/3EdgeConnected.a ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I ${libPath}/ -c ${libSources}
	ar rc cactusCore.a *.o
	ranlib cactusCore.a 
	rm *.o
	mv cactusCore.a ${libPath}/
	cp ${libHeaders} ${libPath}/

${binPath}/cactus_core : cactus_core2.c ${libPath}/cactusCore.a ${libPath}/cactusGraph.a ${libPath}/pinchGraph.a ${libPath}/3EdgeConnected.a ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_core cactus_core2.c ${libPath}/cactusCore.a ${libPath}/cactusGraph.a ${libPath}/pinchGraph.a ${libPath}/3EdgeConnected.a ${libPath}/sonLib.a ${libPath}/cactusLib.a ${libPath}/sonLib.a ${dblibs}

clean : 
	rm -f *.o
	rm -f ${binPath}/cactus_core ${libPath}/cactusCore.a

