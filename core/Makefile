include ../../../include.mk
binPath = ../../../bin
libPath = ../../../lib

cflags += ${tokyoCabinetIncl}

libSources = impl/*.c
libHeaders = inc/*.h

all : ${libPath}/cactusCore.a ${binPath}/cactus_core

${libPath}/cactusCore.a : ${libSources} ${libHeaders} ${libPath}/pinchGraph.a ${libPath}/3EdgeConnected.a ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I inc -I ${libPath}/ -c ${libSources}
	ar rc cactusCore.a *.o
	ranlib cactusCore.a 
	rm *.o
	mv cactusCore.a ${libPath}/
	cp ${libHeaders} ${libPath}/

${binPath}/cactus_core : cactus_core2.c ${libPath}/cactusCore.a ${libPath}/pinchGraph.a ${libPath}/3EdgeConnected.a ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I inc -I${libPath} -o ${binPath}/cactus_core cactus_core2.c ${libPath}/cactusCore.a ${libPath}/pinchGraph.a ${libPath}/3EdgeConnected.a ${libPath}/sonLib.a ${libPath}/cactusLib.a ${libPath}/sonLib.a ${dblibs}

clean : 
	rm -f *.o
	rm -f ${binPath}/cactus_core ${libPath}/cactusCore.a

