include ../../../include.mk
binPath = ../../../bin
libPath = ../../../lib

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

all : ${binPath}/cactus_baseAligner ${binPath}/cactus_baseAlignerTests 
 
clean : 
	rm -f ${binPath}/cactus_baseAligner ${binPath}/cactus_baseAlignerTests 

${binPath}/cactus_baseAligner : cactus_baseAligner.c ${libTests} ${libSources} ${libHeaders} ${libPath}/cactusLib.a ${libPath}/sonLib.a ${libPath}/cuTest.a ${libPath}/cactusCore.a ${libPath}/pinchGraph.a ${libPath}/3EdgeConnected.a
	${cxx} ${cflags} -I inc -I${libPath} -o ${binPath}/cactus_baseAligner cactus_baseAligner.c ${libSources} ${libPath}/cactusCore.a ${libPath}/pinchGraph.a ${libPath}/3EdgeConnected.a ${libPath}/cactusLib.a ${libPath}/sonLib.a ${dblibs}
	
${binPath}/cactus_baseAlignerTests : ${libTests} ${libSources} ${libHeaders} tests/*.h ${libPath}/cactusLib.a ${libPath}/cactusLib.a ${libPath}/sonLib.a ${libPath}/cuTest.a ${libPath}/cactusCore.a ${libPath}/pinchGraph.a ${libPath}/3EdgeConnected.a
	${cxx} ${cflags} -I inc -I${libPath}  -o ${binPath}/cactus_baseAlignerTests ${libTests} ${libSources} ${libPath}/cactusLib.a ${libPath}/sonLib.a ${libPath}/cuTest.a ${dblibs}
