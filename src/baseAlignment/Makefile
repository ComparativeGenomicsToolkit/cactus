rootPath = ../../
include ../../include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

all : ${binPath}/cactus_baseAligner ${binPath}/cactus_baseAlignerTests 
 
clean : 
	rm -f ${binPath}/cactus_baseAligner ${binPath}/cactus_baseAlignerTests 

${binPath}/cactus_baseAligner : cactus_baseAligner.c ${libTests} ${libSources} ${libHeaders} ${libPath}/cactusLib.a ${libPath}/cactusCore.a ${libPath}/pinchGraph.a ${libPath}/3EdgeConnected.a ${basicLibsDependencies}
	${cxx} ${cflags} -I inc -I${libPath} -o ${binPath}/cactus_baseAligner cactus_baseAligner.c ${libSources} ${libPath}/cactusCore.a ${libPath}/pinchGraph.a ${libPath}/3EdgeConnected.a ${libPath}/cactusLib.a ${basicLibs}

${binPath}/cactus_baseAlignerTests : ${libTests} ${libSources} ${libHeaders} tests/*.h ${libPath}/cactusLib.a ${libPath}/cactusCore.a ${libPath}/pinchGraph.a ${libPath}/3EdgeConnected.a ${basicLibsDependencies}
	${cxx} ${tokyoCabinetIncl} -Wall --pedantic -g -O3 -ffast-math -DBEN_DEBUG -lm -I ${sonLibPath} -I inc -I${libPath}  -o ${binPath}/cactus_baseAlignerTests ${libTests} ${libSources} ${libPath}/cactusLib.a ${basicLibs}
#${cflags}