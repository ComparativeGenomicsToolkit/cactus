rootPath = ../
include ../include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

commonBarLibs = ${libPath}/stCaf.a ${libPath}/stCactusGraphs.a ${libPath}/stPinchGraphs.a ${libPath}/cactusLib.a  ${libPath}/3EdgeConnected.a
stBarDependencies =  ${commonBarLibs} ${basicLibsDependencies}
stBarLibs = ${commonBarLibs} ${basicLibs}

all : ${binPath}/cactus_baseAligner ${binPath}/cactus_baseAlignerTests 
 
clean : 
	rm -f ${binPath}/cactus_baseAligner ${binPath}/cactus_baseAlignerTests 

${binPath}/cactus_baseAligner : cactus_baseAligner.c ${libTests} ${libSources} ${libHeaders} ${stBarDependencies}
	${cxx} ${cflags} -I inc -I${libPath} -o ${binPath}/cactus_baseAligner cactus_baseAligner.c ${libSources} ${stBarLibs}

${binPath}/cactus_baseAlignerTests : ${libTests} ${libSources} ${libHeaders} tests/*.h ${libPath}/
	${cxx} ${cflags} -I inc -I${libPath} -Wno-error -o ${binPath}/cactus_baseAlignerTests ${libTests} ${libSources} ${stBarLibs}