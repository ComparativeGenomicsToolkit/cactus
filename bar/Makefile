rootPath = ../
include ../include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

commonBarLibs = ${libPath}/stCaf.a ${sonLibPath}/stPinchesAndCacti.a ${libPath}/cactusLib.a  ${sonLibPath}/3EdgeConnected.a
stBarDependencies =  ${commonBarLibs} ${basicLibsDependencies}
stBarLibs = ${commonBarLibs} ${basicLibs}

all : ${binPath}/cactus_bar ${binPath}/cactus_barTests 
 
clean : 
	rm -f ${binPath}/cactus_bar ${binPath}/cactus_barTests 

${binPath}/cactus_bar : cactus_bar.c ${libSources} ${libHeaders} ${stBarDependencies}
	${cxx} ${cflags} -I inc -I${libPath} -o ${binPath}/cactus_bar cactus_bar.c ${libSources} ${stBarLibs}

${binPath}/cactus_barTests : ${libTests} ${libSources} ${libHeaders} tests/*.h ${stBarDependencies}
	${cxx} ${cflags} -I inc -I${libPath} -Wno-error -o ${binPath}/cactus_barTests ${libTests} ${libSources} ${stBarLibs}