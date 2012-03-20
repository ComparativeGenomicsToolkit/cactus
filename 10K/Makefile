rootPath = ../
include ../include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

all : ${binPath}/cactus_10KGenerator ${binPath}/cactus_10KGeneratorTests 
 
clean : 
	rm -f ${binPath}/cactus_10KGenerator ${binPath}/cactus_10KGeneratorTests 

${binPath}/cactus_10KGenerator : cactus_10KGenerator.c ${libTests} ${libSources} ${libHeaders} ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I inc -I${libPath} -o ${binPath}/cactus_10KGenerator cactus_10KGenerator.c ${libSources} ${libPath}/cactusLib.a ${basicLibs}

${binPath}/cactus_10KGeneratorTests : ${libTests} ${libSources} ${libHeaders} ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I inc -I${libPath} -Wno-error -o ${binPath}/cactus_10KGeneratorTests ${libTests} ${libSources} ${libPath}/cactusLib.a ${basicLibs}