rootPath = ../
include ../include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

all : ${binPath}/cactus_halGenerator ${binPath}/cactus_halGeneratorTests 
 
clean : 
	rm -f ${binPath}/cactus_halGenerator ${binPath}/cactus_halGeneratorTests 

${binPath}/cactus_halGenerator : cactus_halGenerator.c ${libTests} ${libSources} ${libHeaders} ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I inc -I${libPath} -o ${binPath}/cactus_halGenerator cactus_halGenerator.c ${libSources} ${libPath}/cactusLib.a ${basicLibs}

${binPath}/cactus_halGeneratorTests : ${libTests} ${libSources} ${libHeaders} ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I inc -I${libPath} -Wno-error -o ${binPath}/cactus_halGeneratorTests ${libTests} ${libSources} ${libPath}/cactusLib.a ${basicLibs}