rootPath = ../
include ../include.mk

cflags += ${tokyoCabinetIncl}

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

all : ${binPath}/cactus_reference ${binPath}/cactus_addReferenceCoordinates ${binPath}/referenceTests

${binPath}/cactus_reference : cactus_reference.c ${libSources} ${libHeaders} ${libPath}/cactusLib.a  ${basicLibsDependencies}
	${cxx} ${cflags} -I inc -I${libPath} -o ${binPath}/cactus_reference cactus_reference.c ${libSources} ${libPath}/cactusLib.a ${basicLibs}

${binPath}/cactus_addReferenceCoordinates : cactus_addReferenceCoordinates.c ${libSources} ${libHeaders} ${libPath}/cactusLib.a  ${basicLibsDependencies}
	${cxx} ${cflags} -I inc -I${libPath} -o ${binPath}/cactus_addReferenceCoordinates cactus_addReferenceCoordinates.c ${libSources} ${libPath}/cactusLib.a ${basicLibs}

${binPath}/referenceTests : ${libTests} ${libSources} ${libHeaders} ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I inc -I${libPath} -o ${binPath}/referenceTests ${libTests} ${libSources} ${libPath}/cactusLib.a ${basicLibs}

clean : 
	rm -f *.o
	rm -f ${binPath}/cactus_reference ${binPath}/referenceTests ${binPath}/cactus_addReferenceCoordinates
