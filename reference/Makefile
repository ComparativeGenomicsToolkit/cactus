rootPath = ../
include ../include.mk

cflags += ${tokyoCabinetIncl}

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

all : ${libPath}/cactusReference.a ${binPath}/cactus_reference ${binPath}/cactus_addReferenceCoordinates ${binPath}/referenceTests

${libPath}/cactusReference.a : ${libSources} ${libHeaders} ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I inc -I ${libPath}/ -c ${libSources}
	ar rc cactusReference.a *.o
	ranlib cactusReference.a 
	rm *.o
	mv cactusReference.a ${libPath}/
	cp ${libHeaders} ${libPath}/
	
${binPath}/cactus_reference : cactus_reference.c ${libSources} ${libHeaders} ${libPath}/cactusReference.a ${libPath}/cactusMatching.a ${libPath}/cactusLib.a  ${basicLibsDependencies}
	${cxx} ${cflags} -I inc -I${libPath} -o ${binPath}/cactus_reference cactus_reference.c ${libPath}/cactusReference.a ${libPath}/cactusMatching.a ${libPath}/cactusLib.a ${basicLibs}

${binPath}/cactus_addReferenceCoordinates : ${libSources} ${libHeaders} cactus_addReferenceCoordinates.c ${libPath}/cactusReference.a ${libPath}/cactusMatching.a ${libSources} ${libPath}/cactusLib.a  ${basicLibsDependencies}
	${cxx} ${cflags} -I inc -I${libPath} -o ${binPath}/cactus_addReferenceCoordinates cactus_addReferenceCoordinates.c ${libPath}/cactusReference.a ${libPath}/cactusMatching.a ${libPath}/cactusLib.a ${basicLibs}

${binPath}/referenceTests : ${libTests} ${libSources} ${libHeaders} ${libPath}/cactusReference.a ${libPath}/cactusMatching.a ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I inc -I impl -I${libPath} -o ${binPath}/referenceTests ${libTests} ${libPath}/cactusReference.a ${libPath}/cactusMatching.a ${libPath}/cactusLib.a ${basicLibs}

clean : 
	rm -f *.o
	rm -f ${libPath}/cactusReference.a ${binPath}/cactus_reference ${binPath}/referenceTests ${binPath}/cactus_addReferenceCoordinates
