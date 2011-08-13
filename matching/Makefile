rootPath = ../
include ../include.mk

cflags += ${tokyoCabinetIncl}

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

all : ${libPath}/cactusMatching.a ${binPath}/matchingTests

${libPath}/cactusMatching.a : ${libSources} ${libHeaders} ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I inc -I ${libPath}/ -c ${libSources}
	ar rc cactusMatching.a *.o
	ranlib cactusMatching.a 
	rm *.o
	mv cactusMatching.a ${libPath}/
	cp ${libHeaders} ${libPath}/

${binPath}/matchingTests : ${libTests} ${libSources} ${libHeaders} ${libPath}/cactusMatching.a ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I inc -I impl -I${libPath} -o ${binPath}/matchingTests ${libTests} ${libPath}/cactusMatching.a ${libPath}/cactusLib.a ${basicLibs}

clean : 
	rm -f *.o
	rm -f ${libPath}/cactusMatching.a ${binPath}/matchingTests 
