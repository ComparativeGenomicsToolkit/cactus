include ../../../include.mk
binPath = ../../../bin
libPath = ../../../lib

cflags += ${tokyoCabinetIncl}

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

all : ${binPath}/cactus_reference ${binPath}/referenceTests

${binPath}/cactus_reference : cactus_reference.c ${libSources} ${libHeaders} ${libPath}/sonLib.a ${libPath}/cactusLib.a 
	${cxx} ${cflags} -I${libPath} -I inc -o ${binPath}/cactus_reference cactus_reference.c ${libSources} ${libPath}/cactusLib.a ${libPath}/sonLib.a ${dblibs}

${binPath}/referenceTests : ${libTests} ${libSources} ${libHeaders} ${libPath}/sonLib.a ${libPath}/cactusLib.a ${libPath}/cuTest.a
	${cxx} ${cflags} -I inc -I ${libPath} -o ${binPath}/referenceTests ${libTests} ${libSources} ${libPath}/sonLib.a ${libPath}/cactusLib.a ${libPath}/cuTest.a ${dblibs} ${mysqlLibs}

clean : 
	rm -f *.o
	rm -f ${binPath}/cactus_reference
