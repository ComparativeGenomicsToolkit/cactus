include ../../../include.mk
binPath = ../../../bin
libPath = ../../../lib

libSources = cactus*.c
libHeaders = cactus*.h
cflags += ${tokyoCabinetIncl}

all : ${libPath}/cactusLib.a ${binPath}/cactusAPITests

clean : 
	rm -f ${libPath}/cactusLib.a ${libPath}/cactus*.h ${binPath}/cactusAPITests

${libPath}/cactusLib.a : ${libSources} ${libHeaders} ${libPath}/sonLib.a 
	${cxx} ${cflags} -I ${libPath}/ -c ${libSources} 
	ar rc cactusLib.a *.o
	ranlib cactusLib.a
	rm *.o
	mv cactusLib.a ${libPath}/
	cp cactus*.h ${libPath}/
	
${binPath}/cactusAPITests : ${libSources} ${libHeaders} allTests.c ${libPath}/cactusLib.a ${libPath}/sonLib.a ${libPath}/cuTest.a
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactusAPITests allTests.c ${libPath}/cactusLib.a ${libPath}/sonLib.a ${libPath}/cuTest.a ${tokyoCabinetLib}
