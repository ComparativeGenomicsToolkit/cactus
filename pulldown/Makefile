include ../../../include.mk
binPath = ../../../bin
libPath = ../../../lib

testSources = cactus_pulldownTest.c
testSharedSources = cactus_pulldownTestShared.c
libSources = cactus_pulldown.c 
libHeaders = cactus_pulldown.h
cflags += ${tokyoCabinetIncl}

pulldownLib=${libPath}/cactus_pulldown.a
pulldownTest=${binPath}/cactus_pulldownTest
sonLib=${libPath}/sonLib.a
cactusLib=${libPath}/cactusLib.a

all : ${pulldownLib} ${pulldownTest}

${pulldownLib} : ${libSources} ${libHeaders} ${sonLib} ${cactusLib}
	${cxx} ${cflags} -I ${libPath}/ -c ${libSources}
	ar rc cactus_pulldown.a *.o
	ranlib cactus_pulldown.a 
	rm *.o
	mv cactus_pulldown.a ${libPath}/
	cp ${libHeaders} ${libPath}/

${pulldownTest}: ${pulldownLib} ${sonLib} ${cactusLib} ${testSources}
	${cxx} ${cflags} -I ${libPath}/ -o $@ ${testSources} ${pulldownLib} ${sonLib} ${cactusLib} ${dblibs}

clean : 
	rm -f ${libPath}/cactus_pulldown.a 

