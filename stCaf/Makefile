rootPath = ../
include ../include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

all : ${libPath}/stCactusGraphs.a ${binPath}/stCactusGraphTests

${libPath}/stCactusGraphs.a : ${libSources} ${libHeaders} ${basicLibsDependencies}
	${cxx} ${cflags} -I inc -I ${libPath}/ -c ${libSources}
	ar rc stCactusGraphs.a *.o
	ranlib stCactusGraphs.a 
	rm *.o
	mv stCactusGraphs.a ${libPath}/
	cp ${libHeaders} ${libPath}/

${binPath}/stCactusGraphTests : ${libTests} ${libSources} ${libHeaders} ${basicLibsDependencies} ${libPath}/3EdgeConnected.a
	${cxx} ${cflags} -I inc -I impl -I${libPath} -o ${binPath}/stCactusGraphTests ${libTests} ${libSources} ${basicLibs}  ${libPath}/3EdgeConnected.a

clean : 
	rm -f *.o
	rm -f ${libPath}/stCactusGraphs.a ${binPath}/stCactusGraphTests

