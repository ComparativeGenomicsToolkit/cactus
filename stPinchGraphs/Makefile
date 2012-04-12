rootPath = ../
include ../include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

all : ${libPath}/stPinchGraphs.a ${binPath}/stPinchGraphTests

${libPath}/stPinchGraphs.a : ${libSources} ${libHeaders} ${basicLibsDependencies}
	${cxx} ${cflags} -I inc -I ${libPath}/ -c ${libSources}
	ar rc stPinchGraphs.a *.o
	ranlib stPinchGraphs.a 
	rm *.o
	mv stPinchGraphs.a ${libPath}/
	cp ${libHeaders} ${libPath}/

${binPath}/stPinchGraphTests : ${libTests} ${libSources} ${libHeaders} ${basicLibsDependencies}
	${cxx} ${cflags} -I inc -I impl -I${libPath} -o ${binPath}/stPinchGraphTests ${libTests} ${libSources} ${basicLibs} 

clean : 
	rm -f *.o
	rm -f ${libPath}/stPinchGraphs.a ${binPath}/stPinchGraphTests

