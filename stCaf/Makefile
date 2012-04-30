rootPath = ../
include ../include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

commonCafLibs = ${libPath}/cactusBlastAlignment.a ${libPath}/stCactusGraphs.a ${libPath}/stPinchGraphs.a ${libPath}/3EdgeConnected.a ${libPath}/cactusLib.a
cafDependencies =  ${commonCafLibs} ${basicLibsDependencies}
cafLibs = ${commonCafLibs} ${basicLibs}

all : ${libPath}/cactusCaf.a ${binPath}/cafTests

${libPath}/cactusCaf.a : ${libSources} ${libHeaders} ${cafDependencies}
	${cxx} ${cflags} -I inc -I ${libPath}/ -c ${libSources}
	ar rc cactusCaf.a *.o
	ranlib cactusCaf.a 
	rm *.o
	mv cactusCaf.a ${libPath}/
	cp ${libHeaders} ${libPath}/

${binPath}/cafTests : ${libTests} ${libPath}/cactusCaf.a ${cafDependencies}
	${cxx} ${cflags} -I inc -I impl -I${libPath} -o ${binPath}/cafTests ${libTests} ${libSources} ${libPath}/cactusCaf.a ${cafLibs}

clean : 
	rm -f *.o
	rm -f ${libPath}/cactusCaf.a ${binPath}/cafTests

