rootPath = ../../
include ../../include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

all : ${libPath}/3EdgeConnected.a ${binPath}/3EdgeTests

${libPath}/3EdgeConnected.a : ${libSources} ${libHeaders} ${basicLibsDependencies}
	${cxx} ${cflags} -I ${libPath}/ -c ${libSources}
	ar rc 3EdgeConnected.a *.o
	ranlib 3EdgeConnected.a 
	rm *.o
	mv 3EdgeConnected.a ${libPath}/
	cp ${libHeaders} ${libPath}/

clean : 
	rm -f ${libPath}/3EdgeConnected.a ${binPath}/3EdgeTests

${binPath}/3EdgeTests : ${libTests} ${libSources} ${libHeaders} ${libPath}/3EdgeConnected.a ${basicLibsDependencies}
	${cxx} ${cflags} -I inc -I${libPath} -o ${binPath}/3EdgeTests ${libTests} ${libPath}/3EdgeConnected.a ${basicLibs}
