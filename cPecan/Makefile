rootPath = ../
include ../include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c
#${libPath}/stCaf.a

cPecanDependencies =  ${basicLibsDependencies}
cPecanLibs = ${basicLibs}

all : ${libPath}/cPecanLib.a ${binPath}/cactus_realign ${binPath}/cactus_expectationMaximisation 
 
clean : 
	rm -f ${binPath}/cactus_realign ${binPath}/cactus_expectationMaximisation ${binPath}/cPecanLibTests  ${libPath}/cPecanLib.a

${binPath}/cactus_realign : cactus_realign.c ${libPath}/cPecanLib.a ${cPecanDependencies} 
	${cxx} ${cflags} -I inc -I${libPath} -o ${binPath}/cactus_realign cactus_realign.c ${libPath}/cPecanLib.a ${cPecanLibs}
	
${binPath}/cactus_expectationMaximisation : cactus_expectationMaximisation.py
	cp cactus_expectationMaximisation.py ${binPath}/cactus_expectationMaximisation
	chmod +x ${binPath}/cactus_expectationMaximisation

${binPath}/cPecanLibTests : ${libTests} tests/*.h ${libPath}/cPecanLibLib.a ${cPecanDependencies}
	${cxx} ${cflags} -I inc -I${libPath} -Wno-error -o ${binPath}/cPecanLibTests ${libTests} ${libPath}/cPecanLib.a ${cPecanLibs}
	
${libPath}/cPecanLib.a : ${libSources} ${libHeaders} ${stBarDependencies}
	${cxx} ${cflags} -I inc -I ${libPath}/ -c ${libSources} 
	ar rc cPecanLib.a *.o
	ranlib cPecanLib.a 
	rm *.o
	mv cPecanLib.a ${libPath}/
	cp ${libHeaders} ${libPath}/
