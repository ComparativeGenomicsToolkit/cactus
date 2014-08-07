rootPath = ../
include ../include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

commonBarLibs = ${libPath}/stCaf.a ${sonLibPath}/stPinchesAndCacti.a ${libPath}/cactusLib.a ${sonLibPath}/3EdgeConnected.a
stBarDependencies =  ${commonBarLibs} ${basicLibsDependencies}
stBarLibs = ${commonBarLibs} ${basicLibs}

all : ${libPath}/cactusBarLib.a ${binPath}/cactus_bar ${binPath}/cactus_realign ${binPath}/cactus_expectationMaximisation ${binPath}/cactus_barTests 
 
clean : 
	rm -f ${binPath}/cactus_bar ${binPath}/cactus_realign ${binPath}/cactus_expectationMaximisation ${binPath}/cactus_barTests  ${libPath}/cactusBarLib.a

${binPath}/cactus_bar : cactus_bar.c  ${libPath}/cactusBarLib.a ${stBarDependencies} 
	${cxx} ${cflags} -I inc -I${libPath} -o ${binPath}/cactus_bar cactus_bar.c ${libPath}/cactusBarLib.a ${stBarLibs}

${binPath}/cactus_realign : cactus_realign.c ${libPath}/cactusBarLib.a ${stBarDependencies} 
	${cxx} ${cflags} -I inc -I${libPath} -o ${binPath}/cactus_realign cactus_realign.c ${libPath}/cactusBarLib.a ${stBarLibs}
	
${binPath}/cactus_expectationMaximisation : cactus_expectationMaximisation.py
	cp cactus_expectationMaximisation.py ${binPath}/cactus_expectationMaximisation
	chmod +x ${binPath}/cactus_expectationMaximisation

${binPath}/cactus_barTests : ${libTests} tests/*.h ${libPath}/cactusBarLib.a ${stBarDependencies}
	${cxx} ${cflags} -I inc -I${libPath} -Wno-error -o ${binPath}/cactus_barTests ${libTests} ${libPath}/cactusBarLib.a ${stBarLibs}
	
${libPath}/cactusBarLib.a : ${libSources} ${libHeaders} ${stBarDependencies}
	${cxx} ${cflags} -I inc -I ${libPath}/ -c ${libSources} 
	ar rc cactusBarLib.a *.o
	ranlib cactusBarLib.a 
	rm *.o
	mv cactusBarLib.a ${libPath}/
	cp ${libHeaders} ${libPath}/