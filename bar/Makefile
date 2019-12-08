rootPath = ..
include ${rootPath}/include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c
#${libPath}/stCaf.a
commonBarLibs =  ${libPath}/stCaf.a ${sonLibPath}/stPinchesAndCacti.a ${libPath}/cactusLib.a ${sonLibPath}/3EdgeConnected.a ${sonLibPath}/cPecanLib.a  
stBarDependencies =  ${commonBarLibs} ${basicLibsDependencies}
stBarLibs = ${commonBarLibs} ${basicLibs}

all: all_libs all_progs
all_libs: ${libPath}/cactusBarLib.a
all_progs: all_libs
	${MAKE} ${binPath}/cactus_bar ${binPath}/cactus_barTests

clean : 
	rm -f ${binPath}/cactus_barTests  ${libPath}/cactusBarLib.a

${binPath}/cactus_bar : cactus_bar.c  ${libPath}/cactusBarLib.a ${stBarDependencies} 
	${cxx} ${cflags} -I inc -I${libPath} -o ${binPath}/cactus_bar cactus_bar.c ${libPath}/cactusBarLib.a ${stBarLibs}

${binPath}/cactus_barTests : ${libTests} tests/*.h ${libPath}/cactusBarLib.a ${stBarDependencies}
	${cxx} ${cflags} -I inc -I${libPath} -Wno-error -o ${binPath}/cactus_barTests ${libTests} ${libPath}/cactusBarLib.a ${stBarLibs}

${libPath}/cactusBarLib.a : ${libSources} ${libHeaders}
	${cxx} ${cflags} -I inc -I ${libPath}/ -c ${libSources} 
	ar rc cactusBarLib.a *.o
	ranlib cactusBarLib.a 
	rm *.o
	mv cactusBarLib.a ${libPath}/
