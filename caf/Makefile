rootPath = ..
include ${rootPath}/include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

commonCafLibs = ${libPath}/cactusBlastAlignment.a ${sonLibPath}/stPinchesAndCacti.a ${sonLibPath}/3EdgeConnected.a ${libPath}/cactusLib.a
stCafDependencies =  ${commonCafLibs} ${basicLibsDependencies}
stCafLibs = ${commonCafLibs} ${basicLibs}

all: all_libs all_progs
all_libs: ${libPath}/stCaf.a
all_progs: all_libs
	${MAKE} ${binPath}/stCafTests ${binPath}/cactus_caf

${libPath}/stCaf.a : ${libSources} ${libHeaders}
	${cxx} ${cflags} -I inc -I ${libPath}/ -c ${libSources}
	ar rc stCaf.a *.o
	ranlib stCaf.a 
	rm *.o
	mv stCaf.a ${libPath}/

${binPath}/stCafTests : ${libTests} ${libPath}/stCaf.a ${stCafDependencies}
	${cxx} ${cflags} -I inc -I impl -I${libPath} -o ${binPath}/stCafTests ${libTests} ${libSources} ${libPath}/stCaf.a ${stCafLibs} -lpthread

${binPath}/cactus_caf : cactus_caf.c ${libPath}/stCaf.a ${stCafDependencies}
	${cxx} ${cflags} -I inc -I impl -I${libPath} -o ${binPath}/cactus_caf cactus_caf.c ${libSources} ${libPath}/stCaf.a ${stCafLibs} -lpthread

clean : 
	rm -f *.o
	rm -f ${libPath}/stCaf.a ${binPath}/stCafTests ${binPath}/cactus_caf

