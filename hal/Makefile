rootPath = ..
include ${rootPath}/include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

commonHalLibs = ${libPath}/stReference.a ${libPath}/cactusLib.a
stHalDependencies =  ${commonHalLibs} ${basicLibsDependencies}
stHalLibs = ${commonHalLibs} ${basicLibs}

all: all_libs all_progs
all_libs: 
all_progs: all_libs
	${MAKE} ${binPath}/cactus_halGenerator ${binPath}/cactus_halGeneratorTests ${binPath}/cactus_fastaGenerator

clean : 
	rm -f ${binPath}/cactus_halGenerator ${binPath}/cactus_halGeneratorTests 

${binPath}/cactus_halGenerator : cactus_halGenerator.c ${libTests} ${libSources} ${libHeaders} ${stHalDependencies}
	${cxx} ${cflags} -I inc -I${libPath} -o ${binPath}/cactus_halGenerator cactus_halGenerator.c ${libSources} ${stHalLibs}

${binPath}/cactus_fastaGenerator : cactus_fastaGenerator.c ${libTests} ${libSources} ${libHeaders} ${stHalDependencies}
	${cxx} ${cflags} -I inc -I${libPath} -o ${binPath}/cactus_fastaGenerator cactus_fastaGenerator.c ${libSources} ${stHalLibs}

${binPath}/cactus_halGeneratorTests : ${libTests} ${libSources} ${libHeaders} ${stHalDependencies}
	${cxx} ${cflags} -I inc -I${libPath} -Wno-error -o ${binPath}/cactus_halGeneratorTests ${libTests} ${libSources} ${stHalLibs}
