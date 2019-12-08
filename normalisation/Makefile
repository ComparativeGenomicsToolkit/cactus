rootPath = ..
include ${rootPath}/include.mk

libSources = main.c normal.c
libHeaders = normal.h

all: all_libs all_progs
all_libs: 
all_progs: all_libs
	${MAKE} ${binPath}/cactus_normalisation

${binPath}/cactus_normalisation : ${libSources} ${libHeaders} ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_normalisation main.c normal.c ${libPath}/cactusLib.a ${basicLibs}

clean : 
	rm -f ${binPath}/cactus_normalisation
