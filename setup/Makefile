rootPath = ..
include ${rootPath}/include.mk

cflags += ${tokyoCabinetIncl}

all: all_libs all_progs
all_libs: 
all_progs: all_libs
	${MAKE} ${binPath}/cactus_setup  

${binPath}/cactus_setup : cactus_setup.c ${basicLibsDependencies} ${libPath}/cactusLib.a
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_setup cactus_setup.c ${libPath}/cactusLib.a ${basicLibs}

clean : 
	rm -f *.o
	rm -f ${binPath}/cactus_setup 
