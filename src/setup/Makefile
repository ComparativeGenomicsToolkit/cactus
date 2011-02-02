rootPath = ../../
include ../../include.mk

cflags += ${tokyoCabinetIncl}

all : ${binPath}/cactus_setup

${binPath}/cactus_setup : cactus_setup.c ${basicLibsDependencies}
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_setup cactus_setup.c ${libPath}/cactusLib.a ${basicLibs}

clean : 
	rm -f *.o
	rm -f ${binPath}/cactus_setup