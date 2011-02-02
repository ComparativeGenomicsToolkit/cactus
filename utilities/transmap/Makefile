rootPath = ../../
include ../../include.mk

cflags += ${tokyoCabinetIncl}

all : ${libPath}/cactus_transmap.o

${libPath}/cactus_transmap.o : *.c ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I${libPath} -c -o cactus_transmap.o cactus_transmap.c ${libPath}/sonLib.a ${libPath}/cactusLib.a ${basicLibs}

clean :
	rm -f *.o
	rm -f cactus_transmap.o
