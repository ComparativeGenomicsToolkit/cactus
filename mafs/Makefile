rootPath = ../
include ../include.mk

cflags += ${tokyoCabinetIncl}

all : ${binPath}/cactus_recursiveMafGenerator

${binPath}/cactus_recursiveMafGenerator : *.c ${basicLibsDependencies}
	${cxx} ${cflags} -I inc -I${libPath} -o ${binPath}/cactus_recursiveMafGenerator *.c ${libPath}/cactusLib.a ${basicLibs} 

clean : 
	rm -f *.o
	rm -f ${binPath}/cactus_recursiveMafGenerator 

