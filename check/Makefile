rootPath = ../
include ../include.mk

cflags += ${tokyoCabinetIncl}

all :  ${binPath}/cactus_check 

${binPath}/cactus_check : *.c ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_check cactus_check.c ${libPath}/cactusLib.a ${basicLibs}

clean :
	rm -f *.o
	rm -f ${binPath}/cactus_check 
