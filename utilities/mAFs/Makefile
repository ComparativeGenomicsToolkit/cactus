include ../../../../include.mk
binPath = ../../../../bin
libPath = ../../../../lib

cflags += ${tokyoCabinetIncl}

all :  ${binPath}/cactus_MAFGenerator

${binPath}/cactus_MAFGenerator : *.c ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_MAFGenerator cactus_MAFGenerator.c  ${libPath}/sonLib.a ${libPath}/cactusLib.a ${libPath}/sonLib.a ${tokyoCabinetLib}

clean :
	rm -f *.o
	rm -f ${binPath}/cactus_MAFGenerator
