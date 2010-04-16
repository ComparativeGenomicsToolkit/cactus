include ../../../../include.mk
binPath = ../../../../bin
libPath = ../../../../lib

cflags += ${tokyoCabinetIncl}

all :  ${binPath}/cactus_treeStats

${binPath}/cactus_treeStats : *.c *.h ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_treeStats main.c treeStats.c ${libPath}/sonLib.a ${libPath}/cactusLib.a ${libPath}/sonLib.a ${tokyoCabinetLib}

clean :
	rm -f *.o
	rm -f ${binPath}/cactus_treeStats 
