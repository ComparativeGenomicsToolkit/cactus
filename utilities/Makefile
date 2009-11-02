include ../../../include.mk
binPath = ../../../bin
libPath = ../../../lib

cflags += ${tokyoCabinetIncl}

all :  ${binPath}/cactus_treeViewer ${binPath}/cactus_treeStats

${binPath}/cactus_treeViewer : cactus_treeViewer.c ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_treeViewer cactus_treeViewer.c  ${libPath}/sonLib.a ${libPath}/cactusLib.a ${libPath}/sonLib.a ${tokyoCabinetLib}

${binPath}/cactus_treeStats : cactus_treeStats.c ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_treeStats cactus_treeStats.c  ${libPath}/sonLib.a ${libPath}/cactusLib.a ${libPath}/sonLib.a ${tokyoCabinetLib}

clean :
	rm -f *.o
	rm -f ${binPath}/cactus_treeViewer ${binPath}/cactus_treeStats
