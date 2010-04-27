include ../../../../include.mk
binPath = ../../../../bin
libPath = ../../../../lib

cflags += ${tokyoCabinetIncl}

all :  ${binPath}/cactus_referenceViewer

${binPath}/cactus_referenceViewer : *.c *.h  ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_referenceViewer main.c referenceViewer.c graphViz.c ${libPath}/sonLib.a ${libPath}/cactusLib.a ${libPath}/sonLib.a ${tokyoCabinetLib}

clean :
	rm -f *.o
	rm -f ${binPath}/cactus_referenceViewer
