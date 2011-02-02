rootPath = ../../
include ../../include.mk

all :  ${binPath}/cactus_referenceViewer

${binPath}/cactus_referenceViewer : *.c *.h ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_referenceViewer main.c referenceViewer.c graphViz.c ${libPath}/cactusLib.a ${basicLibs}

clean :
	rm -f *.o
	rm -f ${binPath}/cactus_referenceViewer
