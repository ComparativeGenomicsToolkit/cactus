include ../../../../include.mk
binPath = ../../../../bin
libPath = ../../../../lib

cflags += ${tokyoCabinetIncl}

all :  ${binPath}/cactus_adjacencyGraphViewer ${binPath}/cactus_treeViewer 

${binPath}/cactus_adjacencyGraphViewer : *.c ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_adjacencyGraphViewer cactus_adjacencyGraphViewer.c ${libPath}/sonLib.a ${libPath}/cactusLib.a ${libPath}/sonLib.a ${dblibs}

${binPath}/cactus_treeViewer : *.c  ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_treeViewer cactus_treeViewer.c ${libPath}/sonLib.a ${libPath}/cactusLib.a ${libPath}/sonLib.a ${dblibs}

clean :
	rm -f *.o
	rm -f ${binPath}/cactus_adjacencyGraphViewer ${binPath}/cactus_treeViewer 
