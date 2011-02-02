rootPath = ../../
include ../../include.mk

all :  ${binPath}/cactus_adjacencyGraphViewer ${binPath}/cactus_treeViewer 

${binPath}/cactus_adjacencyGraphViewer : *.c ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_adjacencyGraphViewer cactus_adjacencyGraphViewer.c ${libPath}/cactusLib.a ${basicLibs}

${binPath}/cactus_treeViewer : *.c  ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_treeViewer cactus_treeViewer.c ${libPath}/cactusLib.a ${basicLibs}

clean :
	rm -f *.o
	rm -f ${binPath}/cactus_adjacencyGraphViewer ${binPath}/cactus_treeViewer 
