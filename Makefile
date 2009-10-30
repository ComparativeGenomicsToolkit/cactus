include ../../include.mk
binPath = ../../bin
libPath = ../../lib

cflags += ${tokyoCabinetIncl}

all : ${binPath}/cactus_3Edge ${binPath}/cactus_core ${binPath}/cactus_phylogeny ${binPath}/cactus_setup ${binPath}/cactus_aligner.py ${binPath}/cactus_workflow.py ${binPath}/cactus_workflow_getNets ${binPath}/cactus_alignerTestAligner.py ${binPath}/cactus_adjacencyTestAdjacencyBuilder.py utilitiesM

${binPath}/cactus_3Edge : 3_Absorb3edge2x.c ${libPath}/sonLib.a
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_3Edge 3_Absorb3edge2x.c ${libPath}/sonLib.a 

${binPath}/cactus_core : *.c *.h ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_core cactus_core.c pinchGraph.c pinchGraphTest.c pinchGraphManipulation.c cactusGraph.c cactusNetFunctions.c ${libPath}/sonLib.a ${libPath}/cactusLib.a ${tokyoCabinetLib}

${binPath}/cactus_phylogeny : *.c *.h ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_phylogeny cactus_phylogeny.c  ${libPath}/sonLib.a ${libPath}/cactusLib.a ${tokyoCabinetLib}

${binPath}/cactus_setup : cactus_setup.c ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_setup cactus_setup.c ${libPath}/cactusLib.a ${libPath}/sonLib.a ${tokyoCabinetLib}

${binPath}/cactus_aligner.py : cactus_aligner.py cactus_aligner.c ${libPath}/cactusLib.a
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_aligner cactus_aligner.c ${libPath}/cactusLib.a ${libPath}/sonLib.a ${tokyoCabinetLib}
	cp cactus_aligner.py ${binPath}/cactus_aligner.py
	chmod +x ${binPath}/cactus_aligner.py

${binPath}/cactus_workflow.py : cactus_workflow.py
	cp cactus_workflow.py ${binPath}/cactus_workflow.py
	chmod +x ${binPath}/cactus_workflow.py

${binPath}/cactus_workflow_getNets : cactus_workflow_getNets.c ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_workflow_getNets cactus_workflow_getNets.c ${libPath}/cactusLib.a ${libPath}/sonLib.a ${tokyoCabinetLib}

${binPath}/cactus_alignerTestAligner.py : cactus_alignerTestAligner.py
	cp cactus_alignerTestAligner.py ${binPath}/cactus_alignerTestAligner.py
	chmod +x ${binPath}/cactus_alignerTestAligner.py

${binPath}/cactus_adjacencyTestAdjacencyBuilder.py : cactus_adjacencyTestAdjacencyBuilder.py
	cp cactus_adjacencyTestAdjacencyBuilder.py ${binPath}/cactus_adjacencyTestAdjacencyBuilder.py
	chmod +x ${binPath}/cactus_adjacencyTestAdjacencyBuilder.py

utilitiesM :
	#Making cactus utilities
	cd utilities && make all
	
clean : utilities.clean
	rm -f *.o
	rm -f ${binPath}/cactus_core ${binPath}/cactus_phylogeny ${binPath}/cactus_3Edge ${binPath}/cactus_setup ${binPath}/cactus_aligner.py ${binPath}/cactus_aligner ${binPath}/cactus_workflow.py ${binPath}/cactus_alignerTestAligner.py ${binPath}/cactus_adjacencyTestAdjacencyBuilder.py 

utilities.clean :
	#Making cactus utilities
	cd utilities && make clean