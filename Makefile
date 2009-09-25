makefile:

include ../../include.mk
binPath = ../../bin
libPath = ../../lib

all : ${binPath}/cactus_core ${binPath}/cactus_3Edge ${binPath}/cactus_checkReconstructionTree ${binPath}/cactus_setup.py ${binPath}/cactus_aligner.py ${binPath}/cactus_workflow.py ${binPath}/cactus_alignerTestAligner.py ${binPath}/cactus_coreTestTreeBuilder.py ${binPath}/cactus_adjacencyTestAdjacencyBuilder.py ${binPath}/cactus_reconstructionTreeViewer.py ${binPath}/cactus_atomGraphViewer.py

${binPath}/cactus_3Edge : 3_Absorb3edge2x.c ${libPath}/sonLib.a
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_3Edge 3_Absorb3edge2x.c ${libPath}/sonLib.a

${binPath}/cactus_core : *.cc *.c *.h ${libPath}/sonLib.a ${libPath}/sonLibPlus.a ${libPath}/xmlLib.a
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_core cactus_core.c pinchGraph.c pinchGraphTest.c pinchGraphManipulation.c cactusGraph.c net.c netSerialisation.c cactusNetFunctions.c ${libPath}/sonLib.a
	rm *.o

${binPath}/cactus_setup.py : cactus_setup.py cactus_setup_processFasta.c ${libPath}/sonLib.a
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_setup_processFasta cactus_setup_processFasta.c ${libPath}/sonLib.a
	cp cactus_setup.py ${binPath}/cactus_setup.py
	chmod +x ${binPath}/cactus_setup.py
	
${binPath}/cactus_aligner.py : cactus_aligner.py
	cp cactus_aligner.py ${binPath}/cactus_aligner.py
	chmod +x ${binPath}/cactus_aligner.py

${binPath}/cactus_workflow.py : cactus_workflow.py
	cp cactus_workflow.py ${binPath}/cactus_workflow.py
	chmod +x ${binPath}/cactus_workflow.py
	
${binPath}/cactus_alignerTestAligner.py : cactus_alignerTestAligner.py
	cp cactus_alignerTestAligner.py ${binPath}/cactus_alignerTestAligner.py
	chmod +x ${binPath}/cactus_alignerTestAligner.py
	
${binPath}/cactus_coreTestTreeBuilder.py : cactus_coreTestTreeBuilder.py
	cp cactus_coreTestTreeBuilder.py ${binPath}/cactus_coreTestTreeBuilder.py
	chmod +x ${binPath}/cactus_coreTestTreeBuilder.py
	
${binPath}/cactus_adjacencyTestAdjacencyBuilder.py : cactus_adjacencyTestAdjacencyBuilder.py
	cp cactus_adjacencyTestAdjacencyBuilder.py ${binPath}/cactus_adjacencyTestAdjacencyBuilder.py
	chmod +x ${binPath}/cactus_adjacencyTestAdjacencyBuilder.py

${binPath}/cactus_reconstructionTreeViewer.py : cactus_reconstructionTreeViewer.py
	cp cactus_reconstructionTreeViewer.py ${binPath}/cactus_reconstructionTreeViewer.py
	chmod +x ${binPath}/cactus_reconstructionTreeViewer.py

${binPath}/cactus_atomGraphViewer.py : cactus_atomGraphViewer.py
	cp cactus_atomGraphViewer.py ${binPath}/cactus_atomGraphViewer.py
	chmod +x ${binPath}/cactus_atomGraphViewer.py

clean :
	rm -f *.o
	rm -f ${binPath}/cactus_core ${binPath}/cactus_3Edge ${binPath}/cactus_checkReconstructionTree ${binPath}/cactus_setup.py ${binPath}/cactus_setup_processFasta ${binPath}/cactus_aligner.py ${binPath}/cactus_workflow.py ${binPath}/cactus_alignerTestAligner.py ${binPath}/cactus_coreTestTreeBuilder.py ${binPath}/cactus_adjacencyTestAdjacencyBuilder.py ${binPath}/cactus_reconstructionTreeViewer.py ${binPath}/cactus_atomGraphViewer.py
