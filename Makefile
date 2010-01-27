include ../../include.mk
binPath = ../../bin
libPath = ../../lib

cflags += ${tokyoCabinetIncl}

all : ${binPath}/cactus_core ${binPath}/cactus_phylogeny ${binPath}/cactus_setup ${binPath}/cactus_aligner.py ${binPath}/cactus_workflow.py ${binPath}/cactus_workflow_getNets ${binPath}/cactus_alignerTestAligner.py utilitiesM ${binPath}/cactus_baseAligner ${binPath}/cactus_batch.py

${binPath}/cactus_core : *.c *.h ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_core cactus_core2.c 3_Absorb3edge2x.c cactus_core.c pinchGraph.c pinchGraphTest.c pinchGraphManipulation.c cactusGraph.c cactusNetFunctions.c ${libPath}/sonLib.a ${libPath}/cactusLib.a ${tokyoCabinetLib}

${binPath}/cactus_phylogeny : *.c *.h ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_phylogeny cactus_phylogeny.c ${libPath}/sonLib.a ${libPath}/cactusLib.a ${tokyoCabinetLib}

${binPath}/cactus_setup : cactus_setup.c ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_setup cactus_setup.c ${libPath}/cactusLib.a ${libPath}/sonLib.a ${tokyoCabinetLib}

${binPath}/cactus_baseAligner : cactus_baseAligner.cc cactus_core.c cactus_core.h 3_Absorb3edge2x.c 3_Absorb3edge2x.h ${libPath}/pecan2Lib.a ${libPath}/sonLib.a ${libPath}/cactusLib.a  ${libPath}/sonLibPlus.a ${libPath}/xmlLib.a
	${cxx} ${cflags} -I${libPath} -c cactus_core.c 3_Absorb3edge2x.c pinchGraph.c pinchGraphTest.c pinchGraphManipulation.c cactusGraph.c cactusNetFunctions.c 
	${cyy} ${cflags} -I${libPath} -o ${binPath}/cactus_baseAligner cactus_baseAligner.cc cactus_core.o 3_Absorb3edge2x.o pinchGraph.o pinchGraphTest.o pinchGraphManipulation.o cactusGraph.o cactusNetFunctions.o ${libPath}/cactusLib.a ${libPath}/pecan2Lib.a ${libPath}/sonLib.a  ${libPath}/sonLibPlus.a ${libPath}/xmlLib.a ${tokyoCabinetLib}
	rm -f cactus_core.o 3_Absorb3edge2x.o pinchGraph.o pinchGraphTest.o pinchGraphManipulation.o cactusGraph.o cactusNetFunctions.o

${binPath}/cactus_aligner.py : cactus_aligner.py cactus_aligner.c ${libPath}/cactusLib.a
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_aligner cactus_aligner.c ${libPath}/cactusLib.a ${libPath}/sonLib.a ${tokyoCabinetLib}
	cp cactus_aligner.py ${binPath}/cactus_aligner.py
	chmod +x ${binPath}/cactus_aligner.py

${binPath}/cactus_batch.py : cactus_batch* cactus_misc.c cactus_misc.h ${libPath}/cactusLib.a
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_batch_chunkSequences cactus_batch_chunkSequences.c ${libPath}/sonLib.a
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_batch_convertCoordinates cactus_batch_convertCoordinates.c cactus_misc.c ${libPath}/sonLib.a 
	cp cactus_batch.py ${binPath}/cactus_batch.py
	chmod +x ${binPath}/cactus_batch.py

${binPath}/cactus_workflow.py : cactus_workflow.py
	cp cactus_workflow.py ${binPath}/cactus_workflow.py
	chmod +x ${binPath}/cactus_workflow.py

${binPath}/cactus_workflow_getNets : cactus_workflow_getNets.c ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_workflow_getNets cactus_workflow_getNets.c ${libPath}/cactusLib.a ${libPath}/sonLib.a ${tokyoCabinetLib}

${binPath}/cactus_alignerTestAligner.py : cactus_alignerTestAligner.py
	cp cactus_alignerTestAligner.py ${binPath}/cactus_alignerTestAligner.py
	chmod +x ${binPath}/cactus_alignerTestAligner.py

utilitiesM :
	#Making cactus utilities
	cd utilities && ${MAKE} all

clean : utilities.clean
	rm -f *.o
	rm -f ${binPath}/cactus_core ${binPath}/cactus_phylogeny ${binPath}/cactus_setup ${binPath}/cactus_aligner.py ${binPath}/cactus_aligner ${binPath}/cactus_workflow.py ${binPath}/cactus_alignerTestAligner.py ${binPath}/cactus_baseAligner ${binPath}/cactus_batch.py ${binPath}/cactus_batch_chunkSequences ${binPath}/cactus_batch_convertCoordinates

utilities.clean :
	#Making cactus utilities
	cd utilities && ${MAKE} clean
