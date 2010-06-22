include ../../../include.mk
binPath = ../../../bin
libPath = ../../../lib

cpflags += ${tokyoCabinetIncl}

all : ${binPath}/cactus_baseAligner 
#${binPath}/cactus_baseAlignerTests

${binPath}/cactus_baseAligner : cactus_baseAligner.cc ${libPath}/pecan2Lib.a ${libPath}/sonLib.a ${libPath}/cactusLib.a  ${libPath}/sonLibPlus.a ${libPath}/xmlLib.a
	${cyy} ${cpflags} -I${libPath} -o ${binPath}/cactus_baseAligner cactus_baseAligner.cc ${libPath}/cactusCore.a ${libPath}/cactusGraph.a ${libPath}/pinchGraph.a ${libPath}/3EdgeConnected.a ${libPath}/cactusLib.a ${libPath}/pecan2Lib.a ${libPath}/sonLib.a  ${libPath}/sonLibPlus.a ${libPath}/xmlLib.a ${tokyoCabinetLib}
	rm -f cactus_core.o 3_Absorb3edge2x.o pinchGraph.o pinchGraphTest.o pinchGraphManipulation.o cactusGraph.o cactusNetFunctions.o

${binPath}/cactus_baseAlignerTests : allTests.c stPosetAlignment* ${libPath}/sonLib.a ${libPath}/cuTest.a
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_baseAlignerTests allTests.c stPosetAlignmentTest.c stPosetAlignment.c ${libPath}/sonLib.a ${libPath}/cuTest.a 

clean : 
	rm -f *.o
	rm -f ${binPath}/cactus_baseAligner
