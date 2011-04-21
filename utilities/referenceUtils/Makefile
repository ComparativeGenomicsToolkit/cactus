rootPath = ../../
include ../../include.mk

progs = $(notdir $(wildcard))

targets = ${progs:%=${binPath}/%} ${binPath}/cactus_addReferenceSeq ${binPath}/cactus_addReferenceSeqTest ${binPath}/cactus_getReferenceSeq

all: ${targets}

#${libPath}/cactus_addReferenceSeq.o : cactus_addReferenceSeq.h cactus_addReferenceSeq.c
#	cp cactus_addReferenceSeq.h ${libPath}
#	${cxx} ${cflags} -c cactus_addReferenceSeq.c  -I ${libPath} -o ${libPath}/cactus_addReferenceSeq.o

${binPath}/cactus_addReferenceSeq: cactus_addReferenceSeq.c ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_addReferenceSeq cactus_addReferenceSeq.c ${libPath}/cactusLib.a ${basicLibs}

${binPath}/cactus_addReferenceSeqTest: cactus_addReferenceSeqTest.c ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_addReferenceSeqTest cactus_addReferenceSeqTest.c ${libPath}/cactusLib.a ${basicLibs}

${binPath}/cactus_getReferenceSeq: cactus_getReferenceSeq.c ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_getReferenceSeq cactus_getReferenceSeq.c ${libPath}/cactusLib.a ${basicLibs}

#${binPath}/cactus_addReferenceSeqTest: ${libPath}/cactus_addReferenceSeq.o cactus_addReferenceSeqTest.c ${libPath}/cactusLib.a ${basicLibsDependencies}
#	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_addReferenceSeqTest cactus_addReferenceSeqTest.c ${libPath}/cactus_addReferenceSeq.o ${libPath}/cactusLib.a ${basicLibs}

#${binPath}/cactus_getReferenceSeq: ${libPath}/cactus_addReferenceSeq.o cactus_getReferenceSeq.c ${libPath}/cactusLib.a ${basicLibsDependencies}
#	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_getReferenceSeq cactus_getReferenceSeq.c ${libPath}/cactus_addReferenceSeq.o ${libPath}/cactusLib.a ${basicLibs}

clean :
	rm -f *.o
	rm -f ${binPath}/cactus_addReferenceSeqTest
	rm -f ${binPath}/cactus_*ReferenceSeq
	rm -f ${libPath}/cactus_addReferenceSeq.*
