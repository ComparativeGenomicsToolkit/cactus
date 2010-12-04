include ../../../../include.mk
binPath = ../../../../bin
libPath = ../../../../lib

cflags += ${tokyoCabinetIncl} -lssl

progs = $(notdir $(wildcard))

targets = ${progs:%=${binPath}/%} ${binPath}/cactus_addReferenceSeqTest ${binPath}/cactus_getReferenceSeq ${libPath}/cactus_addReferenceSeq.o

all: ${targets}

${libPath}/cactus_addReferenceSeq.o : cactus_addReferenceSeq.h cactus_addReferenceSeq.c
	cp cactus_addReferenceSeq.h ${libPath}
	${cxx} ${cflags} -c cactus_addReferenceSeq.c  -I ${libPath} -o ${libPath}/cactus_addReferenceSeq.o

${binPath}/cactus_addReferenceSeqTest: ${libPath}/cactus_addReferenceSeq.o cactus_addReferenceSeqTest.c ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_addReferenceSeqTest cactus_addReferenceSeqTest.c ${libPath}/cactus_addReferenceSeq.o ${libPath}/cactusLib.a ${libPath}/sonLib.a ${dblibs}

${binPath}/cactus_getReferenceSeq: ${libPath}/cactus_addReferenceSeq.o cactus_getReferenceSeq.c ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_getReferenceSeq cactus_getReferenceSeq.c ${libPath}/cactus_addReferenceSeq.o ${libPath}/cactusLib.a ${libPath}/sonLib.a ${dblibs}

clean :
	rm -f *.o
	rm -f ${binPath}/cactus_addReferenceSeqTest
	rm -f ${binPath}/cactus_getReferenceSeq
	rm -f ${libPath}/cactus_addReferenceSeq.*
