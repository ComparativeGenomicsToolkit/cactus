include ../../../../include.mk
binPath = ../../../../bin
libPath = ../../../../lib
refDir = ../referenceUtils

cflags += ${tokyoCabinetIncl} -lssl

progs = $(notdir $(wildcard))

targets = ${progs:%=${binPath}/%} ${binPath}/cactus_bedGenerator

all: ${targets}

${libPath}/cactus_addReferenceSeq.h: ${refDir}/cactus_addReferenceSeq.h
	cp ${refDir}/cactus_addReferenceSeq.h ${libPath}

${libPath}/cactus_addReferenceSeq.o : ${refDir}/cactus_addReferenceSeq.c
	${cxx} ${cflags} -c ${refDir}/cactus_addReferenceSeq.c  -I ${libPath} -o ${libPath}/cactus_addReferenceSeq.o

${binPath}/cactus_bedGenerator:${libPath}/cactus_addReferenceSeq.h ${libPath}/cactus_addReferenceSeq.o cactus_bedGenerator.c ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I ${libPath} -I ${kentInc} -o ${binPath}/cactus_bedGenerator cactus_bedGenerator.c ${libPath}/cactus_addReferenceSeq.o ${libPath}/cactusLib.a ${libPath}/sonLib.a ${dblibs}

clean :
	rm -f *.o
	rm -f ${binPath}/cactus_bedGenerator
