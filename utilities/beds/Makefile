include ../../../../include.mk
binPath = ../../../../bin
libPath = ../../../../lib

cflags += ${tokyoCabinetIncl} -lssl

progs = $(notdir $(wildcard))

targets = ${progs:%=${binPath}/%} ${binPath}/cactus_bedGenerator

all: ${targets}

${binPath}/cactus_bedGenerator: ${libPath}/cactus_addReferenceSeq.o cactus_bedGenerator.c ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I ${libPath} -I ${kentInc} -o ${binPath}/cactus_bedGenerator cactus_bedGenerator.c ${libPath}/cactus_addReferenceSeq.o ${libPath}/cactusLib.a ${libPath}/sonLib.a ${dblibs}

clean :
	rm -f *.o
	rm -f ${binPath}/cactus_bedGenerator
