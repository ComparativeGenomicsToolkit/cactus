rootPath = ../../
include ../../include.mk
refDir = ../referenceUtils

progs = $(notdir $(wildcard))

targets = ${progs:%=${binPath}/%} ${binPath}/cactus_bedGenerator

all: ${targets}

#${libPath}/cactus_addReferenceSeq.h: ${refDir}/cactus_addReferenceSeq.h
#	cp ${refDir}/cactus_addReferenceSeq.h ${libPath}

#${libPath}/cactus_addReferenceSeq.o : ${refDir}/cactus_addReferenceSeq.c
#	${cxx} ${cflags} -c ${refDir}/cactus_addReferenceSeq.c  -I ${libPath} -o ${libPath}/cactus_addReferenceSeq.o

${binPath}/cactus_bedGenerator : cactus_bedGenerator.c ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_bedGenerator cactus_bedGenerator.c ${libPath}/cactusUtils.a ${libPath}/cactusLib.a ${basicLibs}

clean :
	rm -f *.o
	rm -f ${binPath}/cactus_bedGenerator
