rootPath = ..
include ${rootPath}/include.mk

cflags += ${tokyoCabinetIncl}

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

commonReferenceLibs = ${sonLibPath}/matchingAndOrdering.a ${libPath}/cactusLib.a 
stReferenceDependencies =  ${commonReferenceLibs} ${basicLibsDependencies}
stReferenceLibs = ${commonReferenceLibs} ${basicLibs}

all: all_libs all_progs
all_libs: ${libPath}/stReference.a
all_progs: all_libs
	${MAKE} ${binPath}/cactus_reference ${binPath}/cactus_addReferenceCoordinates ${binPath}/referenceTests ${binPath}/cactus_getReferenceSeq

${binPath}/cactus_reference : cactus_reference.c ${libSources} ${libHeaders} ${stReferenceDependencies}
	${cxx} ${cflags} -I inc -I${libPath} -o ${binPath}/cactus_reference cactus_reference.c ${libSources} ${stReferenceLibs}

${binPath}/cactus_addReferenceCoordinates : cactus_addReferenceCoordinates.c ${libSources} ${libHeaders} ${stReferenceDependencies}
	${cxx} ${cflags} -I inc -I${libPath} -o ${binPath}/cactus_addReferenceCoordinates cactus_addReferenceCoordinates.c ${libSources} ${stReferenceLibs}

${binPath}/cactus_getReferenceSeq: cactus_getReferenceSeq.c ${stReferenceDependencies}
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_getReferenceSeq cactus_getReferenceSeq.c ${stReferenceLibs}

${binPath}/referenceTests : ${libTests} ${libSources} ${libHeaders} ${stReferenceDependencies}
	${cxx} ${cflags} -I inc -I impl -I${libPath} -o ${binPath}/referenceTests ${libTests} ${libSources} ${stReferenceLibs}

${libPath}/stReference.a : ${libSources} ${libHeaders}
	${cxx} ${cflags} -I inc -I ${libPath}/ -c ${libSources}
	ar rc stReference.a *.o
	ranlib stReference.a 
	rm *.o
	mv stReference.a ${libPath}/

clean : 
	rm -f *.o
	rm -f ${libPath}/stReference.a ${binPath}/cactus_reference ${binPath}/referenceTests ${binPath}/cactus_addReferenceCoordinates ${binPath}/cactus_getReferenceSeq
