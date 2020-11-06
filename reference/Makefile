rootPath = ..
include ${rootPath}/include.mk

CFLAGS += ${tokyoCabinetIncl} ${hiredisIncl}

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

commonReferenceLibs = ${sonLibDir}/matchingAndOrdering.a ${LIBDIR}/cactusLib.a 
stReferenceDependencies =  ${commonReferenceLibs} ${LIBDEPENDS}
stReferenceLibs = ${commonReferenceLibs} ${LDLIBS}

all: all_libs all_progs
all_libs: ${LIBDIR}/stReference.a
all_progs: all_libs
	${MAKE} ${BINDIR}/cactus_reference ${BINDIR}/cactus_addReferenceCoordinates ${BINDIR}/referenceTests ${BINDIR}/cactus_getReferenceSeq

${BINDIR}/cactus_reference : cactus_reference.c ${libSources} ${libHeaders} ${stReferenceDependencies}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/cactus_reference cactus_reference.c ${libSources} ${stReferenceLibs} ${LDLIBS}

${BINDIR}/cactus_addReferenceCoordinates : cactus_addReferenceCoordinates.c ${libSources} ${libHeaders} ${stReferenceDependencies}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/cactus_addReferenceCoordinates cactus_addReferenceCoordinates.c ${libSources} ${stReferenceLibs} ${LDLIBS}

${BINDIR}/cactus_getReferenceSeq: cactus_getReferenceSeq.c ${stReferenceDependencies}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/cactus_getReferenceSeq cactus_getReferenceSeq.c ${stReferenceLibs} ${LDLIBS}

${BINDIR}/referenceTests : ${libTests} ${libSources} ${libHeaders} ${stReferenceDependencies}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -I${LIBDIR} -o ${BINDIR}/referenceTests ${libTests} ${libSources} ${stReferenceLibs}

${LIBDIR}/stReference.a : ${libSources} ${libHeaders}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -c ${libSources}
	${AR} rc stReference.a *.o
	${RANLIB} stReference.a 
	mv stReference.a ${LIBDIR}/

clean : 
	rm -f *.o
	rm -f ${LIBDIR}/stReference.a ${BINDIR}/cactus_reference ${BINDIR}/referenceTests ${BINDIR}/cactus_addReferenceCoordinates ${BINDIR}/cactus_getReferenceSeq
