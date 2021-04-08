rootPath = ..
include ${rootPath}/include.mk

CFLAGS += ${hiredisIncl}

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

commonReferenceLibs = ${sonLibDir}/matchingAndOrdering.a ${LIBDIR}/cactusLib.a 
stReferenceDependencies =  ${commonReferenceLibs} ${LIBDEPENDS}
stReferenceLibs = ${commonReferenceLibs} ${LDLIBS}

all: all_libs all_progs
all_libs: ${LIBDIR}/stReference.a
all_progs: all_libs
	${MAKE} ${BINDIR}/referenceTests

${BINDIR}/referenceTests : ${libTests} ${libSources} ${libHeaders} ${stReferenceDependencies}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -I${LIBDIR} -o ${BINDIR}/referenceTests ${libTests} ${libSources} ${stReferenceLibs}

${LIBDIR}/stReference.a : ${libSources} ${libHeaders} ${stReferenceDependencies}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -c ${libSources}
	${AR} rc stReference.a *.o
	${RANLIB} stReference.a 
	mv stReference.a ${LIBDIR}/

clean : 
	rm -f *.o
	rm -f ${LIBDIR}/stReference.a ${BINDIR}/referenceTests

