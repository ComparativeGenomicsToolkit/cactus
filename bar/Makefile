rootPath = ..
include ${rootPath}/include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/adjacencySequencesTest.c tests/allTests.c tests/endAlignerTest.c tests/flowerAlignerTest.c tests/rescueTest.c tests/poaBarTest.c
libRunEndAlignment = tests/runEndAlignment.c

commonBarLibs = ${LIBDIR}/stCaf.a ${LIBDIR}/stPaf.a ${sonLibDir}/stPinchesAndCacti.a ${LIBDIR}/cactusLib.a ${sonLibDir}/cPecanLib.a
# abPOA and minipoa are only used from here, so they go here and not in LIBDEPENDS.  What links
# bar's code is what needs relinking, and that all depends on cactusBarLib.a.  The headers are
# what the submodules publish to include/, which poaBarAligner.c compiles against, and simde is
# there through abpoa.h.  The whole simde directory is copied whenever any of it changes, so the
# directory stands for its headers.
stBarDependencies =  ${commonBarLibs} ${LIBDEPENDS} ${LIBDIR}/libabpoa.a ${minipoaDepends} $(wildcard ${INCLDIR}/*.h ${INCLDIR}/simde)
LDLIBS += ${commonBarLibs} ${sonLibDir}/sonLib.a ${databaseLibs} -lm
# simde (included via abPOA) doesn't compile with --Werror --pedantic
CFLAGS:=$(filter-out --pedantic,$(CFLAGS))

all: all_libs all_progs
all_libs: ${LIBDIR}/cactusBarLib.a
all_progs: all_libs
	${MAKE} ${BINDIR}/cactus_barTests

clean : 
	rm -f ${BINDIR}/cactus_barTests  ${LIBDIR}/cactusBarLib.a *.o

${BINDIR}/cactus_barTests : ${libTests} tests/*.h ${LIBDIR}/cactusBarLib.a ${stBarDependencies}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -Wno-error -o ${BINDIR}/cactus_barTests ${libTests} ${LIBDIR}/cactusBarLib.a ${LDLIBS}

${LIBDIR}/cactusBarLib.a : ${libSources} ${libHeaders} ${stBarDependencies}
# the -Wno-unused-function is required to include abpoa.h with CGL_DEBUG defined
	${CC} ${CPPFLAGS} ${CFLAGS} -c ${libSources} -Wno-unused-function 
	${AR} rc cactusBarLib.a *.o
	${RANLIB} cactusBarLib.a 
	mv cactusBarLib.a ${LIBDIR}/
