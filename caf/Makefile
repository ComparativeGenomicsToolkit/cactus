rootPath = ..
include ${rootPath}/include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

commonCafLibs = ${sonLibDir}/stPinchesAndCacti.a ${sonLibDir}/3EdgeConnected.a ${LIBDIR}/cactusLib.a ${LIBDIR}/stPaf.a
stCafDependencies =  ${commonCafLibs} ${LIBDEPENDS}
stCafLibs = ${commonCafLibs} ${LDLIBS}

all: all_libs all_progs
all_libs: ${LIBDIR}/stCaf.a
all_progs: all_libs
	${MAKE} ${BINDIR}/stCafTests

${LIBDIR}/stCaf.a : ${libSources} ${libHeaders}  ${stCafDependencies}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -c ${libSources}
	${AR} rc stCaf.a *.o
	${RANLIB} stCaf.a 
	mv stCaf.a ${LIBDIR}/

${BINDIR}/stCafTests : ${libTests} ${LIBDIR}/stCaf.a ${stCafDependencies}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/stCafTests ${libTests} ${libSources} ${LIBDIR}/stCaf.a ${stCafLibs} ${LDLIBS}

clean : 
	rm -f *.o
	rm -f ${LIBDIR}/stCaf.a ${BINDIR}/stCafTests

