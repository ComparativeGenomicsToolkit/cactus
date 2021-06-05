rootPath = ..
include ${rootPath}/include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

commonHalLibs = ${LIBDIR}/stReference.a ${LIBDIR}/cactusLib.a
stHalDependencies =  ${commonHalLibs} ${LIBDEPENDS}

LDLIBS := ${commonHalLibs} ${LDLIBS}



all: all_libs all_progs
all_libs: 
all_progs: all_libs
	${MAKE} ${LIBDIR}/stCactusToHal.a ${BINDIR}/cactus_halGeneratorTests

clean : 
	rm -f ${BINDIR}/cactus_halGeneratorTests

${BINDIR}/cactus_halGeneratorTests : ${libTests} ${LIBDIR}/stCactusToHal.a ${stHalDependencies}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -Wno-error -o ${BINDIR}/cactus_halGeneratorTests ${libTests} ${LIBDIR}/stCactusToHal.a ${LDLIBS}

${LIBDIR}/stCactusToHal.a : ${libSources} ${libHeaders}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -c ${libSources}
	${AR} rc stCactusToHal.a *.o
	${RANLIB} stCactusToHal.a
	mv stCactusToHal.a ${LIBDIR}/