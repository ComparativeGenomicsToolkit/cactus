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
	${MAKE} ${BINDIR}/cactus_halGenerator ${BINDIR}/cactus_halGeneratorTests ${BINDIR}/cactus_fastaGenerator

clean : 
	rm -f ${BINDIR}/cactus_halGenerator ${BINDIR}/cactus_halGeneratorTests 

${BINDIR}/cactus_halGenerator : cactus_halGenerator.c ${libTests} ${libSources} ${libHeaders} ${stHalDependencies}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/cactus_halGenerator cactus_halGenerator.c ${libSources} ${LDLIBS}

${BINDIR}/cactus_fastaGenerator : cactus_fastaGenerator.c ${libTests} ${libSources} ${libHeaders} ${stHalDependencies}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/cactus_fastaGenerator cactus_fastaGenerator.c ${libSources} ${LDLIBS}

${BINDIR}/cactus_halGeneratorTests : ${libTests} ${libSources} ${libHeaders} ${stHalDependencies}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -Wno-error -o ${BINDIR}/cactus_halGeneratorTests ${libTests} ${libSources} ${LDLIBS}
