rootPath = ..
include ${rootPath}/include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

stTafDependencies =  ${LIBDEPENDS}
stTafLibs = ${commonTafLibs} ${LDLIBS}

all: all_libs all_progs
all_libs: ${LIBDIR}/stTaf.a
all_progs: all_libs
	${MAKE} ${BINDIR}/stTafTests

${LIBDIR}/stTaf.a : ${libSources} ${libHeaders}  ${stTafDependencies}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -c ${libSources}
	${AR} rc stTaf.a *.o
	${RANLIB} stTaf.a 
	mv stTaf.a ${LIBDIR}/

${BINDIR}/stTafTests : ${libTests} ${LIBDIR}/stTaf.a ${stTafDependencies}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/stTafTests ${libTests} ${libSources} ${LIBDIR}/stTaf.a ${stTafLibs} ${LDLIBS}

clean : 
	rm -f *.o
	rm -f ${LIBDIR}/stTaf.a ${BINDIR}/stTafTests

