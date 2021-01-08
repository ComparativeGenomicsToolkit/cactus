rootPath = ..
include ${rootPath}/include.mk


libSources = impl/cactus*.c
libHeaders = inc/cactus*.h
libInternalHeaders = impl/*.h
libTests = tests/cactus*.c
libTestsHeaders = tests/cactus*.h

CPPFLAGS += -I{tests} -Iimpl

all: all_libs all_progs
all_libs: ${LIBDIR}/cactusLib.a
all_progs: all_libs
	${MAKE} ${BINDIR}/cactusAPITests 

clean : 
	rm -f ${LIBDIR}/cactusLib.a ${LIBDIR}/cactus*.h ${BINDIR}/cactusAPITests *.o

${LIBDIR}/cactusLib.a : ${libSources} ${libHeaders} ${libInternalHeaders}
	${CC} ${CPPFLAGS} ${CFLAGS} -I inc -I ${LIBDIR}/ -c ${libSources}
	${AR} rc cactusLib.a *.o
	${RANLIB} cactusLib.a 
	mv cactusLib.a ${LIBDIR}/

${BINDIR}/cactusAPITests : ${libTests} ${libTestsHeaders} ${libSources} ${libHeaders} ${libInternalHeaders} tests/allTests.c ${LIBDIR}/cactusLib.a ${LIBDEPENDS}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/cactusAPITests tests/allTests.c ${libTests} ${LIBDIR}/cactusLib.a ${LDLIBS}
