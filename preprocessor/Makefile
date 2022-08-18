rootPath = ..
include ${rootPath}/include.mk

CFLAGS += ${hiredisIncl}

all: all_libs all_progs
all_libs: 
all_progs: all_libs
	${MAKE} ${BINDIR}/cactus_analyseAssembly ${BINDIR}/cactus_makeAlphaNumericHeaders.py ${BINDIR}/cactus_filterSmallFastaSequences.py ${BINDIR}/cactus_softmask2hardmask ${BINDIR}/cactus_sanitizeFastaHeaders
	cd lastzRepeatMasking && ${MAKE} all

${BINDIR}/cactus_filterSmallFastaSequences.py : cactus_filterSmallFastaSequences.py
	cp cactus_filterSmallFastaSequences.py ${BINDIR}/cactus_filterSmallFastaSequences.py
	chmod +x ${BINDIR}/cactus_filterSmallFastaSequences.py

${BINDIR}/cactus_makeAlphaNumericHeaders.py : cactus_makeAlphaNumericHeaders.py
	cp cactus_makeAlphaNumericHeaders.py ${BINDIR}/cactus_makeAlphaNumericHeaders.py
	chmod +x ${BINDIR}/cactus_makeAlphaNumericHeaders.py

${BINDIR}/cactus_analyseAssembly : cactus_analyseAssembly.c ${LIBDEPENDS} ${LIBDIR}/cactusLib.a
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/cactus_analyseAssembly cactus_analyseAssembly.c ${LIBDIR}/cactusLib.a ${LDLIBS}

${BINDIR}/cactus_softmask2hardmask : cactus_softmask2hardmask.c ${LIBDEPENDS} ${LIBDIR}/cactusLib.a
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/cactus_softmask2hardmask cactus_softmask2hardmask.c ${LIBDIR}/cactusLib.a ${LDLIBS}

${BINDIR}/cactus_sanitizeFastaHeaders : cactus_sanitizeFastaHeaders.c ${LIBDEPENDS} ${LIBDIR}/cactusLib.a
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/cactus_sanitizeFastaHeaders cactus_sanitizeFastaHeaders.c ${LIBDIR}/cactusLib.a ${LDLIBS}

clean : 
	rm -f *.o
	rm -f ${BINDIR}/cactus_analyseAssembly ${BINDIR}/cactus_softmask2hardmask
	cd lastzRepeatMasking && ${MAKE} clean
