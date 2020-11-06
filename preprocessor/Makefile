rootPath = ..
include ${rootPath}/include.mk

CFLAGS += ${tokyoCabinetIncl} ${hiredisIncl}

all: all_libs all_progs
all_libs: 
all_progs: all_libs
	${MAKE} ${BINDIR}/cactus_analyseAssembly ${BINDIR}/cactus_batch_mergeChunks ${BINDIR}/cactus_makeAlphaNumericHeaders.py ${BINDIR}/cactus_filterSmallFastaSequences.py
	cd lastzRepeatMasking && ${MAKE} all

${BINDIR}/cactus_filterSmallFastaSequences.py : cactus_filterSmallFastaSequences.py
	cp cactus_filterSmallFastaSequences.py ${BINDIR}/cactus_filterSmallFastaSequences.py
	chmod +x ${BINDIR}/cactus_filterSmallFastaSequences.py

${BINDIR}/cactus_makeAlphaNumericHeaders.py : cactus_makeAlphaNumericHeaders.py
	cp cactus_makeAlphaNumericHeaders.py ${BINDIR}/cactus_makeAlphaNumericHeaders.py
	chmod +x ${BINDIR}/cactus_makeAlphaNumericHeaders.py

${BINDIR}/cactus_analyseAssembly : cactus_analyseAssembly.c ${LIBDEPENDS} ${LIBDIR}/cactusLib.a
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/cactus_analyseAssembly cactus_analyseAssembly.c ${LIBDIR}/cactusLib.a ${LDLIBS}

${BINDIR}/cactus_batch_mergeChunks : *.c ${LIBDIR}/cactusLib.a ${LIBDEPENDS}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/cactus_batch_mergeChunks cactus_batch_mergeChunks.c ${LIBDIR}/cactusLib.a ${LDLIBS}

clean : 
	rm -f *.o
	rm -f ${BINDIR}/cactus_analyseAssembly ${BINDIR}/cactus_batch_mergeChunks
	cd lastzRepeatMasking && ${MAKE} clean
