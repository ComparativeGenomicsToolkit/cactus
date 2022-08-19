rootPath = ..
include ${rootPath}/include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

stPafDependencies = ${LIBDEPENDS}
stPafLibs = ${LDLIBS}

all: all_libs all_progs
all_libs: ${LIBDIR}/stPaf.a
all_progs: all_libs ${BINDIR}/stPafTests ${BINDIR}/paf_chain ${BINDIR}/paf_invert ${BINDIR}/paf_tile ${BINDIR}/paf_trim ${BINDIR}/paf_shatter ${BINDIR}/paf_view ${BINDIR}/paf_dechunk ${BINDIR}/paf_dedupe ${BINDIR}/paf_to_bed ${BINDIR}/paf_upconvert

${LIBDIR}/stPaf.a : ${libSources} ${libHeaders}  ${stPafDependencies}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -c ${libSources}
	${AR} rc stPaf.a *.o
	${RANLIB} stPaf.a
	mv stPaf.a ${LIBDIR}/

${BINDIR}/paf_dechunk : paf_dechunk.c ${LIBDEPENDS} ${commonPafLibs} ${libSources} ${libHeaders}
	${CC} ${CPPFLAGS} ${CFLAGS} -o ${BINDIR}/paf_dechunk paf_dechunk.c ${libSources} ${commonPafLibs} ${LDLIBS}

${BINDIR}/paf_chain : paf_chain.c ${LIBDEPENDS} ${commonPafLibs} ${libSources} ${libHeaders}
	${CC} ${CPPFLAGS} ${CFLAGS} -o ${BINDIR}/paf_chain paf_chain.c ${libSources} ${commonPafLibs} ${LDLIBS}

${BINDIR}/paf_invert : paf_invert.c ${LIBDEPENDS} ${commonPafLibs} ${libSources} ${libHeaders}
	${CC} ${CPPFLAGS} ${CFLAGS} -o ${BINDIR}/paf_invert paf_invert.c ${libSources} ${commonPafLibs} ${LDLIBS}

${BINDIR}/paf_tile : paf_tile.c ${LIBDEPENDS} ${commonPafLibs} ${libSources} ${libHeaders}
	${CC} ${CPPFLAGS} ${CFLAGS} -o ${BINDIR}/paf_tile paf_tile.c ${libSources} ${commonPafLibs} ${LDLIBS}

${BINDIR}/paf_trim : paf_trim.c ${LIBDEPENDS} ${commonPafLibs} ${libSources} ${libHeaders}
	${CC} ${CPPFLAGS} ${CFLAGS} -o ${BINDIR}/paf_trim paf_trim.c ${libSources} ${commonPafLibs} ${LDLIBS}

${BINDIR}/paf_shatter : paf_shatter.c ${LIBDEPENDS} ${commonPafLibs} ${libSources} ${libHeaders}
	${CC} ${CPPFLAGS} ${CFLAGS} -o ${BINDIR}/paf_shatter paf_shatter.c ${libSources} ${commonPafLibs} ${LDLIBS}

${BINDIR}/paf_view : paf_view.c ${LIBDEPENDS} ${commonPafLibs} ${libSources} ${libHeaders}
	${CC} ${CPPFLAGS} ${CFLAGS} -o ${BINDIR}/paf_view paf_view.c ${libSources} ${commonPafLibs} ${LDLIBS}

${BINDIR}/paf_dedupe : paf_dedupe.c ${LIBDEPENDS} ${commonPafLibs} ${libSources} ${libHeaders}
	${CC} ${CPPFLAGS} ${CFLAGS} -o ${BINDIR}/paf_dedupe paf_dedupe.c ${libSources} ${commonPafLibs} ${LDLIBS}

${BINDIR}/paf_to_bed : paf_to_bed.c ${LIBDEPENDS} ${commonPafLibs} ${libSources} ${libHeaders}
	${CC} ${CPPFLAGS} ${CFLAGS} -o ${BINDIR}/paf_to_bed paf_to_bed.c ${libSources} ${commonPafLibs} ${LDLIBS}

${BINDIR}/paf_upconvert : paf_upconvert.c ${LIBDEPENDS} ${commonPafLibs} ${libSources} ${libHeaders}
	${CC} ${CPPFLAGS} ${CFLAGS} -o ${BINDIR}/paf_upconvert paf_upconvert.c ${libSources} ${commonPafLibs} ${LDLIBS}

${BINDIR}/stPafTests : ${libTests} ${LIBDIR}/stPaf.a ${stPafDependencies}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/stPafTests ${libTests} ${libSources} ${LIBDIR}/stPaf.a ${stCafLibs} ${LDLIBS}

clean : 
	rm -f *.o
	rm -f ${LIBDIR}/stPaf.a ${BINDIR}/stPafTests ${BINDIR}/paf_chain ${BINDIR}/paf_invert ${BINDIR}/paf_tile ${BINDIR}/paf_trim ${BINDIR}/paf_shatter ${BINDIR}/paf_view ${BINDIR}/paf_dechunk ${BINDIR}/paf_dedupe ${BINDIR}/paf_to_bed ${BINDIR}/paf_upconvert
