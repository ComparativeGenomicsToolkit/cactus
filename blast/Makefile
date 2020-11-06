rootPath = ..
include ${rootPath}/include.mk

CFLAGS += ${tokyoCabinetIncl} ${hiredisIncl}

all: all_libs all_progs
all_libs: 
all_progs: all_libs
	${MAKE} ${BINDIR}/cactus_convertAlignmentsToInternalNames ${BINDIR}/cactus_stripUniqueIDs ${BINDIR}/cactus_blast_convertCoordinates ${BINDIR}/cactus_blast_chunkSequences ${BINDIR}/cactus_blast_chunkFlowerSequences ${BINDIR}/cactus_blast_sortAlignments ${BINDIR}/cactus_calculateMappingQualities ${BINDIR}/cactus_mirrorAndOrientAlignments ${BINDIR}/cactus_splitAlignmentOverlaps ${BINDIR}/cactus_coverage

${BINDIR}/cactus_blast_chunkFlowerSequences : *.c ${LIBDIR}/cactusBlastAlignment.a ${LIBDIR}/cactusLib.a ${LIBDEPENDS}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/cactus_blast_chunkFlowerSequences cactus_blast_chunkFlowerSequences.c ${LIBDIR}/cactusBlastAlignment.a ${LIBDIR}/cactusLib.a ${LDLIBS}

${BINDIR}/cactus_blast_chunkSequences : *.c ${LIBDIR}/cactusBlastAlignment.a ${LIBDIR}/cactusLib.a ${LIBDEPENDS}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/cactus_blast_chunkSequences cactus_blast_chunkSequences.c ${LIBDIR}/cactusBlastAlignment.a ${LIBDIR}/cactusLib.a ${LDLIBS}

${BINDIR}/cactus_blast_convertCoordinates : *.c ${LIBDIR}/cactusBlastAlignment.a ${LIBDIR}/cactusLib.a ${LIBDEPENDS}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/cactus_blast_convertCoordinates cactus_blast_convertCoordinates.c ${LIBDIR}/cactusBlastAlignment.a ${LIBDIR}/cactusLib.a ${LDLIBS}

${BINDIR}/cactus_blast_sortAlignments : cactus_blast_sortAlignments.c ${LIBDIR}/stCaf.a ${LIBDIR}/cactusLib.a ${LIBDEPENDS}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/cactus_blast_sortAlignments cactus_blast_sortAlignments.c ${LIBDIR}/stCaf.a ${LIBDIR}/cactusBlastAlignment.a ${LIBDIR}/cactusLib.a ${LDLIBS}

${BINDIR}/cactus_calculateMappingQualities : cactus_calculateMappingQualities.c ${LIBDIR}/stCaf.a ${LIBDIR}/cactusLib.a ${LIBDEPENDS}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/cactus_calculateMappingQualities cactus_calculateMappingQualities.c ${LIBDIR}/stCaf.a ${LIBDIR}/cactusBlastAlignment.a ${LIBDIR}/cactusLib.a ${LDLIBS}

${BINDIR}/cactus_mirrorAndOrientAlignments : cactus_mirrorAndOrientAlignments.c ${LIBDIR}/stCaf.a ${LIBDIR}/cactusLib.a ${LIBDEPENDS}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/cactus_mirrorAndOrientAlignments cactus_mirrorAndOrientAlignments.c ${LIBDIR}/stCaf.a ${LIBDIR}/cactusBlastAlignment.a ${LIBDIR}/cactusLib.a ${LDLIBS}

${BINDIR}/cactus_splitAlignmentOverlaps : cactus_splitAlignmentOverlaps.c ${LIBDIR}/stCaf.a ${LIBDIR}/cactusLib.a ${LIBDEPENDS}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/cactus_splitAlignmentOverlaps cactus_splitAlignmentOverlaps.c ${LIBDIR}/stCaf.a ${LIBDIR}/cactusBlastAlignment.a ${LIBDIR}/cactusLib.a ${LDLIBS}

${BINDIR}/cactus_coverage : cactus_coverage.c ${LIBDEPENDS}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/cactus_coverage cactus_coverage.c ${LDLIBS}

${BINDIR}/cactus_convertAlignmentsToInternalNames : cactus_convertAlignmentsToInternalNames.c ${LIBDIR}/cactusLib.a
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/cactus_convertAlignmentsToInternalNames cactus_convertAlignmentsToInternalNames.c ${LIBDIR}/cactusLib.a ${LDLIBS}

${BINDIR}/cactus_stripUniqueIDs : cactus_stripUniqueIDs.c ${LIBDIR}/cactusLib.a
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/cactus_stripUniqueIDs cactus_stripUniqueIDs.c ${LIBDIR}/cactusLib.a ${LDLIBS}

clean : 
	rm -f *.o
	rm -f ${LIBDIR}/cactusBlastAlignment.a ${BINDIR}/cactus_blast.py ${BINDIR}/cactus_blast_chunkSequences ${BINDIR}/cactus_blast_sortAlignments ${BINDIR}/cactus_calculateMappingQualities ${BINDIR}/cactus_mirrorAndOrientAlignments ${BINDIR}/cactus_splitAlignmentOverlaps ${BINDIR}/cactus_blast_chunkFlowerSequences ${BINDIR}/cactus_blast_convertCoordinates 
