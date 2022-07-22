rootPath = ..
include ${rootPath}/include.mk

libTests = tests/*.c

all: all_progs
all_libs:
all_progs: ${BINDIR}/stFastaTests ${BINDIR}/fasta_chunk ${BINDIR}/fasta_merge ${BINDIR}/fasta_extract

${BINDIR}/fasta_chunk : ../fasta/fasta_chunk.c ${LIBDEPENDS}
	${CC} ${CPPFLAGS} ${CFLAGS} -o ${BINDIR}/fasta_chunk fasta_chunk.c ${libSources} ${commonPafLibs} ${LDLIBS}

${BINDIR}/fasta_merge : ../fasta/fasta_merge.c ${LIBDEPENDS}
	${CC} ${CPPFLAGS} ${CFLAGS} -o ${BINDIR}/fasta_merge fasta_merge.c ${libSources} ${commonPafLibs} ${LDLIBS}

${BINDIR}/fasta_extract : ../fasta/fasta_extract.c ${LIBDEPENDS}
	${CC} ${CPPFLAGS} ${CFLAGS} -o ${BINDIR}/fasta_extract fasta_extract.c ${libSources} ${commonPafLibs} ${LDLIBS}

${BINDIR}/stFastaTests : ${libTests} ${LIBDEPENDS}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/stFastaTests ${libTests} ${libSources} ${stCafLibs} ${LDLIBS}

clean : 
	rm -f *.o
	rm -f ${BINDIR}/stFastaTests ${BINDIR}/fasta_chunk ${BINDIR}/fasta_merge ${BINDIR}/fasta_extract
