rootPath = ../..
include ${rootPath}/include.mk

CFLAGS += ${tokyoCabinetIncl} ${hiredisIncl}

all :  ${BINDIR}/cactus_fasta_fragments.py ${BINDIR}/cactus_fasta_softmask_intervals.py ${BINDIR}/cactus_covered_intervals

${BINDIR}/cactus_covered_intervals : *.c  ${LIBDEPENDS}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/cactus_covered_intervals cactus_covered_intervals.c ${LDLIBS}

${BINDIR}/cactus_fasta_fragments.py : cactus_fasta_fragments.py
	cp cactus_fasta_fragments.py ${BINDIR}/cactus_fasta_fragments.py
	chmod +x ${BINDIR}/cactus_fasta_fragments.py

${BINDIR}/cactus_fasta_softmask_intervals.py : cactus_fasta_softmask_intervals.py
	cp cactus_fasta_softmask_intervals.py ${BINDIR}/cactus_fasta_softmask_intervals.py
	chmod +x ${BINDIR}/cactus_fasta_softmask_intervals.py

clean : 
	rm -f *.o
	rm -f  ${BINDIR}/cactus_lastzRepeatMask.py ${BINDIR}/cactus_fasta_fragments.py  ${BINDIR}/cactus_fasta_softmask_intervals.py ${BINDIR}/cactus_covered_intervals
