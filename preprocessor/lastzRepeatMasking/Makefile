rootPath = ../..
include ${rootPath}/include.mk

cflags += ${tokyoCabinetIncl}

all :  ${binPath}/cactus_fasta_fragments.py ${binPath}/cactus_fasta_softmask_intervals.py ${binPath}/cactus_covered_intervals

${binPath}/cactus_covered_intervals : *.c  ${basicLibsDependencies}
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_covered_intervals cactus_covered_intervals.c  ${basicLibs}

${binPath}/cactus_fasta_fragments.py : cactus_fasta_fragments.py
	cp cactus_fasta_fragments.py ${binPath}/cactus_fasta_fragments.py
	chmod +x ${binPath}/cactus_fasta_fragments.py

${binPath}/cactus_fasta_softmask_intervals.py : cactus_fasta_softmask_intervals.py
	cp cactus_fasta_softmask_intervals.py ${binPath}/cactus_fasta_softmask_intervals.py
	chmod +x ${binPath}/cactus_fasta_softmask_intervals.py

clean : 
	rm -f *.o
	rm -f  ${binPath}/cactus_lastzRepeatMask.py ${binPath}/cactus_fasta_fragments.py  ${binPath}/cactus_fasta_softmask_intervals.py ${binPath}/cactus_covered_intervals
