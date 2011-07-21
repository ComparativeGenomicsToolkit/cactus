rootPath = ../
include ../include.mk

cflags += ${tokyoCabinetIncl}

all : ${binPath}/cactus_aligner.py ${binPath}/cactus_alignerTestAligner.py ${binPath}/cactus_batch.py

${binPath}/cactus_aligner : *.c ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_aligner cactus_aligner.c ${libPath}/cactusLib.a ${basicLibs}

${binPath}/cactus_aligner.py : ${binPath}/cactus_aligner
	cp cactus_aligner.py ${binPath}/cactus_aligner.py
	chmod +x ${binPath}/cactus_aligner.py

${binPath}/cactus_batch_chunkSequences : *.c ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_batch_chunkSequences cactus_batch_chunkSequences.c ${libPath}/cactusLib.a ${basicLibs}

${binPath}/cactus_batch_mergeChunks : *.c ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_batch_mergeChunks cactus_batch_mergeChunks.c ${libPath}/cactusLib.a ${basicLibs}

${binPath}/cactus_batch_convertCoordinates : *.c ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_batch_convertCoordinates cactus_batch_convertCoordinates.c ${libPath}/cactusLib.a ${basicLibs}

${binPath}/cactus_batch.py : ${binPath}/cactus_batch_convertCoordinates ${binPath}/cactus_batch_chunkSequences  ${binPath}/cactus_batch_mergeChunks
	cp cactus_batch.py ${binPath}/cactus_batch.py
	chmod +x ${binPath}/cactus_batch.py

${binPath}/cactus_alignerTestAligner.py : cactus_alignerTestAligner.py
	cp cactus_alignerTestAligner.py ${binPath}/cactus_alignerTestAligner.py
	chmod +x ${binPath}/cactus_alignerTestAligner.py

${binPath}/cactus_lastzRepeatMask.py : cactus_lastzRepeatMask.py
	cp cactus_lastzRepeatMask.py ${binPath}/cactus_lastzRepeatMask.py
	chmod +x ${binPath}/cactus_lastzRepeatMask.py

clean : 
	rm -f *.o
	rm -f ${binPath}/cactus_aligner.py ${binPath}/cactus_aligner ${binPath}/cactus_alignerTestAligner.py ${binPath}/cactus_batch.py ${binPath}/cactus_batch_chunkSequences ${binPath}/cactus_batch_mergeChunks ${binPath}/cactus_batch_convertCoordinates ${binPath}/cactus_${binPath}/cactus_batch_convertCoordinates