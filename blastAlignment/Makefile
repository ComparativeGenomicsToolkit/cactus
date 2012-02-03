rootPath = ../
include ../include.mk

cflags += ${tokyoCabinetIncl}

libSources = blastAlignmentLib.c
libHeaders = blastAlignmentLib.h

all : ${binPath}/cactus_aligner.py ${binPath}/cactus_alignerTestAligner.py ${binPath}/cactus_batch.py

${libPath}/cactusBlastAlignment.a : ${libSources} ${libHeaders} ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I inc -I ${libPath}/ -c ${libSources}
	ar rc cactusBlastAlignment.a *.o
	ranlib cactusBlastAlignment.a 
	rm *.o
	mv cactusBlastAlignment.a ${libPath}/
	cp ${libHeaders} ${libPath}/

${binPath}/cactus_aligner : *.c ${libPath}/cactusBlastAlignment.a ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_aligner cactus_aligner.c ${libPath}/cactusBlastAlignment.a ${libPath}/cactusLib.a ${basicLibs}

${binPath}/cactus_aligner.py : ${binPath}/cactus_aligner
	cp cactus_aligner.py ${binPath}/cactus_aligner.py
	chmod +x ${binPath}/cactus_aligner.py

${binPath}/cactus_batch_chunkSequences : *.c ${libPath}/cactusBlastAlignment.a ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_batch_chunkSequences cactus_batch_chunkSequences.c ${libPath}/cactusBlastAlignment.a ${libPath}/cactusLib.a ${basicLibs}

${binPath}/cactus_batch_convertCoordinates : *.c ${libPath}/cactusBlastAlignment.a ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_batch_convertCoordinates cactus_batch_convertCoordinates.c ${libPath}/cactusBlastAlignment.a ${libPath}/cactusLib.a ${basicLibs}

${binPath}/cactus_batch.py : ${binPath}/cactus_batch_convertCoordinates ${binPath}/cactus_batch_chunkSequences
	cp cactus_batch.py ${binPath}/cactus_batch.py
	chmod +x ${binPath}/cactus_batch.py

${binPath}/cactus_alignerTestAligner.py : cactus_alignerTestAligner.py
	cp cactus_alignerTestAligner.py ${binPath}/cactus_alignerTestAligner.py
	chmod +x ${binPath}/cactus_alignerTestAligner.py

clean : 
	rm -f *.o
	rm -f ${binPath}/cactus_aligner.py ${binPath}/cactus_aligner ${binPath}/cactus_alignerTestAligner.py ${binPath}/cactus_batch.py ${binPath}/cactus_batch_chunkSequences ${binPath}/cactus_batch_convertCoordinates