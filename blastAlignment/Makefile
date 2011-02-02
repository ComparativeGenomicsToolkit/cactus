rootPath = ../
include ../include.mk

cflags += ${tokyoCabinetIncl}

all : ${binPath}/cactus_aligner.py ${binPath}/cactus_alignerTestAligner.py ${binPath}/cactus_batch.py

${binPath}/cactus_aligner.py : cactus_aligner.py cactus_aligner.c ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_aligner cactus_aligner.c ${libPath}/cactusLib.a ${basicLibs}
	cp cactus_aligner.py ${binPath}/cactus_aligner.py
	chmod +x ${binPath}/cactus_aligner.py

${binPath}/cactus_batch.py : cactus_batch* ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_batch_chunkSequences cactus_batch_chunkSequences.c ${libPath}/cactusLib.a ${basicLibs}
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_batch_convertCoordinates cactus_batch_convertCoordinates.c ${libPath}/cactusLib.a ${basicLibs}
	cp cactus_batch.py ${binPath}/cactus_batch.py
	chmod +x ${binPath}/cactus_batch.py

${binPath}/cactus_alignerTestAligner.py : cactus_alignerTestAligner.py
	cp cactus_alignerTestAligner.py ${binPath}/cactus_alignerTestAligner.py
	chmod +x ${binPath}/cactus_alignerTestAligner.py

clean : 
	rm -f *.o
	rm -f ${binPath}/cactus_aligner.py ${binPath}/cactus_aligner ${binPath}/cactus_alignerTestAligner.py ${binPath}/cactus_batch.py ${binPath}/cactus_batch_chunkSequences ${binPath}/cactus_batch_convertCoordinates