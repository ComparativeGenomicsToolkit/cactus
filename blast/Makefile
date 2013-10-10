rootPath = ../
include ../include.mk

cflags += ${tokyoCabinetIncl}

all : ${binPath}/cactus_blast.py

${binPath}/cactus_blast_chunkFlowerSequences : *.c ${libPath}/cactusBlastAlignment.a ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_blast_chunkFlowerSequences cactus_blast_chunkFlowerSequences.c ${libPath}/cactusBlastAlignment.a ${libPath}/cactusLib.a ${basicLibs}

${binPath}/cactus_blast_chunkSequences : *.c ${libPath}/cactusBlastAlignment.a ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_blast_chunkSequences cactus_blast_chunkSequences.c ${libPath}/cactusBlastAlignment.a ${libPath}/cactusLib.a ${basicLibs}

${binPath}/cactus_blast_convertCoordinates : *.c ${libPath}/cactusBlastAlignment.a ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_blast_convertCoordinates cactus_blast_convertCoordinates.c ${libPath}/cactusBlastAlignment.a ${libPath}/cactusLib.a ${basicLibs}

${binPath}/cactus_realign : cactus_realign.c ${libPath}/cactusBarLib.a ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I inc -I${libPath} -o ${binPath}/cactus_realign cactus_realign.c ${libPath}/cactusBarLib.a ${libPath}/cactusLib.a ${basicLibs}

${binPath}/cactus_blast.py : cactus_blast.py ${binPath}/cactus_blast_convertCoordinates ${binPath}/cactus_blast_chunkSequences ${binPath}/cactus_blast_chunkFlowerSequences 
	cp cactus_blast.py ${binPath}/cactus_blast.py
	chmod +x ${binPath}/cactus_blast.py

clean : 
	rm -f *.o
	rm -f ${libPath}/cactusBlastAlignment.a ${binPath}/cactus_blast.py ${binPath}/cactus_blast_chunkSequences ${binPath}/cactus_blast_chunkFlowerSequences ${binPath}/cactus_blast_convertCoordinates ${binPath}/cactus_realign