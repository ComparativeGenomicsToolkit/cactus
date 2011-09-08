rootPath = ../
include ../include.mk

cflags += ${tokyoCabinetIncl}

all : ${binPath}/cactus_preprocessor.py ${binPath}/cactus_lastzRepeatMask.py ${binPath}/cactus_addFastaHeaderDots.py ${binPath}/cactus_filterSmallFastaSequences.py

${binPath}/cactus_preprocessor.py : cactus_preprocessor.py  ${binPath}/cactus_batch_mergeChunks
	cp cactus_preprocessor.py ${binPath}/cactus_preprocessor.py
	chmod +x ${binPath}/cactus_preprocessor.py

${binPath}/cactus_batch_mergeChunks : *.c ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_batch_mergeChunks cactus_batch_mergeChunks.c ${libPath}/cactusLib.a ${basicLibs}

${binPath}/cactus_lastzRepeatMask.py : cactus_lastzRepeatMask.py
	cp cactus_lastzRepeatMask.py ${binPath}/cactus_lastzRepeatMask.py
	chmod +x ${binPath}/cactus_lastzRepeatMask.py

${binPath}/cactus_addFastaHeaderDots.py : cactus_addFastaHeaderDots.py
	cp cactus_addFastaHeaderDots.py ${binPath}/cactus_addFastaHeaderDots.py
	chmod +x ${binPath}/cactus_addFastaHeaderDots.py
	
${binPath}/cactus_filterSmallFastaSequences.py : cactus_filterSmallFastaSequences.py
	cp cactus_filterSmallFastaSequences.py ${binPath}/cactus_filterSmallFastaSequences.py
	chmod +x ${binPath}/cactus_filterSmallFastaSequences.py
	
clean : 
	rm -f *.o
	rm -f ${binPath}/cactus_preprocessor.py ${binPath}/cactus_lastzRepeatMask.py ${binPath}/cactus_batch_mergeChunks ${binPath}/cactus_addFastaHeaderDots.py ${binPath}/cactus_filterSmallFastaSequences.py