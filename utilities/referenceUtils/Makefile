include ../../../../include.mk
binPath = ../../../../bin
libPath = ../../../../lib

cflags += ${tokyoCabinetIncl}

progs = $(notdir $(wildcard cactus_mafToReferenceSeq.py))
targets = ${progs:%=${binPath}/%}

all :  ${targets} cactus_addReferenceSeq.o ${binPath}/cactus_addReferenceSeqTest ${binPath}/cactus_getReferenceSeq

${binPath}/%: %
	@mkdir -p $(dir $@)
	cp -f $< $@
	chmod 775 $@

cactus_addReferenceSeq.o : cactus_addReferenceSeq.c
	${cxx} ${cflags} -c cactus_addReferenceSeq.c  -I ${libPath}
	mv cactus_addReferenceSeq.o ${libPath} 
	#${cxx} ${cflags} -c cactus_addReferenceSeq.c  -I ${libPath} ${libPath}/cactusLib.a ${libPath}/sonLib.a ${dblibs}

${binPath}/cactus_getReferenceSeq : ${libPath}/cactus_addReferenceSeq.o cactus_getReferenceSeq.c ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_getReferenceSeq cactus_getReferenceSeq.c  ${libPath}/cactus_addReferenceSeq.o ${libPath}/sonLib.a ${libPath}/cactusLib.a ${libPath}/sonLib.a ${dblibs}

${binPath}/cactus_addReferenceSeqTest : ${libPath}/cactus_addReferenceSeq.o cactus_addReferenceSeqTest.c ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_addReferenceSeqTest cactus_addReferenceSeqTest.c  ${libPath}/cactus_addReferenceSeq.o ${libPath}/sonLib.a ${libPath}/cactusLib.a ${libPath}/sonLib.a ${dblibs}

clean :
	rm -f *.o
	rm -f ${libPath}/cactus_addReferenceSeq.o
	rm -f ${binPath}/cactus_addReferenceSeqTest
	rm -f ${progs:%=${binPath}/%}
