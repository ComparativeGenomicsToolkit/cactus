include ../../../../include.mk
binPath = ../../../../bin
libPath = ../../../../lib

cflags += ${tokyoCabinetIncl}

progs = $(notdir $(wildcard cactus_mafToReferenceSeq.py))
targets = ${progs:%=${binPath}/%}

all :  ${targets} ${binPath}/cactus_addReferenceSeq

${binPath}/%: %
	@mkdir -p $(dir $@)
	cp -f $< $@
	chmod 775 $@

${binPath}/cactus_addReferenceSeq : *.c ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_addReferenceSeq cactus_addReferenceSeq.c  ${libPath}/sonLib.a ${libPath}/cactusLib.a ${libPath}/sonLib.a ${dblibs}

clean :
	rm -f *.o
	rm -f ${binPath}/cactus_addReferenceSeq
	rm -f ${progs:%=${binPath}/%}
