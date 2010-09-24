include ../../../../include.mk
binPath = ../../../../bin
libPath = ../../../../lib

cflags += ${tokyoCabinetIncl}

progs = $(notdir $(wildcard cactus_mafToReferenceSeq.py))
targets = ${progs:%=${binPath}/%}

all :  ${targets} ${binPath}/cactus_MAFGenerator

${binPath}/%: %
	@mkdir -p $(dir $@)
	cp -f $< $@
	chmod 775 $@

${binPath}/cactus_MAFGenerator : *.c ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_MAFGenerator cactus_MAFGenerator.c  ${libPath}/sonLib.a ${libPath}/cactusLib.a ${libPath}/sonLib.a ${dblibs}

clean :
	rm -f *.o
	rm -f ${binPath}/cactus_MAFGenerator
	rm -f ${progs:%=${binPath}/%}
