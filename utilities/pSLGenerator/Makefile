include ../../../include.mk
binPath = ../../../bin
libPath = ../../../lib

cflags += ${tokyoCabinetIncl} -lssl

progs = $(notdir $(wildcard eval_mergeMAFComparatorResults.py eval_MAFLength.py eval_maskedFaLength.py ))

targets = ${progs:%=${binPath}/%} ${binPath}/eval_MFAToMAF mAFComparison
ifneq ($(wildcard ${kentLibWeb}),)
  targets += ${binPath}/eval_pslGenerator ${binPath}/eval_pslMapCoor
  # FIXME: tmp hack, both sonlib commonC and kent libraries have logInfo and logDebug functions
  # kent ones are not actually used, but get dragged in.  As a tmp work around, we force
  # the sonlib ones to be brought in first by the linker.
  jkLibHack = -Wl,-ulogDebug ${libPath}/sonLib.a 
endif

all: ${targets}

mAFComparison : mAFComparison/*.c
	cd mAFComparison && make

${binPath}/eval_MFAToMAF : eval_MFAToMAF.c ${libPath}/sonLib.a 
	${cxx} ${cflags} -I ${libPath} ${tokyoCabinetIncl} -o ${binPath}/eval_MFAToMAF eval_MFAToMAF.c  ${libPath}/sonLib.a ${tokyoCabinetLib}

${binPath}/%: %
	@mkdir -p $(dir $@)
	cp -f $< $@
	chmod 775 $@

${binPath}/eval_pslMapCoor : eval_pslMapCoor.c
	${cxx} ${cflags} eval_pslMapCoor.c -lpthread -lz -I ${kentInc} -o ${binPath}/eval_pslMapCoor ${jkLibHack} ${kentLibWeb}

${binPath}/eval_pslGenerator : eval_pslGenerator.c ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I ${libPath} -I ${kentInc} -o ${binPath}/eval_pslGenerator eval_pslGenerator.c ${jkLibHack} ${kentLibWeb} ${libPath}/cactusLib.a ${libPath}/sonLib.a ${tokyoCabinetLib}

clean :
	rm -f *.o
	rm -f ${binPath}/eval_MAFComparator ${binPath}/eval_MFAToMAF
	rm -f ${binPath}/eval_pslGenerator ${binPath}/eval_pslMapCoor
	rm -f ${progs:%=${binPath}/%}
	cd mAFComparison && make clean
