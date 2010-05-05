include ../../../../include.mk
binPath = ../../../../bin
libPath = ../../../../lib

cflags += ${tokyoCabinetIncl} -lssl

progs = $(notdir $(wildcard))

targets = ${progs:%=${binPath}/%} 
ifneq ($(wildcard ${kentLibWeb}),)
  targets += ${binPath}/cactus_bedGenerator
  # FIXME: tmp hack, both sonlib commonC and kent libraries have logInfo and logDebug functions
  # kent ones are not actually used, but get dragged in.  As a tmp work around, we force
  # the sonlib ones to be brought in first by the linker.
  jkLibHack = -Wl,-ulogDebug ${libPath}/sonLib.a 
endif

all: ${targets}

${binPath}/cactus_bedGenerator : cactus_bedGenerator.c ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I ${libPath} -I ${kentInc} -o ${binPath}/cactus_bedGenerator cactus_bedGenerator.c ${jkLibHack} ${kentLibWeb} ${libPath}/cactusLib.a ${libPath}/sonLib.a ${tokyoCabinetLib}

clean :
	rm -f *.o
	rm -f ${binPath}/cactus_bedGenerator
