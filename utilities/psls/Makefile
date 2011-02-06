rootPath = ../../
include ../../include.mk

progs = $(notdir $(wildcard))

targets = ${progs:%=${binPath}/%} 
ifneq ($(wildcard ${kentLibWeb}),)
  targets += ${binPath}/cactus_pslGenerator
  # FIXME: tmp hack, both sonlib commonC and kent libraries have logInfo and logDebug functions
  # kent ones are not actually used, but get dragged in.  As a tmp work around, we force
  # the sonlib ones to be brought in first by the linker.
  jkLibHack = -Wl,-ulogDebug ${basicLibs}
endif

all: ${targets}

${binPath}/cactus_pslGenerator : cactus_pslGenerator.c ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I ${libPath} -I ${kentInc} -o ${binPath}/cactus_pslGenerator cactus_pslGenerator.c ${jkLibHack} ${kentLibWeb} ${libPath}/cactusLib.a ${basicLibs}

clean :
	rm -f *.o
	rm -f ${binPath}/cactus_pslGenerator
