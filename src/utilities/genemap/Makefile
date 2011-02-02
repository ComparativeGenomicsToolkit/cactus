rootPath = ../../../
include ../../../include.mk
libPath = ../../../lib

progs = $(notdir $(wildcard  geneMap.py))

targets = ${progs:%=${binPath}/%}

all: ${targets} ${binPath}/cactus_geneMap

${binPath}/%: %
	@mkdir -p $(dir $@)
	cp -f $< $@
	chmod 775 $@

${binPath}/cactus_geneMap : cactus_geneMap.c ${basicLibsDependencies}
	${cxx} ${cflags} -I ${libPath} cactus_geneMap.c -lpthread -lz -o ${binPath}/cactus_geneMap ${libPath}/cactusLib.a ${basicLibs}

clean :
	rm -f ${binPath}/cactus_geneMap
	rm -f ${progs:%=${binPath}/%}
