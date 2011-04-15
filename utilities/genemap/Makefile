rootPath = ../../
include ../../include.mk

progs = $(notdir $(wildcard  geneMap.py))

targets = ${progs:%=${binPath}/%}

all: ${targets} ${binPath}/cactus_genemap ${binPath}/cactus_genemap2

${binPath}/%: %
	@mkdir -p $(dir $@)
	cp -f $< $@
	chmod 775 $@

${binPath}/cactus_genemap : cactus_genemap.c ${basicLibsDependencies}
	${cxx} ${cflags} -I ${libPath} cactus_genemap.c -lpthread -lz -o ${binPath}/cactus_genemap ${libPath}/cactusUtils.a ${libPath}/cactusLib.a ${basicLibs}

${binPath}/cactus_genemap2 : cactus_genemap2.c ${basicLibsDependencies}
	${cxx} ${cflags} -I ${libPath} cactus_genemap2.c -lpthread -lz -o ${binPath}/cactus_genemap2 ${libPath}/cactusUtils.a ${libPath}/cactusLib.a ${basicLibs}

clean :
	rm -f ${binPath}/cactus_geneMap
	rm -f ${progs:%=${binPath}/%}
