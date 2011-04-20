rootPath = ../../
include ../../include.mk

progs = $(notdir $(wildcard  genemapCactusVsMultiz.py genemapChain.py genemapMerge.py genemapMakeLatexTab.py cactus_runGenemap.py))

targets = ${progs:%=${binPath}/%}

all: ${targets} ${binPath}/cactus_genemapChain ${binPath}/cactus_genemapHomolog

${binPath}/%: %
	@mkdir -p $(dir $@)
	cp -f $< $@
	chmod 775 $@

${binPath}/cactus_genemapChain : cactus_genemapChain.c ${basicLibsDependencies}
	${cxx} ${cflags} -I ${libPath} cactus_genemapChain.c -lpthread -lz -o ${binPath}/cactus_genemapChain ${libPath}/cactusUtils.a ${libPath}/cactusLib.a ${basicLibs}

${binPath}/cactus_genemapHomolog : cactus_genemapHomolog.c ${basicLibsDependencies}
	${cxx} ${cflags} -I ${libPath} cactus_genemapHomolog.c -lpthread -lz -o ${binPath}/cactus_genemapHomolog ${libPath}/cactusUtils.a ${libPath}/cactusLib.a ${basicLibs}

clean :
	rm -f ${binPath}/cactus_genemap*
	rm -f ${progs:%=${binPath}/%}
