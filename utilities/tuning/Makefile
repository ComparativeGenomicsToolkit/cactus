include ../../../../include.mk
binPath = ../../../../bin
libPath = ../../../../lib

progs = $(notdir $(wildcard cactus_tuning.py))

targets = ${progs:%=${binPath}/%}

all: ${targets}

${binPath}/%: %
	@mkdir -p $(dir $@)
	cp -f $< $@
	chmod 775 $@

clean :
	rm -f ${progs:%=${binPath}/%}
