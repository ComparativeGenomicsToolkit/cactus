rootPath = ../../
include ../../include.mk

#progs = $(notdir $(wildcard cactus_tuning.py cactus_tuning_parameters.py cactus_tuning_summary2.py))
progs = $(notdir $(wildcard cactus_tuning.py))

targets = ${progs:%=${binPath}/%}

all: ${targets}

${binPath}/%: %
	@mkdir -p $(dir $@)
	cp -f $< $@
	chmod 775 $@

clean :
	rm -f ${progs:%=${binPath}/%}
