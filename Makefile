# order is important, libraries first
modules = api 3EdgeConnected pinchGraphs cactusGraphs core setup blastAlignment baseAlignment phylogeny faces reference normalisation pipeline check utilities
.PHONY: all %.all clean %.clean
#baseAlignment
#coreModule make file
all : ${modules:%=all.%}

all.%:
	cd $* && make all

clean:  ${modules:%=clean.%}

clean.%:
	cd $* && make clean


