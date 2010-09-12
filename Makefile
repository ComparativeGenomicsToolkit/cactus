# order is important, libraries first
modules = api 3EdgeConnected pinchGraphs core setup blastAlignment baseAlignment normalisation reference phylogeny faces check pipeline utilities pulldown
.PHONY: all %.all clean %.clean
#baseAlignment
#coreModule make file
all : ${modules:%=all.%}

all.%:
	cd $* && make all

clean:  ${modules:%=clean.%}

clean.%:
	cd $* && make clean


