# order is important, libraries first
modules = api threeEdgeConnected pinchGraphs core setup blastAlignment baseAlignment normalisation reference phylogeny faces check pipeline utilities pulldown
.PHONY: all %.all clean %.clean

all : ${modules:%=all.%}

all.%:
	cd $* && make all

clean:  ${modules:%=clean.%}

clean.%:
	cd $* && make clean


