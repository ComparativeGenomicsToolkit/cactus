# order is important, libraries first
modules = externalTools api pinchGraphs core setup blastAlignment baseAlignment normalisation reference phylogeny faces check pipeline 

.PHONY: all %.all clean %.clean

all : ${modules:%=all.%}

all.%:
	cd $* && make all

clean:  ${modules:%=clean.%}
	rm -rf lib/*.h bin/*.dSYM

clean.%:
	cd $* && make clean