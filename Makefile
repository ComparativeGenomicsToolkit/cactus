# order is important, libraries first
modules = externalTools api stPinchGraphs stCactusGraphs pinchGraphs blastAlignment core setup baseAlignment normalisation reference phylogeny faces check pipeline progressive preprocessor hal

.PHONY: all %.all clean %.clean

all : ${modules:%=all.%}

all.%:
	cd $* && make all

clean:  ${modules:%=clean.%}
	rm -rf lib/*.h bin/*.dSYM

clean.%:
	cd $* && make clean
	
test: all
	python allTests.py