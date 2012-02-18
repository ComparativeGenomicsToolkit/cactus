# order is important, libraries first
modules = externalTools api pinchGraphs blastAlignment core setup baseAlignment normalisation matching reference phylogeny faces check pipeline progressive preprocessor mafs

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