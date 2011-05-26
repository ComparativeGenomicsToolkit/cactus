# order is important, libraries first
modules = externalTools api pinchGraphs core setup blastAlignment baseAlignment normalisation reference phylogeny faces check pipeline 

.PHONY: all %.all clean %.clean

all : ${modules:%=all.%}

all.%:
	cd $* && ${MAKE} all

clean:  ${modules:%=clean.%}
	rm -rf lib/*.h bin/*.dSYM

clean.%:
	cd $* && ${MAKE} clean
	
test: all
	python allTests.py

