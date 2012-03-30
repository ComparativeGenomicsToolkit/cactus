# order is important, libraries first
modules = blossom matchGraph
#pulldown
.PHONY: all %.all clean %.clean

all : ${modules:%=all.%}

all.%:
	cd $* && ${MAKE} all

clean:  ${modules:%=clean.%}
	rm -rf lib/*.h bin/*.dSYM

clean.%:
	cd $* && ${MAKE} clean

