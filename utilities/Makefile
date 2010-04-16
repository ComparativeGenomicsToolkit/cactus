# order is important, libraries first
modules = graphVizPlots mAFs stats
.PHONY: all %.all clean %.clean
#coreModule make file
all : ${modules:%=all.%}

all.%:
	cd $* && make all

clean:  ${modules:%=clean.%}

clean.%:
	cd $* && make clean