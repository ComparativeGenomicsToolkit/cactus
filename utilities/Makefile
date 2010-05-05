# order is important, libraries first
modules = graphVizPlots mAFs stats referenceViewer pSLGenerator bEDs
.PHONY: all %.all clean %.clean
#coreModule make file
all : ${modules:%=all.%}

all.%:
	cd $* && make all

clean:  ${modules:%=clean.%}

clean.%:
	cd $* && make clean
