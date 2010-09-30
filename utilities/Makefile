# order is important, libraries first
modules = graphVizPlots mafs stats referenceViewer psls tuning 
#beds  
.PHONY: all %.all clean %.clean
#coreModule make file
all : ${modules:%=all.%}

all.%:
	cd $* && make all

clean:  ${modules:%=clean.%}

clean.%:
	cd $* && make clean
