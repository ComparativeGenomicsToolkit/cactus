# order is important, libraries first
modules = genemap graphVizPlots mafs stats referenceViewer psls tuning beds referenceUtils 
#beds  
.PHONY: all %.all clean %.clean
#coreModule make file
all : ${modules:%=all.%}

all.%:
	cd $* && ${MAKE} all

clean:  ${modules:%=clean.%}

clean.%:
	cd $* && ${MAKE} clean
