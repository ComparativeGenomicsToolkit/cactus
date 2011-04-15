rootPath = ../
include ../include.mk

# order is important, libraries first
modules = genemap graphVizPlots mafs stats referenceViewer psls tuning beds referenceUtils 
#beds  
.PHONY: all %.all clean %.clean
#coreModule make file
all : ${libPath}/cactusUtils.a ${modules:%=all.%}

${libPath}/cactusUtils.a : cactusUtils.h cactusUtils.c ${basicLibsDependencies}
	${cxx} ${cflags} -c cactusUtils.c  -I ${libPath}
	ar rc cactusUtils.a *.o
	ranlib cactusUtils.a
	rm *.o
	mv cactusUtils.a ${libPath}/
	cp cactusUtils.h ${libPath}/

all.%:
	cd $* && make all

clean:  ${modules:%=clean.%} clean.cactusUtils

clean.%:
	cd $* && make clean

clean.cactusUtils:
	rm -f ${libPath}/cactusUtils.a
