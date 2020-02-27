SHELL = /bin/bash

# if include.local.mk exists, include it first to set various options
# it shuld not be checked in
includeLocal = ${rootPath}/include.local.mk
ifneq ($(wildcard ${includeLocal}),)
   include ${includeLocal}
endif

# special handling to get C++ ABI right on UCSC Centos 6 servers
ifeq (${CXX_ABI_DEF},)
ifneq ($(wildcard /etc/redhat-release),)
ifeq ($(shell hostname -d), gi.ucsc.edu)
    export CXX_ABI_DEF = -D_GLIBCXX_USE_CXX11_ABI=0
endif
endif
endif


#Location of sonLib
binPath=${rootPath}/bin/
libPath=${rootPath}/lib/

#Modify this variable to set the location of sonLib
sonLibRootPath ?= ${rootPath}/submodules/sonLib
sonLibPath=${sonLibRootPath}/lib

include ${sonLibRootPath}/include.mk

dataSetsPath=/Users/benedictpaten/Dropbox/Documents/work/myPapers/genomeCactusPaper/dataSets

inclDirs = api/inc bar/inc caf/inc hal/inc reference/inc submodules/sonLib/C/inc \
	blastLib submodules/sonLib/externalTools/cutest submodules/pinchesAndCacti/inc \
	submodules/matchingAndOrdering/inc submodules/cPecan/inc

cflags += ${inclDirs:%=-I${rootPath}/%}
basicLibs = ${sonLibPath}/sonLib.a ${sonLibPath}/cuTest.a ${dblibs}
basicLibsDependencies = ${sonLibPath}/sonLib.a ${sonLibPath}/cuTest.a

kyotoTycoonIncl=-I${rootPath}/include -DHAVE_KYOTO_TYCOON=1
kyotoTycoonLib=-L${rootPath}/lib -Wl,-rpath,${rootPath}/lib -lkyototycoon -lkyotocabinet -lz -lbz2 -lpthread -lm -lstdc++
