SHELL = /bin/bash

##
# Users can set CPPFLAGS, CFLAGS, LIBS to reference external packages.  These
# can be set in environment of config.local.mk.  LDLIBS should not be modified,
# as it is not seen my kyoto configure.
##

# if include.local.mk exists, include it first to set various options
# it should not be checked into git
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
BINDIR = ${rootPath}/bin
LIBDIR = ${rootPath}/lib
INCLDIR = ${rootPath}/include

#Modify this variable to set the location of sonLib
sonLibRootDir ?= ${rootPath}/submodules/sonLib
sonLibDir=${sonLibRootDir}/lib

include ${sonLibRootDir}/include.mk

#Turn asserts back on in spite of sonLib
#https://github.com/ComparativeGenomicsToolkit/cactus/issues/235
CFLAGS += -UNDEBUG

dataSetsPath=/Users/benedictpaten/Dropbox/Documents/work/myPapers/genomeCactusPaper/dataSets

inclDirs = api/inc bar/inc caf/inc hal/inc reference/inc submodules/sonLib/C/inc \
	blastLib submodules/sonLib/externalTools/cutest submodules/pinchesAndCacti/inc \
	submodules/matchingAndOrdering/inc submodules/cPecan/inc

kyotoTycoonIncl=-I${rootPath}/include -DHAVE_KYOTO_TYCOON=1
kyotoTycoonLib=-L${rootPath}/lib -Wl,-rpath,${rootPath}/lib -lkyototycoon -lkyotocabinet -lz -lbz2 -lpthread -lm -lstdc++

CPPFLAGS += ${inclDirs:%=-I${rootPath}/%} -I${LIBDIR} ${kyotoTycoonIncl}

# libraries can't be added until they are build, so add as to LDLIBS until needed
cactusLibs = ${LIBDIR}/stCaf.a ${LIBDIR}/stReference.a ${LIBDIR}/cactusBarLib.a ${LIBDIR}/cactusBlastAlignment.a ${LIBDIR}/cactusLib.a
sonLibLibs = ${sonLibDir}/sonLib.a ${sonLibDir}/cuTest.a

databaseLibs = ${kyotoTycoonLib} ${tokyoCabinetLib}

LDLIBS += ${cactusLibs} ${sonLibLibs} ${databaseLibs} ${LIBS} -lm -labpoa
LIBDEPENDS = ${sonLibDir}/sonLib.a ${sonLibDir}/cuTest.a
