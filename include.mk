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

# Hack to include xml2
ifeq (${CACTUS_LIBXML2_INCLUDE_PATH},)
CACTUS_LIBXML2_INCLUDE_PATH = /usr/include/libxml2
endif
CFLAGS+= -I${CACTUS_LIBXML2_INCLUDE_PATH}

ifndef TARGETOS
  TARGETOS := $(shell uname -s)
endif

# Control variable for jemalloc
# Disable on Mac, disable when address sanitizer is active, disable for legacy arch builds
jemalloc = on
ifeq ($(TARGETOS), Darwin)
	jemalloc = off
endif
ifeq (${CGL_DEBUG},ultra)
	jemalloc = off
endif
ifdef CACTUS_LEGACY_ARCH
	jemalloc = off
endif

# Hack to include openmp on os x after "brew install lomp
ifeq ($(TARGETOS), Darwin)
	CFLAGS+= -Xpreprocessor -fopenmp -lomp
else
	CFLAGS+= -fopenmp
endif

# Hack in ARM support
# Toggle on if "arm" is set, or if uname -m returns aarch64
ifeq ($(shell uname -m || true), aarch64)
	arm=1
endif
ifeq ($(shell arch || true), aarch64)
	arm=1
endif
ifeq ($(shell arch || true), arm64)
	arm=1
endif
ifdef arm
#	flags to build abpoa
	export armv8 = 1
	export aarch64 = 1
#	flags to include simde abpoa in cactus on ARM
	CFLAGS+= -march=armv8-a+simd
else
	ifdef CACTUS_LEGACY_ARCH
		export sse2 = 1
		CFLAGS+= -msse2
	else
#		flags to build abpoa
		export avx2 = 1
#		flags to include simde abpoa in cactus on X86
		CFLAGS+= -mavx2
endif
endif

# flags needed to include simde abpoa in cactus on any architecture
ifdef CACTUS_LEGACY_ARCH
	CFLAGS+= -D__SSE2__ -DUSE_SIMDE -DSIMDE_ENABLE_NATIVE_ALIASES
else
	CFLAGS+= -D__AVX2__ -DUSE_SIMDE -DSIMDE_ENABLE_NATIVE_ALIASES
endif

dataSetsPath=/Users/benedictpaten/Dropbox/Documents/work/myPapers/genomeCactusPaper/dataSets

inclDirs = hal/inc api/inc setup/inc bar/inc caf/inc paf/inc hal/inc reference/inc pipeline/inc submodules/sonLib/C/inc \
	blastLib submodules/sonLib/externalTools/cutest submodules/pinchesAndCacti/inc \
	submodules/matchingAndOrdering/inc submodules/cPecan/inc

CPPFLAGS += ${inclDirs:%=-I${rootPath}/%} -I${LIBDIR} -I${rootPath}/include

# libraries can't be added until they are build, so add as to LDLIBS until needed
cactusLibs = ${LIBDIR}/stCaf.a ${LIBDIR}/stReference.a ${LIBDIR}/cactusBarLib.a ${LIBDIR}/cactusLib.a
sonLibLibs = ${sonLibDir}/sonLib.a ${sonLibDir}/cuTest.a

# Add jemalloc if enabled
ifeq ($(jemalloc),on)
	jemallocLib = -ljemalloc -lm
	jemallocDepends = ${LIBDIR}/libjemalloc.a
endif

# note: the CACTUS_STATIC_LINK_FLAGS below can generally be empty -- it's used by the static builder script only
LDLIBS += ${cactusLibs} ${sonLibLibs} ${LIBS} -L${rootPath}/lib -Wl,-rpath,${rootPath}/lib -labpoa ${jemallocLib} -lz -lbz2 -lpthread -lm -lstdc++ -lm -lxml2 ${CACTUS_STATIC_LINK_FLAGS}
LIBDEPENDS = ${sonLibDir}/sonLib.a ${sonLibDir}/cuTest.a ${jemallocDepends}
