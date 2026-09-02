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
# CPU baseline for the code we generate.  Override for a machine-specific build:
#   CACTUS_ARCH_FLAGS="-march=native" make
#   CACTUS_ARCH_FLAGS="-march=x86-64-v3 -mtune=znver3" make
#
# x86-64-v3 rather than the bare -mavx2 this used to be.  That is not a portability change:
# the release has required AVX2 since 2021, and the shipped v3.3.0 cactus_consolidated proves
# it -- simd_abpoa_align_sequence_to_subgraph and paf_write_to_buffer carry unguarded AVX2,
# so it already SIGILLs on anything pre-Haswell.  v3 is exactly that AVX2 baseline plus the
# rest of the same CPU generation (FMA, BMI1/2, LZCNT, MOVBE, F16C), which every AVX2-capable
# CPU has.  Measured: minigraph 1.3% faster, lastz 2.1%, output byte-identical in both.
#
# Deliberately NOT x86-64-v4.  AVX-512 would emit binaries that will not run on the
# Skylake-class machines we develop and test on, and it is worth little anyway: abPOA is the
# only hand-vectorised code in the tree, and doubling its vector width 128->256 bits bought
# 8.8%, so 256->512 would be a few percent at best, before AVX-512 downclocking.
#
# No -mtune here.  -mtune=skylake measured a further 2.7% on minigraph and emits no new
# instructions, so it costs no portability -- but it is a bet on one microarchitecture that
# may go the other way on AMD, which we have not measured.  Set it per-cluster via the knob.
# Two baselines, and which you get depends on whether the build will be shipped.
#
# PORTABLE is for anything we distribute.  NATIVE is for a plain "make", which is someone
# compiling cactus for the machine in front of them and should use that machine.  Anything
# that ships its output sets CACTUS_PORTABLE_BUILD=1 -- makeBinRelease and the Dockerfile
# both do -- so the portable value lives in exactly one place and cannot drift.
#
# -march=native is not reliably faster: measured here it gained 4.6% on minigraph and lost
# 5.7% on lastz, the loss isolating to -mtune, which native implies.  It is still the right
# default for a machine-specific build; measure on your own hardware if it matters.
#
# NOT native on ARM: Apple's clang spells it -mcpu=native and rejects -march=native.
# NOT native for CACTUS_LEGACY_ARCH either: that build exists precisely to be portable.
ifdef arm
#	flags to build abpoa
	export armv8 = 1
	export aarch64 = 1
#	flags to include simde abpoa in cactus on ARM
	CACTUS_PORTABLE_ARCH_FLAGS = -march=armv8-a+simd
	CACTUS_NATIVE_ARCH_FLAGS = ${CACTUS_PORTABLE_ARCH_FLAGS}
else ifdef CACTUS_LEGACY_ARCH
	export sse2 = 1
	CACTUS_PORTABLE_ARCH_FLAGS = -msse2
	CACTUS_NATIVE_ARCH_FLAGS = ${CACTUS_PORTABLE_ARCH_FLAGS}
else
#	flags to build abpoa
	export avx2 = 1
#	flags to include simde abpoa in cactus on X86
	CACTUS_PORTABLE_ARCH_FLAGS = -march=x86-64-v3
	CACTUS_NATIVE_ARCH_FLAGS = -march=native
endif

ifdef CACTUS_PORTABLE_BUILD
	CACTUS_ARCH_FLAGS ?= ${CACTUS_PORTABLE_ARCH_FLAGS}
else
	CACTUS_ARCH_FLAGS ?= ${CACTUS_NATIVE_ARCH_FLAGS}
endif

# Both, deliberately.  Until now only CFLAGS carried an arch flag, so Red, hal's C++ and
# cactus2hal were built at the plain x86-64 baseline -- confirmed against the shipped
# binaries, where Red's only AVX2 is glibc's runtime-dispatched strcasecmp, never its own
# code.  They are the components with the most headroom here, having had none at all.
CFLAGS   += ${CACTUS_ARCH_FLAGS}
CXXFLAGS += ${CACTUS_ARCH_FLAGS}

# ...but only for what this makefile compiles itself.  CFLAGS and CXXFLAGS are ordinary make
# variables here, not exported, so "cd submodules/x && ${MAKE}" starts a make that re-derives
# them from sonLib's include.mk and sees no arch flag at all.  Only the release build ever
# propagated them, by way of the "static:" target putting them in the environment -- so in a
# plain build sonLib, cPecan, pinchesAndCacti, matchingAndOrdering, paffy, lastz, hal and
# cactus2hal were all still at the bare baseline, which is most of the tree.
#
# Prefix a sub-make with ${archEnv} to hand them over.  Every submodule we do this for appends
# with +=, so its own -O3 and include paths survive; the arch flags simply land in front.
# Deliberately not "export CFLAGS": that would also push --pedantic, -D__AVX2__ -DUSE_SIMDE
# and our libxml2 include path into every submodule, and abPOA already has to filter
# --pedantic back out.  This hands over the baseline and nothing else.
archEnv = CFLAGS="$${CFLAGS} ${CACTUS_ARCH_FLAGS}" CXXFLAGS="$${CXXFLAGS} ${CACTUS_ARCH_FLAGS}"

# ...but CFLAGS specifically must not be handed to a submodule that recursively builds bundled
# third-party C.  Doing so makes CFLAGS environment-origin, and make then exports it onward
# from that submodule's own make, carrying sonLib's -Werror --pedantic down into code that has
# never compiled clean under them: sonLib/quicktree, sonLib/cutest, cPecan/lastz,
# pinchesAndCacti/threeEdgeConnected and matchingAndOrdering/matchGraph are all in that
# position, and CGL_DEBUG=ultra then fails on unrelated fscanf and format warnings -- which is
# how CI broke.  paffy is in the same position at one remove: it builds its own nested
# sonLib, quicktree and all.  Bundled C++ (blossom) is fine, since CXXFLAGS_ultraDbg
# carries no -Werror.
#
# Hand the baseline to those four through CC instead.  A command-line override reaches the
# nested makes through MAKEFLAGS just as effectively, but carries nothing but the arch flag.
# Not CXX, which would clobber sonLib's own clang++/g++ selection.
archEnvExt = CXXFLAGS="$${CXXFLAGS} ${CACTUS_ARCH_FLAGS}"
archCC = CC="$${CC:-cc} ${CACTUS_ARCH_FLAGS}"

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

# jemalloc flags for submodules that are built by cd'ing into their own tree.  jemallocLib
# deliberately carries no -L, because cactus's own link lines pick up -L${LIBDIR} from LDLIBS
# below -- but a submodule that has already cd'd elsewhere cannot resolve a bare -ljemalloc
# against the copy we build into lib/, and ${LIBDIR} is relative to the cactus root so it
# would point at the wrong directory from there.  Empty when jemalloc is off.
jemallocSubLibs = $(if ${jemallocLib},-L$(abspath ${LIBDIR}) ${jemallocLib})

# note: the CACTUS_STATIC_LINK_FLAGS below can generally be empty -- it's used by the static builder script only
LDLIBS += ${cactusLibs} ${sonLibLibs} ${LIBS} -L${rootPath}/lib -Wl,-rpath,${rootPath}/lib -labpoa ${jemallocLib} -lz -lbz2 -lpthread -lm -lstdc++ -lm -lxml2 ${CACTUS_STATIC_LINK_FLAGS}
LIBDEPENDS = ${sonLibDir}/sonLib.a ${sonLibDir}/cuTest.a ${jemallocDepends}
