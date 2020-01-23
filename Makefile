rootPath = .

include ${rootPath}/include.mk

modules = api setup blastLib caf bar blast normalisation hal phylogeny reference faces check pipeline preprocessor hal dbTest

# submodules are in multiple pass to handle dependencies cactus2hal being dependent on
# both cactus and sonLib
submodules1 = sonLib cPecan hal matchingAndOrdering pinchesAndCacti
submodules2 = cactus2hal
submodules = ${submodules1} ${submodules2}

git_commit ?= $(shell git rev-parse HEAD)
dockstore = quay.io/comparative-genomics-toolkit
name = ${dockstore}/cactus
tag = ${git_commit}
runtime_fullpath = $(realpath runtime)

# FIXME: this is also built by setup.py need to sort this out
versionPy = src/cactus/shared/version.py

CWD = ${PWD}

# these must be absolute, as used in submodules.
export sonLibRootPath = ${CWD}/submodules/sonLib
.PHONY: all all.% clean clean.% selfClean suball suball.% subclean.%

##
# Building.  First build submodules, then a pass for libs and a pass for bins
##
all:
	${MAKE} suball1
	${MAKE} all_libs
	${MAKE} all_progs
	${MAKE} suball2

all_libs:
	${MAKE} ${modules:%=all_libs.%}
all_progs: all_libs
	${MAKE} ${modules:%=all_progs.%}

all_libs.%:
	cd $* && ${MAKE} all_libs
all_progs.%:
	cd $* && ${MAKE} all_progs

# special ordering
all_libs.blastLib: all_libs.api

##
# tests, see DEVELOPMENT.md for environment variables controling tests.
##

# under test modules under src/cactus/, split up to allowing run them in parallel.
testModules = \
    bar/cactus_barTest.py \
    blast/blastTest.py \
    blast/cactus_coverageTest.py \
    blast/cactus_realignTest.py \
    blast/mappingQualityRescoringAndFilteringTest.py \
    blast/trimSequencesTest.py \
    faces/cactus_fillAdjacenciesTest.py \
    hal/cactus_halTest.py \
    normalisation/cactus_normalisationTest.py \
    phylogeny/cactus_phylogenyTest.py \
    pipeline/cactus_evolverTest.py \
    pipeline/cactus_workflowTest.py \
    preprocessor/cactus_preprocessorTest.py \
    preprocessor/lastzRepeatMasking/cactus_lastzRepeatMaskTest.py \
    progressive/cactus_progressiveTest.py \
    progressive/multiCactusTreeTest.py \
    progressive/outgroupTest.py \
    progressive/scheduleTest.py \
    reference/cactus_referenceTest.py \
    shared/commonTest.py \
    shared/experimentWrapperTest.py

export PATH := ${CWD}/bin:${PATH}
export PYTHONPATH = ${CWD}/submodules/sonLib/src:${CWD}/submodules:${CWD}/src
pytestOpts = --tb=native --durations=0 -rsx
testLogDir = testlogs

test: ${testModules:%=%_runtest}
test_blast: ${testModules:%=%_runtest}
test_nonblast: ${testModules:%=%_runtest_nonblast}

# run one test and save output
%_runtest: ${versionPy}
	@mkdir -p ${testLogDir}
	${PYTHON} -m pytest ${pytestOpts} src/cactus/$* >& ${testLogDir}/$(basename $(notdir $*)).log

%_runtest_blast: ${versionPy}
	@mkdir -p ${testLogDir}
	${PYTHON} -m pytest ${pytestOpts} src/cactus/$* --suite=blast >& ${testLogDir}/$(basename $(notdir $*)).blast.log

%_runtest_nonblast: ${versionPy}
	@mkdir -p ${testLogDir}
	${PYTHON} -m pytest ${pytestOpts} src/cactus/$* --suite=nonblast >& ${testLogDir}/$(basename $(notdir $*)).nonblast.log

${versionPy}:
	echo "cactus_commit = '${git_commit}'" >$@


##
# clean targets
##
selfClean: ${modules:%=clean.%}
	rm -rf lib/*.h bin/*.dSYM ${versionPy} ${testLogDir}

clean.%:
	cd $* && ${MAKE} clean

clean: selfClean ${submodules:%=subclean.%}

##
# submodules
##
suball1: ${submodules1:%=suball.%}
suball2: ${submodules2:%=suball.%}
suball.sonLib:
	cd submodules/sonLib && ${MAKE}
	mkdir -p bin
	ln -f submodules/sonLib/bin/[a-zA-Z]* bin/

suball.pinchesAndCacti: suball.sonLib
	cd submodules/pinchesAndCacti && ${MAKE}

suball.matchingAndOrdering: suball.sonLib
	cd submodules/matchingAndOrdering && ${MAKE}

suball.cPecan: suball.sonLib
	cd submodules/cPecan && ${MAKE}

suball.cactus2hal: suball.sonLib suball.hal all_libs.api
	cd submodules/cactus2hal && ${MAKE}
	mkdir -p bin
	-ln -f submodules/cactus2hal/bin/* bin/

suball.hal: suball.sonLib
	cd submodules/hal &&  ${MAKE}
	mkdir -p bin
	-ln -f submodules/hal/bin/* bin/
	-ln -f submodules/hal/lib/libHal.a submodules/hal/lib/halLib.a

subclean.%:
	cd submodules/$* && ${MAKE} clean

##
# docker
##
docker: Dockerfile
	-docker rmi -f ${name}:latest
	docker build --network=host -t ${name}:${tag} . --build-arg CACTUS_COMMIT=${git_commit}
	docker tag ${name}:${tag} ${name}:latest

# Requires ~/.dockercfg
push: docker
	docker push ${name}

