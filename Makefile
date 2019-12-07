modules = api setup blastLib caf bar blast normalisation phylogeny reference faces check pipeline preprocessor hal dbTest

# submodules are in two passes, since cactus2hal requires api which requires sonLib
submodules1 = cPecan hal hdf5 matchingAndOrdering pinchesAndCacti sonLib
submodules2 = cactus2hal
submodules = ${submodules1} ${submodules2}

git_commit ?= $(shell git rev-parse HEAD)
dockstore = quay.io/comparative-genomics-toolkit
name = ${dockstore}/cactus
tag = ${git_commit}
runtime_fullpath = $(realpath runtime)

CWD = ${PWD}

export sonLibRootPath = ${CWD}/submodules/sonLib
hdf5Bin = ${CWD}/submodules/hdf5/bin
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
	${MAKE} ${modules:%=all_progs.%s}

all_libs.%:
	cd $* && ${MAKE} all_libs
all_progs.%:
	cd $* && ${MAKE} all_progs



##
# tests
##
test:
	pytest .

test_blast:
	pytest . --suite=blast

test_nonblast:
	pytest . --suite=nonblast

##
# clean targets
##
selfClean: ${modules:%=clean.%}
	rm -rf lib/*.h bin/*.dSYM

clean.%:
	cd $* && ${MAKE} clean

clean: ${modules:%=clean.%} ${submodules:%=subclean.%}

##
# submodules
##
suball1: ${submodules1:%=suball.%}
suball2: ${submodules2:%=suball.%}
suball.sonLib:
	cd submodules/sonLib && ${MAKE}
	mkdir -p bin
	ln -f submodules/sonLib/bin/* bin/

suball.pinchesAndCacti: suball.sonLib
	cd submodules/pinchesAndCacti && ${MAKE}

suball.matchingAndOrdering: suball.sonLib
	cd submodules/matchingAndOrdering && ${MAKE}

suball.cPecan: suball.sonLib
	cd submodules/cPecan && ${MAKE}

suball.cactus2hal: suball.sonLib suball.hal all_libs.api
	cd submodules/cactus2hal && PATH=${hdf5Bin}:$(PATH) ${MAKE}
	mkdir -p bin
	ln -f submodules/cactus2hal/bin/* bin/

suball.hal: suball.hdf5 suball.sonLib
	cd submodules/hal && PATH=${CWD}/submodules/hdf5/bin:$(PATH) ${MAKE}
	mkdir -p bin
	ln -f submodules/hal/bin/* bin/

# hdf5 insists on running configure, so it isn't quick.  We shouldn't be building this anyway.
someHdf5File = submodules/hdf5/bin/yodconfigure 
suball.hdf5: ${someHdf5File}
${someHdf5File}:
	cd submodules/hdf5 && ./configure --prefix=${CWD}/submodules/hdf5 --enable-cxx && CFLAGS=-std=c99 AM_MAKEFLAGS=-e ${MAKE} -j4 -e && ${MAKE} install

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

