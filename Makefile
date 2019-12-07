# order is important, libraries first
modules = api setup blastLib caf bar blast normalisation phylogeny reference faces check pipeline preprocessor hal dbTest

git_commit ?= $(shell git rev-parse HEAD)
dockstore = quay.io/comparative-genomics-toolkit
name = ${dockstore}/cactus
tag = ${git_commit}
runtime_fullpath = $(realpath runtime)

export sonLibRootPath = ${PWD}/submodules/sonLib

libSonLib = ${PWD}/submodules/sonLib/sonLib.a
libPinchesAndCacti = ${PWD}/submodules/sonLib/lib/stPinchesAndCacti.a
libCPecan = ${PWD}/submodules/sonLib/lib/cPecanLib.a
libMatchingAndOrdering = ${PWD}/submodules/sonLib/lib/matchingAndOrdering.a

halAppendCactusSubtree = ${PWD}/submodules/cactus2hal/bin/halAppendCactusSubtree
h5c++ = ${PWD}/submodules/hdf5/bin/h5c++

.PHONY: all all.% clean clean.% selfClean copyToBin

all: deps ${modules:%=all.%} ${halAppendCactusSubtree} copyToBin

all.%:
	cd $* && ${MAKE} all

copyToBin:
	cp -rp submodules/sonLib/bin/* bin/
	cp -rp submodules/hal/bin/* bin/
	cp -rp submodules/cactus2hal/bin/* bin/

deps: ${libSonLib} ${libPinchesAndCacti} ${libMatchingAndOrdering} ${libCPecan} hdf5Rule halRule

selfClean:  ${modules:%=clean.%}
	rm -rf lib/*.h bin/*.dSYM

clean.%:
	cd $* && ${MAKE} clean

test:
	pytest .

test_blast:
	pytest . --suite=blast

test_nonblast:
	pytest . --suite=nonblast

docker: Dockerfile
	-docker rmi -f ${name}:latest
	docker build --network=host -t ${name}:${tag} . --build-arg CACTUS_COMMIT=${git_commit}
	docker tag ${name}:${tag} ${name}:latest

push: docker
	# Requires ~/.dockercfg
	docker push ${name}

${libSonLib}:
	@echo "Building dependency sonLib"
	@cd ${PWD}/submodules/sonLib && ${MAKE}

sonLibRule: ${libSonLib}

${libPinchesAndCacti}: ${libSonLib}
	@echo "Building dependency pinchesAndCacti"
	@cd ${PWD}/submodules/pinchesAndCacti && ${MAKE}

${libMatchingAndOrdering}: ${libSonLib}
	@echo "Building dependency matchingAndOrdering"
	@cd ${PWD}/submodules/matchingAndOrdering && ${MAKE}

${libCPecan}: ${libSonLib}
	@echo "Building dependency cPecan"
	@cd ${PWD}/submodules/cPecan && ${MAKE}

${halAppendCactusSubtree}: ${h5c++} halRule all.api
	@echo "Building dependency cactus2hal"
	@cd ${PWD}/submodules/cactus2hal && PATH=${PWD}/submodules/hdf5/bin:$(PATH) ${MAKE}

hdf5Rule: ${h5c++}

${h5c++}:
	@echo "Building dependency hdf5"
	@cd ${PWD}/submodules/hdf5 && ./configure --prefix=$(PWD)/submodules/hdf5 --enable-cxx && CFLAGS=-std=c99 AM_MAKEFLAGS=-e ${MAKE} -j4 -e && ${MAKE} install

halRule: ${h5c++} ${libSonLib}
	@echo "Building dependency hal"
	@cd ${PWD}/submodules/hal && PATH=${PWD}/submodules/hdf5/bin:$(PATH) ${MAKE}

ucscClean: selfClean
	cd ${PWD}/submodules/sonLib && ${MAKE} clean
	cd ${PWD}/submodules/pinchesAndCacti && ${MAKE} clean
	cd ${PWD}/submodules/matchingAndOrdering && ${MAKE} clean

clean: ucscClean selfClean
	cd ${PWD}/submodules/hdf5 && ${MAKE} clean
	cd ${PWD}/submodules/hal && ${MAKE} clean
	cd ${PWD}/submodules/cactus2hal && ${MAKE} clean
