# order is important, libraries first
modules = api setup blastLib caf bar blast normalisation phylogeny reference faces check pipeline preprocessor hal dbTest repeats

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
	cd $* && make all

copyToBin:
	cp -rp submodules/sonLib/bin/* bin/
	cp -rp submodules/hal/bin/* bin/
	cp -rp submodules/cactus2hal/bin/* bin/

deps: ${libSonLib} ${libPinchesAndCacti} ${libMatchingAndOrdering} ${libCPecan} hdf5Rule halRule

selfClean:  ${modules:%=clean.%}
	rm -rf lib/*.h bin/*.dSYM

clean.%:
	cd $* && make clean

test:
	python allTests.py

build_output: Dockerfile
	mkdir -p ${runtime_fullpath}/tools
	docker build -t cactusbuild:${tag} .
	docker run -v ${runtime_fullpath}/tools:/data cactusbuild:${tag} sh -c 'cp /home/cactus/bin/* /data'
	docker run -v ${runtime_fullpath}/tools:/data cactusbuild:${tag} sh -c 'cp /home/cactus/submodules/sonLib/bin/* /data'
	docker run -v ${runtime_fullpath}/tools:/data cactusbuild:${tag} sh -c 'cp /home/cactus/submodules/cactus2hal/bin/* /data'

docker: build_output runtime/Dockerfile
	-docker rmi -f ${name}:latest
	docker build -t ${name}:${tag} ./runtime/ --build-arg CACTUS_COMMIT=${git_commit}
	docker tag ${name}:${tag} ${name}:latest

push: docker
	# Requires ~/.dockercfg
	docker push ${name}

${libSonLib}:
	@echo "Building dependency sonLib"
	@cd ${PWD}/submodules/sonLib && (output=$$(make 2>&1) || (echo "$$output"; exit 1))

sonLibRule: ${libSonLib}

${libPinchesAndCacti}: ${libSonLib}
	@echo "Building dependency pinchesAndCacti"
	@cd ${PWD}/submodules/pinchesAndCacti && (output=`make 2>&1` || (echo "$$output"; exit 1))

${libMatchingAndOrdering}: ${libSonLib}
	@echo "Building dependency matchingAndOrdering"
	@cd ${PWD}/submodules/matchingAndOrdering && (output=`make 2>&1` || (echo "$$output"; exit 1))

${libCPecan}: ${libSonLib}
	@echo "Building dependency cPecan"
	@cd ${PWD}/submodules/cPecan && (output=`make 2>&1` || (echo "$$output"; exit 1))

${halAppendCactusSubtree}: ${h5c++} halRule all.api
	@echo "Building dependency cactus2hal"
	@cd ${PWD}/submodules/cactus2hal && (PATH=${PWD}/submodules/hdf5/bin:$(PATH) output=`make 2>&1` || (echo "$$output"; exit 1))

hdf5Rule: ${h5c++}

${h5c++}:
	@echo "Building dependency hdf5"
	@cd ${PWD}/submodules/hdf5 && (output=`(./configure --prefix=$(PWD)/submodules/hdf5 --enable-cxx && CFLAGS=-std=c99 make -j4 -e && make install) 2>&1` || (echo "$$output"; exit 1))

halRule: ${h5c++} ${libSonLib}
	@echo "Building dependency hal"
	@cd ${PWD}/submodules/hal && (PATH=${PWD}/submodules/hdf5/bin:$(PATH) output=`make 2>&1` || (echo "$$output"; exit 1))

ucscClean: selfClean
	cd ${PWD}/submodules/sonLib && make clean
	cd ${PWD}/submodules/pinchesAndCacti && make clean
	cd ${PWD}/submodules/matchingAndOrdering && make clean

clean: ucscClean selfClean
	cd ${PWD}/submodules/hdf5 && make clean
	cd ${PWD}/submodules/hal && make clean
	cd ${PWD}/submodules/cactus2hal && make clean
