# order is important, libraries first
modules = api setup blastLib caf bar blast normalisation phylogeny reference faces check pipeline preprocessor hal dbTest

git_commit ?= $(shell git rev-parse HEAD)
dockstore = quay.io/comparative-genomics-toolkit
name = ${dockstore}/cactus
tag = ${git_commit}
runtime_fullpath = $(realpath runtime)

export sonLibRootPath = ${PWD}/submodules/sonLib

libtokyocabinet = ${PWD}/submodules/tokyocabinet/libtokyocabinet.a
libkyototycoon = ${PWD}/submodules/kyototycoon/libkyototycoon.a
libkyotocabinet = ${PWD}/submodules/kyotocabinet/libkyotocabinet.a


export tcPrefix = $(PWD)/submodules/tokyocabinet
export tokyoCabinetIncl = -I ${tcPrefix}/include -DHAVE_TOKYO_CABINET=1
export tokyoCabinetLib = -L${tcPrefix}/lib -Wl,-Bstatic -ltokyocabinet -Wl,-Bdynamic -lz -lpthread -lm

export kcPrefix =$(PWD)/submodules/kyotocabinet
export ttPrefix =$(PWD)/submodules/kyototycoon
export kyotoTycoonIncl = -I${kcPrefix}/include -I${ttPrefix}/include -DHAVE_KYOTO_TYCOON=1 -I$(PWD)/zlib/include 
export kyotoTycoonLib = -L${ttPrefix}/lib -L${kcPrefix}/lib -Wl,-Bstatic -lkyototycoon -lkyotocabinet -Wl,-Bdynamic -lz -lpthread -lm -lstdc++

libSonLib = ${PWD}/submodules/sonLib/sonLib.a
libPinchesAndCacti = ${PWD}/submodules/sonLib/lib/stPinchesAndCacti.a
libCPecan = ${PWD}/submodules/sonLib/lib/cPecanLib.a
libMatchingAndOrdering = ${PWD}/submodules/sonLib/lib/matchingAndOrdering.a

halAppendCactusSubtree = ${PWD}/submodules/cactus2hal/bin/halAppendCactusSubtree
h5c++ = ${PWD}/submodules/hdf5/bin/h5c++

.PHONY: all %.all clean %.clean

all : ${libSonLib} ${libPinchesAndCacti} ${libMatchingAndOrdering} ${libCPecan} ${modules:%=all.%} ${halAppendCactusSubtree}


all.%:
	cd $* && make all

clean:  ${modules:%=clean.%}
	rm -rf lib/*.h bin/*.dSYM

clean.%:
	cd $* && make clean

test: all
	python allTests.py


build_output: Dockerfile
	mkdir -p ${runtime_fullpath}/tools
	docker build -t cactusbuild:${tag} .
	docker run -v ${runtime_fullpath}/tools:/data cactusbuild:${tag} sh -c 'cp /home/cactus/bin/* /data'
	docker run -v ${runtime_fullpath}/tools:/data cactusbuild:${tag} sh -c 'cp /home/cactus/submodules/sonLib/bin/* /data'
	docker run -v ${runtime_fullpath}/tools:/data cactusbuild:${tag} sh -c 'cp /home/cactus/submodules/cactus2hal/bin/* /data'
	docker rmi -f cactusbuild:${tag}

docker: build_output runtime/Dockerfile
	docker rmi -f ${name}:latest
	docker build -t ${name}:${tag} ./runtime/ --build-arg CACTUS_COMMIT=${git_commit}
	docker tag ${name}:${tag} ${name}:latest

push: docker
	# Requires ~/.dockercfg
	docker push ${name}


${libSonLib}: ${libkyotocabinet} ${libtokyocabinet} ${libkyototycoon}
	cd ${PWD}/submodules/sonLib && make

sonLibRule: ${libSonLib}

tokyocabinetRule: ${libtokyocabinet}

kyototycoonRule: ${libkyototycoon}

kyotocabinetRule: ${libkyotocabinet}

${libtokyocabinet}:
	cd ${PWD}/submodules/tokyocabinet && ./configure --prefix=${PWD}/submodules/tokyocabinet --enable-static --disable-shared --disable-bzip && make && make install

${libkyototycoon}:
	cd ${PWD}/submodules/kyototycoon && ./configure --prefix=${PWD}/submodules/kyototycoon --enable-static --disable-shared --with-kc=${PWD}/submodules/kyotocabinet && make && make install

${libkyotocabinet}:
	cd ${PWD}/submodules/kyotocabinet && ./configure --prefix=${PWD}/submodules/kyotocabinet --enable-static --disable-shared && make && make install


${libPinchesAndCacti}:
	cd ${PWD}/submodules/pinchesAndCacti && make

${libMatchingAndOrdering}:
	cd ${PWD}/submodules/matchingAndOrdering && make

${libCPecan}:
	cd ${PWD}/submodules/cPecan && make

${halAppendCactusSubtree}: ${h5c++} halRule
	cd ${PWD}/submodules/cactus2hal && PATH=${PWD}/submodules/hdf5/bin:$(PATH) make

hdf5Rule: ${h5c++}

${h5c++}:
	cd ${PWD}/submodules/hdf5 &&./configure --prefix=$(PWD)/submodules/hdf5 --enable-cxx && CFLAGS=-std=c99 make -e && make install

halRule:
	cd ${PWD}/submodules/hal && PATH=${PWD}/submodules/hdf5/bin:$(PATH) make

ucscClean:
	make clean
	cd ${PWD}/submodules/sonLib && make clean
	cd ${PWD}/submodules/pinchesAndCacti && make clean
	cd ${PWD}/submodules/matchingAndOrdering && make clean

clean: ucscClean
	cd ${PWD}/submodules/kyotocabinet && make clean
	cd ${PWD}/submodules/kyototycoon && make clean
	cd ${PWD}/submodules/tokyocabinet && make clean
	cd ${PWD}/submodules/hdf5 && make clean
	cd ${PWD}/submodules/hal && make clean
	cd ${PWD}/submodules/cactus2hal && make clean
