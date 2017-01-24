# order is important, libraries first
modules = api setup blastLib caf bar blast normalisation phylogeny reference faces check pipeline preprocessor hal dbTest

git_commit ?= $(shell git rev-parse HEAD)
dockstore = quay.io/comparative-genomics-toolkit
name = ${dockstore}/cactus
tag = ${git_commit}
runtime_fullpath = $(realpath runtime)

.PHONY: all %.all clean %.clean

all : ${modules:%=all.%}


all.%:
	cd $* && make all

clean:  ${modules:%=clean.%}
	rm -rf lib/*.h bin/*.dSYM

clean.%:
	cd $* && make clean

test: all
	python allTests.py

docker:
	cd containers && make build

docker_push: docker
	cd containers && make push


build_output: Dockerfile
	mkdir -p ${runtime_fullpath}/tools
	docker build -t cactusbuild:${tag} .
	docker run -v ${runtime_fullpath}/tools:/data cactusbuild:${tag} sh -c 'cp /home/cactus/bin/* /data'
	docker run -v ${runtime_fullpath}/tools:/data cactusbuild:${tag} sh -c 'cp /home/cactus/submodules/sonLib/bin/* /data'
	docker run -v ${runtime_fullpath}/tools:/data cactusbuild:${tag} sh -c 'cp /home/cactus/submodules/cactus2hal/bin/* /data'

build: build_output runtime/Dockerfile
	docker build -t ${name}:${tag} ./runtime/ --build-arg CACTUS_COMMIT=${git_commit}
	docker tag ${name}:${tag} ${name}:latest

push: build
	# Requires ~/.dockercfg
	docker push ${name}:${tag}

