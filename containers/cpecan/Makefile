include ../include.mk

name = ${dockstore}/cpecan
tag = ${git_commit}
runtime_fullpath = $(realpath runtime)


build_output: build/Dockerfile
	mkdir -p ${runtime_fullpath}/tools
	docker build -t cpecanbuild:${tag} ./build
	docker run -v ${runtime_fullpath}/tools/:/data cpecanbuild:${tag} sh -c 'cp /home/sonLib/bin/* /data'

build: build_output runtime/Dockerfile
	docker build -t ${name}:${tag} ./runtime

push: build
	# Requires ~/.dockercfg
	docker push ${name}:${tag}

test: build
	python test.py

