include ../include.mk

name = ${dockstore}/cactus
tag = ${git_commit}
runtime_fullpath = $(realpath runtime)

build_output: build/Dockerfile
	mkdir -p ${runtime_fullpath}/tools
	docker build -t cactusbuild:${tag} build/ --build-arg CACTUS_COMMIT=${git_commit}
	docker run -v ${runtime_fullpath}/tools:/data cactusbuild:${tag} sh -c 'cp /home/cactus/bin/* /data'
	docker run -v ${runtime_fullpath}/tools:/data cactusbuild:${tag} sh -c 'cp /home/sonLib/bin/* /data'
	docker run -v ${runtime_fullpath}/tools:/data cactusbuild:${tag} sh -c 'cp /home/cactus2hal/bin/* /data'

build: build_output runtime/Dockerfile
	docker build -t ${name}:${tag} ./runtime/

push: build
	# Requires ~/.dockercfg
	docker push ${name}:${tag}

test: build
	python test.py
