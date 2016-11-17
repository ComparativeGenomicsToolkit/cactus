include ../include.mk

name = ${dockstore}/cactus
tag = ${git_commit}

build: Dockerfile
	docker build -t ${name}:${tag} . --build-arg CACTUS_COMMIT=${git_commit}

push: build
	# Requires ~/.dockercfg
	docker push ${name}:${tag}

test: build
	python test.py
