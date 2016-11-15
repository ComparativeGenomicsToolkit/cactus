include ../include.mk

name = ${dockstore}/ktremotemgr
tag = ${git_commit}

build: Dockerfile
	docker build -t ${name}:${tag} .
	docker tag ${name}:${tag} ${name}:latest

push: build
	# Requires ~/.dockercfg
	docker push ${name}:${tag}

test: build
	python test.py
