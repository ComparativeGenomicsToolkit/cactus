include ../include.mk

name = ${dockstore}/cpecan
tag = ${git_commit}


build: Dockerfile
	docker build -t ${name}:${tag} .

push: build
	# Requires ~/.dockercfg
	docker push ${name}:${tag}

test: build
	python test.py

