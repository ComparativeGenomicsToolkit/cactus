git_commit ?= $(shell git rev-parse HEAD)
name = quay.io/adderan/cpecan
tag = ${git_commit}


build: Dockerfile
	docker build -t ${name}:${tag} .
	docker tag ${name}:${tag} ${name}:latest

push: build
	# Requires ~/.dockercfg
	docker push ${name}:${tag}

test: build
	python test.py

