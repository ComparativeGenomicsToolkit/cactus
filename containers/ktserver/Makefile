git_commit ?= $(shell git rev-parse HEAD)
name = quay.io/adderan/ktserver
tag = ${git_commit}


build: Dockerfile
	docker build -t ${name}:${tag} .
	docker tag ${name}:${tag} ${name}:latest

push: build
	# Requires ~/.dockercfg
	docker push ${name}:${tag}
	docker push ${name}:latest

test: build
	python test.py
