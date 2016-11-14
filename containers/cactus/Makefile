build_tool = runtime-container.DONE
git_commit ?= $(shell git rev-parse HEAD)
name = quay.io/adderan/cactus
tag = ${git_commit}

build: ${build_tool}

${build_tool}: Dockerfile
	docker build -t ${name}:${tag} . --build-arg CACTUS_COMMIT=${git_commit}
	docker tag ${name}:${tag} ${name}:latest
	touch ${build_tool}

push: build
	# Requires ~/.dockercfg
	docker push ${name}:${tag}

test: build
	python test.py

clean:
	-rm ${build_tool}
