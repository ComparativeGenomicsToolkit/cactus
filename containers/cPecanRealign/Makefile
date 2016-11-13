build_output = cPecanRealign
build_tool = runtime-container.DONE
runtime_fullpath = $(realpath .)
git_commit ?= $(shell git log --pretty=oneline -n 1 -- ../cPecanRealign | cut -f1 -d " ")
name = quay.io/adderan/cpecan-realign
tag = 1.0--${git_commit}

build: ../cPecan/Dockerfile Dockerfile
	cd ../cPecan && make
	docker run -v ${runtime_fullpath}:/data cpecanbuild cp /home/sonLib/bin/${build_output} /data
	docker build -t ${name}:${tag} .
	docker tag ${name}:${tag} ${name}:latest
	docker rmi -f cpecanbuild
	touch ${build_tool}


push: build
	# Requires ~/.dockercfg
	docker push ${name}:${tag}
	docker push ${name}:latest

test: build
	python test.py

clean:
	-rm ${build_tool}
	-rm ${build_output}

all: ${build_tool}
