# Definitions
build_output = runtime/lastz
runtime_fullpath = $(realpath runtime)
build_tool = runtime-container.DONE
git_commit ?= $(shell git log --pretty=oneline -n 1 -- ../lastz | cut -f1 -d " ")
name = quay.io/adderan/lastz
tag = 1.0.0--${git_commit}

# Steps
build: ${build_output} ${build_tool}

${build_output}: build/Dockerfile
	cd build && docker build -t lastzbuild .
	docker run -v ${runtime_fullpath}:/data lastzbuild cp /home/cPecan/externalTools/lastz-distrib-1.03.54/src/lastz /data

${build_tool}: ${build_output} runtime/Dockerfile
	cd runtime && docker build -t ${name}:${tag} .
	docker tag ${name}:${tag} ${name}:latest
	#docker rmi lastzbuild
	touch ${build_tool}

push: build
	docker push ${name}:${tag}
	docker push ${name}:latest

test: build
	python test.py

clean:
	-rm ${build_tool}
	-rm ${build_output}
