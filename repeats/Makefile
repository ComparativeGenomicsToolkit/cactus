rootPath=../
include ../include.mk

libs=../../TwoPaCo/build/graphconstructor/libgraphconstructor.a

includes=-I ../../TwoPaCo/src/common -I ../../TwoPaCo/src/graphconstructor -I ../../TwoPaCo/src/graphdump -I ../../TwoPaCo/src/naiveconstructor

${binPath}/cactus_repeats_seedProbabilities: cactus_repeats_seedProbabilities.cpp ${libs}
	${cpp} -std=c++11 -g -Wall ${includes} -o ${binPath}/cactus_repeats_seedProbabilities cactus_repeats_seedProbabilities.cpp ${libs} -ltbb


all: ${binPath}/cactus_repeats_seedProbabilities

clean:
	rm ${binPath}/cactus_repeats_seedProbabilities

