rootPath = ../
include ../include.mk

all: ${binPath}/cactus_coveredSeeds ${binPath}/cactus_adjustSamplingRates

clean:
	rm ${binPath}/cactus_coveredSeeds ${binPath}/cactus_adjustSamplingRates

${binPath}/cactus_coveredSeeds: cactus_coveredSeeds.py
	cp cactus_coveredSeeds.py ../bin/cactus_coveredSeeds
	chmod +x ../bin/cactus_coveredSeeds

${binPath}/cactus_adjustSamplingRates: cactus_adjustSamplingRates.py
	cp cactus_adjustSamplingRates.py ../bin/cactus_adjustSamplingRates
	chmod +x ../bin/cactus_adjustSamplingRates
