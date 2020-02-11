#Location of sonLib
binPath=${rootPath}/bin/
libPath=${rootPath}/lib/

#Modify this variable to set the location of sonLib
sonLibRootPath?=${PWD}/../sonLib
sonLibPath=${sonLibRootPath}/lib

include ${sonLibRootPath}/include.mk

dataSetsPath=/Users/benedictpaten/Dropbox/Documents/work/myPapers/genomeCactusPaper/dataSets

inclDirs = api/inc bar/inc caf/inc hal/inc reference/inc submodules/sonLib/C/inc \
	blastLib submodules/sonLib/externalTools/cutest submodules/pinchesAndCacti/inc \
	submodules/matchingAndOrdering/inc submodules/cPecan/inc

cflags += ${inclDirs:%=-I${rootPath}/%}
basicLibs = ${sonLibPath}/sonLib.a ${sonLibPath}/cuTest.a ${dblibs}
basicLibsDependencies = ${sonLibPath}/sonLib.a ${sonLibPath}/cuTest.a 
