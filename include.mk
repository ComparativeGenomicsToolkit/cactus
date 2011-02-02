#Location of sonLib
binPath=${rootPath}bin
libPath=${rootPath}lib
#Modify this variable to set the location of sonLib
sonLibRootPath=${rootPath}../sonLib 
sonLibPath=${sonLibRootPath}/lib

include  ${sonLibRootPath}/include.mk

dataSetsPath=/Users/benedictpaten/Dropbox/Documents/work/myPapers/genomeCactusPaper/dataSets

cflags += -I ${sonLibPath}
dblibs = ${tokyoCabinetLib} ${mysqlLibs} ${pgsqlLibs}
basicLibs = ${sonLibPath}/sonLib.a ${sonLibPath}/cutest.a ${dblibs}
basicLibsDependencies = ${sonLibPath}/sonLib.a ${sonLibPath}/cutest.a 

# optional kent library stuff
KENTDIR=/hive/groups/recon/local/kent/src
ifneq ($(wildcard ${KENTDIR}),)
    kentInc = ${KENTDIR}/inc
    kentLib = ${KENTDIR}/lib
    kentLibWeb = ${kentLib}/${MACH}/jkweb.a
endif
