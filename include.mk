#Location of sonLib
binPath=${rootPath}bin/
libPath=${rootPath}lib/

#Modify this variable to set the location of sonLib
sonLibRootPath=${rootPath}/../sonLib
sonLibPath=${sonLibRootPath}/lib

include  ${sonLibRootPath}/include.mk

dataSetsPath=/Users/benedictpaten/Dropbox/Documents/work/myPapers/genomeCactusPaper/dataSets

cflags += -I ${sonLibPath} ${tokyoCabinetIncl} ${kyotoTycoonIncl} ${mysqlIncl} ${pgsqlIncl}
basicLibs = ${sonLibPath}/sonLib.a ${sonLibPath}/cuTest.a ${dblibs}
basicLibsDependencies = ${sonLibPath}/sonLib.a ${sonLibPath}/cuTest.a 
