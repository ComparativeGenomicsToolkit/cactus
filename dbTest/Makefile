rootPath = ../
include ../include.mk

host=localhost
port=1978
databaseDir=./testDb
databaseOptions=-ls -tout 200000 -th 64 -bgsi 100000000 #opts=ls#bnum=30m#msiz=50g#ktopts=p
minRecordSize=0
maxRecordSize=100
keysPerJob=100
totalJobs=2

all : ${binPath}/dbTestScript 

${binPath}/dbTestScript  : ${basicLibsDependencies} *.c
	${cxx} ${cflags} -I inc -I${libPath} -Wno-error -o ${binPath}/dbTestScript *.c ${basicLibs}
	
clean :
	rm -rf ${binPath}/dbTestScript

test :
	ktserver -log ./log.txt -host ${host} -port ${port} ${databaseOptions} &
	rm -rf ${databaseDir} ./jobTree
	python dbTestScript.py --logLevel DEBUG --jobTree ./jobTree --stats --keysPerJob ${keysPerJob} --totalJobs ${totalJobs} --host ${host} --port ${port} --databaseDir ${databaseDir} --databaseOption "${databaseOptions}" --minRecordSize ${minRecordSize} --maxRecordSize ${maxRecordSize}
	jobTreeStats --jobTree ./jobTree --outputFile ./jobTreeStatsTest.xml
	rm -rf ${databaseDir} ./jobTree
	ktremotemgr report -host ${host} -port ${port}
	ps ax | grep 'ktserver' | cut -f1 -d' ' | xargs kill