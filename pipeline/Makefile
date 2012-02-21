rootPath = ../
include ../include.mk

all : ${binPath}/cactus_workflow.py ${binPath}/cactus_workflow_getFlowers ${binPath}/cactus_workflow_extendFlowers ${binPath}/cactus_workflow_flowerStats

${binPath}/cactus_workflow.py : cactus_workflow.py
	cp cactus_workflow.py ${binPath}/cactus_workflow.py
	chmod +x ${binPath}/cactus_workflow.py

${binPath}/cactus_workflow_getFlowers : *.c *.h ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_workflow_getFlowers cactus_workflow_getFlowers.c flowerWriter.c ${libPath}/cactusLib.a ${basicLibs}

${binPath}/cactus_workflow_extendFlowers : *.c *.h ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_workflow_extendFlowers cactus_workflow_extendFlowers.c flowerWriter.c ${libPath}/cactusLib.a ${basicLibs}

${binPath}/cactus_workflow_flowerStats : *.c *.h ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_workflow_flowerStats cactus_workflow_flowerStats.c ${libPath}/cactusLib.a ${basicLibs}

clean :  
	rm -f *.o
	rm -f ${binPath}/cactus_workflow.py ${binPath}/cactus_workflow_getFlowers ${binPath}/cactus_workflow_extendFlowers ${binPath}/cactus_workflow_flowerStats
