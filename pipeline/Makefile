include ../../../include.mk
binPath = ../../../bin
libPath = ../../../lib

cflags += ${tokyoCabinetIncl}

all : ${binPath}/cactus_workflow.py ${binPath}/cactus_workflow_getNets ${binPath}/cactus_workflow_extendNets ${binPath}/cactus_workflow_getUniqueName 

${binPath}/cactus_workflow.py : cactus_workflow.py
	cp cactus_workflow.py ${binPath}/cactus_workflow.py
	chmod +x ${binPath}/cactus_workflow.py

${binPath}/cactus_workflow_getNets : cactus_workflow_getNets.c ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_workflow_getNets cactus_workflow_getNets.c ${libPath}/cactusLib.a ${libPath}/sonLib.a ${tokyoCabinetLib}

${binPath}/cactus_workflow_extendNets : cactus_workflow_extendNets.c ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_workflow_extendNets cactus_workflow_extendNets.c ${libPath}/cactusLib.a ${libPath}/sonLib.a ${tokyoCabinetLib}

${binPath}/cactus_workflow_getUniqueName : cactus_workflow_getUniqueName.c ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_workflow_getUniqueName cactus_workflow_getUniqueName.c ${libPath}/cactusLib.a ${libPath}/sonLib.a ${tokyoCabinetLib}

clean : 
	rm -f *.o
	rm -f ${binPath}/cactus_workflow.py ${binPath}/cactus_workflow_getNets ${binPath}/cactus_workflow_extendNets ${binPath}/cactus_workflow_getUniqueName 

