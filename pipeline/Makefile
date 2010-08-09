include ../../../include.mk
binPath = ../../../bin
libPath = ../../../lib

cflags += ${tokyoCabinetIncl}

all : ${binPath}/cactus_workflow.py ${binPath}/cactus_workflow_getFlowers ${binPath}/cactus_workflow_extendFlowers ${binPath}/cactus_workflow_getUniqueName ${binPath}/cactus_workflow_makeTerminalNormal

${binPath}/cactus_workflow.py : cactus_workflow.py
	cp cactus_workflow.py ${binPath}/cactus_workflow.py
	chmod +x ${binPath}/cactus_workflow.py

${binPath}/cactus_workflow_getFlowers : cactus_workflow_getFlowers.c ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_workflow_getFlowers cactus_workflow_getFlowers.c ${libPath}/cactusLib.a ${libPath}/sonLib.a ${tokyoCabinetLib}

${binPath}/cactus_workflow_extendFlowers : cactus_workflow_extendFlowers.c ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_workflow_extendFlowers cactus_workflow_extendFlowers.c ${libPath}/cactusLib.a ${libPath}/sonLib.a ${tokyoCabinetLib}

${binPath}/cactus_workflow_getUniqueName : cactus_workflow_getUniqueName.c ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_workflow_getUniqueName cactus_workflow_getUniqueName.c ${libPath}/cactusLib.a ${libPath}/sonLib.a ${tokyoCabinetLib}

${binPath}/cactus_workflow_makeTerminalNormal : cactus_workflow_makeTerminalNormal.c ${libPath}/sonLib.a ${libPath}/cactusLib.a
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_workflow_makeTerminalNormal cactus_workflow_makeTerminalNormal.c ${libPath}/cactusLib.a ${libPath}/sonLib.a ${tokyoCabinetLib}

clean : 
	rm -f *.o
	rm -f ${binPath}/cactus_workflow.py ${binPath}/cactus_workflow_getFlowers ${binPath}/cactus_workflow_extendFlowers ${binPath}/cactus_workflow_getUniqueName ${binPath}/cactus_workflow_makeTerminalNormal

