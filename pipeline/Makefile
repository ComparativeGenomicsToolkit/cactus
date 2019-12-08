rootPath = ..
include ${rootPath}/include.mk

all: all_libs all_progs
all_libs: 
all_progs: all_libs
	${MAKE} ${binPath}/cactus_workflow_getFlowers ${binPath}/cactus_workflow_extendFlowers ${binPath}/cactus_workflow_flowerStats ${binPath}/cactus_workflow_convertAlignmentCoordinates ${binPath}/cactus_secondaryDatabase ${binPath}/docker_test_script

${binPath}/cactus_workflow_getFlowers : *.c *.h ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_workflow_getFlowers cactus_workflow_getFlowers.c ${libPath}/cactusLib.a ${basicLibs}

${binPath}/cactus_workflow_extendFlowers : *.c *.h ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_workflow_extendFlowers cactus_workflow_extendFlowers.c ${libPath}/cactusLib.a ${basicLibs}

${binPath}/cactus_workflow_flowerStats : *.c *.h ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_workflow_flowerStats cactus_workflow_flowerStats.c ${libPath}/cactusLib.a ${basicLibs}

${binPath}/cactus_workflow_convertAlignmentCoordinates : *.c *.h ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_workflow_convertAlignmentCoordinates cactus_workflow_convertAlignmentCoordinates.c ${libPath}/cactusLib.a ${basicLibs}

${binPath}/cactus_secondaryDatabase : *.c *.h ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I${libPath} -o ${binPath}/cactus_secondaryDatabase cactus_secondaryDatabase.c ${libPath}/cactusLib.a ${basicLibs}

${binPath}/docker_test_script : docker_test_script.py
	cp docker_test_script.py ${binPath}/docker_test_script
	chmod +x ${binPath}/docker_test_script

clean :  
	rm -f *.o
	rm -f ${binPath}/cactus_workflow.py ${binPath}/cactus_workflow_getFlowers ${binPath}/cactus_workflow_extendFlowers ${binPath}/cactus_workflow_flowerStats ${binPath}/cactus_workflow_convertAlignmentCoordinates ${binPath}/cactus_secondaryDatabase ${binPath}/docker_test_script
