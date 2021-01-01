rootPath = ..
include ${rootPath}/include.mk

# Hack to include libxml2
CFLAGS+= -I/usr/include/libxml2 -lxml2 ${LDFLAGS}

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

all: all_libs all_progs
all_libs: 
all_progs: all_libs
	${MAKE} ${BINDIR}/stPipelineTests ${BINDIR}/cactus_consolidated ${BINDIR}/cactus_workflow_getFlowers ${BINDIR}/cactus_workflow_extendFlowers ${BINDIR}/cactus_workflow_flowerStats ${BINDIR}/cactus_workflow_convertAlignmentCoordinates ${BINDIR}/cactus_secondaryDatabase ${BINDIR}/docker_test_script

${BINDIR}/stPipelineTests : ${libTests} ${libSources} ${libHeaders} ${LIBDIR}/cactusLib.a ${LIBDEPENDS}
	${CC} ${CPPFLAGS} ${CFLAGS} -o ${BINDIR}/stPipelineTests ${libTests} ${libSources} ${LIBDIR}/cactusLib.a ${LDLIBS}

${BINDIR}/cactus_consolidated : ${libSources} ${libHeaders} *.c *.h ${LIBDIR}/cactusLib.a ${LIBDEPENDS}
	${CC} ${CPPFLAGS} ${CFLAGS} -o ${BINDIR}/cactus_consolidated cactus_consolidated.c ${libSources} ${LIBDIR}/cactusLib.a ${LDLIBS}

${BINDIR}/cactus_workflow_getFlowers : *.c *.h ${LIBDIR}/cactusLib.a ${LIBDEPENDS}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/cactus_workflow_getFlowers cactus_workflow_getFlowers.c ${LIBDIR}/cactusLib.a ${LDLIBS}

${BINDIR}/cactus_workflow_extendFlowers : *.c *.h ${LIBDIR}/cactusLib.a ${LIBDEPENDS}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/cactus_workflow_extendFlowers cactus_workflow_extendFlowers.c ${LIBDIR}/cactusLib.a ${LDLIBS}

${BINDIR}/cactus_workflow_flowerStats : *.c *.h ${LIBDIR}/cactusLib.a ${LIBDEPENDS}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/cactus_workflow_flowerStats cactus_workflow_flowerStats.c ${LIBDIR}/cactusLib.a ${LDLIBS}

${BINDIR}/cactus_workflow_convertAlignmentCoordinates : *.c *.h ${LIBDIR}/cactusLib.a ${LIBDEPENDS}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/cactus_workflow_convertAlignmentCoordinates cactus_workflow_convertAlignmentCoordinates.c ${LIBDIR}/cactusLib.a ${LDLIBS}

${BINDIR}/cactus_secondaryDatabase : *.c *.h ${LIBDIR}/cactusLib.a ${LIBDEPENDS}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/cactus_secondaryDatabase cactus_secondaryDatabase.c ${LIBDIR}/cactusLib.a ${LDLIBS}

${BINDIR}/docker_test_script : docker_test_script.py
	cp docker_test_script.py ${BINDIR}/docker_test_script
	chmod +x ${BINDIR}/docker_test_script

clean :  
	rm -f *.o
	rm -f ${BINDIR}/cactus_workflow.py ${BINDIR}/cactus_workflow_getFlowers ${BINDIR}/cactus_workflow_extendFlowers ${BINDIR}/cactus_workflow_flowerStats ${BINDIR}/cactus_workflow_convertAlignmentCoordinates ${BINDIR}/cactus_secondaryDatabase ${BINDIR}/docker_test_script
