rootPath = ..
include ${rootPath}/include.mk

libTests = tests/*.c
libSources = impl/*.c
libHeaders = inc/*.h

commonCafLibs = ${LIBDIR}/stCactusToHal.a ${LIBDIR}/stReference.a ${LIBDIR}/cactusBarLib.a ${LIBDIR}/stCaf.a ${LIBDIR}/stCactusSetup.a ${sonLibDir}/stPinchesAndCacti.a ${sonLibDir}/matchingAndOrdering.a ${sonLibDir}/3EdgeConnected.a ${sonLibDir}/cPecanLib.a  ${LIBDIR}/cactusLib.a ${LIBDIR}/stPaf.a
# simde (included via abPOA) doesn't compile with --Werror --pedantic
CFLAGS:=$(filter-out --pedantic,$(CFLAGS))

all: all_libs all_progs
all_libs: 
all_progs: all_libs
	${MAKE} ${BINDIR}/stPipelineTests ${BINDIR}/cactus_consolidated ${BINDIR}/docker_test_script

${BINDIR}/stPipelineTests : ${libTests} ${LIBDIR}/cactusLib.a ${LIBDEPENDS}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/stPipelineTests ${libTests} ${LIBDIR}/cactusLib.a ${LDLIBS}

${BINDIR}/cactus_consolidated : cactus_consolidated.c ${LIBDEPENDS} ${commonCafLibs} ${libSources} ${libHeaders}
# the -Wno-unused-function is required to include abpoa.h with CGL_DEBUG defined
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/cactus_consolidated cactus_consolidated.c ${libSources} ${commonCafLibs} ${LDLIBS} -Wno-unused-function

${BINDIR}/docker_test_script : docker_test_script.py
	cp docker_test_script.py ${BINDIR}/docker_test_script
	chmod +x ${BINDIR}/docker_test_script

clean :  
	rm -f *.o
	rm -f ${BINDIR}/cactus_workflow.py ${BINDIR}/cactus_consolidated ${BINDIR}/docker_test_script
