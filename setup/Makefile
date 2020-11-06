rootPath = ..
include ${rootPath}/include.mk

CFLAGS += ${tokyoCabinetIncl} ${hiredisIncl}

all: all_libs all_progs
all_libs: 
all_progs: all_libs
	${MAKE} ${BINDIR}/cactus_setup  

${BINDIR}/cactus_setup : cactus_setup.c ${LIBDEPENDS} ${LIBDIR}/cactusLib.a
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/cactus_setup cactus_setup.c ${LIBDIR}/cactusLib.a ${LDLIBS}

clean : 
	rm -f *.o
	rm -f ${BINDIR}/cactus_setup 
