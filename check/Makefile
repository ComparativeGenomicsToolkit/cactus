rootPath = ..
include ${rootPath}/include.mk

CFLAGS += ${tokyoCabinetIncl} ${hiredisIncl}

all: all_libs all_progs
all_libs: 
all_progs: all_libs
	${MAKE} ${BINDIR}/cactus_check 

${BINDIR}/cactus_check : *.c ${LIBDIR}/cactusLib.a ${LIBDEPENDS}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/cactus_check cactus_check.c ${LIBDIR}/cactusLib.a ${LDLIBS}

clean :
	rm -f *.o
	rm -f ${BINDIR}/cactus_check 
