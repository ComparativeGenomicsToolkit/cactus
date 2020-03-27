rootPath = ../..
include ${rootPath}/include.mk

treelibSrc = src/*.c
treelibInc = inc
treelibObj = *.o

CFLAGS += -DNDEBUG -I${treelibInc}

all: libtree.a

libtree.a: ${treelibSrc}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -c ${treelibSrc}
	${AR} rc libtree.a ${treelibObj}
	${RANLIB} libtree.a
	rm -f ${treelibObj}

clean:
	rm -f libtree.a ${treelibObj}
