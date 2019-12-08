rootPath = ../..
include ${rootPath}/include.mk

treelibSrc = src/*.c
treelibInc = inc
treelibObj = *.o

cflags += -DNDEBUG -I${treelibInc}

all: libtree.a

libtree.a: ${treelibSrc}
	${cxx} ${cflags} -c ${treelibSrc}
	ar rc libtree.a ${treelibObj}
	ranlib libtree.a
	rm -f ${treelibObj}

clean:
	rm -f libtree.a ${treelibObj}
