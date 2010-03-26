makefile:

include ../../../include.mk
binPath = ../../../bin
libPath = ../../../lib

libSources = 3_Absorb3edge2x.c
libHeaders = 3_Absorb3edge2x.h

all : ${libPath}/3EdgeConnected.a 

${libPath}/3EdgeConnected.a : ${libSources} ${libHeaders} ${libPath}/sonLib.a
	${cxx} ${cflags} -I ${libPath}/ -c ${libSources}
	ar rc 3EdgeConnected.a *.o
	ranlib 3EdgeConnected.a 
	rm *.o
	mv 3EdgeConnected.a ${libPath}/
	cp ${libHeaders} ${libPath}/

clean : 
	rm -f ${libPath}/3EdgeConnected.a 

