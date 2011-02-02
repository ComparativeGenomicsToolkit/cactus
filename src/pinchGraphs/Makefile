rootPath = ../../
include ../../include.mk

libSources = pinchGraph.c pinchGraphTest.c
libHeaders = pinchGraph.h
cflags += ${tokyoCabinetIncl}

all : ${libPath}/pinchGraph.a 

${libPath}/pinchGraph.a : ${libSources} ${libHeaders} ${basicLibsDependencies}
	${cxx} ${cflags} -I ${libPath}/ -c ${libSources}
	ar rc pinchGraph.a *.o
	ranlib pinchGraph.a 
	rm *.o
	mv pinchGraph.a ${libPath}/
	cp ${libHeaders} ${libPath}/

clean : 
	rm -f ${libPath}/pinchGraph.a

