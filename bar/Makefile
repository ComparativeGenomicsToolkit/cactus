rootPath = ../
include ../include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c
#${libPath}/stCaf.a
commonBarLibs =  ${libPath}/cactusLib.a 
stBarDependencies =  ${commonBarLibs} ${basicLibsDependencies}
stBarLibs = ${commonBarLibs} ${basicLibs}

all : ${libPath}/cactusBarLib.a 
#${binPath}/cactus_bar ${binPath}/cactus_barTests 
 
clean : 
	rm -f ${libPath}/cactusBarLib.a

${libPath}/cactusBarLib.a : ${libSources} ${libHeaders} ${stBarDependencies}
	${cxx} ${cflags} -I inc -I ${libPath}/ -c ${libSources} 
	ar rc cactusBarLib.a *.o
	ranlib cactusBarLib.a 
	rm *.o
	mv cactusBarLib.a ${libPath}/
	cp ${libHeaders} ${libPath}/