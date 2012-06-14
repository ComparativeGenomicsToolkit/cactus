rootPath = ../
include ../include.mk

all : ${binPath}/dbTestScript 

${binPath}/dbTestScript  : ${basicLibsDependencies} *.c
	${cxx} ${cflags} -I inc -I${libPath} -Wno-error -o ${binPath}/dbTestScript *.c ${basicLibs}
	
${binPath}/dbTestScript.py : *.py
	mv dbTestScript.py ${binPath}/dbTestScript.py
	
clean :
	rm -rf ${binPath}/dbTestScript