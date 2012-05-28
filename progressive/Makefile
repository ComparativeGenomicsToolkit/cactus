rootPath = ../
include ../include.mk

all : ${binPath}/cactus_progressive.py ${binPath}/cactus_createMultiCactusProject.py ${binPath}/cactus_resetProjectPaths.py 

${binPath}/cactus_progressive.py : cactus_progressive.py
	cp cactus_progressive.py ${binPath}/cactus_progressive.py
	chmod +x ${binPath}/cactus_progressive.py
	
${binPath}/cactus_createMultiCactusProject.py : cactus_createMultiCactusProject.py
	cp cactus_createMultiCactusProject.py ${binPath}/cactus_createMultiCactusProject.py
	chmod +x ${binPath}/cactus_createMultiCactusProject.py
	
${binPath}/cactus_resetProjectPaths.py : cactus_resetProjectPaths.py
	cp cactus_resetProjectPaths.py ${binPath}/cactus_resetProjectPaths.py
	chmod +x ${binPath}/cactus_resetProjectPaths.py

clean :  
	rm -f *.o
	rm -f ${binPath}/cactus_progressive.py 
	rm -f ${binPath}/cactus_createMultiCactusProject.py
	rm -f ${binPath}/cactus_resetProjectPaths.py
 