rootPath = ../
include ../include.mk

all : ${binPath}/cactus_progressive.py 

${binPath}/cactus_progressive.py : cactus_progressive.py
	cp cactus_progressive.py ${binPath}/cactus_progressive.py
	chmod +x ${binPath}/cactus_progressive.py

clean :  
	rm -f *.o
	rm -f ${binPath}/cactus_progressive.py 
 