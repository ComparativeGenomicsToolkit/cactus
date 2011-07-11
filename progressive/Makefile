rootPath = ../
include ../include.mk

all : ${binPath}/cactus_progressive.py ${binPath}/cactus_prepKtserver.py 

${binPath}/cactus_progressive.py : cactus_progressive.py
	cp cactus_progressive.py ${binPath}/cactus_progressive.py
	chmod +x ${binPath}/cactus_progressive.py
	
${binPath}/cactus_prepKtserver.py : cactus_prepKtserver.py
	cp cactus_prepKtserver.py ${binPath}/cactus_prepKtserver.py
	chmod +x ${binPath}/cactus_prepKtserver.py

clean :  
	rm -f *.o
	rm -f ${binPath}/cactus_progressive.py 
	rm -f ${binPath}/cactus_prepKtserver.py
 