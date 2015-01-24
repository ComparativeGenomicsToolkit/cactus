rootPath = ../
include ../include.mk

all: 
	cd lastz-distrib-1.03.54 && make install
	mv ${binPath}/lastz ${binPath}/cactus_lastz
	mv ${binPath}/lastz_D ${binPath}/cactus_lastz_D
	cd ../last-543 && make install
	


clean: 
	cd lastz-distrib-1.03.54 && make clean 
	rm -rf ${binPath}/cactus_lastz ${binPath}/cactus_lastz_D
 
