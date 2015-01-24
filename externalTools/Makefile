rootPath = ../
include ../include.mk

all: 
	cd lastz-distrib-1.03.54 && make install
	mv ${binPath}/lastz ${binPath}/cactus_lastz
	mv ${binPath}/lastz_D ${binPath}/cactus_lastz_D
	cd last-534 && make install
	mv ${binPath}/last-534/lastdb ${binPath}/cactus_lastdb
	mv ${binPath}/last-534/lastal ${binPath}/cactus_lastal
	


clean: 
	cd lastz-distrib-1.03.54 && make clean 
	rm -rf ${binPath}/cactus_lastz ${binPath}/cactus_lastz_D
	cd last-534 && make clean
	rm -rf ${binPath}/cactus_lastdb ${binPath}/cactus_lastal
 
