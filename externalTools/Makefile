rootPath = ../
include ../include.mk

all : ${binPath}/cactus_lastz

${binPath}/cactus_lastz :
	cd lastz-distrib-1.03.46 && make install
	mv ${binPath}/lastz ${binPath}/cactus_lastz
	
${binPath}/cactus_lastz_D : ${binPath}/cactus_lastz
	mv ${binPath}/lastz_D ${binPath}/cactus_lastz_D

clean: 
	cd lastz-distrib-1.03.46 && make clean 
	rm -rf ${binPath}/cactus_lastz ${binPath}/cactus_lastz_D
