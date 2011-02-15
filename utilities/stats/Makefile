rootPath = ../../
include ../../include.mk

all : ${libPath}/cactusTreeStats.a  ${binPath}/cactus_treeStats ${binPath}/cactus_treeStatsToLatexTables.py

${binPath}/cactus_treeStats : *.c *.h ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_treeStats main.c treeStats.c ${libPath}/cactusLib.a ${basicLibs}

${binPath}/cactus_treeStatsToLatexTables.py : cactus_treeStatsToLatexTables.py
	cp cactus_treeStatsToLatexTables.py ${binPath}/cactus_treeStatsToLatexTables.py
	chmod +x ${binPath}/cactus_treeStatsToLatexTables.py


${libPath}/cactusTreeStats.a : treeStats.c treeStats.h ${basicLibsDependencies}
	${cxx} ${cflags} -I ${libPath}/ -c treeStats.c
	ar rc cactusTreeStats.a *.o
	ranlib cactusTreeStats.a 
	rm *.o
	mv cactusTreeStats.a ${libPath}/
	cp treeStats.h ${libPath}/

clean :
	rm -f *.o
	rm -f ${binPath}/cactus_treeStats ${binPath}/cactus_treeStatsToLatexTables.py ${libPath}/cactusTreeStats.a
