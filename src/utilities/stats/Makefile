rootPath = ../../../
include ../../../include.mk

all :  ${binPath}/cactus_treeStats ${binPath}/cactus_treeStatsPlotter.py ${binPath}/cactus_treeStatsScatterPlots.py ${binPath}/cactus_treeStatsToLatexTables.py

${binPath}/cactus_treeStats : *.c *.h ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_treeStats main.c treeStats.c ${libPath}/cactusLib.a ${basicLibs}

${binPath}/cactus_treeStatsPlotter.py  : cactus_treeStatsPlotter.py
	cp cactus_treeStatsPlotter.py ${binPath}/cactus_treeStatsPlotter.py
	chmod +x ${binPath}/cactus_treeStatsPlotter.py

${binPath}/cactus_treeStatsScatterPlots.py : cactus_treeStatsScatterPlots.py 
	cp cactus_treeStatsScatterPlots.py ${binPath}/cactus_treeStatsScatterPlots.py
	chmod +x ${binPath}/cactus_treeStatsScatterPlots.py

${binPath}/cactus_treeStatsToLatexTables.py : cactus_treeStatsToLatexTables.py
	cp cactus_treeStatsToLatexTables.py ${binPath}/cactus_treeStatsToLatexTables.py
	chmod +x ${binPath}/cactus_treeStatsToLatexTables.py

clean :
	rm -f *.o
	rm -f ${binPath}/cactus_treeStats ${binPath}/cactus_treeStatsPlotter.py ${binPath}/cactus_treeStatsScatterPlots.py ${binPath}/cactus_treeStatsToLatexTables.py
