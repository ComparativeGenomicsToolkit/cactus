rootPath = ../../
include ../../include.mk
refDir = ../referenceUtils

progs = $(notdir $(wildcard cactus_mafToReferenceSeq.py))
targets = ${progs:%=${binPath}/%}

#all :  ${targets} ${binPath}/cactus_MAFGenerator ${binPath}/cactus_augmentedMaf 
all :  ${targets} ${binPath}/cactus_MAFGenerator ${libPath}/cactusMafs.a ${binPath}/cactus_augmentedMaf

${binPath}/%: %
	@mkdir -p $(dir $@)
	cp -f $< $@
	chmod 775 $@

#${libPath}/cactus_addReferenceSeq.h: ${refDir}/cactus_addReferenceSeq.h
#	cp ${refDir}/cactus_addReferenceSeq.h ${libPath}

#${libPath}/cactus_addReferenceSeq.o : ${refDir}/cactus_addReferenceSeq.c
#	${cxx} ${cflags} -c ${refDir}/cactus_addReferenceSeq.c  -I ${libPath} -o ${libPath}/cactus_addReferenceSeq.o

#${binPath}/cactus_MAFGenerator : ${libPath}/cactus_addReferenceSeq.h ${libPath}/cactus_addReferenceSeq.o cactus_MAFGenerator.c mafs.c ${libPath}/cactusLib.a ${basicLibsDependencies}
#	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_MAFGenerator cactus_MAFGenerator.c mafs.c  ${libPath}/cactus_addReferenceSeq.o ${libPath}/cactusLib.a ${basicLibs}

#mycommon.o : mycommon.h mycommon.c
	#cp mycommon.h ${libPath}
	#${cxx} ${cflags} -c mycommon.c  -I ${libPath} -o mycommon.o

${binPath}/cactus_MAFGenerator : cactus_MAFGenerator.c mafs.c ${libPath}/cactusLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_MAFGenerator cactus_MAFGenerator.c mafs.c ${libPath}/cactusUtils.a ${libPath}/cactusLib.a ${basicLibs}

${binPath}/cactus_augmentedMaf : ${libPath}/cactusUtils.h ${libPath}/cactusUtils.a cactus_augmentedMaf.c ${basicLibsDependencies}
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/cactus_augmentedMaf cactus_augmentedMaf.c ${libPath}/cactusUtils.a ${libPath}/cactusLib.a ${basicLibs}

clean :
	rm -f *.o
	rm -f ${binPath}/cactus_MAFGenerator ${libPath}/cactusMafs.a
	rm -f ${progs:%=${binPath}/%}
	
${libPath}/cactusMafs.a : mafs.c cactusMafs.h ${basicLibsDependencies}
	${cxx} ${cflags} -I ${libPath}/ -c mafs.c
	ar rc cactusMafs.a *.o 
	ranlib cactusMafs.a 
	rm *.o
	mv cactusMafs.a ${libPath}/
	cp cactusMafs.h ${libPath}/ 
