include ../make-include.mak
include version.mak

CC=gcc

# default targets

default: lastz lastz_D

special: lastz_32

builds = lastz lastz_D lastz_32

#---------
# program build
#
# normally creates two versions:
#	lastz:    standard lastz (integer scoring)
#	lastz_D:  lastz with double-float scoring
#
# and can also create this special version:
#	lastz_32: standard lastz with 32-bit sequence indexing
#
#---------
# build-time options (presently only one)
#
# This will enable certain options within the compiled program.  To build with
# one of these option, you would do "make lastz <option_name>=ON".  The result
# will be an executable named lastz, i.e. mostly indistinguishable from the
# executable without this option.  The only place the executable will report
# the option is in response to "lastz --version".
#
# Note that you should do a "make clean" before doing a build with one of these
# options; otherwise there is no guarantee that all modules will get built with
# the same option settings.
#
#	allowBackToBackGaps		gapped_extend.c is modified to allow the opening of
#							.. a delete right after an insert, or vice versa
#
#---------
# Notes re optimization flags:
#
#	On a 2GHz intel core duo iMac:
#
# 		O3 is a definite improvement over no optimization, improving many of
#		the most-used routines down to as low as 60 to 70% of unoptimized run
#		time.
#
#		However, using -funroll-loops actually slowed things down a little.
#
#---------

definedForAll = -Wall -Wextra -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE
flagsFor32    = -Dmax_sequence_index=32 -Dmax_malloc_index=40 -Ddiag_hash_size=4194304

allowBackToBackGaps ?= 0
ifeq ($(allowBackToBackGaps), ON)
	definedForAll += -DallowBackToBackGaps
endif


VERSION_FLAGS= \
	-DVERSION_MAJOR="\"${VERSION_MAJOR}"\" \
	-DVERSION_MINOR="\"${VERSION_MINOR}"\" \
	-DVERSION_SUBMINOR="\"${VERSION_SUBMINOR}"\" \
	-DREVISION_DATE="\"${REVISION_DATE}"\" \
	-DSUBVERSION_REV="\"${SUBVERSION_REV}"\"


CFLAGS = -O3 ${definedForAll} ${VERSION_FLAGS}


srcFiles = lastz infer_scores \
           seeds pos_table quantum seed_search diag_hash                  \
           chain gapped_extend tweener masking                            \
           segment edit_script                                            \
           identity_dist coverage_dist continuity_dist                    \
           output gfa lav axt maf cigar sam genpaf text_align align_diffs \
           utilities dna_utilities sequences capsule

incFiles = lastz.h infer_scores.h \
           seeds.h pos_table.h quantum.h seed_search.h diag_hash.h                            \
           chain.h gapped_extend.h tweener.h masking.h                                        \
           segment.h edit_script.h                                                            \
           identity_dist.h coverage_dist.h continuity_dist.h                                  \
           output.h gfa.h lav.h axt.h maf.h sam.h cigar.h genpaf.h text_align.h align_diffs.h \
           utilities.h dna_utilities.h sequences.h capsule.h

%.o: %.c version.mak ${incFiles}
	${CC} -c ${CFLAGS} -Dscore_type=\'I\' $< -o $@

%_D.o: %.c version.mak ${incFiles}
	${CC} -c ${CFLAGS} -Dscore_type=\'D\' $< -o $@

%_32.o: %.c version.mak ${incFiles}
	${CC} -c ${CFLAGS} ${flagsFor32} $< -o $@


lastz: $(foreach part,${srcFiles},${part}.o)
	${CC} $(foreach part,${srcFiles},${part}.o) -lm -o $@

lastz_D: $(foreach part,${srcFiles},${part}_D.o)
	${CC} $(foreach part,${srcFiles},${part}_D.o) -lm -o $@

lastz_32: $(foreach part,${srcFiles},${part}_32.o)
	${CC} $(foreach part,${srcFiles},${part}_32.o) -lm -o $@

# cleanup

clean: cleano clean_builds clean_test

cleano:
	rm -f *.o

clean_builds:
	rm -f lastz
	rm -f lastz_D
	rm -f lastz_32

# installation;  change installDir to suit your needs (in ../make-include.mak)

install: lastz lastz_D
	${INSTALL} -d      ${installDir}
	${INSTALL} lastz   ${installDir}
	${INSTALL} lastz_D ${installDir}

install_lastz: lastz
	${INSTALL} -d    ${installDir}
	${INSTALL} lastz ${installDir}

install_D: lastz_D
	${INSTALL} -d      ${installDir}
	${INSTALL} lastz_D ${installDir}

install_32: lastz_32
	${INSTALL} -d       ${installDir}
	${INSTALL} lastz_32 ${installDir}

#---------
# testing
#
# A small test to give some comfort level that the program has built properly,
# or that changes you've made to the source code haven't broken it.  If the
# test succeeds, there will be no output from the diff.
#---------

clean_test:
	rm -f ../test_results/base_test*.*

test: lastz
	@rm -f ../test_results/base_test.default.lav
	@./lastz                                     \
	  ../test_data/pseudocat.fa                  \
	  ../test_data/pseudopig.fa                  \
	  | sed "s/\"lastz\.[^ ]* //g"               \
	  > ../test_results/base_test.default.lav
	@diff                                        \
	  ../test_data/base_test.default.lav         \
	  ../test_results/base_test.default.lav

