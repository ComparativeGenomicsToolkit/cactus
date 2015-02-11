CXX = g++
CC  = gcc

CXXFLAGS = -O3 -Wall -Wextra -Wcast-qual -Wswitch-enum -Wundef	\
-Wcast-align -Wno-long-long -ansi -pedantic
# -Wconversion
# -fomit-frame-pointer ?

CFLAGS = -Wall -O2

DBOBJ = Alphabet.o MultiSequence.o CyclicSubsetSeed.o	\
SubsetSuffixArray.o LastdbArguments.o io.o fileMap.o	\
SubsetSuffixArraySort.o MultiSequenceQual.o lastdb.o

ALOBJ = Alphabet.o MultiSequence.o CyclicSubsetSeed.o			\
SubsetSuffixArray.o LastalArguments.o io.o fileMap.o ScoreMatrix.o	\
DiagonalTable.o SegmentPair.o Alignment.o GappedXdropAligner.o		\
SegmentPairPot.o AlignmentPot.o GeneralizedAffineGapCosts.o		\
Centroid.o LambdaCalculator.o TwoQualityScoreMatrix.o			\
OneQualityScoreMatrix.o QualityPssmMaker.o GeneticCode.o		\
gaplessXdrop.o gaplessPssmXdrop.o gaplessTwoQualityXdrop.o		\
SubsetSuffixArraySearch.o AlignmentWrite.o MultiSequenceQual.o		\
GappedXdropAlignerPssm.o GappedXdropAligner2qual.o			\
GappedXdropAligner3frame.o lastal.o CA_code/lambda_calculator.o

EXOBJ = Alphabet.o LastexArguments.o io.o ScoreMatrix.o lastex.o	\
gumbel_params/mcf_local_alignment_evaluer.o				\
gumbel_params/njn_dynprogprob.o gumbel_params/njn_dynprogproblim.o	\
gumbel_params/njn_dynprogprobproto.o gumbel_params/njn_ioutil.o		\
gumbel_params/njn_localmaxstat.o					\
gumbel_params/njn_localmaxstatmatrix.o					\
gumbel_params/njn_localmaxstatutil.o gumbel_params/njn_matrix.o		\
gumbel_params/random_gen.o gumbel_params/sls_alp.o			\
gumbel_params/sls_alp_data.o gumbel_params/sls_alp_regression.o		\
gumbel_params/sls_alp_sim.o gumbel_params/sls_pvalues.o

SPOBJ = Alphabet.o MultiSequence.o fileMap.o split/cbrc_linalg.o	\
split/last-split.o split/cbrc_split_aligner.o split/last-split-main.o	\
split/cbrc_unsplit_alignment.o

PPOBJ = last-pair-probs.o last-pair-probs-main.o io.o

MBOBJ = last-merge-batches.o

ALL = lastdb lastal lastex last-split last-merge-batches last-pair-probs

all: $(ALL)

lastdb: $(DBOBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(DBOBJ)

lastal: $(ALOBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(ALOBJ)

lastex: $(EXOBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(EXOBJ)

last-split: $(SPOBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(SPOBJ)

last-pair-probs: $(PPOBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(PPOBJ)

last-merge-batches: $(MBOBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(MBOBJ)

.SUFFIXES:
.SUFFIXES: .o .c .cc .cpp

.c.o:
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<

.cc.o:
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Igumbel_params -I. -c -o $@ $<

.cpp.o:
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Igumbel_params -c -o $@ $<

clean:
	rm -f $(ALL) *.o */*.o

VERSION = \"`hg id -n`\"
VERSION2 = \"`grep latesttagdistance ../.hg_archival.txt | cut -d' ' -f2`\"
VERSION3 = \"UNKNOWN\"

version.hh: FORCE
	if test -e ../.hg ; \
	then echo $(VERSION) | cmp -s $@ - || echo $(VERSION) > $@ ; \
	elif test -e ../.hg_archival.txt ; \
	then echo $(VERSION2) | cmp -s $@ - || echo $(VERSION2) > $@ ; \
	else test -e $@ || echo $(VERSION3) > $@ ; \
	fi

FORCE:

depend:
	sed '/[m][v]/q' makefile > m
	$(CXX) -MM *.cc >> m
	$(CC) -MM CA_code/lambda_calculator.c | sed 's|.*:|CA_code/&|' >> m
	$(CXX) -MM -Igumbel_params */*.cpp | sed 's|.*:|gumbel_params/&|' >> m
	$(CXX) -MM -I. split/*.cc | sed 's|.*:|split/&|' >> m
	mv m makefile
Alignment.o: Alignment.cc Alignment.hh ScoreMatrixRow.hh SegmentPair.hh \
 Alphabet.hh Centroid.hh GappedXdropAligner.hh \
 GeneralizedAffineGapCosts.hh OneQualityScoreMatrix.hh GeneticCode.hh \
 TwoQualityScoreMatrix.hh
AlignmentPot.o: AlignmentPot.cc AlignmentPot.hh Alignment.hh \
 ScoreMatrixRow.hh SegmentPair.hh
AlignmentWrite.o: AlignmentWrite.cc Alignment.hh ScoreMatrixRow.hh \
 SegmentPair.hh GeneticCode.hh MultiSequence.hh VectorOrMmap.hh Mmap.hh \
 fileMap.hh stringify.hh Alphabet.hh
Alphabet.o: Alphabet.cc Alphabet.hh
Centroid.o: Centroid.cc Centroid.hh GappedXdropAligner.hh \
 ScoreMatrixRow.hh GeneralizedAffineGapCosts.hh SegmentPair.hh \
 OneQualityScoreMatrix.hh GappedXdropAlignerInl.hh
CyclicSubsetSeed.o: CyclicSubsetSeed.cc CyclicSubsetSeed.hh io.hh \
 stringify.hh
DiagonalTable.o: DiagonalTable.cc DiagonalTable.hh
GappedXdropAligner.o: GappedXdropAligner.cc GappedXdropAligner.hh \
 ScoreMatrixRow.hh GappedXdropAlignerInl.hh
GappedXdropAligner2qual.o: GappedXdropAligner2qual.cc \
 GappedXdropAligner.hh ScoreMatrixRow.hh GappedXdropAlignerInl.hh \
 TwoQualityScoreMatrix.hh
GappedXdropAligner3frame.o: GappedXdropAligner3frame.cc \
 GappedXdropAligner.hh ScoreMatrixRow.hh GappedXdropAlignerInl.hh
GappedXdropAligner3framePssm.o: GappedXdropAligner3framePssm.cc \
 GappedXdropAligner.hh ScoreMatrixRow.hh GappedXdropAlignerInl.hh
GappedXdropAlignerPssm.o: GappedXdropAlignerPssm.cc GappedXdropAligner.hh \
 ScoreMatrixRow.hh GappedXdropAlignerInl.hh
GeneralizedAffineGapCosts.o: GeneralizedAffineGapCosts.cc \
 GeneralizedAffineGapCosts.hh
GeneticCode.o: GeneticCode.cc GeneticCode.hh Alphabet.hh
LambdaCalculator.o: LambdaCalculator.cc LambdaCalculator.hh \
 CA_code/lambda_calculator.h
LastalArguments.o: LastalArguments.cc LastalArguments.hh \
 SequenceFormat.hh stringify.hh
LastdbArguments.o: LastdbArguments.cc LastdbArguments.hh \
 SequenceFormat.hh stringify.hh
LastexArguments.o: LastexArguments.cc LastexArguments.hh stringify.hh
MultiSequence.o: MultiSequence.cc MultiSequence.hh ScoreMatrixRow.hh \
 VectorOrMmap.hh Mmap.hh fileMap.hh stringify.hh io.hh
MultiSequenceQual.o: MultiSequenceQual.cc MultiSequence.hh \
 ScoreMatrixRow.hh VectorOrMmap.hh Mmap.hh fileMap.hh stringify.hh
OneQualityScoreMatrix.o: OneQualityScoreMatrix.cc \
 OneQualityScoreMatrix.hh ScoreMatrixRow.hh qualityScoreUtil.hh \
 stringify.hh
QualityPssmMaker.o: QualityPssmMaker.cc QualityPssmMaker.hh \
 ScoreMatrixRow.hh qualityScoreUtil.hh stringify.hh
ScoreMatrix.o: ScoreMatrix.cc ScoreMatrix.hh io.hh
SegmentPair.o: SegmentPair.cc SegmentPair.hh
SegmentPairPot.o: SegmentPairPot.cc SegmentPairPot.hh SegmentPair.hh
SubsetSuffixArray.o: SubsetSuffixArray.cc SubsetSuffixArray.hh \
 CyclicSubsetSeed.hh VectorOrMmap.hh Mmap.hh fileMap.hh stringify.hh \
 io.hh
SubsetSuffixArraySearch.o: SubsetSuffixArraySearch.cc \
 SubsetSuffixArray.hh CyclicSubsetSeed.hh VectorOrMmap.hh Mmap.hh \
 fileMap.hh stringify.hh
SubsetSuffixArraySort.o: SubsetSuffixArraySort.cc SubsetSuffixArray.hh \
 CyclicSubsetSeed.hh VectorOrMmap.hh Mmap.hh fileMap.hh stringify.hh
TwoQualityScoreMatrix.o: TwoQualityScoreMatrix.cc \
 TwoQualityScoreMatrix.hh ScoreMatrixRow.hh qualityScoreUtil.hh \
 stringify.hh
fileMap.o: fileMap.cc fileMap.hh stringify.hh
gaplessPssmXdrop.o: gaplessPssmXdrop.cc gaplessPssmXdrop.hh \
 ScoreMatrixRow.hh
gaplessTwoQualityXdrop.o: gaplessTwoQualityXdrop.cc \
 gaplessTwoQualityXdrop.hh TwoQualityScoreMatrix.hh ScoreMatrixRow.hh
gaplessXdrop.o: gaplessXdrop.cc gaplessXdrop.hh ScoreMatrixRow.hh
io.o: io.cc io.hh
last-pair-probs-main.o: last-pair-probs-main.cc last-pair-probs.hh \
 stringify.hh
last-pair-probs.o: last-pair-probs.cc last-pair-probs.hh io.hh \
 stringify.hh
lastal.o: lastal.cc LastalArguments.hh SequenceFormat.hh \
 QualityPssmMaker.hh ScoreMatrixRow.hh OneQualityScoreMatrix.hh \
 TwoQualityScoreMatrix.hh qualityScoreUtil.hh stringify.hh \
 LambdaCalculator.hh GeneticCode.hh SubsetSuffixArray.hh \
 CyclicSubsetSeed.hh VectorOrMmap.hh Mmap.hh fileMap.hh Centroid.hh \
 GappedXdropAligner.hh GeneralizedAffineGapCosts.hh SegmentPair.hh \
 AlignmentPot.hh Alignment.hh SegmentPairPot.hh ScoreMatrix.hh \
 Alphabet.hh MultiSequence.hh DiagonalTable.hh gaplessXdrop.hh \
 gaplessPssmXdrop.hh gaplessTwoQualityXdrop.hh io.hh version.hh
lastdb.o: lastdb.cc LastdbArguments.hh SequenceFormat.hh \
 SubsetSuffixArray.hh CyclicSubsetSeed.hh VectorOrMmap.hh Mmap.hh \
 fileMap.hh stringify.hh Alphabet.hh MultiSequence.hh ScoreMatrixRow.hh \
 io.hh qualityScoreUtil.hh version.hh
lastex.o: lastex.cc LastexArguments.hh ScoreMatrix.hh Alphabet.hh io.hh \
 gumbel_params/mcf_local_alignment_evaluer.hpp \
 gumbel_params/sls_pvalues.hpp gumbel_params/sls_normal_distr_array.hpp \
 version.hh
CA_code/lambda_calculator.o: CA_code/lambda_calculator.c CA_code/nrutil.h \
 CA_code/nrutil.c CA_code/ludcmp.c CA_code/lubksb.c \
 CA_code/lambda_calculator.h
gumbel_params/mcf_local_alignment_evaluer.o: \
 gumbel_params/mcf_local_alignment_evaluer.cpp \
 gumbel_params/mcf_local_alignment_evaluer.hpp \
 gumbel_params/sls_pvalues.hpp gumbel_params/corelib/ncbistl.hpp \
 gumbel_params/common/ncbi_export.h gumbel_params/corelib/ncbitype.h \
 gumbel_params/corelib/ncbi_limits.h \
 gumbel_params/sls_normal_distr_array.hpp \
 gumbel_params/njn_localmaxstatmatrix.hpp \
 gumbel_params/njn_localmaxstat.hpp gumbel_params/sls_alp_sim.hpp \
 gumbel_params/sls_alp_data.hpp gumbel_params/util/random_gen.hpp \
 gumbel_params/corelib/ncbistd.hpp gumbel_params/corelib/ncbidbg.hpp \
 gumbel_params/sls_alp_regression.hpp gumbel_params/sls_alp.hpp
gumbel_params/njn_dynprogprob.o: gumbel_params/njn_dynprogprob.cpp \
 gumbel_params/ncbi_pch.hpp gumbel_params/corelib/ncbitype.h \
 gumbel_params/corelib/ncbi_limits.h gumbel_params/njn_dynprogprob.hpp \
 gumbel_params/corelib/ncbidbg.hpp gumbel_params/njn_dynprogprobproto.hpp \
 gumbel_params/corelib/ncbistl.hpp gumbel_params/common/ncbi_export.h \
 gumbel_params/njn_memutil.hpp gumbel_params/njn_ioutil.hpp
gumbel_params/njn_dynprogproblim.o: gumbel_params/njn_dynprogproblim.cpp \
 gumbel_params/ncbi_pch.hpp gumbel_params/njn_dynprogproblim.hpp \
 gumbel_params/corelib/ncbitype.h gumbel_params/corelib/ncbi_limits.h \
 gumbel_params/njn_dynprogprob.hpp gumbel_params/corelib/ncbidbg.hpp \
 gumbel_params/njn_dynprogprobproto.hpp gumbel_params/corelib/ncbistl.hpp \
 gumbel_params/common/ncbi_export.h gumbel_params/njn_memutil.hpp \
 gumbel_params/njn_ioutil.hpp
gumbel_params/njn_dynprogprobproto.o: gumbel_params/njn_dynprogprobproto.cpp \
 gumbel_params/ncbi_pch.hpp gumbel_params/njn_dynprogprobproto.hpp \
 gumbel_params/corelib/ncbistl.hpp gumbel_params/common/ncbi_export.h \
 gumbel_params/corelib/ncbitype.h
gumbel_params/njn_ioutil.o: gumbel_params/njn_ioutil.cpp gumbel_params/ncbi_pch.hpp \
 gumbel_params/corelib/ncbistre.hpp gumbel_params/njn_ioutil.hpp \
 gumbel_params/corelib/ncbistl.hpp gumbel_params/common/ncbi_export.h
gumbel_params/njn_localmaxstat.o: gumbel_params/njn_localmaxstat.cpp \
 gumbel_params/ncbi_pch.hpp gumbel_params/corelib/ncbistl.hpp \
 gumbel_params/common/ncbi_export.h gumbel_params/sls_alp_data.hpp \
 gumbel_params/util/random_gen.hpp gumbel_params/corelib/ncbistd.hpp \
 gumbel_params/corelib/ncbidbg.hpp gumbel_params/corelib/ncbitype.h \
 gumbel_params/sls_alp_regression.hpp gumbel_params/corelib/ncbi_limits.h \
 gumbel_params/njn_localmaxstat.hpp gumbel_params/njn_memutil.hpp \
 gumbel_params/njn_ioutil.hpp gumbel_params/njn_dynprogproblim.hpp \
 gumbel_params/njn_dynprogprob.hpp gumbel_params/njn_dynprogprobproto.hpp \
 gumbel_params/njn_function.hpp gumbel_params/njn_doubletype.hpp \
 gumbel_params/njn_integer.hpp gumbel_params/njn_localmaxstatutil.hpp \
 gumbel_params/njn_matrix.hpp gumbel_params/njn_approx.hpp \
 gumbel_params/njn_vector.hpp
gumbel_params/njn_localmaxstatmatrix.o: gumbel_params/njn_localmaxstatmatrix.cpp \
 gumbel_params/ncbi_pch.hpp gumbel_params/corelib/ncbistl.hpp \
 gumbel_params/common/ncbi_export.h \
 gumbel_params/njn_localmaxstatmatrix.hpp \
 gumbel_params/njn_localmaxstat.hpp gumbel_params/corelib/ncbitype.h \
 gumbel_params/njn_localmaxstatutil.hpp \
 gumbel_params/corelib/ncbi_limits.h gumbel_params/njn_matrix.hpp \
 gumbel_params/njn_approx.hpp gumbel_params/njn_doubletype.hpp \
 gumbel_params/njn_ioutil.hpp gumbel_params/njn_vector.hpp \
 gumbel_params/njn_memutil.hpp
gumbel_params/njn_localmaxstatutil.o: gumbel_params/njn_localmaxstatutil.cpp \
 gumbel_params/ncbi_pch.hpp gumbel_params/njn_localmaxstatutil.hpp \
 gumbel_params/corelib/ncbistl.hpp gumbel_params/common/ncbi_export.h \
 gumbel_params/corelib/ncbitype.h gumbel_params/corelib/ncbi_limits.h \
 gumbel_params/njn_matrix.hpp gumbel_params/njn_approx.hpp \
 gumbel_params/njn_doubletype.hpp gumbel_params/njn_ioutil.hpp \
 gumbel_params/njn_vector.hpp gumbel_params/njn_dynprogproblim.hpp \
 gumbel_params/njn_dynprogprob.hpp gumbel_params/corelib/ncbidbg.hpp \
 gumbel_params/njn_dynprogprobproto.hpp gumbel_params/njn_integer.hpp \
 gumbel_params/njn_memutil.hpp gumbel_params/njn_root.hpp \
 gumbel_params/njn_function.hpp gumbel_params/sls_alp_data.hpp \
 gumbel_params/util/random_gen.hpp gumbel_params/corelib/ncbistd.hpp \
 gumbel_params/sls_alp_regression.hpp
gumbel_params/njn_matrix.o: gumbel_params/njn_matrix.cpp gumbel_params/ncbi_pch.hpp \
 gumbel_params/njn_matrix.hpp gumbel_params/corelib/ncbitype.h \
 gumbel_params/njn_approx.hpp gumbel_params/corelib/ncbistl.hpp \
 gumbel_params/common/ncbi_export.h gumbel_params/njn_doubletype.hpp \
 gumbel_params/njn_ioutil.hpp gumbel_params/njn_vector.hpp
gumbel_params/random_gen.o: gumbel_params/random_gen.cpp gumbel_params/ncbi_pch.hpp \
 gumbel_params/util/random_gen.hpp gumbel_params/corelib/ncbistd.hpp \
 gumbel_params/corelib/ncbistl.hpp gumbel_params/common/ncbi_export.h \
 gumbel_params/corelib/ncbidbg.hpp gumbel_params/corelib/ncbitype.h
gumbel_params/sls_alp.o: gumbel_params/sls_alp.cpp gumbel_params/ncbi_pch.hpp \
 gumbel_params/sls_alp.hpp gumbel_params/corelib/ncbistl.hpp \
 gumbel_params/common/ncbi_export.h gumbel_params/sls_alp_data.hpp \
 gumbel_params/util/random_gen.hpp gumbel_params/corelib/ncbistd.hpp \
 gumbel_params/corelib/ncbidbg.hpp gumbel_params/corelib/ncbitype.h \
 gumbel_params/sls_alp_regression.hpp gumbel_params/corelib/ncbi_limits.h
gumbel_params/sls_alp_data.o: gumbel_params/sls_alp_data.cpp gumbel_params/ncbi_pch.hpp \
 gumbel_params/sls_alp_data.hpp gumbel_params/corelib/ncbistl.hpp \
 gumbel_params/common/ncbi_export.h gumbel_params/util/random_gen.hpp \
 gumbel_params/corelib/ncbistd.hpp gumbel_params/corelib/ncbidbg.hpp \
 gumbel_params/corelib/ncbitype.h gumbel_params/sls_alp_regression.hpp \
 gumbel_params/corelib/ncbi_limits.h
gumbel_params/sls_alp_regression.o: gumbel_params/sls_alp_regression.cpp \
 gumbel_params/ncbi_pch.hpp gumbel_params/sls_alp_regression.hpp \
 gumbel_params/corelib/ncbistl.hpp gumbel_params/common/ncbi_export.h \
 gumbel_params/corelib/ncbitype.h gumbel_params/corelib/ncbi_limits.h \
 gumbel_params/sls_alp_data.hpp gumbel_params/util/random_gen.hpp \
 gumbel_params/corelib/ncbistd.hpp gumbel_params/corelib/ncbidbg.hpp
gumbel_params/sls_alp_sim.o: gumbel_params/sls_alp_sim.cpp gumbel_params/ncbi_pch.hpp \
 gumbel_params/sls_alp_sim.hpp gumbel_params/corelib/ncbistl.hpp \
 gumbel_params/common/ncbi_export.h gumbel_params/sls_alp_data.hpp \
 gumbel_params/util/random_gen.hpp gumbel_params/corelib/ncbistd.hpp \
 gumbel_params/corelib/ncbidbg.hpp gumbel_params/corelib/ncbitype.h \
 gumbel_params/sls_alp_regression.hpp gumbel_params/corelib/ncbi_limits.h \
 gumbel_params/sls_alp.hpp
gumbel_params/sls_pvalues.o: gumbel_params/sls_pvalues.cpp gumbel_params/ncbi_pch.hpp \
 gumbel_params/sls_pvalues.hpp gumbel_params/corelib/ncbistl.hpp \
 gumbel_params/common/ncbi_export.h gumbel_params/corelib/ncbitype.h \
 gumbel_params/corelib/ncbi_limits.h \
 gumbel_params/sls_normal_distr_array.hpp gumbel_params/sls_alp_data.hpp \
 gumbel_params/util/random_gen.hpp gumbel_params/corelib/ncbistd.hpp \
 gumbel_params/corelib/ncbidbg.hpp gumbel_params/sls_alp_regression.hpp
split/cbrc_linalg.o: split/cbrc_linalg.cc split/cbrc_linalg.hh
split/cbrc_split_aligner.o: split/cbrc_split_aligner.cc \
 split/cbrc_split_aligner.hh split/cbrc_unsplit_alignment.hh \
 split/cbrc_int_exponentiator.hh Alphabet.hh MultiSequence.hh \
 ScoreMatrixRow.hh VectorOrMmap.hh Mmap.hh fileMap.hh stringify.hh \
 split/cbrc_linalg.hh
split/cbrc_unsplit_alignment.o: split/cbrc_unsplit_alignment.cc \
 split/cbrc_unsplit_alignment.hh
split/last-split-main.o: split/last-split-main.cc split/last-split.hh \
 stringify.hh version.hh
split/last-split.o: split/last-split.cc split/last-split.hh \
 split/cbrc_split_aligner.hh split/cbrc_unsplit_alignment.hh \
 split/cbrc_int_exponentiator.hh Alphabet.hh MultiSequence.hh \
 ScoreMatrixRow.hh VectorOrMmap.hh Mmap.hh fileMap.hh stringify.hh
