// Copyright 2012 Risa Kawaguchi
// Copyright 2013, 2014 Martin C. Frith

//
// 1. UnsplitAlignment -> container of maf file information
//

#ifndef CBRC_UNSPLIT_ALIGNMENT_HH
#define CBRC_UNSPLIT_ALIGNMENT_HH

#include <string>
#include <vector>

namespace cbrc {

typedef std::vector<std::string>::iterator       StringIt;
typedef std::vector<std::string>::const_iterator StringCi;

class UnsplitAlignment {
public:
    StringIt linesBeg;
    StringIt linesEnd;
    const char *qname;
    unsigned qstart;
    unsigned qend;
    char qstrand;
    unsigned rstart;
    unsigned rend;
    const char *rname;
    const char *ralign;
    const char *qalign;
    const char *qQual;
    UnsplitAlignment(){}
    UnsplitAlignment(StringIt linesBegIn, StringIt linesEndIn)
      : linesBeg(linesBegIn), linesEnd(linesEndIn) { init(); }
    void init();
};

void flipMafStrands(StringIt linesBeg, StringIt linesEnd);

std::vector<std::string> mafSlice(StringCi linesBeg, StringCi linesEnd,
				  unsigned alnBeg, unsigned alnEnd);

void mafSliceBeg(const char* rAln, const char* qAln,
		 unsigned qBeg, unsigned& qSliceBeg, unsigned& alnBeg);

void mafSliceEnd(const char* rAln, const char* qAln,
		 unsigned qEnd, unsigned& qSliceEnd, unsigned& alnEnd);

void printMaf(const std::vector<std::string>& maf);

std::string pLineFromProbs(const std::vector<double>& p);

}

#endif
