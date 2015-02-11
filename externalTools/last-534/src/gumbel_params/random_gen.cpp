/*  $Id: random_gen.cpp 103491 2007-05-04 17:18:18Z kazimird $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * Authors: Clifford Clausen, Denis Vakatov, Jim Ostell, Jonathan Kans,
 *          Greg Schuler
 * Contact: Clifford Clausen
 *
 * File Description:
 *   CRandom is a lagged Fibonacci (LFG) random number generator (RNG)
 *   with lags 33 and 13, modulus 2^31, and operation '+'. It is a slightly
 *   modified version of Nlm_Random() found in the NCBI C toolkit.
 *   It generates uniform random numbers between 0 and 2^31 - 1 (inclusive).
 *
 *   CRandom has been tested using the Diehard RNG test package
 *   developed by George Marsaglia, Prof., Dept of Statistics, Florida
 *   State University. CRandom in particular, and LFG type RNGs in general, 
 *   cannot pass all of the Diehard RNG tests. Specifically, it fails the
 *   "Birthday" test as do other LFG RNGs. CRandom performs as well as
 *   other LFG RNGS, as provided in the Diehard test package. The LFG
 *   class of RNGs was chosen as the RNG for the NCBI C++ Toolkit as it
 *   provides the best tradeoff between time to generate a random number
 *   and performance on tests for randomness.
 *   
 *   For a download of Diehard software and documentation
 *   see http://stat.fsu.edu/~geo/diehard.html. 
 *
 *   For further information also see
 *   http://random.mat.sbg.ac.at/,  
 *   http://csep1.phy.ornl.gov/rn/rn.html, and
 *   http://www.agner.org/random/.
 *
 *   Some relevant papers are:
 *   1. Hellekalek, P.: "Inversive pseudorandom number generators: concepts
 *   results, and links", In Alexopoulos, C and Kang, K, and Lilegdon, WR,
 *   and Goldsman, D, editor(s), Proceedings of the 1995 Winter Simulation
 *   Conference, pp 255-262, 1995.
 *   2. Leeb, H: "Random Numbers for Computer Simulation", Master's thesis,
 *   University of Salzburg, 1995.
 *   3. Marsaglia, G., "A Current View of Random Number Generators",
 *   Proceedings of 16th Symposium on the Interface, Atlanta, 1984, Elsevier
 *   Press.
 *   4. Marsaglia, G. "Monkey Tests for Random Number Generators", Computers
 *   & Mathematics with Applications, Vol 9, pp. 1-10, 1993.
 *
 *   For a list of other published papers, see
 *   http://random.mat.sbg.ac.at/literature and
 *   http://www.evensen.org/marsaglia/.
 *
 *   class CRandom:: 
 */

#include <ncbi_pch.hpp>
#include <util/random_gen.hpp>


BEGIN_NCBI_SCOPE


const size_t CRandom::kStateSize = sizeof(CRandom::sm_State)
    / sizeof(CRandom::sm_State[0]);

const size_t kStateOffset = 12;

const CRandom::TValue CRandom::sm_State[kStateSize] = {
    0xd53f1852,  0xdfc78b83,  0x4f256096,  0xe643df7,
    0x82c359bf,  0xc7794dfa,  0xd5e9ffaa,  0x2c8cb64a,
    0x2f07b334,  0xad5a7eb5,  0x96dc0cde,  0x6fc24589,
    0xa5853646,  0xe71576e2,  0xdae30df,   0xb09ce711,
    0x5e56ef87,  0x4b4b0082,  0x6f4f340e,  0xc5bb17e8,
    0xd788d765,  0x67498087,  0x9d7aba26,  0x261351d4,
    0x411ee7ea,  0x393a263,   0x2c5a5835,  0xc115fcd8,
    0x25e9132c,  0xd0c6e906,  0xc2bc5b2d,  0x6c065c98,
    0x6e37bd55
};


CRandom::CRandom(void)
{
    Reset();
}


CRandom::CRandom(TValue seed)
{
    SetSeed(seed);
}


void CRandom::Reset(void)
{
    _ASSERT(sizeof(sm_State) / sizeof(sm_State[0]) == kStateSize);
    _ASSERT(kStateOffset < kStateSize);

    for (size_t i = 0;  i < kStateSize;  ++i) {
        m_State[i] = sm_State[i];
    }

    m_RJ = &m_State[kStateOffset];
    m_RK = &m_State[kStateSize - 1];
}


void CRandom::SetSeed(TValue seed)
{
    _ASSERT(kStateOffset < kStateSize);

    m_State[0] = m_Seed = seed;

    // linear congruential initializer
    for (size_t i = 1;  i < kStateSize;  ++i) {
        m_State[i] = 1103515245 * m_State[i-1] + 12345;
    }

    m_RJ = &m_State[kStateOffset];
    m_RK = &m_State[kStateSize - 1];

    for (size_t i = 0;  i < 10 * kStateSize;  ++i) {
        GetRand();
    }
}


END_NCBI_SCOPE
