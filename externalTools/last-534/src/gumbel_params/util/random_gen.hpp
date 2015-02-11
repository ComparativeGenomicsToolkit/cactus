#ifndef RANDOM_GEN__HPP
#define RANDOM_GEN__HPP

/*  $Id: random_gen.hpp 188609 2010-04-13 12:56:05Z ivanov $
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
 *   CRandom implements a lagged Fibonacci (LFG) random number generator (RNG)
 *   with lags 33 and 13, modulus 2^31, and operation '+'. It is a slightly
 *   modified version of Nlm_Random() found in the NCBI C toolkit.
 *   It generates uniform random numbers between 0 and 2^31 - 1 (inclusive).
 *
 *   More details and literature refs are provided in "random_gen.cpp".
 */

#include <corelib/ncbistd.hpp>


/** @addtogroup RandomGen
 *
 * @{
 */


BEGIN_NCBI_SCOPE


/////////////////////////////////////////////////////////////////////////////
//  CRandom::
//

class NCBI_XUTIL_EXPORT CRandom
{
public:
    // Type of the generated integer value and/or the seed value
    typedef Uint4 TValue;

    // Constructors
    CRandom(void);
    CRandom(TValue seed);

    // Initialize and Seed the random number generator
    void   SetSeed(TValue seed);
    TValue GetSeed(void);

    // Reset random number generator to initial startup condition
    void Reset(void);

    // Get the next random number in the interval [0..GetMax()] (inclusive)
    TValue GetRand(void);
    TValue GetRand(TValue min_value, TValue max_value); 

    // The max. value GetRand() returns
    static TValue GetMax(void);

private:
    // Static array used to initialize "m_State", and its size
    static const TValue sm_State[33];
    static const size_t kStateSize;

    // Instance data members
    TValue  m_State[sizeof(sm_State) / sizeof(sm_State[0])];
    TValue* m_RJ;
    TValue* m_RK;
    TValue  m_Seed;

private:
    // prevent copying
    CRandom(const CRandom&);
    CRandom& operator=(const CRandom&);
};


/* @} */


/////////////////////////////////////////////////////////////////////////////
//  IMPLEMENTATION of INLINE functions
/////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////
//  CRandom::
//

inline CRandom::TValue CRandom::GetSeed(void)
{
    return m_Seed;
}


inline CRandom::TValue CRandom::GetRand(void)
{
    register TValue r;

    r = *m_RK;
    r += *(m_RJ--);
    *(m_RK--) = r;

    if (m_RK < m_State) {
        m_RK = &m_State[sizeof(m_State) / sizeof(m_State[0]) - 1];
    }
    else if (m_RJ < m_State) {
        m_RJ = &m_State[sizeof(m_State) / sizeof(m_State[0]) - 1];
    }

    return (r >> 1) & 0x7fffffff;  // discard the least-random bit
}

inline CRandom::TValue CRandom::GetRand(TValue min_value, TValue max_value)
{
  return min_value + (GetRand() % (max_value - min_value + 1));
}


inline CRandom::TValue CRandom::GetMax(void)
{
    return 0x7fffffff;
}


END_NCBI_SCOPE

#endif  /* RANDOM_GEN__HPP */
