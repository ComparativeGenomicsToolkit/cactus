#ifndef ALGO_BLAST_GUMBEL_PARAMS__INCLUDED_NJN_DYNPROGPROB
#define ALGO_BLAST_GUMBEL_PARAMS__INCLUDED_NJN_DYNPROGPROB

/* $Id: njn_dynprogprob.hpp 183505 2010-02-18 16:10:58Z boratyng $
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's offical duties as a United States Government employee and
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
* ===========================================================================*/

/*****************************************************************************

File name: njn_dynprogprob.hpp

Author: John Spouge

Contents: 

******************************************************************************/

#include <corelib/ncbitype.h>
#include <corelib/ncbi_limits.h>
#include <corelib/ncbidbg.hpp>

#include "njn_dynprogprobproto.hpp"

BEGIN_NCBI_SCOPE
BEGIN_SCOPE(blast)

BEGIN_SCOPE(Njn)

    class DynProgProb : public DynProgProbProto {

        // DynProgProb performs updates for probabilities in a dynamic programming computation.

        //    The object expands storage as necessary to hold all probabilities.
        //
        // The object behaves as follows:
        //
        // Default:
        //    (1) The initial state of the dynamic programming computation is 0 with probability 1.0.
        //
        // Behavior:
        //    (2) If input_ is the computation's input, it replaces oldValue_ with ValueFct (oldValue_, input_)
        //    (3) The dynamic programming function can be reset with setValueFct (...).
        //    (4) The probability for the input can be reset with setInput (...).
        //    (5) The probability of input_ = [0, dimInputProb_) is inputProb_ [input_].
        //    (6) getProb (Int4 i_) returns the probability corresponding to the Int4 value i_.

        public:

        static const size_t VALUE_BEGIN;
        static const size_t ARRAY_CAPACITY;

        inline DynProgProb ( // range for Int4 values = [valueLower_, valueUpper_) 
        ValueFct *valueFct_ = 0, // function for updating dynamic programming values
        size_t dimInputProb_ = 0, 
        const double *inputProb_ = 0, // array of input states : d_inputProb_p [0...dimInputProb - 1]
        // The following behave like arguments to clear ().
        Int4 valueLower_ = 0, // lower Int4 value corresponding to the "probability" array
        Int4 valueUpper_ = 0, // one beyond present upper Int4 value corresponding to the "probability" array
        const double *prob_ = 0) // "probabilities" prob [valueLower_, valueUpper_) corresponding to the Int4s
        // default prob_ == 0 assigns prob_ [0] = 1.0
        // if (valueLower_ == 0 && valueUpper_ == 0) prob_ [-(ARRAY_CAPACITY / 2) + 1...ARRAY_CAPACITY / 2]
        //    otherwise
        //                                           prob_ [valueLower_...valueUpper_) 
            : d_step (0), d_arrayCapacity (0), d_valueBegin (0), 
            d_valueLower (0), d_valueUpper (0),
            d_valueFct (0), d_dimInputProb (0), d_inputProb_p (0)
        {
            d_array_p [0] = d_array_p [1] = 0;
            setValueFct (valueFct_);
            setInput (dimInputProb_, inputProb_);

            clear (valueLower_, valueUpper_, prob_);
        }

        inline DynProgProb (const DynProgProb &dynProgProb_)
            : d_step (0), d_arrayCapacity (0), d_valueBegin (0), 
            d_valueLower (0), d_valueUpper (0),
            d_valueFct (0), d_dimInputProb (0), d_inputProb_p (0)
        {
            copy (dynProgProb_);
        }

        virtual inline ~DynProgProb () 
        {
            free ();
            freeInput ();
        }

        virtual inline operator bool () // ? is the object ready for computation ?
        const {
            return getArrayCapacity () != 0 && 
                d_valueFct && d_dimInputProb != 0 && d_inputProb_p;
        }

        virtual inline DynProgProb &operator= (const DynProgProb &dynProgProb_)
        {
            if (this != &dynProgProb_) copy (dynProgProb_);
            return *this;
        }

        virtual inline void copy (const DynProgProb &dynProgProb_)
        {
            copy (dynProgProb_.getStep (), 
            dynProgProb_.getArray (), dynProgProb_.getArrayCapacity (), 
            dynProgProb_.getValueBegin (), dynProgProb_.getValueLower (), dynProgProb_.getValueUpper (),
            dynProgProb_.getValueFct (), dynProgProb_.getDimInputProb (), dynProgProb_.getInputProb ());
        }

        virtual void copy (
            size_t step_, // current index : starts at 0 
            const double *const *array_, // two corresponding arrays of probabilities 
        size_t arrayCapacity_, // present capacity of the array
        Int4 valueBegin_ = 0, // lower limit for Int4 values in the array (an offset)
        Int4 valueLower_ = 0, // present lower Int4 value in the array
        Int4 valueUpper_ = 0, // one beyond present upper Int4 value in the array
        ValueFct *valueFct_ = 0, // function for updating dynamic programming values
        size_t dimInputProb_ = 0, 
        const double *inputProb_ = 0); // array of input states : d_inputProb_p [0...dimInputProb - 1]

        virtual void clear ( // restarts the computation
        Int4 valueLower_, // lower Int4 value corresponding to the "probability" array
        Int4 valueUpper_ = 0, // one beyond present upper Int4 value corresponding to the "probability" array
        const double *prob_ = 0); // "probabilities" prob_ [valueLower_, valueUpper_) corresponding to the Int4s
        // default prob_ == 0 assigns prob_ [0] = 1.0
        // assumes prob_ [valueLower_, valueUpper_) 

        virtual inline void clear () {clear (0);}

        virtual inline void setValueFct (ValueFct *valueFct_) // function for updating dynamic programming values
        {
            d_valueFct = valueFct_;
        }

        virtual void setInput (
        size_t dimInputProb_, 
        const double *inputProb_); // array of input states : d_inputProb_p [0...dimInputProb - 1]

        virtual void update (); // updates dynamic prog probs 
        // assert (getValueFct ());
        // assert (getDimInputProb ());
        // assert (getInputProb ());

        virtual inline double getProb (Int4 value_) const // probability value
        {
            _ASSERT (getArray ());
            _ASSERT (getArray () [getStep () % 2]);
            if (value_ < getValueBegin ()) return 0.0;
            if (getValueEnd () <= value_) return 0.0;
            return getArray () [getStep () % 2][getArrayPos (value_)];
        }

        virtual inline size_t getStep () const {return d_step;} // current index : starts at 0 

        virtual inline const double *const *getArray () const {return d_array_p;} // two corresponding arrays of probabilities d_array_p [0,1][0...d_arrayCapacity - 1]
        virtual inline size_t getArrayCapacity () const {return d_arrayCapacity;} // # (different values)
        virtual inline Int4 getValueBegin () const {return d_valueBegin;} // lower limit for Int4 values in the array (an offset)
        virtual inline Int4 getValueLower () const {return d_valueLower;} // present lower Int4 value in the array
        virtual inline Int4 getValueUpper () const {return d_valueUpper;} // one beyond present upper Int4 value in the array

        virtual inline ValueFct *getValueFct ()  const {return d_valueFct;} // function for updating dynamic programming values
        virtual inline size_t getDimInputProb ()  const {return d_dimInputProb;}
        virtual inline const double *getInputProb ()  const {return d_inputProb_p;} // array of input states : d_inputProb_p [0...dimInputProb - 1]

        private:

        size_t d_step; // current index for time-step : starts at 0 
            double *d_array_p [2]; // two corresponding arrays of probabilities d_array_p [0,1][0...d_arrayCapacity - 1]
        // d_array_p [0...1][0...d_valueBound [1] - d_valueBound [0] - 1]
        size_t d_arrayCapacity; // present capacity of the array
        Int4 d_valueBegin; // lower limit for Int4 values in the array (an offset)
        Int4 d_valueLower; // present lower Int4 value in the array
        Int4 d_valueUpper; // one beyond present upper Int4 value in the array

        // parameters for update (which might be constant throughout the calculation)
        ValueFct *d_valueFct; // function for updating dynamic programming values
        size_t d_dimInputProb; 
        double *d_inputProb_p; // array of input states : d_inputProb_p [0...dimInputProb - 1]

        virtual void initInput (size_t dimInputProb_); // array of input states : d_inputProb_p [0...dimInputProb - 1]
        virtual void freeInput ();

        virtual inline Int4 getValue (size_t arrayPos_) const { // value corresponding to array position 
            return static_cast <Int4> (arrayPos_) + getValueBegin ();
        }
        
        virtual void init (size_t arrayCapacity_); // range for d_array_p [0,1][0...arrayCapacity_ - 1]
        virtual void free ();

        protected:

        virtual void clear (
        Int4 valueBegin_, // lower limit for Int4 values in the array (an offset)
        size_t arrayCapacity_); // new array capacity 

        virtual inline Int4 getArrayPos (Int4 value_) const // offset for array position containing value_
        { // no range-checking
            return value_ - getValueBegin ();
        }

        virtual inline Int4 getValueEnd () const // one beyond largest possible Int4 value in present range
        {
            return getValue (getArrayCapacity ());
        }

        void reserve (size_t arrayCapacity_); // new array capacity
        // increases capacity of and copies d_array_p, while updating other variables
        
        virtual void setValueBegin (Int4 valueBegin_); // lowest possible Int4 value in the array

        virtual inline size_t &lgetStep () {return d_step;} // current index : starts at 0 
        virtual inline double **lgetArray () {return d_array_p;} // two corresponding arrays of probabilities d_array_p [0,1][0...d_arrayCapacity - 1]
        virtual inline size_t &lgetArrayCapacity () {return d_arrayCapacity;} // # (different values)
        virtual inline Int4 &lgetValueBegin () {return d_valueBegin;} // lower limit for Int4 values in the array (an offset)
        virtual inline Int4 &lgetValueLower () {return d_valueLower;} // present lower Int4 value in the array
        virtual inline Int4 &lgetValueUpper () {return d_valueUpper;} // one beyond present upper Int4 value in the array
    };

END_SCOPE(Njn)

END_SCOPE(blast)
END_NCBI_SCOPE

#endif //! ALGO_BLAST_GUMBEL_PARAMS__INCLUDED_NJN_DYNPROGPROB
