// Copyright 2008 Michiaki Hamada

#include "LambdaCalculator.hh"
#include <vector>
#include <cassert>

namespace lambda
{
extern "C" {
#include "CA_code/lambda_calculator.h"
}
}

namespace cbrc{

void LambdaCalculator::setBad(){
  lambda_ = -1;
  letterProbs1_.clear();
  letterProbs2_.clear();
}

void LambdaCalculator::calculate( const int matrix[MAT][MAT], int alphSize ){
  assert( alphSize < MAT );

  // We need to pass the parameters as 1-based pointers, hence the +1s
  // and -1s.

  std::vector< double > x( alphSize * alphSize + 1 );
  std::vector< const double* > y( alphSize + 1 );

  for( int i = 0; i < alphSize; ++i ){
    y[ i+1 ] = &x[ i * alphSize ];
    for( int j = 0; j < alphSize; ++j ){
      x[ i * alphSize + j + 1 ] = matrix[i][j];
    }
  }

  letterProbs1_.resize(alphSize);
  letterProbs2_.resize(alphSize);

  lambda_ = lambda::calculate_lambda( &y[0], alphSize,
                                      &letterProbs1_[0] - 1,
                                      &letterProbs2_[0] - 1 );

  if( lambda_ < 0 ) setBad();
}

}  // end namespace
