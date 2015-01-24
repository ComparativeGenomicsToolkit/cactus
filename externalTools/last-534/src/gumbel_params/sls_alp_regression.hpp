/* $Id: sls_alp_regression.hpp 183505 2010-02-18 16:10:58Z boratyng $
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

File name: sls_alp_regression.hpp

Author: Sergey Sheetlin

Contents: Regression methods

******************************************************************************/

#ifndef ALGO_BLAST_GUMBEL_PARAMS__INCLUDED_SLS_ALP_REGRESSION
#define ALGO_BLAST_GUMBEL_PARAMS__INCLUDED_SLS_ALP_REGRESSION

#include <complex>
#include <iostream>
#include <map>
#include <vector>
#include <fstream>
#include <float.h>
#include <algorithm>

#include <corelib/ncbistl.hpp>
#include <corelib/ncbitype.h>
#include <corelib/ncbi_limits.h>

#include "sls_alp_data.hpp"


BEGIN_NCBI_SCOPE
BEGIN_SCOPE(blast)


BEGIN_SCOPE(Sls)

        typedef double function_type(double x_,void* func_number_);


        class alp_reg{

        

        public:


                        alp_reg(//constructor
                                );


                        ~alp_reg();//destructor

                        static void find_tetta_general(
                        function_type *func_,
                        void* func_pointer_,
                        double a_,//[a,b] is the interval for search of equation solution
                        double b_,
                        Int4 n_partition_,
                        double eps_,
                        std::vector<double> &res_);

                        static double find_single_tetta_general(
                        function_type *func_,
                        void* func_pointer_,
                        double a_,//[a,b] is the interval for search of equation solution
                        double b_,
                        double eps_);

                        static void correction_of_errors(
                        double *errors_,
                        Int4 number_of_elements_);


                        static void robust_regression_sum_with_cut_LSM(
                        Int4 min_length_,
                        Int4 number_of_elements_,
                        double *values_,
                        double *errors_,
                        bool cut_left_tail_,
                        bool cut_right_tail_,
                        double y_,
                        double &beta0_,
                        double &beta1_,
                        double &beta0_error_,
                        double &beta1_error_,
                        Int4 &k1_opt_,
                        Int4 &k2_opt_,
                        bool &res_was_calculated_);

                        static double function_for_robust_regression_sum_with_cut_LSM(
                        double *values_,
                        double *errors_,
                        Int4 number_of_elements_,
                        Int4 k_start_,
                        double c_,
                        double &beta0_,
                        double &beta1_,
                        double &beta0_error_,
                        double &beta1_error_,
                        bool &res_was_calculated_);

                        static void robust_regression_sum_with_cut_LSM_beta1_is_defined(
                        Int4 min_length_,
                        Int4 number_of_elements_,
                        double *values_,
                        double *errors_,
                        bool cut_left_tail_,
                        bool cut_right_tail_,
                        double y_,
                        double &beta0_,
                        double beta1_,
                        double &beta0_error_,
                        double beta1_error_,
                        Int4 &k1_opt_,
                        Int4 &k2_opt_,
                        bool &res_was_calculated_);

                        static double function_for_robust_regression_sum_with_cut_LSM_beta1_is_defined(
                        double *values_,
                        double *errors_,
                        Int4 number_of_elements_,
                        Int4 k_start_,
                        double c_,
                        double &beta0_,
                        double beta1_,
                        double &beta0_error_,
                        double beta1_error_,
                        bool &res_was_calculated_);

                        static double error_of_the_lg(//lg(v1_)
                        double v1_,
                        double v1_error_);

                        static double error_of_the_sqrt(//sqrt(v1_)
                        double v1_,
                        double v1_error_);

                        static double error_of_the_ratio(//v1_/v2_
                        double v1_,
                        double v1_error_,
                        double v2_,
                        double v2_error_);

                        static double error_of_the_product(//v1_*v2_
                        double v1_,
                        double v1_error_,
                        double v2_,
                        double v2_error_);

                        static double error_of_the_sum(//v1_+v2_
                        double v1_,
                        double v1_error_,
                        double v2_,
                        double v2_error_);


                        inline static double sqrt_for_errors(
                                double x_)
                                {
                                        if(x_<=0)
                                        {
                                                return 0.0;
                                        }
                                        else
                                        {
                                                return sqrt(x_);
                                        };
                                };


                        static double median(
                        Int4 dim_,
                        double *array_);

                        static double robust_sum(
                        double *values,
                        Int4 dim,
                        Int4 N_points,
                        bool *&remove_flag);


        };

END_SCOPE(Sls)

END_SCOPE(blast)
END_NCBI_SCOPE

#endif //! ALGO_BLAST_GUMBEL_PARAMS__INCLUDED_SLS_ALP_REGRESSION
