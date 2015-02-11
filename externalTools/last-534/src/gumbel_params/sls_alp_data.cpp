/* $Id: sls_alp_data.cpp 189387 2010-04-21 17:55:34Z boratyng $
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

File name: sls_alp_data.cpp

Author: Sergey Sheetlin

Contents: Input data for the ascending ladder points simulation

******************************************************************************/


#include <ncbi_pch.hpp>

#include <ncbi_pch.hpp>
#include "sls_alp_data.hpp"

USING_NCBI_SCOPE;
USING_SCOPE(blast);
USING_SCOPE(Sls);


alp_data::alp_data(//constructor
                   Int4 open_,
                   Int4 epen_,
                   double eps_lambda_,
                   double eps_K_,
                   const Int4 *const *scoreMatrix,
                   Int4 numResidues,
                   const vector<double> &seq1ResidueProbs,
                   const vector<double> &seq2ResidueProbs,
                   double max_time_,
                   double max_mem_,
                   Int4 rand_)
{


        bool ee_error_flag=false;
        error ee_error("",0);

        d_smatr=NULL;
        d_RR1=NULL;
        d_RR1_sum=NULL;
        d_RR1_sum_elements=NULL;

        d_RR2=NULL;
        d_RR2_sum=NULL;
        d_RR2_sum_elements=NULL;

        d_is=NULL;
        d_r_i_dot=NULL;
        d_r_dot_j=NULL;

        d_rand_all=NULL;

        d_rand_object=NULL;




        try
        {
        try
        {
                
                d_sentinels_flag=false;


                d_memory_size_in_MB=0;

                #ifndef NCBI_OS_MSWIN

                #else
                        _CrtMemCheckpoint( &d_s1 );
                #endif


                Int4 number_of_AA_RR1;
                Int4 number_of_AA_RR2;

                

                Int4 i,j;
                d_number_of_AA_smatr=numResidues;

                if(d_number_of_AA_smatr<=0)
                {
                        throw error("Error - number of letters in the scoring matrix file must be greater than 0\n",3);
                };

                get_memory_for_matrix(d_number_of_AA_smatr,d_smatr);


                for(i=0;i<d_number_of_AA_smatr;i++)
                {
                        for(j=0;j<d_number_of_AA_smatr;j++)
                        {
                                d_smatr[i][j]=scoreMatrix[i][j];
                        };
                };

                d_smatr_symmetric_flag=false;


                read_RR(
                seq1ResidueProbs,
                d_RR1,
                d_RR1_sum,
                d_RR1_sum_elements,
                number_of_AA_RR1);


                read_RR(
                seq2ResidueProbs,
                d_RR2,
                d_RR2_sum,
                d_RR2_sum_elements,
                number_of_AA_RR2);


                if(number_of_AA_RR1==d_number_of_AA_smatr)
                {
                        d_number_of_AA=d_number_of_AA_smatr;
                }
                else
                {
                        throw error("Number of letters is different for the scoring matrix and probabilities array\n",3);
                };

                if(number_of_AA_RR2!=d_number_of_AA_smatr)
                {
                        throw error("Number of letters is different for the scoring matrix and probabilities array\n",3);
                };


                d_open=open_+epen_;

                d_epen=epen_;

                d_max_time=max_time_;

                d_max_mem=max_mem_;

                d_eps_lambda=eps_lambda_;

                d_eps_K=eps_K_;

                d_out_file_name="test.out";
                d_minimum_realizations_number=40;

                d_rand_all=new struct_for_randomization;
                alp_data::assert_mem(d_rand_all);
                d_memory_size_in_MB+=sizeof(struct_for_randomization)/mb_bytes;


                //randomization
                Uint4 random_factor=rand_;
                d_rand_flag=false;

#if 0
                CRef<CGumbelParamsRandDiagnostics> AdvancedParams_tmp
                    = rand_params_;

                if(AdvancedParams_tmp.Empty())
                {
                        random_factor=(Uint4)time(NULL);
                        #ifndef NCBI_OS_MSWIN //UNIX program
                                struct timeval tv;
                                struct timezone tz;
                                gettimeofday(&tv, &tz);
                                random_factor+=tv.tv_usec*10000000;
                        #else
                                struct _timeb timebuffer;
                                char *timeline;
                                _ftime( &timebuffer );
                                timeline = ctime( & ( timebuffer.time ) );
                                random_factor+=timebuffer.millitm*10000000;
                        #endif

                        d_rand_flag=false;

                }
                else
                {
                        d_rand_flag=true;
                        if(d_rand_flag)
                        {
                            random_factor=AdvancedParams_tmp->GetRandomSeed();

                            Int4 size=AdvancedParams_tmp->GetFirstStagePrelimReNumbers().size();
                                d_rand_all->d_first_stage_preliminary_realizations_numbers_ALP.resize(size);

                                Int4 i;
                                for(i=0;i<size;i++)
                                {
                                    d_rand_all->d_first_stage_preliminary_realizations_numbers_ALP[i]=AdvancedParams_tmp->GetFirstStagePrelimReNumbers()[i];
                                };


                                size=AdvancedParams_tmp->GetPrelimReNumbers().size();
                                d_rand_all->d_preliminary_realizations_numbers_ALP.resize(size);
                                for(i=0;i<size;i++)
                                {
                                    d_rand_all->d_preliminary_realizations_numbers_ALP[i]=AdvancedParams_tmp->GetPrelimReNumbers()[i];
                                };


                                size=AdvancedParams_tmp->GetPrelimReNumbersKilling().size();
                                d_rand_all->d_preliminary_realizations_numbers_killing.resize(size);
                                for(i=0;i<size;i++)
                                {
                                    d_rand_all->d_preliminary_realizations_numbers_killing[i]=AdvancedParams_tmp->GetPrelimReNumbersKilling()[i];
                                };



                                d_rand_all->d_total_realizations_number_with_ALP=AdvancedParams_tmp->GetTotalReNumber();
                                d_rand_all->d_total_realizations_number_with_killing=AdvancedParams_tmp->GetTotalReNumberKilling();


                        };
                };
#endif


                d_random_factor=random_factor;

                d_rand_object=new CRandom;
                d_rand_object->SetSeed(d_random_factor);

                


                d_is=new importance_sampling(
                        this,
                        d_open,
                        d_epen,
                        d_number_of_AA,
                        d_smatr,
                        d_RR1,
                        d_RR2);

                alp_data::assert_mem(d_is);

                d_memory_size_in_MB+=sizeof(d_is)/mb_bytes;

                d_r_i_dot=new double[d_number_of_AA];
                alp_data::assert_mem(d_r_i_dot);
                d_r_dot_j=new double[d_number_of_AA];
                alp_data::assert_mem(d_r_dot_j);
                Int4 k;
                for(k=0;k<d_number_of_AA;k++)
                {
                        d_r_i_dot[k]=0;
                        if(d_RR1[k]!=0)
                        {
                                Int4 i;
                                for(i=0;i<d_number_of_AA;i++)
                                {
                                        if(d_RR2[i]!=0)
                                        {
                                                d_r_i_dot[k]+=d_is->d_exp_s[k][i]*d_RR2[i];
                                        };
                                };
                        };
                };

                for(k=0;k<d_number_of_AA;k++)
                {
                        d_r_dot_j[k]=0;
                        if(d_RR2[k]!=0)
                        {
                                Int4 i;
                                for(i=0;i<d_number_of_AA;i++)
                                {
                                        if(d_RR1[i]!=0)
                                        {
                                                d_r_dot_j[k]+=d_is->d_exp_s[i][k]*d_RR1[i];
                                        };
                                };
                        };
                };


                d_memory_size_in_MB+=(double)(sizeof(double)*d_number_of_AA*2.0)/mb_bytes;

                double tmp_size1=kMax_I4;

                double tmp_size=Tmin((double)(tmp_size1),
                        (

                        (double)mb_bytes*d_max_mem/(double)d_minimum_realizations_number
                        )
                        /(
                        (double)(sizeof(double)*12)+(double)(sizeof(Int4)*17)
                        )
                        );

                d_dim1_tmp=(Int4)tmp_size;
                d_dim2_tmp=(Int4)tmp_size;
        }
        catch (error er)
        {
                ee_error_flag=true;
                ee_error=er;                
        };
        }
        catch (...)
        { 
                ee_error_flag=true;
                ee_error=error("Internal error in the program\n",4);
        };


        if(ee_error_flag)
        {
                this->~alp_data();
                throw error(ee_error.st,ee_error.error_code);
        };

};

alp_data::alp_data(//constructor
Int4 rand_,//randomization number
Int4 open_,//gap opening penalty
Int4 epen_,//gap extension penalty
string smatr_file_name_,//scoring matrix file name
string RR1_file_name_,//probabilities1 file name
string RR2_file_name_,//probabilities2 file name
double max_time_,//maximum allowed calculation time in seconds
double max_mem_,//maximum allowed memory usage in MB
double eps_lambda_,//relative error for lambda calculation
double eps_K_,//relative error for K calculation
string out_file_name_)//output file name
{

        ifstream frand;
        bool ee_error_flag=false;
        error ee_error("",0);

        d_smatr=NULL;
        d_RR1=NULL;
        d_RR1_sum=NULL;
        d_RR1_sum_elements=NULL;

        d_RR2=NULL;
        d_RR2_sum=NULL;
        d_RR2_sum_elements=NULL;

        d_is=NULL;
        d_r_i_dot=NULL;
        d_r_dot_j=NULL;

        d_rand_all=NULL;



        try
        {
        try
        {
                d_sentinels_flag=false;


                d_memory_size_in_MB=0;

                #ifndef NCBI_OS_MSWIN //UNIX program

                #else
                        _CrtMemCheckpoint( &d_s1 );
                #endif


                Int4 number_of_AA_RR1;
                Int4 number_of_AA_RR2;

                read_smatr(
                smatr_file_name_,
                d_smatr,
                d_number_of_AA_smatr);
                
                

                read_RR(
                RR1_file_name_,
                d_RR1,
                d_RR1_sum,
                d_RR1_sum_elements,
                number_of_AA_RR1);


                read_RR(
                RR2_file_name_,
                d_RR2,
                d_RR2_sum,
                d_RR2_sum_elements,
                number_of_AA_RR2);


                if(number_of_AA_RR1==d_number_of_AA_smatr)
                {
                        d_number_of_AA=d_number_of_AA_smatr;
                }
                else
                {
                        throw error("Number of letters is different in the files "+smatr_file_name_+" and "+RR1_file_name_+"\n",3);
                };

                if(number_of_AA_RR2!=d_number_of_AA_smatr)
                {
                        throw error("Number of letters is different in the files "+smatr_file_name_+" and "+RR2_file_name_+"\n",3);
                };

                Int4 t;
                for(t=0;t<number_of_AA_RR1;t++)
                {
                        if(d_RR1[t]!=d_RR2[t])
                        {
                                d_smatr_symmetric_flag=false;
                                break;
                        };
                };


                check_out_file(out_file_name_);

                d_open=open_+epen_;
                d_epen=epen_;
                d_max_time=max_time_;
                d_max_mem=max_mem_;
                d_eps_lambda=eps_lambda_;
                d_eps_K=eps_K_;
                d_out_file_name=out_file_name_;
                d_minimum_realizations_number=40;

                d_rand_all=new struct_for_randomization;
                alp_data::assert_mem(d_rand_all);
                d_memory_size_in_MB+=sizeof(struct_for_randomization)/mb_bytes;

                //randomization
                Uint4 random_factor=rand_;


                if(random_factor<0)
                {
                        random_factor=(Uint4)time(NULL);
                        #ifndef NCBI_OS_MSWIN //UNIX program
                                struct timeval tv;
                                struct timezone tz;
                                gettimeofday(&tv, &tz);
                                random_factor+=tv.tv_usec*10000000;
                        #else
                                struct _timeb timebuffer;
                                char *timeline;
                                _ftime( &timebuffer );
                                timeline = ctime( & ( timebuffer.time ) );
                                random_factor+=timebuffer.millitm*10000000;
                        #endif

                        d_rand_flag=false;

                }
                else
                {
                        d_rand_flag=true;
                        if(d_rand_flag)
                        {
                                string rand_st="rand_"+alp_data::long_to_string(random_factor)+".out";
                                frand.open(rand_st.data(),ios::in);
                                if(!frand)
                                {
                                        d_rand_flag=false;
                                }
                                else
                                {

                                        Int4 i,size;


                                        
                                        frand>>d_rand_all->d_random_factor;


                                        if((Int4)random_factor!=d_rand_all->d_random_factor)
                                        {
                                                throw error("Unexpected error in randomization seed\n",3);
                                        };



                                        frand>>size;
                                        for(i=0;i<size;i++)
                                        {
                                                Int4 tmp;
                                                frand>>tmp;
                                                d_rand_all->d_first_stage_preliminary_realizations_numbers_ALP.push_back(tmp);
                                        };

                                        frand>>size;
                                        for(i=0;i<size;i++)
                                        {
                                                Int4 tmp;
                                                frand>>tmp;
                                                d_rand_all->d_preliminary_realizations_numbers_ALP.push_back(tmp);
                                        };

                                        frand>>size;
                                        for(i=0;i<size;i++)
                                        {
                                                Int4 tmp;
                                                frand>>tmp;
                                                d_rand_all->d_preliminary_realizations_numbers_killing.push_back(tmp);
                                        };


                                        frand>>d_rand_all->d_total_realizations_number_with_ALP;
                                        frand>>d_rand_all->d_total_realizations_number_with_killing;

                                        frand.close();
                                };
                        };
                };


                d_random_factor=random_factor;

                d_rand_object=new CRandom;
                d_rand_object->SetSeed(d_random_factor);

                


                d_is=new importance_sampling(
                        this,
                        d_open,
                        d_epen,
                        d_number_of_AA,
                        d_smatr,
                        d_RR1,
                        d_RR2);

                alp_data::assert_mem(d_is);

                d_memory_size_in_MB+=sizeof(d_is)/mb_bytes;

                d_r_i_dot=new double[d_number_of_AA];
                alp_data::assert_mem(d_r_i_dot);
                d_r_dot_j=new double[d_number_of_AA];
                alp_data::assert_mem(d_r_dot_j);
                Int4 k;
                for(k=0;k<d_number_of_AA;k++)
                {
                        d_r_i_dot[k]=0;
                        if(d_RR1[k]!=0)
                        {
                                Int4 i;
                                for(i=0;i<d_number_of_AA;i++)
                                {
                                        if(d_RR2[i]!=0)
                                        {
                                                d_r_i_dot[k]+=d_is->d_exp_s[k][i]*d_RR2[i];
                                        };
                                };
                        };
                };

                for(k=0;k<d_number_of_AA;k++)
                {
                        d_r_dot_j[k]=0;
                        if(d_RR2[k]!=0)
                        {
                                Int4 i;
                                for(i=0;i<d_number_of_AA;i++)
                                {
                                        if(d_RR1[i]!=0)
                                        {
                                                d_r_dot_j[k]+=d_is->d_exp_s[i][k]*d_RR1[i];
                                        };
                                };
                        };
                };


                d_memory_size_in_MB+=(double)(sizeof(double)*d_number_of_AA*2.0)/mb_bytes;

                double tmp_size1=kMax_I4;

                double tmp_size=Tmin((double)(tmp_size1),
                        (

                        (double)mb_bytes*d_max_mem/(double)d_minimum_realizations_number
                        )
                        /(
                        (double)(sizeof(double)*12)+(double)(sizeof(Int4)*17)
                        )
                        );

                d_dim1_tmp=(Int4)tmp_size;
                d_dim2_tmp=(Int4)tmp_size;
        }
        catch (error er)
        {
                ee_error_flag=true;
                ee_error=er;                
        };
        }
        catch (...)
        { 
                ee_error_flag=true;
                ee_error=error("Internal error in the program\n",4);
        };

        if(frand.is_open())
        {
                frand.close();
        };

        if(ee_error_flag)
        {
                this->~alp_data();
                throw error(ee_error.st,ee_error.error_code);
        };

};

Int4 alp_data::random_long(
double value_,
Int4 dim_)
{
        if(value_<0||value_>1.0||dim_<=0)
        {
                throw error("Unexpected error",4);
        };

        if(dim_==1)
        {
                return 0;
        };

        Int4 tmp=(Int4)floor(value_*(double)dim_);
        tmp=Tmin(tmp,dim_-1);
        return tmp;
};


alp_data::~alp_data()//destructor
{
        delete d_rand_object;

        delete[]d_RR1;d_RR1=NULL;
        delete[]d_RR1_sum;d_RR1_sum=NULL;
        delete[]d_RR1_sum_elements;d_RR1_sum_elements=NULL;

        delete[]d_RR2;d_RR2=NULL;
        delete[]d_RR2_sum;d_RR2_sum=NULL;
        delete[]d_RR2_sum_elements;d_RR2_sum_elements=NULL;


        d_memory_size_in_MB-=(double)(2.0*sizeof(double)+sizeof(Int4))*(double)d_number_of_AA/mb_bytes;

        delete_memory_for_matrix(d_number_of_AA_smatr,d_smatr);

        delete d_is;d_is=NULL;

        d_memory_size_in_MB-=sizeof(d_is)/mb_bytes;

        delete[]d_r_i_dot;d_r_i_dot=NULL;
        delete[]d_r_dot_j;d_r_dot_j=NULL;
        d_memory_size_in_MB-=(double)(sizeof(double)*d_number_of_AA*2.0)/mb_bytes;

        delete d_rand_all;d_rand_all=NULL;
        d_memory_size_in_MB-=sizeof(struct_for_randomization)/mb_bytes;


};

void alp_data::check_out_file(
        string out_file_name_)
{
        bool ee_error_flag=false;
        error ee_error("",0);
        ifstream f;
        char *str_ch=NULL;

        try
        {
        try
        {
                f.open(out_file_name_.data(),ios::in);
                if(!f)
                {
                        return;
                };

                bool symmetric_case_flag;
                
                string str;
                getline(f,str);
                str_ch=new char[str.length()+1];
                if(!str_ch)
                {
                        throw error("Memory allocation error\n",41);
                };

                Int4 k;
                for(k=0;k<(Int4)str.length();k++)
                {
                        str_ch[k]=str[k];
                };
                str_ch[str.length()]='\0';


                char str_for_test0[]="number of realizations with killing";
                char *test_flag0= strstr(str_ch,str_for_test0);

                if(!test_flag0)
                {
                        throw error("The output file "+out_file_name_+" exists and does not have correct format;\nplease delete the file and rerun the program\n",3);
                };

                char str_for_test[]="0.5*";

                char*test_flag= strstr(str_ch,str_for_test);
                if(test_flag)
                {
                        symmetric_case_flag=true;
                }
                else
                {
                        symmetric_case_flag=false;
                };


                

                if(symmetric_case_flag)
                {
                        if(!d_smatr_symmetric_flag)
                        {
                                throw error("The output file "+out_file_name_+" exists and corresponds to symmetric case; \ncurrent calculation uses non-symmetric parameters;\nplease define another output file name\n",3);
                        };
                };

                if(!symmetric_case_flag)
                {
                        if(d_smatr_symmetric_flag)
                        {
                                throw error("The output file "+out_file_name_+" exists and corresponds to non-symmetric case; \ncurrent calculation uses symmetric parameters;\nplease define another output file name\n",3);
                        };
                };

                f.close();
        }
        catch (error er)
        {
                ee_error_flag=true;
                ee_error=er;                
        };
        }
        catch (...)
        { 
                ee_error_flag=true;
                ee_error=error("Internal error in the program\n",4);
        };

        delete[]str_ch;str_ch=NULL;

        if(f.is_open())
        {
                f.close();
        };

        if(ee_error_flag)
        {
                throw error(ee_error.st,ee_error.error_code);
        };

};


double alp_data::get_allocated_memory_in_MB()
{

        #ifndef NCBI_OS_MSWIN //UNIX program

                return 0;

        #else
                _CrtMemCheckpoint( &d_s2 );

                _CrtMemDifference( &d_s3, &d_s1, &d_s2);

                double total=0;
                int use;
                for (use = 0; use < _MAX_BLOCKS; use++)
                {
                        total+=d_s3.lSizes[use];
                }

                total/=(double)1048576;
                return total;

        #endif

};

// Kludge: limit optimization by 32-bit ICC 10.x to avoid undefined
// references to __svml_exp2.
#if defined(NCBI_COMPILER_ICC)  &&  NCBI_PLATFORM_BITS == 32 \
    &&  NCBI_COMPILER_VERSION >= 1000  &&  NCBI_COMPILER_VERSION < 1100 \
    &&  defined(__OPTIMIZE__)
#  pragma optimization_level 1
#endif
double importance_sampling::lambda_equation(double x_,void* func_number_)
{
        data_for_lambda_equation *data=(data_for_lambda_equation*)func_number_;
        Int4 d_number_of_AA=data->d_number_of_AA;
        Int4** d_smatr=data->d_smatr;
        double *d_RR1=data->d_RR1;
        double *d_RR2=data->d_RR2;

        double res=0;
        Int4 i,j;

        for(i=0;i<d_number_of_AA;i++)
        {
                for(j=0;j<d_number_of_AA;j++)
                {
                        res+=d_RR1[i]*d_RR2[j]*exp(x_*d_smatr[i][j]);
                };
        };

        return res-1.0;
};

void alp_data::read_smatr(
string smatr_file_name_,
Int4 **&smatr_,
Int4 &number_of_AA_smatr_)
{
        bool ee_error_flag=false;
        error ee_error("",0);
        ifstream f;

        try
        {
        try
        {

                Int4 i,j;
                f.open(smatr_file_name_.data(),ios::in);
                if(!f)
                {
                        throw error("Error - file "+smatr_file_name_+" is not found\n",3);
                };

                f>>number_of_AA_smatr_;

                if(number_of_AA_smatr_<=0)
                {
                        throw error("Error - number of letters in the scoring matrix file must be greater than 0\n",3);
                };

                get_memory_for_matrix(number_of_AA_smatr_,smatr_);


                for(i=0;i<number_of_AA_smatr_;i++)
                {
                        for(j=0;j<number_of_AA_smatr_;j++)
                        {
                                f>>smatr_[i][j];
                        };
                };

                f.close();

                bool flag=true;
                for(i=0;i<number_of_AA_smatr_;i++)
                {
                        for(j=0;j<i;j++)
                        {
                                if(smatr_[i][j]!=smatr_[j][i])
                                {
                                        flag=false;
                                };
                        };
                };

                d_smatr_symmetric_flag=flag;

                d_smatr_symmetric_flag=false;

        }
        catch (error er)
        {
                ee_error_flag=true;
                ee_error=er;                
        };
        }
        catch (...)
        { 
                ee_error_flag=true;
                ee_error=error("Internal error in the program\n",4);
        };

        //memory release
        if(f.is_open())
        {
                f.close();
        };

        if(ee_error_flag)
        {
                throw error(ee_error.st,ee_error.error_code);
        };

};

void alp_data::read_RR(
string RR_file_name_,
double *&RR_,
double *&RR_sum_,
Int4 *&RR_sum_elements_,
Int4 &number_of_AA_RR_)
{
        bool ee_error_flag=false;
        error ee_error("",0);
        ifstream f;

        try
        {
        try
        {

                Int4 i;
                f.open(RR_file_name_.data(),ios::in);
                if(!f)
                {
                        throw error("Error - file "+RR_file_name_+" is not found\n",3);
                };

                f>>number_of_AA_RR_;

                if(number_of_AA_RR_<=0)
                {
                        throw error("Error - number of letters in the probabilities file must be greater than 0\n",3);
                };
                
                RR_=new double[number_of_AA_RR_];
                assert_mem(RR_);

                RR_sum_=new double[number_of_AA_RR_];
                assert_mem(RR_sum_);

                RR_sum_elements_=new Int4 [number_of_AA_RR_];
                assert_mem(RR_sum_elements_);

                d_memory_size_in_MB+=(double)(2.0*sizeof(double)+sizeof(Int4))*(double)number_of_AA_RR_/mb_bytes;


                for(i=0;i<number_of_AA_RR_;i++)
                {
                        f>>RR_[i];

                        if(RR_[i]<0)
                        {
                                throw error("Error - input letter's probability number "+long_to_string(i+1)+" is negative\n",3);
                        };

                        if(RR_[i]>1.0)
                        {
                                throw error("Error - input letter's probability number "+long_to_string(i+1)+" is greater than 1.0\n",3);
                        };


                        if(i!=0)
                        {
                                RR_sum_[i]=RR_sum_[i-1]+RR_[i];
                        }
                        else
                        {
                                RR_sum_[i]=RR_[i];
                        };
                        RR_sum_elements_[i]=i;
                };

                if(fabs(RR_sum_[number_of_AA_RR_-1]-1.0)>0.000000000001)
                {
                        //cout<<"Warning: sum of probabilities in the file "<<RR_file_name_<<" is not equal to 1\n\n";
                };


                f.close();
        }
        catch (error er)
        {
                ee_error_flag=true;
                ee_error=er;                
        };
        }
        catch (...)
        { 
                ee_error_flag=true;
                ee_error=error("Internal error in the program\n",4);
        };

        //memory release
        if(f.is_open())
        {
                f.close();
        };

        if(ee_error_flag)
        {
                throw error(ee_error.st,ee_error.error_code);
        };

};

void alp_data::read_RR(
const vector<double> &vector_,
double *&RR_,
double *&RR_sum_,
Int4 *&RR_sum_elements_,
Int4 &number_of_AA_RR_)
{
        bool ee_error_flag=false;
        error ee_error("",0);

        try
        {
        try
        {

                Int4 i;

                number_of_AA_RR_=vector_.size();

                if(number_of_AA_RR_<=0)
                {
                        throw error("Error - number of letters in the probabilities file must be greater than 0\n",3);
                };
                
                RR_=new double[number_of_AA_RR_];
                assert_mem(RR_);

                RR_sum_=new double[number_of_AA_RR_];
                assert_mem(RR_sum_);

                RR_sum_elements_=new Int4 [number_of_AA_RR_];
                assert_mem(RR_sum_elements_);

                d_memory_size_in_MB+=(double)(2.0*sizeof(double)+sizeof(Int4))*(double)number_of_AA_RR_/mb_bytes;


                for(i=0;i<number_of_AA_RR_;i++)
                {
                        RR_[i]=vector_[i];

                        if(RR_[i]<0)
                        {
                                throw error("Error - input letter's probability number "+long_to_string(i+1)+" is negative\n",3);
                        };

                        if(RR_[i]>1.0)
                        {
                                throw error("Error - input letter's probability number "+long_to_string(i+1)+" is greater than 1.0\n",3);
                        };


                        if(i!=0)
                        {
                                RR_sum_[i]=RR_sum_[i-1]+RR_[i];
                        }
                        else
                        {
                                RR_sum_[i]=RR_[i];
                        };
                        RR_sum_elements_[i]=i;
                };

                if(fabs(RR_sum_[number_of_AA_RR_-1]-1.0)>0.000000000001)
                {
                        //cout<<"Warning: sum of probabilities in the file "<<RR_file_name_<<" is not equal to 1\n\n";
                };


        }
        catch (error er)
        {
                ee_error_flag=true;
                ee_error=er;                
        };
        }
        catch (...)
        { 
                ee_error_flag=true;
                ee_error=error("Internal error in the program\n",4);
        };

        if(ee_error_flag)
        {
                throw error(ee_error.st,ee_error.error_code);
        };

};


string alp_data::long_to_string(//convert interer ot string
Int4 number_)
{
        string res_="";
        string tmp_string;
        if(number_>0)
        {
                tmp_string="";
        }
        else
        {
                if(number_==0)
                {
                        tmp_string="";
                }
                else
                {
                        tmp_string="-";
                };
        };
        number_=abs(number_);
        do{
                Int4 reminder=number_%10;
                number_=(number_-reminder)/10;
                res_=digit_to_string(reminder)+res_;
                if (number_==0)
                {
                        break;
                };
        }
        while (true);

        return tmp_string+res_;
};

char alp_data::digit_to_string(//convert interer ot string
Int4 digit_)
{
        switch(digit_)
        {
        case 0:return '0';
        case 1:return '1';
        case 2:return '2';
        case 3:return '3';
        case 4:return '4';
        case 5:return '5';
        case 6:return '6';
        case 7:return '7';
        case 8:return '8';
        case 9:return '9';
        default:return '?';
        };
};




void alp_data::assert_mem(void *pointer_)
{
        if(!pointer_)
        {
                throw error("Memory allocation error\n",41);
        };
};

double alp_data::round(//returns nearest integer to x_
const double &x_)
{
        double x_floor=floor(x_);
        double x_ceil=ceil(x_);
        if(fabs(x_-x_floor)<0.5)
        {
                return x_floor;
        };
        return x_ceil;
};



importance_sampling::importance_sampling(
alp_data *alp_data_,
Int4 open_,
Int4 epen_,
Int4 number_of_AA_,
Int4 **smatr_,
double *RR1_,
double *RR2_)
{
        d_elements=NULL;
        d_elements_values=NULL;

        d_exp_s=NULL;


        d_alp_data=alp_data_;
        if(!d_alp_data)
        {
                throw error("Unexpected error",4);
        };

        bool ee_error_flag=false;
        error ee_error("",0);

        try
        {
        try
        {



                {

                        //calculation of the importance sampling theta

                        data_for_lambda_equation tmp_ptr;
                        tmp_ptr.d_number_of_AA=number_of_AA_;
                        tmp_ptr.d_RR1=RR1_;
                        tmp_ptr.d_RR2=RR2_;
                        tmp_ptr.d_smatr=smatr_;

                        //calculate maximum of smatr_ elements
                        Int4 smatr_max=smatr_[0][0];
                        Int4 smatr_max_i=0;
                        Int4 smatr_max_j=0;
                        Int4 smatr_min=smatr_[0][0];

                        Int4 smatr_pos_max=kMin_I4;
                        Int4 smatr_neg_min=kMax_I4;

                        double eps=0.00001;
                        double threshold=DBL_MIN*10.0;

                        double aver_score=0;
                        Int4 i,j;
                        for(i=0;i<number_of_AA_;i++)
                        {
                                for(j=0;j<number_of_AA_;j++)
                                {
                                        if(RR1_[i]*RR2_[j]<=threshold)
                                        {
                                                continue;
                                        };

                                        aver_score+=RR1_[i]*RR2_[j]*smatr_[i][j];

                                        if(smatr_max<smatr_[i][j])
                                        {
                                                smatr_max=smatr_[i][j];
                                                smatr_max_i=i;
                                                smatr_max_j=j;
                                        };
                                        smatr_min=alp_data::Tmin(smatr_min,smatr_[i][j]);
                                        

                                        if(smatr_[i][j]>0)
                                        {
                                                smatr_pos_max=alp_data::Tmax(smatr_pos_max,smatr_[i][j]);
                                        };

                                        if(smatr_[i][j]<0)
                                        {
                                                smatr_neg_min=alp_data::Tmin(smatr_neg_min,smatr_[i][j]);
                                        };

                                };
                        };

                        if(aver_score>=-threshold)
                        {
                                throw error("Error - sum[i,j] RR1[i]*RR2[j]*smatr[i][j]>=0; the program cannot continue the calculation\n",3);
                        };

                        if(smatr_max<=0)
                        {
                                throw error("Error - at least one element of the scoring matrix must be positive\n",3);
                        };

                        

                        double a=eps;

                        while(importance_sampling::lambda_equation(a,(void*)(&tmp_ptr))>0)
                        {
                                a/=2.0;

                                if(a<threshold*100.0)
                                {
                                        throw error("Error - the input parameters correspond to non-logarithmic regime\n",3);
                                };
                        };

                        if(a<threshold*100.0)
                        {
                                throw error("Error - the input parameters define the regime which is too close to the critical regime\n",3);
                        };

                        eps=a/10.0;


                        double tmp_pr=RR1_[smatr_max_i]*RR2_[smatr_max_j];
                        double b=(log(1+10*eps)-log(tmp_pr))/(double)smatr_max;

                        
                        Int4 n_partition=2;
                        std::vector<double> res_lambda;
                        
                        
                        alp_reg::find_tetta_general(
                        importance_sampling::lambda_equation,
                        (void*)(&tmp_ptr),
                        a,
                        b,
                        n_partition,
                        eps,
                        res_lambda);

                        sort(res_lambda.begin(),res_lambda.end());

                        if(res_lambda.size()==0)
                        {
                                throw error("Error - the program is not able to find the ungapped lambda\n",3);
                        };

                        d_lambda=res_lambda[res_lambda.size()-1];
                        d_ungap_lambda=d_lambda;

                        //cout<<"\nUngapped lambda is "<<d_ungap_lambda<<endl;

                        d_lambda*=1.07;
                };


                
                d_is_number_of_AA=number_of_AA_;

                d_elements=new q_elem[number_of_AA_*number_of_AA_];
                alp_data::assert_mem(d_elements);

                d_elements_values=new double[number_of_AA_*number_of_AA_];
                alp_data::assert_mem(d_elements_values);



                d_alp_data->get_memory_for_matrix(d_is_number_of_AA,d_exp_s);

                Int4 ind=0;
                double sum=0;
                Int4 a,b;
                for(a=0;a<number_of_AA_;a++)
                {
                        for(b=0;b<number_of_AA_;b++)
                        {
                                d_exp_s[a][b]=exp(d_lambda*smatr_[a][b]);
                                d_elements_values[ind]=RR1_[a]*RR2_[b]*d_exp_s[a][b];
                                sum+=d_elements_values[ind];
                                ind++;
                        };
                };


                for(a=0;a<number_of_AA_;a++)
                {
                        for(b=0;b<number_of_AA_;b++)
                        {
                                d_exp_s[a][b]/=sum;
                        };
                };


                for(ind=0;ind<number_of_AA_*number_of_AA_;ind++)
                {
                        d_elements_values[ind]/=sum;
                };

                
                for(ind=1;ind<number_of_AA_*number_of_AA_;ind++)
                {
                        d_elements_values[ind]=d_elements_values[ind-1]+d_elements_values[ind];
                };

                
                ind=0;
                for(a=0;a<number_of_AA_;a++)
                {
                        for(b=0;b<number_of_AA_;b++)
                        {
                                q_elem elem_tmp;

                                elem_tmp.d_a=a;
                                elem_tmp.d_b=b;

                                d_elements[ind]=elem_tmp;
                                d_elements_values[ind]=d_elements_values[ind];

                                ind++;

                        };
                };



                d_mu=exp(-fabs(d_lambda)*open_);
                d_nu=exp(-fabs(d_lambda)*epen_);

                double tmp=1+d_mu-d_nu;

                d_eta=(1-d_nu)*(1-d_nu)/(tmp*tmp);
                d_mu_SI=1-d_nu;
                d_mu_IS=d_mu*(1-d_nu)/(tmp*tmp);
                d_mu_DS=d_mu/tmp;
                d_mu_SD=(1-d_nu)*(1-d_nu)/tmp;
                d_mu_ID=d_mu*(1-d_nu)/tmp;


                d_for_D[0]=d_nu;                                d_for_D_states[0]='D';
                d_for_D[1]=d_for_D[0]+d_mu_SD;        d_for_D_states[1]='S';
                d_for_D[2]=d_for_D[1]+d_mu_ID;        d_for_D_states[2]='I';

                d_for_I[0]=d_nu;                                d_for_I_states[0]='I';
                d_for_I[1]=d_for_I[0]+d_mu_SI;        d_for_I_states[1]='S';

                d_for_S[0]=d_eta;                                d_for_S_states[0]='S';
                d_for_S[1]=d_for_S[0]+d_mu_DS;        d_for_S_states[1]='D';
                d_for_S[2]=d_for_S[1]+d_mu_IS;        d_for_S_states[2]='I';

                d_alp_data->d_memory_size_in_MB+=sizeof(double)*number_of_AA_/mb_bytes;
                d_alp_data->d_memory_size_in_MB+=sizeof(q_elem)*number_of_AA_/mb_bytes;
        }
        catch (error er)
        {
                ee_error_flag=true;
                ee_error=er;                
        };
        }
        catch (...)
        { 
                ee_error_flag=true;
                ee_error=error("Internal error in the program\n",4);
        };

        //memory release

        if(ee_error_flag)
        {
                this->~importance_sampling();
                throw error(ee_error.st,ee_error.error_code);
        };

};

importance_sampling::~importance_sampling()
{
        delete []d_elements;d_elements=NULL;
        delete []d_elements_values;d_elements_values=NULL;

        if(d_alp_data)
        {
                d_alp_data->delete_memory_for_matrix(d_is_number_of_AA,d_exp_s);
                d_alp_data->d_memory_size_in_MB-=sizeof(double)*d_is_number_of_AA/mb_bytes;
                d_alp_data->d_memory_size_in_MB-=sizeof(q_elem)*d_is_number_of_AA/mb_bytes;
        };

};

void alp_data::get_current_time(
double &seconds_)
{
#ifndef NCBI_OS_MSWIN //UNIX program
         struct timeval tv;
     struct timezone tz;
     time_t t;

     gettimeofday(&tv, &tz);
     t = tv.tv_sec;
     localtime(&t);

     seconds_=(double)(t)+(double)(tv.tv_usec) * 0.000001;

#else

   struct _timeb timebuffer;

   _ftime( &timebuffer );

        seconds_=timebuffer.time+(double)(timebuffer.millitm)/1000.0;

#endif
};


