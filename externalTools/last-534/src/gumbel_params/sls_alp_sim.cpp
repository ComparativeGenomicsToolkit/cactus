/* $Id: sls_alp_sim.cpp 189337 2010-04-21 13:14:53Z boratyng $
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

File name: sls_alp_sim.cpp

Author: Sergey Sheetlin

Contents: Calculation of Gumbel parameters

******************************************************************************/


#include <ncbi_pch.hpp>

#include "sls_alp_sim.hpp"

USING_NCBI_SCOPE;
USING_SCOPE(blast);
USING_SCOPE(Sls);


alp_sim::alp_sim(//constructor
alp_data *alp_data_)
{

        
        d_alp_obj=NULL;


        d_lambda_tmp=NULL;
        d_lambda_tmp_errors=NULL;

        d_C_tmp=NULL;
        d_C_tmp_errors=NULL;


        d_alp_data=alp_data_;
        if(!d_alp_data)
        {
                throw error("Unexpected error",4);
        };
        bool ee_error_flag=false;
        error ee_error("",0);
        ofstream frand;

        try
        {
        try
        {


                d_alp_obj=new array_positive<alp*> (d_alp_data);
                alp_data::assert_mem(d_alp_obj);
                d_n_alp_obj=0;

                d_alp_data->d_memory_size_in_MB+=(double)(sizeof(array_positive<alp*>))/mb_bytes;


                double memory_before1=d_alp_data->d_memory_size_in_MB;
                double time_before1;
                alp_data::get_current_time(time_before1);


                d_alp_data->d_time_before1=time_before1;

                if(!d_alp_data->d_rand_flag)
                {
                        d_alp_data->d_rand_all->d_random_factor=d_alp_data->d_random_factor;
                };


                //trial simulation for estimation time and memory 
                Int4 M_min;
                Int4 nalp;

                double time_after_tmp;

                d_lambda_tmp=new array_positive<double>(d_alp_data);
                alp_data::assert_mem(d_lambda_tmp);

                d_lambda_tmp_errors=new array_positive<double>(d_alp_data);
                alp_data::assert_mem(d_lambda_tmp_errors);

                d_C_tmp=new array_positive<double>(d_alp_data);
                alp_data::assert_mem(d_C_tmp);

                d_C_tmp_errors=new array_positive<double>(d_alp_data);
                alp_data::assert_mem(d_C_tmp_errors);


                
                Int4 maximum_number_of_realizations_for_preliminary_simulation=1000;

                bool loop_break_flag;
                Int4 rand_i=-1;

                bool lambda_accuracy_flag=true;
                Int4 sim_number=1;
                do{

                        Int4 nalp_lambda;
                        bool C_calculation=false;
                        
                        Int4 number_tmp;

                        bool check_time_flag=true;

                        if(d_alp_data->d_rand_flag)
                        {
                                check_time_flag=false;
                                rand_i++;
                                number_tmp=d_alp_data->d_rand_all->d_first_stage_preliminary_realizations_numbers_ALP[rand_i];
                        }
                        else
                        {
                                number_tmp=alp_data::Tmin(maximum_number_of_realizations_for_preliminary_simulation-1,d_n_alp_obj+sim_number*d_alp_data->d_minimum_realizations_number-1);
                        };


                get_minimal_simulation(
                        0,
                        number_tmp,
                        M_min,
                        nalp,
                        nalp_lambda,
                        C_calculation,
                        check_time_flag);

                if(!d_alp_data->d_rand_flag)
                {
                        d_alp_data->d_rand_all->d_first_stage_preliminary_realizations_numbers_ALP.push_back(number_tmp);
                };

                

                sim_number*=2;

                if(d_lambda_tmp->d_elem[nalp]>=0)
                {
                        if(d_lambda_tmp_errors->d_elem[nalp]/d_lambda_tmp->d_elem[nalp]<alp_data_->d_eps_lambda)
                        {
                                lambda_accuracy_flag=false;
                        };
                };

                alp_data::get_current_time(time_after_tmp);

                
                if(number_tmp>=maximum_number_of_realizations_for_preliminary_simulation-1)
                {
                        if(!d_alp_data->d_rand_flag)
                        {
                                break;
                        };
                };
                
                if(d_alp_data->d_rand_flag)
                {
                        loop_break_flag=(rand_i>=(Int4)(d_alp_data->d_rand_all->d_first_stage_preliminary_realizations_numbers_ALP.size()-1));
                }
                else
                {
                loop_break_flag=!(maximum_number_of_realizations_for_preliminary_simulation>d_n_alp_obj-1&&
                        lambda_accuracy_flag&&((time_after_tmp-time_before1)<=0||((time_after_tmp-time_before1)<0.01*d_alp_data->d_max_time&&d_alp_data->d_memory_size_in_MB<d_alp_data->d_max_mem)));
                };

                }
                while(!loop_break_flag);


                

                double memory_after1=d_alp_data->d_memory_size_in_MB;
                double time_after1;
                alp_data::get_current_time(time_after1);

                if(memory_after1<=memory_before1)
                {
                        delete d_alp_obj;d_alp_obj=NULL;
                        throw error("Unexpected error\n",4);
                };

                Int4 limit_by_memory=(Int4)alp_data::Tmin((double)kMax_I4-1,alp_data::round(d_alp_data->d_max_mem/2.0/(memory_after1-memory_before1)*d_n_alp_obj));

                if(d_alp_data->d_memory_size_in_MB>d_alp_data->d_max_mem)
                {
                        //cout<<"\nWarning: memory limit is reached.\nTo get the requested accuracy please\nincrease maximum allowed amount of memory if possible\n\n";
                };


                Int4 limit_by_time;
                if(time_after1<=time_before1)
                {
                        limit_by_time=kMax_I4;
                }
                else
                {
                        limit_by_time=(Int4)alp_data::Tmin((double)kMax_I4-1,alp_data::round(d_alp_data->d_max_time/3.0/4.0/(time_after1-time_before1)*d_n_alp_obj));
                };


                Int4 realizations_number2=alp_data::Tmin((Int4)alp_data::round(limit_by_time)-1,limit_by_memory-1);
                
                realizations_number2=alp_data::Tmin(maximum_number_of_realizations_for_preliminary_simulation-(Int4)1,realizations_number2);

                realizations_number2=alp_data::Tmax(d_n_alp_obj-(Int4)1,realizations_number2);
                


                delete d_lambda_tmp;d_lambda_tmp=NULL;
                delete d_lambda_tmp_errors;d_lambda_tmp_errors=NULL;

                delete d_C_tmp;d_C_tmp=NULL;
                delete d_C_tmp_errors;d_C_tmp_errors=NULL;


                //trial simulation for estimation M_min and nalp
                d_lambda_tmp=new array_positive<double>(d_alp_data);
                alp_data::assert_mem(d_lambda_tmp);

                d_lambda_tmp_errors=new array_positive<double>(d_alp_data);
                alp_data::assert_mem(d_lambda_tmp_errors);

                d_C_tmp=new array_positive<double>(d_alp_data);
                alp_data::assert_mem(d_C_tmp);

                d_C_tmp_errors=new array_positive<double>(d_alp_data);
                alp_data::assert_mem(d_C_tmp_errors);



                Int4 nalp_lambda;
                bool C_calculation=false;

                


                double lambda;
                double eps_K=d_alp_data->d_eps_K;
                double K_C;
                double K_C_error;
                Int4 level;
                Int4 diff_opt;


                rand_i=-1;
                loop_break_flag=false;

                double time_before_ALP;
                double time_during_ALP;
                Int4 number_of_realizations_with_ALP=alp_data::Tmin(realizations_number2,d_n_alp_obj-(Int4)1+d_alp_data->d_minimum_realizations_number);
                Int4 number_of_realizations_with_ALP_pred;

                alp_data::get_current_time(time_before_ALP);
                do{

                        bool check_time_flag=true;

                        if(d_alp_data->d_rand_flag)
                        {
                                check_time_flag=false;
                                rand_i++;
                                number_of_realizations_with_ALP=d_alp_data->d_rand_all->d_preliminary_realizations_numbers_ALP[rand_i];
                        };

                get_minimal_simulation(
                        0,
                        number_of_realizations_with_ALP,
                        M_min,
                        nalp,
                        nalp_lambda,
                        C_calculation,
                        check_time_flag);

                if(!d_alp_data->d_rand_flag)
                {
                        d_alp_data->d_rand_all->d_preliminary_realizations_numbers_ALP.push_back(number_of_realizations_with_ALP);
                };




                lambda=d_lambda_tmp->d_elem[nalp];

                double tmp_lambda=2.0;

                if(d_lambda_tmp->d_elem[nalp]>0)
                {
                        tmp_lambda=((d_lambda_tmp_errors->d_elem[nalp]/d_lambda_tmp->d_elem[nalp])/d_alp_data->d_eps_lambda);
                };
                

                        
                        

                number_of_realizations_with_ALP_pred=number_of_realizations_with_ALP;



                alp_data::get_current_time(time_during_ALP);

                if(d_alp_data->d_rand_flag)
                {
                        loop_break_flag=(rand_i>=(Int4)(d_alp_data->d_rand_all->d_preliminary_realizations_numbers_ALP.size()-1));
                        if(loop_break_flag)
                        {
                                break;
                        };
                };


                if(time_during_ALP-time_before1>=d_alp_data->d_max_time*0.25||number_of_realizations_with_ALP>=realizations_number2||tmp_lambda<=1.0)
                {
                        if(!d_alp_data->d_rand_flag)
                        {
                                break;
                        };
                };

                if(time_during_ALP<=time_before_ALP)
                {
                        number_of_realizations_with_ALP=alp_data::Tmin(realizations_number2,number_of_realizations_with_ALP+d_alp_data->d_minimum_realizations_number);
                }
                else
                {

                        Int4 max_number_of_realizations=(Int4)floor(number_of_realizations_with_ALP*(d_alp_data->d_max_time*0.35-(time_before_ALP-time_before1))/(time_during_ALP-time_before_ALP));
                        number_of_realizations_with_ALP=alp_data::Tmin(realizations_number2,(Int4)floor(0.5*number_of_realizations_with_ALP+0.5*max_number_of_realizations));
                        if(number_of_realizations_with_ALP>=max_number_of_realizations)
                        {
                                number_of_realizations_with_ALP=d_alp_data->Tmin(realizations_number2,number_of_realizations_with_ALP+d_alp_data->d_minimum_realizations_number);
                        };

                        if((double)(number_of_realizations_with_ALP-number_of_realizations_with_ALP_pred)/(double)number_of_realizations_with_ALP_pred<0.005)
                        {
                                number_of_realizations_with_ALP=number_of_realizations_with_ALP_pred;
                                if(!d_alp_data->d_rand_flag)
                                {
                                        break;
                                };

                        };
                };



                }
                while(!loop_break_flag);

                realizations_number2=number_of_realizations_with_ALP;
                Int4 realizations_number2_lambda=number_of_realizations_with_ALP;

                rand_i=-1;
                loop_break_flag=false;

                double time_before_kill;
                double time_during_kill;
                Int4 number_of_realizations_with_killing=d_alp_data->Tmin(realizations_number2,d_alp_data->d_minimum_realizations_number-(Int4)1);
                Int4 number_of_realizations_with_killing_pred;

                alp_data::get_current_time(time_before_kill);
                do{

                        bool check_time_flag=false;

                        if(d_alp_data->d_rand_flag)
                        {
                                check_time_flag=false;
                                rand_i++;
                                number_of_realizations_with_killing=d_alp_data->d_rand_all->d_preliminary_realizations_numbers_killing[rand_i];
                        };

                        
                        kill(
                                check_time_flag,
                                0,
                                number_of_realizations_with_killing,
                                M_min,
                                lambda,
                                eps_K,
                                K_C,
                                K_C_error,
                                level,
                                diff_opt);

                        if(!d_alp_data->d_rand_flag)
                        {
                                d_alp_data->d_rand_all->d_preliminary_realizations_numbers_killing.push_back(number_of_realizations_with_killing);
                        };

        

                        number_of_realizations_with_killing_pred=number_of_realizations_with_killing;


                alp_data::get_current_time(time_during_kill);

                double tmp_K=2.0;

                if(K_C>0)
                {
                        tmp_K=((K_C_error/K_C)/d_alp_data->d_eps_K);
                };

                if(d_alp_data->d_rand_flag)
                {
                        loop_break_flag=(rand_i>=(Int4)(d_alp_data->d_rand_all->d_preliminary_realizations_numbers_killing.size()-1));
                        if(loop_break_flag)
                        {
                                break;
                        };
                };


                if(time_during_kill-time_before1>=d_alp_data->d_max_time||number_of_realizations_with_killing>=realizations_number2||tmp_K<=1.0)
                {
                        if(!d_alp_data->d_rand_flag)
                        {
                                break;
                        };
                };



                if(time_during_kill<=time_before_kill)
                {
                        number_of_realizations_with_killing=d_alp_data->Tmin(realizations_number2,number_of_realizations_with_killing+d_alp_data->d_minimum_realizations_number);
                }
                else
                {

                        Int4 max_number_of_realizations=(Int4)floor(number_of_realizations_with_killing*(d_alp_data->d_max_time-(time_before_kill-time_before1))/(time_during_kill-time_before_kill));
                        number_of_realizations_with_killing=d_alp_data->Tmin(realizations_number2,(Int4)floor(0.5*number_of_realizations_with_killing+0.5*max_number_of_realizations));
                        if(number_of_realizations_with_killing>=max_number_of_realizations)
                        {
                                number_of_realizations_with_killing=d_alp_data->Tmin(realizations_number2,number_of_realizations_with_killing+d_alp_data->d_minimum_realizations_number);
                        };

                        if((double)(number_of_realizations_with_killing-number_of_realizations_with_killing_pred)/(double)number_of_realizations_with_killing_pred<0.005)
                        {
                                number_of_realizations_with_killing=number_of_realizations_with_killing_pred;
                                if(!d_alp_data->d_rand_flag)
                                {
                                        break;
                                };
                        };
                };


                }
                while(!loop_break_flag);

                double memory_after2=d_alp_data->d_memory_size_in_MB;


                Int4 k;
                for(k=0;k<=number_of_realizations_with_killing;k++)
                {
                        d_alp_obj->d_elem[k]->partially_release_memory();                
                };

                        
                realizations_number2=number_of_realizations_with_killing;
                Int4 realizations_number2_K=number_of_realizations_with_killing;



                d_alp_data->d_realizations_number2=realizations_number2_lambda;


                if(K_C<=0)
                {
                        throw error("The program cannot estimate the parameters; please repeat the calculation",2);
                };

                double tmp=((K_C_error/K_C)/d_alp_data->d_eps_K);
                tmp=alp_data::Tmin(ceil((realizations_number2_K+1)*tmp*tmp),(double)kMax_I4);
                Int4 realizations_number_killing=(Int4)tmp;
                tmp=((d_lambda_tmp_errors->d_elem[nalp]/d_lambda_tmp->d_elem[nalp])/d_alp_data->d_eps_lambda);
                tmp=alp_data::Tmin(ceil((realizations_number2_lambda+1)*tmp*tmp),(double)kMax_I4);
                Int4 realizations_number_lambda=(Int4)tmp;

                double time_after2;
                alp_data::get_current_time(time_after2);
                time_after2=alp_data::Tmax(time_after2,time_after_tmp);

                delete d_lambda_tmp;d_lambda_tmp=NULL;
                delete d_lambda_tmp_errors;d_lambda_tmp_errors=NULL;

                delete d_C_tmp;d_C_tmp=NULL;
                delete d_C_tmp_errors;d_C_tmp_errors=NULL;

                


                //main simulation

                Int4 j;
                j=1;
                Int4 kill_j=0;




                bool kill_flag=(realizations_number_killing>realizations_number2_K+1+j);
                bool lambda_flag=(realizations_number_lambda>realizations_number2_lambda+1+j);
                Int4 nalp_for_simulation=nalp;

                if(d_alp_data->d_rand_flag)
                {
                        lambda_flag=(d_alp_data->d_rand_all->d_total_realizations_number_with_ALP>realizations_number2_K+j-1);
                        kill_flag=(d_alp_data->d_rand_all->d_total_realizations_number_with_killing>realizations_number2_K+j-1);
                };




                if(kill_flag||lambda_flag)
                {

                        

                        Int4 step_for_time=1;

                        Int4 number_of_unsuccesful_objects=0;

                        
                        
                        while(d_n_alp_obj<kMax_I4-realizations_number2_lambda-j)
                        {
                                kill_flag=(realizations_number_killing>realizations_number2_K+j);
                                lambda_flag=(realizations_number_lambda>realizations_number2_lambda+j);

                                if(d_alp_data->d_rand_flag)
                                {
                                        lambda_flag=(d_alp_data->d_rand_all->d_total_realizations_number_with_ALP>realizations_number2_K+j-1);
                                        kill_flag=(d_alp_data->d_rand_all->d_total_realizations_number_with_killing>realizations_number2_K+j-1);
                                };

                                if(!(kill_flag||lambda_flag))
                                {
                                        if(!d_alp_data->d_rand_flag)
                                        {
                                                //cout<<"\nThe parameters were calculated with the required accuracy\n\n";
                                                break;
                                        }
                                        else
                                        {
                                                break;
                                        };
                                };



                                

                                if(!kill_flag)
                                {
                                        nalp_for_simulation=alp_data::Tmin(nalp_lambda,nalp);
                                };

                                bool sucess_flag=false;
                                


                                if(realizations_number2_K+j<=realizations_number2_lambda)
                                {
                        
                                }
                                else
                                {
                                        d_alp_obj->set_elem(realizations_number2_K+j,NULL);
                                        d_n_alp_obj++;
                                };

                                alp *&obj=d_alp_obj->d_elem[realizations_number2_K+j];



                                try
                                {
                                        double eps_tmp;

                                while(!sucess_flag)
                                {

                                        
                                        bool check_time_flag=true;

                                        if(d_alp_data->d_rand_flag)
                                        {
                                                check_time_flag=false;
                                        };

                                                get_single_realization(
                                                check_time_flag,
                                                M_min,
                                                nalp_for_simulation,
                                                kill_flag,
                                                level,
                                                diff_opt,
                                                obj,
                                                sucess_flag,
                                                eps_tmp);

                                                
                                                


                                        if(!sucess_flag)
                                        {
                                                number_of_unsuccesful_objects++;
                                                
                                                if(number_of_unsuccesful_objects>5+j*eps_tmp)
                                                {
                                                        if(realizations_number2_K+j>realizations_number2_lambda)
                                                        {
                                                                d_n_alp_obj--;
                                                        };
                                                        throw error("The program cannot calculate the parameters\n",3);
                                                };
                                        };
                                };

                                }
                                catch (error_for_single_realization er)
                                {
                                        if(d_alp_data->d_rand_flag)
                                        {
                                                throw error("Unexpected error in ramdomization\n",4);
                                        };

                                        //cout<<"\nWarning: time limit is reached.\nTo get the requested accuracy please\nincrease the allowed calculation time if possible\n\n";
                                        if(realizations_number2_K+j>realizations_number2_lambda)
                                        {
                                                delete obj;obj=NULL;
                                                d_n_alp_obj--;
                                        };
                                        break;
                                };

                                

                                

                                if(realizations_number2_K+j>realizations_number2_lambda)
                                {

                                        if(kill_flag)
                                        {
                                                kill_j=j;
                                        };

                                };


                                d_alp_obj->d_elem[realizations_number2_K+j]->partially_release_memory();

                                j++;

                                

                                if(d_alp_data->d_memory_size_in_MB>d_alp_data->d_max_mem)
                                {
                                        if(!d_alp_data->d_rand_flag)
                                        {
                                                //cout<<"\nWarning: memory limit is reached.\nTo get the requested accuracy please\nincrease the allowed amount of memory if possible\n\n";
                                                break;
                                        };
                                };
                                

                                if(j%step_for_time==0)
                                {
                                        double time_after3;
                                        alp_data::get_current_time(time_after3);

                                        if((time_after3-time_before1)>d_alp_data->d_max_time)
                                        {
                                                if(!d_alp_data->d_rand_flag)
                                                {
                                                        //cout<<"\nWarning: time limit is reached.\nTo get the requested accuracy please\nincrease the allowed calculation time if possible\n\n";
                                                        break;
                                                };
                                        };
                                };
                        };
                }
                else
                {
                        //cout<<"\nThe parameters were calculated with required accuracy\n\n";
                };
                

                Int4 final_realizations_number_killing=kill_j+realizations_number2_K+1;
                Int4 final_realizations_number_lambda=alp_data::Tmax(realizations_number2_lambda+1,j+realizations_number2_K);
                d_n_alp_obj=final_realizations_number_lambda;


                double memory_after3=d_alp_data->d_memory_size_in_MB;

                double time_after100;
                alp_data::get_current_time(time_after100);
                //cout<<"\nActual calculation time is "<<time_after100-time_before1<<" seconds\n";
                //cout<<"Actual memory usage is "<<d_alp_data->Tmax(memory_after2,memory_after3)<<" MBs\n\n";

                if(!d_alp_data->d_rand_flag)
                {
                        d_alp_data->d_rand_all->d_total_realizations_number_with_ALP=final_realizations_number_lambda-1;
                        d_alp_data->d_rand_all->d_total_realizations_number_with_killing=final_realizations_number_killing-1;
                };

        

                bool inside_simulation_flag;

                this->m_CalcTime=time_after100-time_before1;

                output_main_parameters2m_new(
                time_after100-time_before1,
                nalp,
                nalp_for_simulation,
                level,
                M_min,
                inside_simulation_flag,
                final_realizations_number_lambda,
                final_realizations_number_killing);


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
        if(frand.is_open())
        {
                frand.close();
        };

        if(ee_error_flag)
        {
                delete d_lambda_tmp;d_lambda_tmp=NULL;
                delete d_lambda_tmp_errors;d_lambda_tmp_errors=NULL;

                delete d_C_tmp;d_C_tmp=NULL;
                delete d_C_tmp_errors;d_C_tmp_errors=NULL;

                this->~alp_sim();
                throw error(ee_error.st,ee_error.error_code);
        };

};


void alp_sim::randomize_realizations(
Int4 final_realizations_number_lambda_,
Int4 final_realizations_number_killing_)
{
        randomize_realizations_ind(0,final_realizations_number_killing_-1);
        randomize_realizations_ind(final_realizations_number_killing_,final_realizations_number_lambda_-1);

};

void alp_sim::randomize_realizations_ind(
Int4 ind1_,
Int4 ind2_)
{
        bool ee_error_flag=false;
        error ee_error("",0);
        alp**array_ind=NULL;
        Int4 *perm=NULL;

        try
        {
        try
        {

                if(ind1_>=ind2_)
                {
                        return;
                };

                if(ind2_>d_n_alp_obj-1)
                {
                        throw error("Unexpected error",4);
                };


                Int4 total_number=ind2_-ind1_+1;

                array_ind=new alp*[total_number];
                alp_data::assert_mem(array_ind);




                
                perm=new Int4 [total_number];
                alp_data::assert_mem(perm);

                generate_random_permulation(perm,total_number);

                Int4 i;
                for(i=0;i<total_number;i++)
                {
                        array_ind[i]=d_alp_obj->d_elem[ind1_+perm[i]];
                };

                for(i=0;i<total_number;i++)
                {
                        d_alp_obj->d_elem[ind1_+i]=array_ind[i];
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

        //memory release
        delete[]array_ind;array_ind=NULL;
        delete[]perm;perm=NULL;


        if(ee_error_flag)
        {
                throw error(ee_error.st,ee_error.error_code);
        };


};

void alp_sim::generate_random_permulation(
Int4 *perm_,
Int4 dim_)
{
        Int4 i;
        for(i=0;i<dim_;i++)
        {
                perm_[i]=i;
        };

        for(i=0;i<dim_-1;i++)
        {
                Int4 ind_swap=i+alp_data::random_long(d_alp_data->ran2(),dim_-i);
                Int4 tmp=perm_[ind_swap];
                perm_[ind_swap]=perm_[i];
                perm_[i]=tmp;
        };
};


void alp_sim::output_main_parameters2m_new(
double time_,
Int4 nalp,
Int4 nalp_for_lambda_simulation,
Int4 level,
Int4 M_min_,
bool &inside_simulation_flag,
Int4 final_realizations_number_lambda_,
Int4 final_realizations_number_killing_)
{
        bool ee_error_flag=false;
        error ee_error("",0);

        try
        {
        try
        {

                double lambda;
                double lambda_error;
                double test_difference;
                double test_difference_error;
                double C;
                double C_error;
                double C2;
                double C2_error;
                double C4;
                double C4_error;
                double K_C;
                double K_C_error;
                double a_I;
                double a_I_error;
                double a_J;
                double a_J_error;
                double sigma;
                double sigma_error;
                double alpha_I;
                double alpha_I_error;
                double alpha_J;
                double alpha_J_error;
                double K;
                double K_error;

                bool flag=false;
                Int4 number_of_trials=0;
                Int4 number_of_trials_threshold=4;

                do{

                calculate_main_parameters2m(
                final_realizations_number_lambda_,
                final_realizations_number_killing_,
                nalp,
                nalp_for_lambda_simulation,
                level,
                inside_simulation_flag,
                lambda,
                lambda_error,
                test_difference,
                test_difference_error,
                C,
                C_error,
                C2,
                C2_error,
                C4,
                C4_error,
                K_C,
                K_C_error,
                a_I,
                a_I_error,
                a_J,
                a_J_error,
                sigma,
                sigma_error,
                alpha_I,
                alpha_I_error,
                alpha_J,
                alpha_J_error,
                K,
                K_error,
                flag);


                number_of_trials++;

                if(!flag)
                {
                        randomize_realizations(
                        final_realizations_number_lambda_,
                        final_realizations_number_killing_);
                        //cout<<"Randomization attempt\t"<<number_of_trials<<endl;
                };
                }
                while(!flag&&number_of_trials<=number_of_trials_threshold);

                if(!flag)
                {
                        throw error("Error - please run the program once again\n",2);
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

        //memory release

        if(ee_error_flag)
        {
                throw error(ee_error.st,ee_error.error_code);
        };
};



double alp_sim::round_doulbe(
double val_,
Int4 digits_)
{
        Int4 i;
        for(i=0;i<digits_;i++)
        {
                val_*=10;
        };
        val_=alp_data::round(val_);
        for(i=0;i<digits_;i++)
        {
                val_/=10.0;
        };

        return val_;
};

double alp_sim::relative_error_in_percents(
double val_,
double val_error_)
{
        if(val_==0)
        {
                return DBL_MAX;
        };

        return fabs(round_doulbe(val_error_/val_*100.0,1));

};


Int4 alp_sim::get_number_of_subsimulations(
Int4 number_of_realizations_)
{
        Int4 max_number_of_subsimulations=20;
        Int4 min_number_of_subsimulations=3;

        Int4 min_number_of_realizations_for_subsimulation=2;

        if(number_of_realizations_<min_number_of_realizations_for_subsimulation*min_number_of_subsimulations)
        {
                throw error("Please repeat the simulation or increase calculation time\n",1);
        };

        Int4 res_subsimulations=(Int4)ceil(sqrt((double)number_of_realizations_));
        res_subsimulations=alp_data::Tmin(res_subsimulations,max_number_of_subsimulations);
        res_subsimulations=alp_data::Tmax(res_subsimulations,min_number_of_subsimulations);

        return res_subsimulations;
};



void alp_sim::calculate_main_parameters2m(
Int4 final_realizations_number_lambda_,
Int4 final_realizations_number_killing_,
Int4 nalp,
Int4 nalp_for_lambda_simulation,
Int4 level,
bool &inside_simulation_flag,
double &lambda,
double &lambda_error,
double &test_difference,
double &test_difference_error,
double &C,
double &C_error,
double &C2,
double &C2_error,
double &C4,
double &C4_error,
double &K_C,
double &K_C_error,
double &a_I,
double &a_I_error,
double &a_J,
double &a_J_error,
double &sigma,
double &sigma_error,
double &alpha_I,
double &alpha_I_error,
double &alpha_J,
double &alpha_J_error,
double &K,
double &K_error,
bool &flag_)
{

        bool ee_error_flag=false;
        error ee_error("",0);

        Int4 *d_mult_realizations=NULL;
        Int4 *d_mult_K_realizations=NULL;

        double *lambda_mult=NULL;
        double *lambda_mult_error=NULL;

        double *C_mult=NULL;
        double *C_mult_error=NULL;

        double *a_I_mult=NULL;
        double *a_I_mult_error=NULL;

        double *a_J_mult=NULL;
        double *a_J_mult_error=NULL;

        double *sigma_mult=NULL;
        double *sigma_mult_error=NULL;

        double *alpha_I_mult=NULL;
        double *alpha_I_mult_error=NULL;

        double *alpha_J_mult=NULL;
        double *alpha_J_mult_error=NULL;

        double *K_C_mult=NULL;
        double *K_C_mult_error=NULL;

        double *K_mult=NULL;
        double *K_mult_error=NULL;

        void **alp_distr=NULL;
        void **alp_distr_errors=NULL;

        void ***alp_mult_distr=NULL;
        void ***alp_mult_distr_errors=NULL;


        try
        {
        try
        {

                flag_=false;


                if(final_realizations_number_killing_>final_realizations_number_lambda_)
                {
                        throw error("Unexpected error\n",4);
                };

                Int4 mult_number_lambda=get_number_of_subsimulations(d_n_alp_obj);
                Int4 mult_number_K=get_number_of_subsimulations(final_realizations_number_killing_);

                d_mult_number=alp_data::Tmin(mult_number_lambda,mult_number_K);

                double mult_number_double_lambda=mult_number_lambda;
                double mult_number_double_K=mult_number_lambda;



                Int4 j;


                d_mult_realizations=new Int4[d_mult_number+1];
                alp_data::assert_mem(d_mult_realizations);



                d_mult_K_realizations=new Int4[d_mult_number+1];
                alp_data::assert_mem(d_mult_K_realizations);

                lambda_mult=new double[d_mult_number+1];
                alp_data::assert_mem(lambda_mult);
                lambda_mult_error=new double[d_mult_number+1];
                alp_data::assert_mem(lambda_mult_error);

                C_mult=new double[d_mult_number+1];
                alp_data::assert_mem(C_mult);
                C_mult_error=new double[d_mult_number+1];
                alp_data::assert_mem(C_mult_error);

                a_I_mult=new double[d_mult_number+1];
                alp_data::assert_mem(a_I_mult);
                a_I_mult_error=new double[d_mult_number+1];
                alp_data::assert_mem(a_I_mult_error);

                a_J_mult=new double[d_mult_number+1];
                alp_data::assert_mem(a_J_mult);
                a_J_mult_error=new double[d_mult_number+1];
                alp_data::assert_mem(a_J_mult_error);

                sigma_mult=new double[d_mult_number+1];
                alp_data::assert_mem(sigma_mult);
                sigma_mult_error=new double[d_mult_number+1];
                alp_data::assert_mem(sigma_mult_error);

                alpha_I_mult=new double[d_mult_number+1];
                alp_data::assert_mem(alpha_I_mult);
                alpha_I_mult_error=new double[d_mult_number+1];
                alp_data::assert_mem(alpha_I_mult_error);

                alpha_J_mult=new double[d_mult_number+1];
                alp_data::assert_mem(alpha_J_mult);
                alpha_J_mult_error=new double[d_mult_number+1];
                alp_data::assert_mem(alpha_J_mult_error);

                K_C_mult=new double[d_mult_number+1];
                alp_data::assert_mem(K_C_mult);
                K_C_mult_error=new double[d_mult_number+1];
                alp_data::assert_mem(K_C_mult_error);

                K_mult=new double[d_mult_number+1];
                alp_data::assert_mem(K_mult);
                K_mult_error=new double[d_mult_number+1];
                alp_data::assert_mem(K_mult_error);


                double lambda_mult2=0;
                double C_mult2=0;
                double K_C_mult2=0;
                double a_I_mult2=0;
                double a_J_mult2=0;
                double sigma_mult2=0;
                double alpha_I_mult2=0;
                double alpha_J_mult2=0;
                double K_mult2=0;

                double lambda_mult2_error=0;
                double C_mult2_error=0;
                double K_C_mult2_error=0;
                double a_I_mult2_error=0;
                double a_J_mult2_error=0;
                double sigma_mult2_error=0;
                double alpha_I_mult2_error=0;
                double alpha_J_mult2_error=0;
                double K_mult2_error=0;




                
                for(j=0;j<=nalp_for_lambda_simulation;j++)
                {
                        get_and_allocate_alp_distribution(
                        0,
                        d_n_alp_obj-1,
                        alp_distr,
                        alp_distr_errors,
                        j);
                };


                alp_mult_distr=new void **[d_mult_number+1];
                alp_data::assert_mem(alp_mult_distr);

                for(j=0;j<=d_mult_number;j++)
                {
                        alp_mult_distr[j]=NULL;
                };

                alp_mult_distr_errors=new void **[d_mult_number+1];
                alp_data::assert_mem(alp_mult_distr_errors);
                for(j=0;j<=d_mult_number;j++)
                {
                        alp_mult_distr_errors[j]=NULL;
                };

                alp_mult_distr[0]=alp_distr;
                alp_mult_distr_errors[0]=alp_distr_errors;

                Int4 real_number=(Int4)floor((double)final_realizations_number_lambda_/(double)d_mult_number);

                d_mult_realizations[0]=final_realizations_number_lambda_;

                Int4 k;
                for(k=1;k<=d_mult_number;k++)
                {
                        d_mult_realizations[k]=real_number;
                };


                Int4 nr_tmp=0;
                for(k=1;k<=d_mult_number;k++)
                {
                        nr_tmp+=d_mult_realizations[k];
                        Int4 j;
                        for(j=0;j<=nalp_for_lambda_simulation;j++)
                        {
                                get_and_allocate_alp_distribution(
                                nr_tmp-d_mult_realizations[k],
                                nr_tmp-1,
                                alp_mult_distr[k],
                                alp_mult_distr_errors[k],
                                j);
                        };
                };

                nr_tmp=0;
                for(k=1;k<=d_mult_number;k++)
                {
                        nr_tmp+=d_mult_realizations[k];
                        Int4 nalp_tmp;

                        double test_difference;
                        double test_difference_error;



                        calculate_lambda(
                        false,
                        nalp_for_lambda_simulation,
                        nalp_tmp,
                        inside_simulation_flag,
                        nr_tmp-d_mult_realizations[k],
                        nr_tmp-1,
                        alp_mult_distr[k],
                        alp_mult_distr_errors[k],
                        lambda_mult[k],
                        lambda_mult_error[k],
                        test_difference,
                        test_difference_error);


                        if(!inside_simulation_flag)
                        {
                                goto label1;
                        };


                        lambda_mult2+=lambda_mult[k];
                        lambda_mult2_error+=lambda_mult[k]*lambda_mult[k];
                };



                Int4 nalp_tmp;

                calculate_lambda(
                false,
                nalp_for_lambda_simulation,
                nalp_tmp,
                inside_simulation_flag,
                0,
                final_realizations_number_lambda_-1,
                alp_distr,
                alp_distr_errors,
                lambda,
                lambda_error,
                test_difference,
                test_difference_error);

                if(!inside_simulation_flag)
                {
                        throw error("Error - please run the program once again\n",2);
                };

                lambda_mult[0]=lambda;
                lambda_mult_error[0]=lambda_error;

                

                nr_tmp=0;
                for(k=1;k<=d_mult_number;k++)
                {

                        nr_tmp+=d_mult_realizations[k];
                        calculate_C(
                        0,
                        nalp_for_lambda_simulation,
                        nr_tmp-d_mult_realizations[k],
                        nr_tmp-1,
                        alp_mult_distr[k],
                        alp_mult_distr_errors[k],
                        lambda_mult[k],
                        lambda_mult_error[k],
                        C_mult[k],
                        C_mult_error[k]);

                        C_mult2+=C_mult[k];
                        C_mult2_error+=C_mult[k]*C_mult[k];


                };
                

                calculate_C(
                0,
                nalp_for_lambda_simulation,
                0,
                final_realizations_number_lambda_-1,
                alp_distr,
                alp_distr_errors,
                lambda,
                lambda_error,
                C,
                C_error);


                C_mult[0]=C;
                C_mult_error[0]=C_error;



                nr_tmp=0;
                for(k=1;k<=d_mult_number;k++)
                {

                        nr_tmp+=d_mult_realizations[k];
                        calculate_FSC(
                        nalp_for_lambda_simulation,
                        nr_tmp-d_mult_realizations[k],
                        nr_tmp-1,
                        alp_mult_distr[k],
                        alp_mult_distr_errors[k],
                        lambda_mult[k],
                        lambda_mult_error[k],
                        a_I_mult[k],
                        a_I_mult_error[k],
                        a_J_mult[k],
                        a_J_mult_error[k],
                        sigma_mult[k],
                        sigma_mult_error[k],
                        alpha_I_mult[k],
                        alpha_I_mult_error[k],
                        alpha_J_mult[k],
                        alpha_J_mult_error[k]);

                        a_I_mult2+=a_I_mult[k];
                        a_I_mult2_error+=a_I_mult[k]*a_I_mult[k];

                        a_J_mult2+=a_J_mult[k];
                        a_J_mult2_error+=a_J_mult[k]*a_J_mult[k];

                        sigma_mult2+=sigma_mult[k];
                        sigma_mult2_error+=sigma_mult[k]*sigma_mult[k];

                        alpha_I_mult2+=alpha_I_mult[k];
                        alpha_I_mult2_error+=alpha_I_mult[k]*alpha_I_mult[k];

                        alpha_J_mult2+=alpha_J_mult[k];
                        alpha_J_mult2_error+=alpha_J_mult[k]*alpha_J_mult[k];


                };


                calculate_FSC(
                nalp_for_lambda_simulation,
                0,
                final_realizations_number_lambda_-1,
                alp_distr,
                alp_distr_errors,
                lambda,
                lambda_error,
                a_I,
                a_I_error,
                a_J,
                a_J_error,
                sigma,
                sigma_error,
                alpha_I,
                alpha_I_error,
                alpha_J,
                alpha_J_error);
                

                a_I_mult[0]=a_I;
                a_I_mult_error[0]=a_I_error;
                a_J_mult[0]=a_J;
                a_J_mult_error[0]=a_J_error;
                sigma_mult[0]=sigma;
                sigma_mult_error[0]=sigma_error;
                alpha_I_mult[0]=alpha_I;
                alpha_I_mult_error[0]=alpha_I_error;
                alpha_J_mult[0]=alpha_J;
                alpha_J_mult_error[0]=alpha_J_error;
                

                real_number=(Int4)floor((double)final_realizations_number_killing_/(double)d_mult_number);


                d_mult_K_realizations[0]=final_realizations_number_killing_;


                for(k=1;k<=d_mult_number;k++)
                {
                        d_mult_K_realizations[k]=real_number;
                };


                nr_tmp=0;
                for(k=1;k<=d_mult_number;k++)
                {
                        nr_tmp+=d_mult_K_realizations[k];

                        Int4 recommended_level;
                        Int4 diff_opt;



                        check_K_criterion_during_killing(
                        nr_tmp-d_mult_K_realizations[k],
                        nr_tmp-1,
                        lambda_mult[k],
                        d_alp_data->d_eps_K,
                        level,
                        recommended_level,
                        diff_opt,
                        K_C_mult[k],
                        K_C_mult_error[k]);


                        K_mult[k]=C_mult[k]*K_C_mult[k];
                        K_mult_error[k]=alp_reg::error_of_the_product(
                        C_mult[k],
                        C_mult_error[k],
                        K_C_mult[k],
                        K_C_mult_error[k]);

                        K_C_mult2+=K_C_mult[k];
                        K_C_mult2_error+=K_C_mult[k]*K_C_mult[k];

                        K_mult2+=K_mult[k];
                        K_mult2_error+=K_mult[k]*K_mult[k];

                };


                Int4 recommended_level;
                Int4 diff_opt;



                check_K_criterion_during_killing(
                0,
                final_realizations_number_killing_-1,
                lambda,
                d_alp_data->d_eps_K,
                level,
                recommended_level,
                diff_opt,
                K_C,
                K_C_error);



                

                K=C*K_C;
                K_error=alp_reg::error_of_the_product(
                C,
                C_error,
                K_C,
                K_C_error);

                K_C_mult[0]=K_C;
                K_C_mult_error[0]=K_C_error;

                K_mult[0]=K;
                K_mult_error[0]=K_error;



                lambda_mult2/=d_mult_number;
                C_mult2/=d_mult_number;
                K_C_mult2/=d_mult_number;
                a_I_mult2/=d_mult_number;
                a_J_mult2/=d_mult_number;
                sigma_mult2/=d_mult_number;
                alpha_I_mult2/=d_mult_number;
                alpha_J_mult2/=d_mult_number;
                K_mult2/=d_mult_number;

                lambda_mult2_error/=d_mult_number;
                C_mult2_error/=d_mult_number;
                K_C_mult2_error/=d_mult_number;
                a_I_mult2_error/=d_mult_number;
                a_J_mult2_error/=d_mult_number;
                sigma_mult2_error/=d_mult_number;
                alpha_I_mult2_error/=d_mult_number;
                alpha_J_mult2_error/=d_mult_number;
                K_mult2_error/=d_mult_number;


                mult_number_double_lambda=(double)final_realizations_number_lambda_/(double)real_number;
                mult_number_double_K=(double)final_realizations_number_killing_/(double)real_number;


                lambda_mult2_error=alp_reg::sqrt_for_errors(lambda_mult2_error-lambda_mult2*lambda_mult2)/sqrt((double)mult_number_double_lambda);
                C_mult2_error=alp_reg::sqrt_for_errors(C_mult2_error-C_mult2*C_mult2)/sqrt((double)mult_number_double_lambda);
                K_C_mult2_error=alp_reg::sqrt_for_errors(K_C_mult2_error-K_C_mult2*K_C_mult2)/sqrt((double)mult_number_double_K);
                a_I_mult2_error=alp_reg::sqrt_for_errors(a_I_mult2_error-a_I_mult2*a_I_mult2)/sqrt((double)mult_number_double_lambda);
                a_J_mult2_error=alp_reg::sqrt_for_errors(a_J_mult2_error-a_J_mult2*a_J_mult2)/sqrt((double)mult_number_double_lambda);
                sigma_mult2_error=alp_reg::sqrt_for_errors(sigma_mult2_error-sigma_mult2*sigma_mult2)/sqrt((double)mult_number_double_lambda);
                alpha_I_mult2_error=alp_reg::sqrt_for_errors(alpha_I_mult2_error-alpha_I_mult2*alpha_I_mult2)/sqrt((double)mult_number_double_lambda);
                alpha_J_mult2_error=alp_reg::sqrt_for_errors(alpha_J_mult2_error-alpha_J_mult2*alpha_J_mult2)/sqrt((double)mult_number_double_lambda);
                K_mult2_error=alp_reg::sqrt_for_errors(K_mult2_error-K_mult2*K_mult2)/sqrt((double)alp_data::Tmin(mult_number_double_lambda,mult_number_double_K));


                error_in_calculate_main_parameters2m(
                lambda,
                lambda_error,
                lambda_mult2,
                lambda_mult2_error);

                error_in_calculate_main_parameters2m(
                C,
                C_error,
                C_mult2,
                C_mult2_error);

                error_in_calculate_main_parameters2m(
                K_C,
                K_C_error,
                K_C_mult2,
                K_C_mult2_error);

                error_in_calculate_main_parameters2m(
                a_I,
                a_I_error,
                a_I_mult2,
                a_I_mult2_error);

                error_in_calculate_main_parameters2m(
                a_J,
                a_J_error,
                a_J_mult2,
                a_J_mult2_error);


                error_in_calculate_main_parameters2m(
                sigma,
                sigma_error,
                sigma_mult2,
                sigma_mult2_error);

                error_in_calculate_main_parameters2m(
                alpha_I,
                alpha_I_error,
                alpha_I_mult2,
                alpha_I_mult2_error);

                error_in_calculate_main_parameters2m(
                alpha_J,
                alpha_J_error,
                alpha_J_mult2,
                alpha_J_mult2_error);

                error_in_calculate_main_parameters2m(
                K,
                K_error,
                K_mult2,
                K_mult2_error);

                flag_=true;


                this->m_AI=a_I;
                this->m_AIError=a_I_error;
                this->m_AJ=a_J;
                this->m_AJError=a_J_error;
                this->m_Sigma=sigma;
                this->m_SigmaError=sigma_error;
                this->m_C=C;
                this->m_CError=C_error;
                this->m_K=K;
                this->m_KError=K_error;
                this->m_Lambda=lambda;
                this->m_LambdaError=lambda_error;

                this->m_AlphaI=alpha_I;
                this->m_AlphaIError=alpha_I_error;
                this->m_AlphaJ=alpha_J;
                this->m_AlphaJError=alpha_J_error;

                this->m_AISbs.resize(d_mult_number);
                this->m_AJSbs.resize(d_mult_number);
                this->m_SigmaSbs.resize(d_mult_number);
                this->m_CSbs.resize(d_mult_number);
                this->m_KSbs.resize(d_mult_number);
                this->m_LambdaSbs.resize(d_mult_number);

                this->m_AlphaISbs.resize(d_mult_number);
                this->m_AlphaJSbs.resize(d_mult_number);


                for(k=1;k<=d_mult_number;k++)
                {
                        this->m_AISbs[k-1]=a_I_mult[k];
                        this->m_AJSbs[k-1]=a_J_mult[k];
                        this->m_SigmaSbs[k-1]=sigma_mult[k];
                        this->m_CSbs[k-1]=C_mult[k];
                        this->m_KSbs[k-1]=K_mult[k];
                        this->m_LambdaSbs[k-1]=lambda_mult[k];

                        this->m_AlphaISbs[k-1]=alpha_I_mult[k];
                        this->m_AlphaJSbs[k-1]=alpha_J_mult[k];
                };



        label1:;



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
        if(alp_distr)
        {
                Int4 j;
                for(j=1;j<=nalp_for_lambda_simulation;j++)
                {
                        delete (array_positive<double>*)alp_distr[j];alp_distr[j]=NULL;
                };

                delete[]alp_distr;alp_distr=NULL;
        };



        if(alp_distr_errors)
        {
                Int4 j;
                for(j=1;j<=nalp_for_lambda_simulation;j++)
                {
                        delete (array_positive<double>*)alp_distr_errors[j];alp_distr_errors[j]=NULL;
                };

                delete[]alp_distr_errors;alp_distr_errors=NULL;
        };


        if(alp_mult_distr)
        {
                Int4 k,j;
                for(k=1;k<=d_mult_number;k++)
                {
                        if(alp_mult_distr[k])
                        {
                                for(j=1;j<=nalp_for_lambda_simulation;j++)
                                {
                                        delete (array_positive<double>*)alp_mult_distr[k][j];alp_mult_distr[k][j]=NULL;
                                };

                                delete[]alp_mult_distr[k];alp_mult_distr[k]=NULL;
                        };
                };

                delete[]alp_mult_distr;alp_mult_distr=NULL;
        };


        if(alp_mult_distr_errors)
        {
                Int4 k,j;
                for(k=1;k<=d_mult_number;k++)
                {
                        if(alp_mult_distr_errors[k])
                        {
                                for(j=1;j<=nalp_for_lambda_simulation;j++)
                                {
                                        delete (array_positive<double>*)alp_mult_distr_errors[k][j];alp_mult_distr_errors[k][j]=NULL;
                                };

                                delete[]alp_mult_distr_errors[k];alp_mult_distr_errors[k]=NULL;
                        };
                };

                delete[]alp_mult_distr_errors;alp_mult_distr_errors=NULL;
        };


        delete[]d_mult_realizations;d_mult_realizations=NULL;
        delete[]d_mult_K_realizations;d_mult_K_realizations=NULL;

        delete[]lambda_mult;lambda_mult=NULL;
        delete[]lambda_mult_error;lambda_mult_error=NULL;

        delete[]C_mult;C_mult=NULL;
        delete[]C_mult_error;C_mult_error=NULL;

        delete[]a_I_mult;a_I_mult=NULL;
        delete[]a_I_mult_error;a_I_mult_error=NULL;
        delete[]a_J_mult;a_J_mult=NULL;
        delete[]a_J_mult_error;a_J_mult_error=NULL;
        delete[]sigma_mult;sigma_mult=NULL;
        delete[]sigma_mult_error;sigma_mult_error=NULL;
        delete[]alpha_I_mult;alpha_I_mult=NULL;
        delete[]alpha_I_mult_error;alpha_I_mult_error=NULL;
        delete[]alpha_J_mult;alpha_J_mult=NULL;
        delete[]alpha_J_mult_error;alpha_J_mult_error=NULL;

        delete[]K_C_mult;K_C_mult=NULL;
        delete[]K_C_mult_error;K_C_mult_error=NULL;

        delete[]K_mult;K_mult=NULL;
        delete[]K_mult_error;K_mult_error=NULL;

        if(ee_error_flag)
        {
                throw error(ee_error.st,ee_error.error_code);
        };

};

void alp_sim::error_in_calculate_main_parameters2m(
double C,
double &C_error,
double C_mult2,
double C_mult2_error)
{
        if(C!=0&&C_mult2!=0)
        {
                C_error=fabs(C*C_mult2_error/C_mult2);
        }
        else
        {
                C_error=C_mult2_error;
        };
};

alp_sim::~alp_sim()//destructor
{
        Int4 i;
        for(i=0;i<d_n_alp_obj;i++)
        {
                delete d_alp_obj->d_elem[i];d_alp_obj->d_elem[i]=NULL;
        };


        if(d_alp_data)
        {
                d_alp_data->d_memory_size_in_MB-=sizeof(alp)*d_n_alp_obj/mb_bytes;
        };

        delete d_alp_obj;d_alp_obj=NULL;
        if(d_alp_data)
        {
                d_alp_data->d_memory_size_in_MB-=(double)(sizeof(array_positive<alp*>))/mb_bytes;
        };


};

void alp_sim::kill(
bool check_time_,
Int4 ind1_,
Int4 ind2_,
Int4 M_min_,
double lambda_,
double eps_K_,
double &K_C_,
double &K_C_error_,
Int4 &level_,
Int4 &diff_opt_)
{

        bool flag=false;
        Int4 current_level=(Int4)floor(M_min_*0.5);
        Int4 recommended_level;
        Int4 number_of_unsucesful_objects=0;

        
        Int4 i;
        for(i=ind1_;i<=ind2_;i++)
        {
        
                alp* &alp_obj_tmp=d_alp_obj->d_elem[i];
                if(i-ind1_+1>alp_obj_tmp->d_alp_data->d_minimum_realizations_number)
                {
                        alp_obj_tmp->d_check_time_flag=check_time_;
                        alp_obj_tmp->d_time_error_flag=check_time_;
                };
        };

        do{
                Int4 i;
                for(i=ind1_;i<=ind2_;i++)
                {
                        alp* &alp_obj_tmp=d_alp_obj->d_elem[i];
                        bool flag=false;
                        while(!flag)
                        {
                                alp_obj_tmp->d_sentinels_flag=false;
                                alp_obj_tmp->kill_upto_level(M_min_,current_level);
                                if(!alp_obj_tmp->d_success)
                                {
                                        number_of_unsucesful_objects++;
                                        if(number_of_unsucesful_objects>5+(ind2_-ind1_+1)*
                                                d_alp_obj->d_alp_data->d_eps_K 
                                                )
                                        {
                                                throw error("The program cannot estimate the parameters.\nPlease try to increase the allowed amount of memory or change parameters of the scoring system to ensure the logarithmic regime of the alignment score\n",1);
                                        };
                                        delete alp_obj_tmp;alp_obj_tmp=NULL;
                                        alp_obj_tmp=new alp(d_alp_data);
                                        alp_data::assert_mem(alp_obj_tmp);

                                        if(i-ind1_+1>alp_obj_tmp->d_alp_data->d_minimum_realizations_number)
                                        {
                                                alp_obj_tmp->d_check_time_flag=check_time_;
                                                alp_obj_tmp->d_time_error_flag=check_time_;
                                        };

                                        bool flag=false;
                                        while(!flag)
                                        {
                                                alp_obj_tmp->simulate_alp_upto_the_given_level(M_min_);
                                                if(!alp_obj_tmp->d_success)
                                                {
                                                        number_of_unsucesful_objects++;
                                                        if(number_of_unsucesful_objects>5+(ind2_-ind1_+1)*
                                                                d_alp_obj->d_alp_data->d_eps_K 
                                                                )
                                                        {
                                                                throw error("The program cannot estimate the parameters.\nPlease try to increase the allowed amount of memory or change parameters of the scoring system to ensure the logarithmic regime of the alignment score\n",1);
                                                        };
                                                };
                                                flag=alp_obj_tmp->d_success;
                                        };

                                };
                                flag=alp_obj_tmp->d_success;
                        };
                };


                flag=check_K_criterion_during_killing(
                ind1_,
                ind2_,
                lambda_,
                eps_K_,
                current_level,
                recommended_level,
                diff_opt_,
                K_C_,
                K_C_error_);

                current_level=recommended_level;

        }
        while(!flag);

        level_=current_level;

};

void alp_sim::get_single_realization(
bool check_time_,
Int4 M_min_,
Int4 nalp_,
bool killing_flag_,
Int4 level_,
Int4 diff_opt_,
alp *&obj_,
bool &sucess_flag_,
double &d_eps_)
{
        if(!obj_)
        {
                obj_=new alp(d_alp_data);
                alp_data::assert_mem(obj_);
                d_alp_data->d_memory_size_in_MB+=sizeof(alp)/mb_bytes;
        };
        obj_->d_single_realiztion_calculation_flag=true;
        obj_->d_check_time_flag=check_time_;

        d_eps_=d_alp_data->Tmin(d_alp_data->d_eps_K,d_alp_data->d_eps_lambda);

        

        alp*&obj=obj_;

        obj->d_diff_opt=diff_opt_;

        obj->d_sentinels_flag=d_alp_data->d_sentinels_flag;

        sucess_flag_=true;

        while(obj->d_nalp<nalp_)
        {
                obj->simulate_next_alp();
                if(!obj->d_success)
                {
                        sucess_flag_=false;
                        delete obj_;obj_=NULL;
                        d_eps_=d_alp_data->d_eps_lambda;
                        d_alp_data->d_memory_size_in_MB-=sizeof(alp)/mb_bytes;
                        return;
                };
        };


        if(killing_flag_)
        {
                obj->kill_upto_level(M_min_,level_);                
                if(!obj->d_success)
                {
                        sucess_flag_=false;
                        delete obj_;obj_=NULL;
                        d_eps_=d_alp_data->d_eps_K;
                        d_alp_data->d_memory_size_in_MB-=sizeof(alp)/mb_bytes;
                        return;
                };
        };

};


void alp_sim::get_minimal_simulation(
Int4 ind1_,
Int4 ind2_,
Int4 &M_min_,
Int4 &nalp_,
Int4 &nalp_lambda_,
bool C_calculation_,
bool check_time_flag_)
{
        bool ee_error_flag=false;
        error ee_error("",0);

        Int4 &alp_number=nalp_;
        alp_number=0;

        void **alp_distr=NULL;
        void **alp_distr_errors=NULL;


        try
        {
        try
        {


                if(d_n_alp_obj<ind1_||d_n_alp_obj-1>ind2_)
                {
                        throw error("Unexpected error\n",4);
                };

                Int4 nalp_lambda_equilibration=-1;

                //create the objects
                Int4 i;
                for(i=d_n_alp_obj;i<=ind2_;i++)
                {
                        d_alp_obj->set_elem(i,NULL);
                        d_n_alp_obj=i+1;


                        alp *&alp_obj_tmp=d_alp_obj->d_elem[i];

                        alp_obj_tmp=new alp(d_alp_data);
                        alp_data::assert_mem(alp_obj_tmp);

                        d_alp_data->d_memory_size_in_MB+=sizeof(alp)/mb_bytes;
                        

                        alp_obj_tmp->d_check_time_flag=check_time_flag_;
                        alp_obj_tmp->d_time_error_flag=check_time_flag_;

                };
        

                bool M_min_flag=false;
                bool nalp_flag=false;


                bool criterion_flag=false;
                Int4 number_of_fails=0;
                Int4 number_of_fails_threshold=5;

                do{

                        Int4 number_of_unsuccessful_objects=0;
                        Int4 i;
                        for(i=ind1_;i<=ind2_;i++)
                        {
                                alp* &alp_obj_tmp=d_alp_obj->d_elem[i];

                                alp_obj_tmp->d_check_time_flag=check_time_flag_;
                                alp_obj_tmp->d_time_error_flag=check_time_flag_;


                                if(alp_obj_tmp->d_nalp<alp_number+1)
                                {

                                        alp_obj_tmp->simulate_alp_upto_the_given_number(alp_number+1);

                                        if(!alp_obj_tmp->d_success)
                                        {
                                                number_of_unsuccessful_objects++;
                                                delete alp_obj_tmp;
                                                alp_obj_tmp=NULL;

                                                if(number_of_unsuccessful_objects>5+(ind2_-d_n_alp_obj+1)*
                                                        d_alp_obj->d_alp_data->d_eps_lambda*(alp_number+1)
                                                        )
                                                {
                                                        
                                                        throw error("The parameters cannot be correctly calculated for the given accuracy, calculation time and memory usage\n",1);
                                                };


                                                bool success2=false;
                                                while(!success2)
                                                {
                                                        alp_obj_tmp=new alp(d_alp_data);
                                                        alp_data::assert_mem(alp_obj_tmp);


                                                        Int4 j;
                                                        for(j=0;j<=alp_number;j++)
                                                        {
                                                                alp_obj_tmp->simulate_alp_upto_the_given_number(j+1);
                                                        };

                                                        success2=alp_obj_tmp->d_success;

                                                        

                                                        if(!success2)
                                                        {
                                                                delete alp_obj_tmp;
                                                                alp_obj_tmp=NULL;
                                                                number_of_unsuccessful_objects++;
                                                                if(number_of_unsuccessful_objects>5+(ind2_-d_n_alp_obj+1)*
                                                                        d_alp_obj->d_alp_data->d_eps_lambda*(alp_number+1)
                                                                        )
                                                                {
                                                                        
                                                                        throw error("The parameters cannot be correctly calculated for the given accuracy, calculation time and memory usage\n",1);
                                                                };
                                                        };

                                                };
                                        };

                                };

                        };

                        alp_number++;
                        


                        bool inside_simulation_flag=false;

                        criterion_flag=the_criterion(
                                alp_number,
                                nalp_lambda_,
                                0,
                                ind2_,
                                alp_distr,
                                alp_distr_errors,
                                M_min_,
                                M_min_flag,
                                nalp_flag,
                                inside_simulation_flag,
                                C_calculation_);

                        
                        if(nalp_lambda_equilibration==-1&&nalp_flag)
                        {
                                nalp_lambda_equilibration=alp_number;
                        };
                        


                        if(!inside_simulation_flag)
                        {
                                number_of_fails++;

                                Int4 i;
                                if(alp_distr)
                                {
                                        for(i=1;i<=alp_number;i++)
                                        {
                                                delete (array_positive<double>*)alp_distr[i];alp_distr[i]=NULL;
                                        };

                                        delete[]alp_distr;alp_distr=NULL;
                                };

                                if(alp_distr_errors)
                                {
                                        for(i=1;i<=alp_number;i++)
                                        {
                                                delete (array_positive<double>*)alp_distr_errors[i];alp_distr_errors[i]=NULL;
                                        };

                                        delete[]alp_distr_errors;alp_distr_errors=NULL;
                                };


                                M_min_flag=false;
                                nalp_flag=false;


                                alp_distr=NULL;
                                alp_distr_errors=NULL;

                                alp_number=0;

                                criterion_flag=false;
                                


                                for(i=ind1_;i<=ind2_;i++)
                                {
                                        alp* &alp_obj_tmp=d_alp_obj->d_elem[i];
                                        delete alp_obj_tmp;alp_obj_tmp=NULL;

                                };

                                if(number_of_fails>number_of_fails_threshold)
                                {
                                        throw error("The program is not able to calculate the parameters accurately.\nPlease try to increase the allowed calculation time and memory limit\n",1);
                                };

                                for(i=ind1_;i<=ind2_;i++)
                                {
                                        alp *alp_obj_tmp;
                                        alp_obj_tmp=new alp(d_alp_data);
                                        alp_data::assert_mem(alp_obj_tmp);
                                        d_alp_obj->set_elem(i,alp_obj_tmp);


                                        alp_obj_tmp->d_check_time_flag=check_time_flag_;
                                        alp_obj_tmp->d_time_error_flag=check_time_flag_;

                                };

                        };

                }
                while(!criterion_flag);

                
                nalp_lambda_=alp_data::Tmax(nalp_lambda_equilibration,nalp_lambda_);
                nalp_lambda_=alp_data::Tmin(nalp_lambda_,nalp_);

                nalp_lambda_=nalp_;

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

        if(alp_distr)
        {
                Int4 i;
                for(i=1;i<=nalp_;i++)
                {
                        delete (array_positive<double>*)alp_distr[i];alp_distr[i]=NULL;
                };

                delete[]alp_distr;alp_distr=NULL;
        };

        if(alp_distr_errors)
        {
                Int4 i;
                for(i=1;i<=nalp_;i++)
                {
                        delete (array_positive<double>*)alp_distr_errors[i];alp_distr_errors[i]=NULL;
                };

                delete[]alp_distr_errors;alp_distr_errors=NULL;
        };

        if(ee_error_flag)
        {
                throw error(ee_error.st,ee_error.error_code);
        };

};

bool alp_sim::the_criterion(//criteria of stopping of the simulating ALP
//if the function returns true then calculates optimal M_min and ALP number
//sets the flags M_min_flag_ and nalp_flag_ checking the corresponding condition
Int4 upto_nalp_,
Int4 &nalp_for_lambda_simulation_,
Int4 ind1_,
Int4 ind2_,
void **&alp_distr,
void **&alp_distr_errors,
Int4 &M_min_,
bool &M_min_flag_,
bool &nalp_flag_,
bool &inside_simulation_flag_,
bool C_calculation_)
{

        nalp_flag_=false;
        M_min_flag_=false;

        if(ind1_>ind2_)
        {
                throw error("Unexpected error\n",4);
        };



        double lambda;
        double lambda_error;

        double test_difference;
        double test_difference_error;

        Int4 nalp=upto_nalp_;

        if(nalp<1)
        {
                throw error("Unexpected error\n",4);
        };


        get_and_allocate_alp_distribution(
        ind1_,
        ind2_,
        alp_distr,
        alp_distr_errors,
        nalp);


        calculate_lambda(
        true,
        upto_nalp_,
        nalp_for_lambda_simulation_,
        inside_simulation_flag_,
        ind1_,
        ind2_,
        alp_distr,
        alp_distr_errors,
        lambda,
        lambda_error,
        test_difference,
        test_difference_error);

        if(!inside_simulation_flag_)
        {
                return false;
        };


        d_lambda_tmp->set_elem(upto_nalp_,lambda);
        d_lambda_tmp_errors->set_elem(upto_nalp_,lambda_error);


        if(C_calculation_)
        {
                double C;
                double C_error;


                calculate_C(
                0,
                upto_nalp_,
                ind1_,
                ind2_,
                alp_distr,
                alp_distr_errors,
                lambda,
                lambda_error,
                C,
                C_error);


                d_C_tmp->set_elem(upto_nalp_,C);
                d_C_tmp_errors->set_elem(upto_nalp_,C_error);

        };




        if(nalp>=1)
        {

                if(test_difference<=test_difference_error)
                {
                        nalp_flag_=true;
                        M_min_flag_=check_K_criterion(
                        nalp,
                        ind1_,
                        ind2_,
                        lambda,
                        d_alp_data->d_eps_K,
                        M_min_);

                        return M_min_flag_;
                };
        };


        return false;
};

double alp_sim::lambda_exp(
Int4 &i_,
double *&exp_array_)
{
        if(exp_array_[i_]==-1)
        {
                throw error("The program is not able to calculate the parameters; rescaling penalties and scoring matrix might help\n",3);
        };

        return exp_array_[i_];
};

void alp_sim::calculate_FSC(
Int4 nalp_,
Int4 ind1_,
Int4 ind2_,
void **alp_distr,
void **alp_distr_errors,
double lambda_,
double lambda_error_,
double &a_I_,
double &a_I_error_,
double &a_J_,
double &a_J_error_,
double &sigma_,
double &sigma_error_,
double &alpha_I_,
double &alpha_I_error_,
double &alpha_J_,
double &alpha_J_error_)
{
        bool ee_error_flag=false;
        error ee_error("",0);

        double *exp_array=NULL;

        double *delta_E=NULL;
        double *delta_E_error=NULL;

        double *delta_E_E=NULL;
        double *delta_E_E_error=NULL;


        double *delta_I=NULL;
        double *delta_I_error=NULL;

        double *delta_J=NULL;
        double *delta_J_error=NULL;

        double *delta_I_I=NULL;
        double *delta_I_I_error=NULL;

        double *delta_I_J=NULL;
        double *delta_I_J_error=NULL;

        double *delta_J_J=NULL;
        double *delta_J_J_error=NULL;

        double *cov_J_J=NULL;
        double *cov_J_J_error=NULL;

        double *cov_I_J=NULL;
        double *cov_I_J_error=NULL;

        double *cov_I_I=NULL;
        double *cov_I_I_error=NULL;

        double *cov_E_E=NULL;
        double *cov_E_E_error=NULL;


        try
        {
        try
        {


                if(nalp_<1)
                {
                        throw error("Unexpected error\n",4);
                };

                array_positive<double>* tmp=((array_positive<double>*)alp_distr[nalp_]);
                Int4 dim=tmp->d_dim;

                exp_array=new double[dim+1];
                alp_data::assert_mem(exp_array);


                Int4 i;
                for(i=0;i<=dim;i++)
                {
                        double tmp=(double)i*lambda_;
                        if(tmp<dbl_max_log)
                        {
                                exp_array[i]=exp(tmp);
                        }
                        else
                        {
                                exp_array[i]=-1;
                        };
                };

                


                delta_E=new double[nalp_];
                alp_data::assert_mem(delta_E);
                delta_E_error=new double[nalp_];
                alp_data::assert_mem(delta_E_error);

                delta_E_E=new double[nalp_];
                alp_data::assert_mem(delta_E_E);
                delta_E_E_error=new double[nalp_];
                alp_data::assert_mem(delta_E_E_error);

                cov_E_E=new double[nalp_];
                alp_data::assert_mem(cov_E_E);
                cov_E_E_error=new double[nalp_];
                alp_data::assert_mem(cov_E_E_error);


                delta_I=new double[nalp_];
                alp_data::assert_mem(delta_I);
                delta_I_error=new double[nalp_];
                alp_data::assert_mem(delta_I_error);

                delta_J=new double[nalp_];
                alp_data::assert_mem(delta_J);
                delta_J_error=new double[nalp_];
                alp_data::assert_mem(delta_J_error);

                delta_I_I=new double[nalp_];
                alp_data::assert_mem(delta_I_I);
                delta_I_I_error=new double[nalp_];
                alp_data::assert_mem(delta_I_I_error);

                delta_I_J=new double[nalp_];
                alp_data::assert_mem(delta_I_J);
                delta_I_J_error=new double[nalp_];
                alp_data::assert_mem(delta_I_J_error);

                delta_J_J=new double[nalp_];
                alp_data::assert_mem(delta_J_J);
                delta_J_J_error=new double[nalp_];
                alp_data::assert_mem(delta_J_J_error);

                cov_J_J=new double[nalp_];
                alp_data::assert_mem(cov_J_J);
                cov_J_J_error=new double[nalp_];
                alp_data::assert_mem(cov_J_J_error);

                cov_I_J=new double[nalp_];
                alp_data::assert_mem(cov_I_J);
                cov_I_J_error=new double[nalp_];
                alp_data::assert_mem(cov_I_J_error);

                cov_I_I=new double[nalp_];
                alp_data::assert_mem(cov_I_I);
                cov_I_I_error=new double[nalp_];
                alp_data::assert_mem(cov_I_I_error);

                Int4 j;
                for(j=0;j<nalp_;j++)
                {
                        delta_E[j]=0.0;
                        delta_E_error[j]=0.0;

                        delta_E_E[j]=0.0;
                        delta_E_E_error[j]=0.0;

                        delta_I[j]=0.0;
                        delta_I_error[j]=0.0;
                        delta_J[j]=0.0;
                        delta_J_error[j]=0.0;

                        delta_I_I[j]=0.0;

                        delta_I_I_error[j]=0.0;
                        delta_I_J[j]=0.0;
                        delta_I_J_error[j]=0.0;
                        delta_J_J[j]=0.0;
                        delta_J_J_error[j]=0.0;
                };

                
                for(i=ind1_;i<=ind2_;i++)
                {
                        alp* alp_obj_tmp=d_alp_obj->d_elem[i];

                        Int4 j;
                        for(j=1;j<=nalp_;j++)
                        {
                                Int4 j_1=j-1;

                                Int4 &E_j_1=alp_obj_tmp->d_alp->d_elem[j_1];
                                Int4 &E_j=alp_obj_tmp->d_alp->d_elem[j];
                                double &weight_j=alp_obj_tmp->d_alp_weights->d_elem[j];

                                Int4 &I_j_1=alp_obj_tmp->d_H_I->d_elem[j_1];
                                Int4 &I_j=alp_obj_tmp->d_H_I->d_elem[j];

                                Int4 &J_j_1=alp_obj_tmp->d_H_J->d_elem[j_1];
                                Int4 &J_j=alp_obj_tmp->d_H_J->d_elem[j];

                                double delta_I_tmp=(I_j-I_j_1)*lambda_exp(E_j,exp_array)*weight_j;
                                double delta_J_tmp=(J_j-J_j_1)*lambda_exp(E_j,exp_array)*weight_j;
                                double delta_E_tmp=(E_j-E_j_1)*lambda_exp(E_j,exp_array)*weight_j;
                                double delta_E_E_tmp=(E_j-E_j_1)*(E_j-E_j_1)*lambda_exp(E_j,exp_array)*weight_j;

                                
                                double delta_I_I_tmp=delta_I_tmp*(I_j-I_j_1);
                                double delta_J_J_tmp=delta_J_tmp*(J_j-J_j_1);
                                double delta_I_J_tmp=delta_I_tmp*(J_j-J_j_1);

        


                                delta_E[j_1]+=delta_E_tmp;
                                delta_E_error[j_1]+=delta_E_tmp*delta_E_tmp;

                                delta_E_E[j_1]+=delta_E_E_tmp;
                                delta_E_E_error[j_1]+=delta_E_E_tmp*delta_E_E_tmp;

                                delta_I[j_1]+=delta_I_tmp;
                                delta_I_error[j_1]+=delta_I_tmp*delta_I_tmp;
                                delta_J[j_1]+=delta_J_tmp;
                                delta_J_error[j_1]+=delta_J_tmp*delta_J_tmp;

                                delta_I_I[j_1]+=delta_I_I_tmp;
                                delta_I_I_error[j_1]+=delta_I_I_tmp*delta_I_I_tmp;

                                delta_I_J[j_1]+=delta_I_J_tmp;
                                delta_I_J_error[j_1]+=delta_I_J_tmp*delta_I_J_tmp;

                                delta_J_J[j_1]+=delta_J_J_tmp;
                                delta_J_J_error[j_1]+=delta_J_J_tmp*delta_J_J_tmp;
                                
                        };
                };


                double ind_diff=(double)(ind2_-ind1_+1);
                for(j=0;j<nalp_;j++)
                {
                        delta_E[j]/=ind_diff;
                        delta_E_error[j]/=ind_diff;
                        delta_E_error[j]-=delta_E[j]*delta_E[j];
                        delta_E_error[j]/=ind_diff;
                        delta_E_error[j]=alp_reg::sqrt_for_errors(delta_E_error[j]);

                        delta_E_E[j]/=ind_diff;
                        delta_E_E_error[j]/=ind_diff;
                        delta_E_E_error[j]-=delta_E_E[j]*delta_E_E[j];
                        delta_E_E_error[j]/=ind_diff;


                        delta_I[j]/=ind_diff;
                        delta_I_error[j]/=ind_diff;
                        delta_I_error[j]-=delta_I[j]*delta_I[j];
                        delta_I_error[j]/=ind_diff;
                        delta_I_error[j]=alp_reg::sqrt_for_errors(delta_I_error[j]);

                        delta_J[j]/=ind_diff;
                        delta_J_error[j]/=ind_diff;
                        delta_J_error[j]-=delta_J[j]*delta_J[j];
                        delta_J_error[j]/=ind_diff;
                        delta_J_error[j]=alp_reg::sqrt_for_errors(delta_J_error[j]);

                        delta_I_J[j]/=ind_diff;
                        delta_I_J_error[j]/=ind_diff;
                        delta_I_J_error[j]-=delta_I_J[j]*delta_I_J[j];
                        delta_I_J_error[j]/=ind_diff;


                        delta_I_I[j]/=ind_diff;
                        delta_I_I_error[j]/=ind_diff;
                        delta_I_I_error[j]-=delta_I_I[j]*delta_I_I[j];
                        delta_I_I_error[j]/=ind_diff;


                        delta_J_J[j]/=ind_diff;
                        delta_J_J_error[j]/=ind_diff;
                        delta_J_J_error[j]-=delta_J_J[j]*delta_J_J[j];
                        delta_J_J_error[j]/=ind_diff;


                        cov_I_J[j]=delta_I_J[j]-delta_I[j]*delta_J[j];
                        cov_I_I[j]=delta_I_I[j]-delta_I[j]*delta_I[j];
                        cov_J_J[j]=delta_J_J[j]-delta_J[j]*delta_J[j];

                        cov_E_E[j]=delta_E_E[j]-delta_E[j]*delta_E[j];

                        cov_I_J_error[j]=alp_reg::error_of_the_product(delta_I[j],delta_I_error[j],delta_J[j],delta_J_error[j]);
                        cov_I_J_error[j]=sqrt(delta_I_J_error[j]+cov_I_J_error[j]*cov_I_J_error[j]);

                        cov_I_I_error[j]=alp_reg::error_of_the_product(delta_I[j],delta_I_error[j],delta_I[j],delta_I_error[j]);
                        cov_I_I_error[j]=sqrt(delta_I_I_error[j]+cov_I_I_error[j]*cov_I_I_error[j]);

                        cov_J_J_error[j]=alp_reg::error_of_the_product(delta_J[j],delta_J_error[j],delta_J[j],delta_J_error[j]);
                        cov_J_J_error[j]=sqrt(delta_J_J_error[j]+cov_J_J_error[j]*cov_J_J_error[j]);

                        cov_E_E_error[j]=alp_reg::error_of_the_product(delta_E[j],delta_E_error[j],delta_E[j],delta_E_error[j]);
                        cov_E_E_error[j]=sqrt(delta_E_E_error[j]+cov_E_E_error[j]*cov_E_E_error[j]);

                };


                //regression

                double beta1=0;
                double beta1_error=0;

                Int4 number_of_elements=nalp_;

                bool cut_left_tail=true;
                bool cut_right_tail=false;

                double y=2;

                Int4 k1_opt;
                Int4 k2_opt;


                double delta_I_aver;
                double delta_I_aver_error;



                bool res_was_calculated;
                
                alp_reg::robust_regression_sum_with_cut_LSM_beta1_is_defined(
                0,
                number_of_elements,
                delta_I,
                delta_I_error,
                cut_left_tail,
                cut_right_tail,
                y,
                delta_I_aver,
                beta1,
                delta_I_aver_error,
                beta1_error,
                k1_opt,
                k2_opt,
                res_was_calculated);

                if(!res_was_calculated)
                {
                        throw error("The program cannot estimate the parameters; please repeat the calculation2\n",2);
                };



                double delta_J_aver;
                double delta_J_aver_error;

                
                
                alp_reg::robust_regression_sum_with_cut_LSM_beta1_is_defined(
                0,
                number_of_elements,
                delta_J,
                delta_J_error,
                cut_left_tail,
                cut_right_tail,
                y,
                delta_J_aver,
                beta1,
                delta_J_aver_error,
                beta1_error,
                k1_opt,
                k2_opt,
                res_was_calculated);

                if(!res_was_calculated)
                {
                        throw error("The program cannot estimate the parameters; please repeat the calculation3\n",2);
                };


                double delta_E_aver;
                double delta_E_aver_error;
                
                alp_reg::robust_regression_sum_with_cut_LSM_beta1_is_defined(
                0,
                number_of_elements,
                delta_E,
                delta_E_error,
                cut_left_tail,
                cut_right_tail,
                y,
                delta_E_aver,
                beta1,
                delta_E_aver_error,
                beta1_error,
                k1_opt,
                k2_opt,
                res_was_calculated);

                if(!res_was_calculated)
                {
                        throw error("The program cannot estimate the parameters; please repeat the calculation4\n",2);
                };

                double cov_I_I_aver;
                double cov_I_I_aver_error;

                double cov_I_J_aver;
                double cov_I_J_aver_error;

                double cov_J_J_aver;
                double cov_J_J_aver_error;

                double cov_E_E_aver;
                double cov_E_E_aver_error;


                alp_reg::robust_regression_sum_with_cut_LSM_beta1_is_defined(
                0,
                number_of_elements,
                cov_I_J,
                cov_I_J_error,
                cut_left_tail,
                cut_right_tail,
                y,
                cov_I_J_aver,
                beta1,
                cov_I_J_aver_error,
                beta1_error,
                k1_opt,
                k2_opt,
                res_was_calculated);

                if(!res_was_calculated)
                {
                        throw error("The program cannot estimate the parameters; please repeat the calculation5\n",2);
                };

                alp_reg::robust_regression_sum_with_cut_LSM_beta1_is_defined(
                0,
                number_of_elements,
                cov_I_I,
                cov_I_I_error,
                cut_left_tail,
                cut_right_tail,
                y,
                cov_I_I_aver,
                beta1,
                cov_I_I_aver_error,
                beta1_error,
                k1_opt,
                k2_opt,
                res_was_calculated);

                if(!res_was_calculated)
                {
                        throw error("The program cannot estimate the parameters; please repeat the calculation6\n",2);
                };

                alp_reg::robust_regression_sum_with_cut_LSM_beta1_is_defined(
                0,
                number_of_elements,
                cov_J_J,
                cov_J_J_error,
                cut_left_tail,
                cut_right_tail,
                y,
                cov_J_J_aver,
                beta1,
                cov_J_J_aver_error,
                beta1_error,
                k1_opt,
                k2_opt,
                res_was_calculated);

                if(!res_was_calculated)
                {
                        throw error("The program cannot estimate the parameters; please repeat the calculation7\n",2);
                };

                alp_reg::robust_regression_sum_with_cut_LSM_beta1_is_defined(
                0,
                number_of_elements,
                cov_E_E,
                cov_E_E_error,
                cut_left_tail,
                cut_right_tail,
                y,
                cov_E_E_aver,
                beta1,
                cov_E_E_aver_error,
                beta1_error,
                k1_opt,
                k2_opt,
                res_was_calculated);

                if(!res_was_calculated)
                {
                        throw error("The program cannot estimate the parameters; please repeat the calculation7_E\n",2);
                };



                if(delta_E_aver<=0)
                {
                        throw error("FSC parameters estimations failed\n",2);
                };



                a_I_=delta_I_aver/delta_E_aver;
                a_I_error_=alp_reg::error_of_the_ratio(delta_I_aver,delta_I_aver_error,delta_E_aver,delta_E_aver_error);
                a_J_=delta_J_aver/delta_E_aver;
                a_J_error_=alp_reg::error_of_the_ratio(delta_J_aver,delta_J_aver_error,delta_E_aver,delta_E_aver_error);
                
                alpha_calculation(
                delta_I_aver,
                delta_I_aver_error,
                delta_J_aver,
                delta_J_aver_error,
                delta_E_aver,
                delta_E_aver_error,
                cov_E_E_aver,
                cov_E_E_aver_error,
                cov_I_J_aver,
                cov_I_J_aver_error,
                sigma_,
                sigma_error_);


                alpha_calculation(
                delta_I_aver,
                delta_I_aver_error,
                delta_I_aver,
                delta_I_aver_error,
                delta_E_aver,
                delta_E_aver_error,
                cov_E_E_aver,
                cov_E_E_aver_error,
                cov_I_I_aver,
                cov_I_I_aver_error,
                alpha_I_,
                alpha_I_error_);



                alpha_calculation(
                delta_J_aver,
                delta_J_aver_error,
                delta_J_aver,
                delta_J_aver_error,
                delta_E_aver,
                delta_E_aver_error,
                cov_E_E_aver,
                cov_E_E_aver_error,
                cov_J_J_aver,
                cov_J_J_aver_error,
                alpha_J_,
                alpha_J_error_);


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
        delete[]exp_array;exp_array=NULL;

        delete[]delta_E;delta_E=NULL;
        delete[]delta_E_error;delta_E_error=NULL;
        delete[]delta_I;delta_I=NULL;
        delete[]delta_I_error;delta_I_error=NULL;
        delete[]delta_J;delta_J=NULL;
        delete[]delta_J_error;delta_J_error=NULL;

        delete[]delta_I_J;delta_I_J=NULL;
        delete[]delta_I_J_error;delta_I_J_error=NULL;
        delete[]delta_J_J;delta_J_J=NULL;
        delete[]delta_J_J_error;delta_J_J_error=NULL;
        delete[]delta_I_I;delta_I_I=NULL;
        delete[]delta_I_I_error;delta_I_I_error=NULL;

        delete[]cov_I_J;cov_I_J=NULL;
        delete[]cov_I_J_error;cov_I_J_error=NULL;
        delete[]cov_J_J;cov_J_J=NULL;
        delete[]cov_J_J_error;cov_J_J_error=NULL;
        delete[]cov_I_I;cov_I_I=NULL;
        delete[]cov_I_I_error;cov_I_I_error=NULL;

        if(ee_error_flag)
        {
                throw error(ee_error.st,ee_error.error_code);
        };

};

void alp_sim::alpha_calculation(
double delta_I_aver_,
double delta_I_aver_error_,
double delta_J_aver_,
double delta_J_aver_error_,
double delta_E_aver_,
double delta_E_aver_error_,
double cov_E_E_aver_,
double cov_E_E_aver_error_,
double cov_I_J_aver_,
double cov_I_J_aver_error_,
double &alpha_,
double &alpha_error_)
{
        double nom1_1=delta_I_aver_*delta_J_aver_;
        double nom2_2=delta_E_aver_*delta_E_aver_;

        double den=nom2_2*delta_E_aver_;

        double nom1=nom1_1*cov_E_E_aver_;
        double nom2=nom2_2*cov_I_J_aver_;

        alpha_=(nom1+nom2)/den;

        
        double nom1_alpha_error=alp_reg::error_of_the_product(delta_I_aver_,delta_I_aver_error_,delta_J_aver_,delta_J_aver_error_);
        nom1_alpha_error=alp_reg::error_of_the_product(nom1_1,nom1_alpha_error,cov_E_E_aver_,cov_E_E_aver_error_);

        
        double nom2_alpha_error_2=alp_reg::error_of_the_product(delta_E_aver_,delta_E_aver_error_,delta_E_aver_,delta_E_aver_error_);
        double nom2_alpha_error=alp_reg::error_of_the_product(nom2_2,nom2_alpha_error_2,cov_I_J_aver_,cov_I_J_aver_error_);

        
        double den_alpha_error=alp_reg::error_of_the_product(nom2_2,nom2_alpha_error_2,delta_E_aver_,delta_E_aver_error_);

        double nom_alpha_error=alp_reg::error_of_the_sum(nom1,nom1_alpha_error,nom2,nom2_alpha_error);

        alpha_error_=alp_reg::error_of_the_ratio(nom1+nom2,nom_alpha_error,den,den_alpha_error);

};



void alp_sim::calculate_C(
Int4 starting_point,
Int4 nalp_,
Int4 ind1_,
Int4 ind2_,
void **alp_distr,
void **alp_distr_errors,
double lambda_,
double lambda_error_,
double &C_,
double &C_error_)
{
        bool ee_error_flag=false;
        error ee_error("",0);

        double *P=NULL;
        double *P_errors=NULL;
        double *values_P_ratio=NULL;
        double *errors_P_ratio=NULL;

        double *E=NULL;
        double *E_errors=NULL;

        double *E_T_beta=NULL;
        double *E_T_beta_errors=NULL;


        try
        {
        try
        {

                Int4 total_number_of_ALP=nalp_;

                if(total_number_of_ALP<1)
                {
                        throw error("Unexpected error\n",4);
                };


                //1)P(beta=infinity)
                Int4 j;


                P=new double[total_number_of_ALP+1];
                alp_data::assert_mem(P);
                P_errors=new double[total_number_of_ALP+1];
                alp_data::assert_mem(P_errors);

                P[0]=1.0;
                P_errors[0]=0.0;

                
                for(j=1;j<=total_number_of_ALP;j++)
                {
                        array_positive<double>* tmp=((array_positive<double>*)alp_distr[j]);
                        array_positive<double>* tmp_errors=((array_positive<double>*)alp_distr_errors[j]);

                        P[j]=0;
                        P_errors[j]=0;
                        Int4 i;
                        for(i=0;i<=tmp->d_dim;i++)
                        {
                                P[j]+=tmp->d_elem[i];
                                P_errors[j]+=tmp_errors->d_elem[i];
                        };

                        P_errors[j]=alp_reg::sqrt_for_errors(P_errors[j]);
                };

                

                values_P_ratio=new double[total_number_of_ALP];
                alp_data::assert_mem(values_P_ratio);
                errors_P_ratio=new double[total_number_of_ALP];
                alp_data::assert_mem(errors_P_ratio);

                

                for(j=0;j<total_number_of_ALP;j++)
                {
                        values_P_ratio[j]=P[j+1]/P[j];
                        errors_P_ratio[j]=alp_reg::error_of_the_ratio(P[j+1],P_errors[j+1],P[j],P_errors[j]);
                };



                double beta1=0;
                double beta1_error=0;

                Int4 number_of_elements=total_number_of_ALP;

                bool cut_left_tail=true;
                bool cut_right_tail=false;

                double y=2;

                Int4 k1_opt;
                Int4 k2_opt;


                double P_beta_inf;
                double P_beta_inf_error=0;

                bool res_was_calculated;
                
                alp_reg::robust_regression_sum_with_cut_LSM_beta1_is_defined(
                0,
                number_of_elements-starting_point,
                values_P_ratio+starting_point,
                errors_P_ratio+starting_point,
                cut_left_tail,
                cut_right_tail,
                y,
                P_beta_inf,
                beta1,
                P_beta_inf_error,
                beta1_error,
                k1_opt,
                k2_opt,
                res_was_calculated);


                

                if(!res_was_calculated)
                {
                        throw error("The program cannot estimate the parameters; please repeat the calculation\n",2);
                };

                P_beta_inf=1-P_beta_inf;

                
                //2)E(exp(lambda*T_beta)) and E(T_beta*exp(lambda*T_beta))
                E=new double[total_number_of_ALP+1];
                alp_data::assert_mem(E);
                E_errors=new double[total_number_of_ALP+1];
                alp_data::assert_mem(E_errors);

                E_T_beta=new double[total_number_of_ALP+1];
                alp_data::assert_mem(E_T_beta);
                E_T_beta_errors=new double[total_number_of_ALP+1];
                alp_data::assert_mem(E_T_beta);


                E[0]=1;
                E_T_beta[0]=0;

                E_errors[0]=0;
                E_T_beta_errors[0]=0;

                


                for(j=1;j<=total_number_of_ALP;j++)
                {
                        array_positive<double>* tmp=((array_positive<double>*)alp_distr[j]);
                        array_positive<double>* tmp_errors=((array_positive<double>*)alp_distr_errors[j]);

                        E[j]=0;
                        E_T_beta[j]=0;

                        E_errors[j]=0;
                        E_T_beta_errors[j]=0;

                        Int4 i;
                        for(i=0;i<=tmp->d_dim;i++)
                        {
                                double tmp_double=exp(lambda_*(double)i);
                                E[j]+=tmp_double*tmp->d_elem[i];
                                E_errors[j]+=tmp_double*tmp_double*tmp_errors->d_elem[i];

                                tmp_double=(double)i*exp(lambda_*(double)i);
                                E_T_beta[j]+=tmp_double*tmp->d_elem[i];
                                E_T_beta_errors[j]+=tmp_double*tmp_double*tmp_errors->d_elem[i];
                        };

                        E_errors[j]=alp_reg::sqrt_for_errors(E_errors[j]);
                        E_T_beta_errors[j]=alp_reg::sqrt_for_errors(E_T_beta_errors[j]);

                };


                double E_aver;
                double E_aver_error;

                double E_T_beta_diff_aver;
                double E_T_beta_diff_aver_error;


                if(total_number_of_ALP==1)
                {
                        E_aver=E[1];
                        E_aver_error=E_errors[1];

                        E_T_beta_diff_aver=E_T_beta[1]-E_T_beta[0];
                        E_T_beta_diff_aver_error=E_T_beta_errors[1];

                }
                else
                {
                        Int4 number_of_elements=total_number_of_ALP;

                        bool cut_left_tail=true;
                        bool cut_right_tail=false;


                        double beta0;
                        double beta1=0;
                        double beta0_error;
                        double beta1_error=0;

                        bool res_was_calculated;

                        alp_reg::robust_regression_sum_with_cut_LSM_beta1_is_defined(
                        0,
                        number_of_elements-starting_point,
                        E+1+starting_point,
                        E_errors+1+starting_point,
                        cut_left_tail,
                        cut_right_tail,
                        y,
                        E_aver,
                        beta1,
                        E_aver_error,
                        beta1_error,
                        k1_opt,
                        k2_opt,
                        res_was_calculated);


                        if(!res_was_calculated)
                        {
                                throw error("The program cannot estimate the parameters; please repeat the calculation\n",2);
                        };


                        number_of_elements=total_number_of_ALP;


                        alp_reg::robust_regression_sum_with_cut_LSM(
                        0,
                        number_of_elements-starting_point,
                        E_T_beta+1+starting_point,
                        E_T_beta_errors+1+starting_point,
                        cut_left_tail,
                        cut_right_tail,
                        y,
                        beta0,
                        beta1,
                        beta0_error,
                        beta1_error,
                        k1_opt,
                        k2_opt,
                        res_was_calculated);

                        

                        if(!res_was_calculated)
                        {
                                throw error("The program cannot estimate the parameters; please repeat the calculation\n",2);
                        };


                        E_T_beta_diff_aver=beta1;
                        E_T_beta_diff_aver_error=beta1_error;

                        
                        
                };

                double exp_lambda_error=exp(-lambda_)*lambda_error_;
                double exp_lambda=(1-exp(-lambda_));

                
                double den_error=alp_reg::error_of_the_product(E_T_beta_diff_aver,E_T_beta_diff_aver_error,exp_lambda,exp_lambda_error);
                double den=(1-exp(-lambda_))*E_T_beta_diff_aver;

                double E_aver_sqr_error=alp_reg::error_of_the_product(E_aver,E_aver_error,E_aver,E_aver_error);
                double E_aver_sqr=E_aver*E_aver;

                double nom_error=alp_reg::error_of_the_product(P_beta_inf,P_beta_inf_error,E_aver_sqr,E_aver_sqr_error);
                double nom=P_beta_inf*E_aver_sqr;

                C_error_=alp_reg::error_of_the_ratio(nom,nom_error,den,den_error);
                C_=nom/den;




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

        delete[]values_P_ratio;values_P_ratio=NULL;
        delete[]errors_P_ratio;errors_P_ratio=NULL;


        delete[]P;P=NULL;
        delete[]P_errors;P_errors=NULL;

        delete[]E;E=NULL;
        delete[]E_T_beta;E_T_beta=NULL;
        delete[]E_errors;E_errors=NULL;
        delete[]E_T_beta_errors;E_T_beta_errors=NULL;

        if(ee_error_flag)
        {
                throw error(ee_error.st,ee_error.error_code);
        };

};

void alp_sim::get_and_allocate_alp_distribution(
Int4 ind1_,
Int4 ind2_,
void **&alp_distr,
void **&alp_distr_errors,
Int4 nalp)
{
        if(nalp<=0)
        {
                if(nalp<0)
                {
                        throw error("Unexpected error\n",4);
                };

                alp_distr=NULL;
                alp_distr_errors=NULL;

                return;
        };

        void **alp_distr_tmp=NULL;
        void **alp_distr_errors_tmp=NULL;

        Int4 allocation_dim=nalp;
        Int4 allocation_dim_tmp=nalp+1;

        bool ee_error_flag=false;
        error ee_error("",0);

        try
        {
        try
        {


                Int4 i;
                alp_distr_tmp=new void*[nalp+1];
                alp_data::assert_mem(alp_distr_tmp);

                alp_distr_errors_tmp=new void*[nalp+1];
                alp_data::assert_mem(alp_distr_errors_tmp);

                for(i=0;i<=nalp;i++)
                {
                        alp_distr_tmp[i]=NULL;
                        alp_distr_errors_tmp[i]=NULL;
                };


                for(i=1;i<=nalp-1;i++)
                {
                        alp_distr_tmp[i]=alp_distr[i];
                        alp_distr_errors_tmp[i]=alp_distr_errors[i];
                };

                delete[]alp_distr;alp_distr=NULL;
                delete[]alp_distr_errors;alp_distr_errors=NULL;

                alp_distr=alp_distr_tmp;alp_distr_tmp=NULL;
                alp_distr_errors=alp_distr_errors_tmp;alp_distr_errors_tmp=NULL;

                allocation_dim=nalp+1;

                

                alp_distr[nalp]=new array_positive<double>(d_alp_data);
                alp_data::assert_mem(alp_distr[nalp]);

                alp_distr_errors[nalp]=new array_positive<double>(d_alp_data);
                alp_data::assert_mem(alp_distr_errors[nalp]);



                for(i=ind1_;i<=ind2_;i++)
                {
                        alp* &alp_obj_tmp=d_alp_obj->d_elem[i];
                        Int4 k=nalp;
                        Int4 &alp_tmp=alp_obj_tmp->d_alp->d_elem[k];
                        double &weight_tmp=alp_obj_tmp->d_alp_weights->d_elem[k];


                        ((array_positive<double>*)alp_distr[k])->increase_elem_by_x(alp_tmp,weight_tmp);
                        ((array_positive<double>*)alp_distr_errors[k])->increase_elem_by_x(alp_tmp,weight_tmp*weight_tmp);
                };

                double ind_diff=(double)(ind2_-ind1_+1);
                Int4 k=nalp;
                array_positive<double>* tmp=((array_positive<double>*)alp_distr[k]);
                array_positive<double>* tmp_errors=((array_positive<double>*)alp_distr_errors[k]);
                Int4 j;
                for(j=0;j<=tmp->d_dim;j++)
                {
                        tmp->d_elem[j]/=ind_diff;
                        tmp_errors->d_elem[j]/=ind_diff;
                        tmp_errors->d_elem[j]-=tmp->d_elem[j]*tmp->d_elem[j];
                        tmp_errors->d_elem[j]/=ind_diff;
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

        //memory release

        if(ee_error_flag)
        {
                Int4 i;
                if(alp_distr_tmp)
                {
                        for(i=0;i<=allocation_dim_tmp;i++)
                        {
                                delete (array_positive<double>*)alp_distr_tmp[i];alp_distr_tmp[i]=NULL;
                        };
                        delete []alp_distr_tmp;alp_distr_tmp=NULL;
                };

                if(alp_distr_errors_tmp)
                {
                        for(i=0;i<=allocation_dim_tmp;i++)
                        {
                                delete (array_positive<double>*)alp_distr_errors_tmp[i];alp_distr_errors_tmp[i]=NULL;
                        };
                        delete []alp_distr_errors_tmp;alp_distr_errors_tmp=NULL;
                };

                if(alp_distr)
                {
                        for(i=0;i<=allocation_dim;i++)
                        {
                                delete (array_positive<double>*)alp_distr[i];alp_distr[i]=NULL;
                        };
                        delete []alp_distr;alp_distr=NULL;
                };

                if(alp_distr_errors)
                {
                        for(i=0;i<=allocation_dim;i++)
                        {
                                delete (array_positive<double>*)alp_distr_errors[i];alp_distr_errors[i]=NULL;
                        };
                        delete []alp_distr_errors;alp_distr_errors=NULL;
                };


                throw error(ee_error.st,ee_error.error_code);
        };

};

bool alp_sim::check_K_criterion(
Int4 nalp_,
Int4 ind1_,
Int4 ind2_,
double lambda_,
double eps_K_,
Int4 &M_min_)
{
        if(nalp_<=0)
        {
                throw error("Unexpected error\n",4);
        };
        bool ee_error_flag=false;
        error ee_error("",0);
        array_positive<double>* diff=NULL;

        try
        {
        try
        {

                diff=new array_positive<double>(d_alp_data);
                alp_data::assert_mem(diff);

                double sum_of_weights=0;
                double M_aver=0;

                Int4 i;
                for(i=ind1_;i<=ind2_;i++)
                {
                        alp* &alp_obj_tmp=d_alp_obj->d_elem[i];
                        Int4 &alp_tmp=alp_obj_tmp->d_alp->d_elem[nalp_];
                        double &weight_tmp=alp_obj_tmp->d_alp_weights->d_elem[nalp_];
                        sum_of_weights+=weight_tmp;
                        M_aver+=alp_tmp*weight_tmp;

                        array<Int4> *cells_counts=alp_obj_tmp->d_cells_counts;

                        Int4 k;
                        for(k=cells_counts->d_ind0;k<=alp_data::Tmin(alp_tmp,cells_counts->d_dim_plus_d_ind0);k++)
                        {
                                diff->increase_elem_by_x(alp_tmp-k,cells_counts->d_elem[k-cells_counts->d_ind0]*weight_tmp);
                        };
                };



                double den=0;
                for(i=0;i<=diff->d_dim;i++)
                {
                        den+=exp(-lambda_*(double)i)*diff->d_elem[i];
                };


                if(den<=0||sum_of_weights<=0)
                {
                        throw error("The program is not able to estimate the parameters\n",2);
                };


                M_aver/=sum_of_weights;


                double delta_val=den*eps_K_*(1-exp(-lambda_));

                Int4 diff_opt=1;;
                for(i=diff->d_dim;i>=0;i--)
                {
                        if(exp(-lambda_*(double)i)*diff->d_elem[i]>delta_val)
                        {
                                diff_opt=i+1;
                                break;
                        };
                };

                


                M_min_=(Int4)alp_data::round(M_aver);


                delete diff;diff=NULL;
                if(M_aver<diff_opt)
                {
                        return false;
                };
                
                return true;
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

        delete diff;diff=NULL;

        if(ee_error_flag)
        {
                throw error(ee_error.st,ee_error.error_code);
        };
        return false;

};

bool alp_sim::check_K_criterion_during_killing(
Int4 ind1_,
Int4 ind2_,
double lambda_,
double eps_K_,
Int4 current_level_,
Int4 &recommended_level_,
Int4 &diff_opt_,
double &K_C_,
double &K_C_error_)
{
        if(ind1_>ind2_)
        {
                throw error("Unexpected error\n",4);
        };

        array_positive<double>* diff=NULL;
        array_positive<double>* diff_error=NULL;

        bool ee_error_flag=false;
        error ee_error("",0);

        try
        {
        try
        {

                diff=new array_positive<double>(d_alp_data);
                alp_data::assert_mem(diff);

                diff_error=new array_positive<double>(d_alp_data);
                alp_data::assert_mem(diff_error);


                double sum_of_weights=0;
                double sum_of_weights_error=0;

                double M_aver=0;

                Int4 i;
                for(i=ind1_;i<=ind2_;i++)
                {
                        alp* &alp_obj_tmp=d_alp_obj->d_elem[i];
                        Int4 &alp_tmp=alp_obj_tmp->d_M;
                        double &weight_tmp=alp_obj_tmp->d_alp_weights->d_elem[alp_obj_tmp->d_nalp_killing];
                        sum_of_weights+=weight_tmp;
                        sum_of_weights_error+=weight_tmp*weight_tmp;
                        M_aver+=alp_tmp*weight_tmp;


                        array<Int4> *cells_counts=alp_obj_tmp->d_cells_counts;

                        Int4 k;
                        for(k=cells_counts->d_ind0;k<=alp_data::Tmin(alp_tmp,cells_counts->d_dim_plus_d_ind0);k++)
                        {
                                double tmp=cells_counts->d_elem[k-cells_counts->d_ind0]*weight_tmp;
                                diff->increase_elem_by_x(alp_tmp-k,tmp);
                                diff_error->increase_elem_by_x(alp_tmp-k,tmp*tmp);
                        };
                };



                double tmp2=(double)(ind2_-ind1_+1);

                sum_of_weights/=tmp2;
                sum_of_weights_error/=tmp2;
                sum_of_weights_error-=sum_of_weights*sum_of_weights;
                sum_of_weights_error/=tmp2;
                sum_of_weights_error=alp_reg::sqrt_for_errors(sum_of_weights_error);


                
                for(i=0;i<=diff->d_dim;i++)
                {
                        diff->d_elem[i]/=tmp2;
                        diff_error->d_elem[i]/=tmp2;
                        diff_error->d_elem[i]-=diff->d_elem[i]*diff->d_elem[i];
                        diff_error->d_elem[i]/=tmp2;
                };



                double den=0;
                double den_error=0;
                for(i=0;i<=diff->d_dim;i++)
                {
                        double tmp=exp(-lambda_*(double)i);
                        den+=tmp*diff->d_elem[i];
                        den_error+=tmp*tmp*diff_error->d_elem[i];

                };



                den_error=alp_reg::sqrt_for_errors(den_error);


                if(den<=0||sum_of_weights<=0)
                {
                        throw error("The program is not able to estimate the parameters\n",2);
                };

                K_C_=sum_of_weights/den;
                K_C_error_=alp_reg::error_of_the_ratio(sum_of_weights,sum_of_weights_error,den,den_error);


                M_aver/=tmp2;
                M_aver/=sum_of_weights;


                double delta_val=den*eps_K_*(1-exp(-lambda_));

                Int4 diff_opt=1;;
                for(i=diff->d_dim;i>=0;i--)
                {
                        if(exp(-lambda_*(double)i)*diff->d_elem[i]>delta_val)
                        {
                                diff_opt=i+1;
                                break;
                        };
                };

                delete diff;diff=NULL;
                delete diff_error;diff_error=NULL;

                if(M_aver-diff_opt<current_level_)
                {
                        recommended_level_=(Int4)floor(M_aver-diff_opt*1.1);
                        diff_opt_=(Int4)ceil(M_aver-recommended_level_);
                        return false;
                };
                recommended_level_=current_level_;
                diff_opt_=(Int4)ceil(M_aver-recommended_level_);
                return true;
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
                delete diff;diff=NULL;
                delete diff_error;diff_error=NULL;

                throw error(ee_error.st,ee_error.error_code);
        };
        return false;

};



void alp_sim::calculate_lambda(
bool check_the_criteria_,
Int4 nalp_,
Int4 &nalp_thr_,
bool &inside_simulation_flag_,
Int4 ind1_,
Int4 ind2_,
void **alp_distr,
void **alp_distr_errors,
double &lambda_,
double &lambda_error_,
double &test_difference_,
double &test_difference_error_)
{
        Int4 nalp=nalp_;

        if(nalp<=0)
        {
                throw error("Unexpected error\n",4);
        };

        

        struct_for_lambda_calculation tmp_struct;
        tmp_struct.d_alp_distr=alp_distr;
        tmp_struct.d_alp_distr_errors=alp_distr_errors;
        tmp_struct.d_nalp=nalp;
        tmp_struct.d_calculate_alp_number=false;


        function_type *func=function_for_lambda_calculation;
        void* func_pointer=&tmp_struct;
        double a=0;
        double b=d_alp_data->d_is->d_lambda*3;
        Int4 n_partition=30;
        double eps=1e-10;
        std::vector<double> res;



        alp_reg::find_tetta_general(
        func,
        func_pointer,
        a,//[a,b] is the interval for search of equation solution
        b,
        n_partition,
        eps,
        res);

        


        inside_simulation_flag_=true;
        if(res.size()==0)
        {
                inside_simulation_flag_=false;
                return;
        };

        

        lambda_=get_root(res,d_alp_data->d_is->d_lambda);
        

        tmp_struct.d_calculate_alp_number=true;
        double f1=func(lambda_,func_pointer);
        nalp_thr_=tmp_struct.d_alp_number;
        tmp_struct.d_calculate_alp_number=false;

        double slope_error=tmp_struct.d_f_error;


        double sum1=tmp_struct.d_last_sum;
        double sum1_error=tmp_struct.d_last_sum_error;

        double delta_lambda=lambda_/100.0;
        double f2=func(lambda_+delta_lambda,func_pointer);
        
        

        if(delta_lambda==0||f1==f2)
        {
                lambda_error_=0.0;
        }
        else
        {
                double derivative=(f2-f1)/delta_lambda;
                lambda_error_=fabs(slope_error/derivative);
        };
        
        if(!check_the_criteria_)
        {
                return;
        };

        if(nalp>1)
        {
                func(d_lambda_tmp->d_elem[nalp-1],func_pointer);
        }
        else
        {
                func(d_alp_data->d_is->d_ungap_lambda,func_pointer);
        };

        
        
        double sum2=tmp_struct.d_last_sum;
        double sum2_error=tmp_struct.d_last_sum_error;

        double max_sum=alp_data::Tmax(fabs(sum1),fabs(sum2));

        if(max_sum!=0)
        {
                test_difference_=fabs((sum1-sum2)/max_sum);
                test_difference_error_=0.5*(sum1_error+sum2_error)/max_sum;
        }
        else
        {
                test_difference_=-1;
                test_difference_error_=0;
        };




        


};

double alp_sim::get_root(
const std::vector<double> &res_tmp_,
double point_)
{
        if(res_tmp_.size()==0)
        {
                throw error("Error in alp_sim::get_root - the equation does not have roots\n",2);
                return 0;
        };

        Int4 i;
        Int4 p=0;
        double d1=fabs(point_-res_tmp_[0]);
        for(i=1;i<(Int4)res_tmp_.size();i++)
        {
                double d2=fabs(point_-res_tmp_[i]);
                if(d2<d1)
                {
                        p=i;
                        d1=d2;
                };
        };

        return res_tmp_[p];
};



double alp_sim::function_for_lambda_calculation(
double lambda_,
void * data_)
{

        double *expect=NULL;
        double *expect_errors=NULL;

        bool ee_error_flag=false;
        error ee_error("",0);

        try
        {
        try
        {

                struct_for_lambda_calculation *tmp_struct=(struct_for_lambda_calculation *)data_;
                void **alp_distr=tmp_struct->d_alp_distr;
                void **alp_distr_errors=tmp_struct->d_alp_distr_errors;
                Int4 nalp=tmp_struct->d_nalp;

                expect=new double[nalp];
                alp_data::assert_mem(expect);
                expect_errors=new double[nalp];
                alp_data::assert_mem(expect_errors);

                if(nalp<1)
                {
                        throw error("Unexpected error\n",4);
                };



                Int4 k;
                for(k=1;k<=nalp;k++)
                {
                        array_positive<double>* tmp=((array_positive<double>*)alp_distr[k]);
                        array_positive<double>* tmp_errors=((array_positive<double>*)alp_distr_errors[k]);

                        double val=0;
                        double val_error=0;

                        Int4 j;
                        for(j=0;j<=tmp->d_dim;j++)
                        {
                                double exp_tmp=exp(lambda_*j);
                                // MCF: added "if"s to avoid overflow
                                // errors if exp_tmp is huge.
                                if(tmp->d_elem[j]>0)
                                  val+=exp_tmp*tmp->d_elem[j];
                                if(tmp_errors->d_elem[j]>0)
                                  val_error+=exp_tmp*exp_tmp*tmp_errors->d_elem[j];
                                /*
                                  std::cerr << "! " << k << " " << nalp << " "
                                  << j << " " << tmp->d_dim << " "
                                  << lambda_ << " " << exp_tmp << " "
                                  << tmp->d_elem[j] << " "
                                  << tmp_errors->d_elem[j] << " "
                                  << val << " " << val_error << "\n";
                                */
                        };
                        val_error=alp_reg::sqrt_for_errors(val_error);
                        expect[k-1]=val;
                        expect_errors[k-1]=val_error;


                };

                tmp_struct->d_last_sum=expect[nalp-1];
                tmp_struct->d_last_sum_error=expect_errors[nalp-1];

                if(tmp_struct->d_calculate_alp_number)
                {
                        double tmp=0.0;
                        Int4 k;
                        for(k=0;k<nalp;k++)
                        {
                                if(expect_errors[k]!=0)
                                {
                                        tmp+=1.0/(expect_errors[k]*expect_errors[k]);
                                };

                        };

                        Int4 tmp_alp=nalp;
                        double tmp1=0.0;
                        for(k=nalp-1;k>=0;k--)
                        {
                                if(expect_errors[k]!=0)
                                {
                                        tmp1+=1.0/(expect_errors[k]*expect_errors[k]);
                                };
                                if(tmp1>0.2*tmp)
                                {
                                        tmp_alp=k+1;
                                        break;
                                };
                        };

                        tmp_struct->d_alp_number=tmp_alp;
                };

                if(nalp==1)
                {
                        double tmp=expect[0]-1.0;
                        tmp_struct->d_f_error=expect_errors[0];
                        delete[]expect;expect=NULL;
                        delete[]expect_errors;expect_errors=NULL;
                        return tmp;
                };


                Int4 min_length=0;
                Int4 number_of_elements=nalp;
                bool cut_left_tail=true;
                bool cut_right_tail=false;
                double y=2;
                double beta0;
                double beta1;
                double beta0_error;
                double beta1_error;
                Int4 k1_opt;
                Int4 k2_opt;

                bool res_was_calculated;

                alp_reg::robust_regression_sum_with_cut_LSM(
                min_length,
                number_of_elements,
                expect,
                expect_errors,
                cut_left_tail,
                cut_right_tail,
                y,
                beta0,
                beta1,
                beta0_error,
                beta1_error,
                k1_opt,
                k2_opt,
                res_was_calculated);

                if(!res_was_calculated)
                {
                        throw error("The program cannot estimate the parameters; please repeat the calculation\n",2);
                };


                delete[]expect;expect=NULL;
                delete[]expect_errors;expect_errors=NULL;

                tmp_struct->d_f_error=beta1_error;
                return beta1;

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
                delete[]expect;expect=NULL;
                delete[]expect_errors;expect_errors=NULL;

                throw error(ee_error.st,ee_error.error_code);
        };

        return 0.0;

};

