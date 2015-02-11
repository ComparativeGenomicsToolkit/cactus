/* $Id: sls_pvalues.cpp 189337 2010-04-21 13:14:53Z boratyng $
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

File name: sls_pvalues.cpp

Author: Sergey Sheetlin

Contents: Calculation of P-values using precalculated Gumbel parameters

******************************************************************************/

#include <ncbi_pch.hpp>

#include <cassert>

#include "sls_pvalues.hpp"
#include "sls_alp_data.hpp"


USING_NCBI_SCOPE;
USING_SCOPE(blast);
USING_SCOPE(Sls);

const bool read_Sbs_par_flag=true;



double pvalues::error_of_the_sum(//v1_+v2_
double v1_,
double v1_error_,
double v2_,
double v2_error_)
{
        if(v1_error_>=1e100||v2_error_>=1e100)
        {
                return 1e100;
        };

        return sqrt(v1_error_*v1_error_+v2_error_*v2_error_);
};

double pvalues::error_of_the_product(//v1_*v2_
double v1_,
double v1_error_,
double v2_,
double v2_error_)
{
        if(v1_error_>=1e100||v2_error_>=1e100)
        {
                return 1e100;
        };

        double a1=(v1_+v1_error_)*(v2_+v2_error_);
        double a2=(v1_-v1_error_)*(v2_+v2_error_);
        double a3=(v1_+v1_error_)*(v2_-v2_error_);
        double a4=(v1_-v1_error_)*(v2_-v2_error_);

        double a=v1_*v2_;

        return alp_data::Tmax(fabs(a1-a),fabs(a2-a),fabs(a3-a),fabs(a4-a));

};


double pvalues::error_of_the_sqrt(//sqrt(v1_)
double v1_,
double v1_error_)
{
        if(v1_error_>=1e100||v1_<0)
        {
                return 1e100;
        };

        double s=sqrt(v1_);
        double s1=sqrt(alp_data::Tmax(0.0,v1_-v1_error_));
        double s2=sqrt(alp_data::Tmax(0.0,v1_+v1_error_));

        return alp_data::Tmax(fabs(s-s1),fabs(s-s2));
};

double pvalues::error_of_the_ratio(//v1_/v2_
double v1_,
double v1_error_,
double v2_,
double v2_error_)
{
        if(v1_error_>=1e100||v2_error_>=1e100)
        {
                return 1e100;
        };


        if(v2_==0)
        {
                return 1e100;
        };

        if(v1_==0&&v1_error_==0)
        {
                return 0.0;
        };

        double a=v1_/v2_;


        if(((v2_+v2_error_)*v2_<=0))
        {
                double a3=(v1_+v1_error_)/(v2_-v2_error_);
                double a4=(v1_-v1_error_)/(v2_-v2_error_);
                return alp_data::Tmax(fabs(a-a3),fabs(a-a4));
        };

        if(((v2_-v2_error_)*v2_<=0))
        {
                double a1=(v1_+v1_error_)/(v2_+v2_error_);
                double a2=(v1_-v1_error_)/(v2_+v2_error_);
                return alp_data::Tmax(fabs(a-a1),fabs(a-a2));
        };


        double a1=(v1_+v1_error_)/(v2_+v2_error_);
        double a2=(v1_-v1_error_)/(v2_+v2_error_);
        double a3=(v1_+v1_error_)/(v2_-v2_error_);
        double a4=(v1_-v1_error_)/(v2_-v2_error_);

        return alp_data::Tmax(fabs(a-a1),fabs(a-a2),fabs(a-a3),fabs(a-a4));
};




double pvalues::one_minus_exp_function(
double y_)
{
        if(fabs(y_)>1e-8)
        {
                return 1.0-exp(y_);
        }
        else
        {
                return -(y_+y_*y_/2.0+y_*y_*y_/6.0+y_*y_*y_*y_/24.0);
        };
};

double pvalues::normal_probability(
double x_,
double eps_)
{

        double pi=3.1415926535897932384626433832795;
        if(x_==0)
        {
                return 0.5;
        };


        eps_=alp_data::Tmin(1.0,eps_);

        double x_max=10*eps_+sqrt(alp_data::Tmax(0.0,-2*log(eps_)));


        if(x_>=x_max)
        {
                double x=x_/sqrt(2.0);
                return 1-0.5*exp(-x*x)/(x*sqrt(pi))*(1-1.0/(2*x*2*x));
        };

        if(x_<=-x_max)
        {
                double x=x_/sqrt(2.0);
                return 0.5*exp(-x*x)/(-x*sqrt(pi))*(1-1.0/(2*x*2*x));
        };


        double const_val=1/sqrt(2.0*pi);

        


        Int4 N=(Int4)alp_data::round(fabs(x_)/eps_)+1;
        double h=x_/(double)N;



        double res=0;
        Int4 i;
        for(i=0;i<=N;i++)
        {
                double y=h*i;
                double tmp=exp(-0.5*y*y);
                if(i==0||i==N)
                {
                        res+=0.5*tmp;
                }
                else
                {
                        res+=tmp;
                };
        };

        res*=h;

        return 0.5+const_val*(res);
};

double pvalues::normal_probability(
double a_,
double b_,
double h_,
Int4 N_,
double *p_,
double x_,
double eps_)
{
        if(x_<a_||x_>b_)
        {
                return normal_probability(x_,eps_);
        };

        Int4 x_n=(Int4)floor((x_-a_)/h_);
        x_n=alp_data::Tmin(N_-1,x_n);
        assert(x_n >= 0);
        return p_[x_n]+(p_[x_n+1]-p_[x_n])*(x_-(h_*x_n+a_))/h_;
};

double pvalues::ln_one_minus_val(
double val_)
{
        if(val_>1e-8)
        {
                return log(1-val_);
        };

        return -val_-val_*val_/2.0-val_*val_*val_/3.0;
};



void pvalues::get_appr_tail_prob_with_cov(
set_of_parameters &par_,
bool blast_,
double y_,
double m_,
double n_,

double &P_,
double &P_error_,

double &area_,

double a_normal_,
double b_normal_,
double h_normal_,
Int4 N_normal_,
double *p_normal_,

bool &area_is_1_flag_)
{


        double lambda_=par_.lambda;
        double lambda_error_=par_.lambda_error;
        double k_=par_.K;
        double k_error_=par_.K_error;

        double ai_hat_=par_.a_I;
        double ai_hat_error_=par_.a_I_error;
        double bi_hat_=2.0*par_.G*(par_.gapless_a-par_.a_I); 
        double bi_hat_error_=2.0*par_.G*error_of_the_sum(par_.gapless_a,par_.gapless_a_error,par_.a_I,par_.a_I_error); 
        double alphai_hat_=par_.alpha_I;
        double alphai_hat_error_=par_.alpha_I_error;
        double betai_hat_=2.0*par_.G*(par_.gapless_alpha-par_.alpha_I); 
        double betai_hat_error_=2.0*par_.G*error_of_the_sum(par_.gapless_alpha,par_.gapless_alpha_error,par_.alpha_I,par_.alpha_I_error); 

        double aj_hat_=par_.a_J;
        double aj_hat_error_=par_.a_J_error;
        double bj_hat_=2.0*par_.G*(par_.gapless_a-par_.a_J); 
        double bj_hat_error_=2.0*par_.G*error_of_the_sum(par_.gapless_a,par_.gapless_a_error,par_.a_J,par_.a_J_error); 
        double alphaj_hat_=par_.alpha_J;
        double alphaj_hat_error_=par_.alpha_J_error;
        double betaj_hat_=2.0*par_.G*(par_.gapless_alpha-par_.alpha_J); 
        double betaj_hat_error_=2.0*par_.G*error_of_the_sum(par_.gapless_alpha,par_.gapless_alpha_error,par_.alpha_J,par_.alpha_J_error); 

        double sigma_hat_=par_.sigma;
        double sigma_hat_error_=par_.sigma_error;
        double tau_hat_=2.0*par_.G*(par_.gapless_alpha-par_.sigma);
        double tau_hat_error_=2.0*par_.G*error_of_the_sum(par_.gapless_alpha,par_.gapless_alpha_error,par_.sigma,par_.sigma_error);
 

        bool where_it_is_works_flag=false;
        bool negative_flag=false;


        if(blast_)
        {
                alphai_hat_=0;
                alphai_hat_error_=0;
                betai_hat_=0;
                betai_hat_error_=0;

                alphaj_hat_=0;
                alphaj_hat_error_=0;
                betaj_hat_=0;
                betaj_hat_error_=0;

                sigma_hat_=0;
                sigma_hat_error_=0;
                tau_hat_=0;
                tau_hat_error_=0;
        };

        double const_val=1/sqrt(2.0*3.1415926535897932384626433832795);
        double eps=0.000001;



        double m_li_y_error;
        double m_li_y;

        bool flag_Mi=true;

        if(flag_Mi||blast_)
        {
                double tmp=alp_data::Tmax(0.0,(ai_hat_*y_+bi_hat_));
                if(where_it_is_works_flag)
                {
                        if(ai_hat_*y_+bi_hat_<0)
                        {
                                negative_flag=true;
                        };
                };
                m_li_y_error=error_of_the_sum(y_*ai_hat_,fabs(y_)*ai_hat_error_,bi_hat_,bi_hat_error_);
                m_li_y=m_-tmp;
        };
        
        double vi_y_error;
        double vi_y;

        bool flag_Mii=true;

        if(flag_Mii||blast_)
        {
                vi_y_error=error_of_the_sum(y_*alphai_hat_,fabs(y_)*alphai_hat_error_,betai_hat_,betai_hat_error_);
                vi_y=alp_data::Tmax(0.0,alphai_hat_*y_+betai_hat_);
                if(where_it_is_works_flag)
                {
                        if(alphai_hat_*y_+betai_hat_<0)
                        {
                                negative_flag=true;
                        };
                };
        };

        double sqrt_vi_y_error=error_of_the_sqrt(vi_y,vi_y_error);
        double sqrt_vi_y=sqrt(vi_y);

        double m_F;
        double m_F_error;

        if(sqrt_vi_y==0.0||blast_)
        {
                m_F=1e100;
                m_F_error=0.0;
        }
        else
        {
                m_F_error=error_of_the_ratio(m_li_y,m_li_y,sqrt_vi_y,sqrt_vi_y_error);
                m_F=m_li_y/sqrt_vi_y;
        };


        double P_m_F=normal_probability(a_normal_,b_normal_,h_normal_,N_normal_,p_normal_,m_F,eps);
        double P_m_F_error=const_val*exp(-0.5*m_F*m_F)*m_F_error;

        double E_m_F=-const_val*exp(-0.5*m_F*m_F);
        double E_m_F_error=fabs(-E_m_F*m_F)*m_F_error;

        double m_li_y_P_m_F_error=error_of_the_product(m_li_y,m_li_y_error,P_m_F,P_m_F_error);
        double m_li_y_P_m_F=m_li_y*P_m_F;

        double sqrt_vi_y_E_m_F_error=error_of_the_product(sqrt_vi_y,sqrt_vi_y_error,E_m_F,E_m_F_error);
        double sqrt_vi_y_E_m_F=sqrt_vi_y*E_m_F;

        double p1_error=error_of_the_sum(m_li_y_P_m_F,m_li_y_P_m_F_error,sqrt_vi_y_E_m_F,sqrt_vi_y_E_m_F_error);
        double p1=m_li_y_P_m_F-sqrt_vi_y_E_m_F;


        double n_lj_y_error;
        double n_lj_y;

        bool flag_Mj=true;

        if(flag_Mj||blast_)
        {
                double tmp=alp_data::Tmax(0.0,(aj_hat_*y_+bj_hat_));
                n_lj_y_error=error_of_the_sum(y_*aj_hat_,fabs(y_)*aj_hat_error_,bj_hat_,bj_hat_error_);
                n_lj_y=n_-tmp;

                if(where_it_is_works_flag)
                {
                        if(aj_hat_*y_+bj_hat_<0)
                        {
                                negative_flag=true;
                        };
                };

        };

        double vj_y_error;
        double vj_y;

        bool flag_Mjj=true;

        if(flag_Mjj||blast_)
        {
                vj_y_error=error_of_the_sum(y_*alphaj_hat_,fabs(y_)*alphaj_hat_error_,betaj_hat_,betaj_hat_error_);
                vj_y=alp_data::Tmax(0.0,alphaj_hat_*y_+betaj_hat_);

                if(where_it_is_works_flag)
                {
                        if(alphaj_hat_*y_+betaj_hat_<0)
                        {
                                negative_flag=true;
                        };
                };

        };

        

        double sqrt_vj_y_error=error_of_the_sqrt(vj_y,vj_y_error);
        double sqrt_vj_y=sqrt(vj_y);

        double n_F;
        double n_F_error;

        if(sqrt_vj_y==0.0||blast_)
        {
                n_F=1e100;
                n_F_error=0.0;
        }
        else
        {
                n_F_error=error_of_the_ratio(n_lj_y,n_lj_y,sqrt_vj_y,sqrt_vj_y_error);
                n_F=n_lj_y/sqrt_vj_y;
        };

        double P_n_F=normal_probability(a_normal_,b_normal_,h_normal_,N_normal_,p_normal_,n_F,eps);
        double P_n_F_error=const_val*exp(-0.5*n_F*n_F)*n_F_error;

        double E_n_F=-const_val*exp(-0.5*n_F*n_F);
        double E_n_F_error=fabs(-E_n_F*n_F)*n_F_error;

        double n_lj_y_P_n_F_error=error_of_the_product(n_lj_y,n_lj_y_error,P_n_F,P_n_F_error);
        double n_lj_y_P_n_F=n_lj_y*P_n_F;

        double sqrt_vj_y_E_n_F_error=error_of_the_product(sqrt_vj_y,sqrt_vj_y_error,E_n_F,E_n_F_error);
        double sqrt_vj_y_E_n_F=sqrt_vj_y*E_n_F;

        double p2_error=error_of_the_sum(n_lj_y_P_n_F,n_lj_y_P_n_F_error,sqrt_vj_y_E_n_F,sqrt_vj_y_E_n_F_error);
        double p2=n_lj_y_P_n_F-sqrt_vj_y_E_n_F;




        double c_y_error;
        double c_y;

        bool flag_Mij=true;

        if(flag_Mij||blast_)
        {
                c_y_error=error_of_the_sum(sigma_hat_*y_,sigma_hat_error_*y_,tau_hat_,tau_hat_error_);
                c_y=alp_data::Tmax(0.0,sigma_hat_*y_+tau_hat_);

                if(where_it_is_works_flag)
                {
                        if(sigma_hat_*y_+tau_hat_<0)
                        {
                                negative_flag=true;
                        };
                };

        };

        double P_m_F_P_n_F_error=error_of_the_product(P_m_F,P_m_F_error,P_n_F,P_n_F_error);
        double P_m_F_P_n_F=P_m_F*P_n_F;

        double c_y_P_m_F_P_n_F_error=error_of_the_product(c_y,c_y_error,P_m_F_P_n_F,P_m_F_P_n_F_error);
        double c_y_P_m_F_P_n_F=c_y*P_m_F_P_n_F;

        double p1_p2_error=error_of_the_product(p1,p1_error,p2,p2_error);
        double p1_p2=p1*p2;

        p1_p2=alp_data::Tmax(0.0,p1_p2);


        double area_error=error_of_the_sum(p1_p2,p1_p2_error,c_y_P_m_F_P_n_F,c_y_P_m_F_P_n_F_error);
        double area=p1_p2+c_y_P_m_F_P_n_F;




        if(!blast_)
        {
                area=alp_data::Tmax(area,1.0);
        }
        else
        {
                if(area<=1.0)
                {
                        area_is_1_flag_=true;
                };

                if(area_is_1_flag_)
                {
                        area=1.0;
                };
        };

        if(negative_flag&&where_it_is_works_flag)
        {
                area=0;
        };



        double exp_lambda_y_error=fabs(lambda_error_*y_*exp(-lambda_*y_));
        double exp_lambda_y=exp(-lambda_*y_);

        double area_k_error=error_of_the_product(area,area_error,k_,k_error_);
        double area_k=area*k_;

        double area_k_exp_lambda_y_error=error_of_the_product(area_k,area_k_error,exp_lambda_y,exp_lambda_y_error);
        double area_k_exp_lambda_y=-area_k*exp_lambda_y;

        P_error_=exp(area_k_exp_lambda_y)*area_k_exp_lambda_y_error;

        P_=one_minus_exp_function(area_k_exp_lambda_y);
//        P_=1-exp(-k_*area*exp(-lambda_*y_));

        area_=area;

        

        
};



void pvalues::get_appr_tail_prob_with_cov_without_errors(
set_of_parameters &par_,
bool blast_,
double y_,
double m_,
double n_,

double &P_,
double &P_error_,

double &area_,

double a_normal_,
double b_normal_,
double h_normal_,
Int4 N_normal_,
double *p_normal_,

bool &area_is_1_flag_)
{



        double lambda_=par_.lambda;
        double k_=par_.K;

        double ai_hat_=par_.a_I;
        double bi_hat_=2.0*par_.G*(par_.gapless_a-par_.a_I);
        double alphai_hat_=par_.alpha_I;
        double betai_hat_=2.0*par_.G*(par_.gapless_alpha-par_.alpha_I);

        double aj_hat_=par_.a_J;
        double bj_hat_=2.0*par_.G*(par_.gapless_a-par_.a_J);
        double alphaj_hat_=par_.alpha_J;
        double betaj_hat_=2.0*par_.G*(par_.gapless_alpha-par_.alpha_J);

        double sigma_hat_=par_.sigma;
        double tau_hat_=2.0*par_.G*(par_.gapless_alpha-par_.sigma);


        bool where_it_is_works_flag=false;
        bool negative_flag=false;


        if(blast_)
        {
                alphai_hat_=0;
                betai_hat_=0;

                alphaj_hat_=0;
                betaj_hat_=0;

                sigma_hat_=0;
                tau_hat_=0;
        };

        double const_val=1/sqrt(2.0*3.1415926535897932384626433832795);
        double eps=0.000001;


        double m_li_y;

        bool flag_Mi=true;

        if(flag_Mi||blast_)
        {
                double tmp=alp_data::Tmax(0.0,(ai_hat_*y_+bi_hat_));
                if(where_it_is_works_flag)
                {
                        if(ai_hat_*y_+bi_hat_<0)
                        {
                                negative_flag=true;
                        };
                };
                m_li_y=m_-tmp;
        };
        
        double vi_y;

        bool flag_Mii=true;

        if(flag_Mii||blast_)
        {
                vi_y=alp_data::Tmax(0.0,alphai_hat_*y_+betai_hat_);
                if(where_it_is_works_flag)
                {
                        if(alphai_hat_*y_+betai_hat_<0)
                        {
                                negative_flag=true;
                        };
                };
        };

        double sqrt_vi_y=sqrt(vi_y);

        double m_F;

        if(sqrt_vi_y==0.0||blast_)
        {
                m_F=1e100;
        }
        else
        {
                m_F=m_li_y/sqrt_vi_y;
        };


        double P_m_F=normal_probability(a_normal_,b_normal_,h_normal_,N_normal_,p_normal_,m_F,eps);

        double E_m_F=-const_val*exp(-0.5*m_F*m_F);

        double m_li_y_P_m_F=m_li_y*P_m_F;

        double sqrt_vi_y_E_m_F=sqrt_vi_y*E_m_F;

        double p1=m_li_y_P_m_F-sqrt_vi_y_E_m_F;


        double n_lj_y;

        bool flag_Mj=true;

        if(flag_Mj||blast_)
        {
                double tmp=alp_data::Tmax(0.0,(aj_hat_*y_+bj_hat_));
                n_lj_y=n_-tmp;

                if(where_it_is_works_flag)
                {
                        if(aj_hat_*y_+bj_hat_<0)
                        {
                                negative_flag=true;
                        };
                };

        };

        double vj_y;

        bool flag_Mjj=true;

        if(flag_Mjj||blast_)
        {
                vj_y=alp_data::Tmax(0.0,alphaj_hat_*y_+betaj_hat_);

                if(where_it_is_works_flag)
                {
                        if(alphaj_hat_*y_+betaj_hat_<0)
                        {
                                negative_flag=true;
                        };
                };

        };

        

        double sqrt_vj_y=sqrt(vj_y);

        double n_F;

        if(sqrt_vj_y==0.0||blast_)
        {
                n_F=1e100;
        }
        else
        {
                n_F=n_lj_y/sqrt_vj_y;
        };

        double P_n_F=normal_probability(a_normal_,b_normal_,h_normal_,N_normal_,p_normal_,n_F,eps);

        double E_n_F=-const_val*exp(-0.5*n_F*n_F);

        double n_lj_y_P_n_F=n_lj_y*P_n_F;

        double sqrt_vj_y_E_n_F=sqrt_vj_y*E_n_F;

        double p2=n_lj_y_P_n_F-sqrt_vj_y_E_n_F;




        double c_y;

        bool flag_Mij=true;

        if(flag_Mij||blast_)
        {
                c_y=alp_data::Tmax(0.0,sigma_hat_*y_+tau_hat_);

                if(where_it_is_works_flag)
                {
                        if(sigma_hat_*y_+tau_hat_<0)
                        {
                                negative_flag=true;
                        };
                };

        };

        double P_m_F_P_n_F=P_m_F*P_n_F;

        double c_y_P_m_F_P_n_F=c_y*P_m_F_P_n_F;

        double p1_p2=p1*p2;

        p1_p2=alp_data::Tmax(0.0,p1_p2);


        double area=p1_p2+c_y_P_m_F_P_n_F;




        if(!blast_)
        {
                area=alp_data::Tmax(area,1.0);
        }
        else
        {
                if(area<=1.0)
                {
                        area_is_1_flag_=true;
                };

                if(area_is_1_flag_)
                {
                        area=1.0;
                };
        };

        if(negative_flag&&where_it_is_works_flag)
        {
                area=0;
        };



        double exp_lambda_y=exp(-lambda_*y_);

        double area_k=area*k_;

        double area_k_exp_lambda_y=-area_k*exp_lambda_y;

        P_=one_minus_exp_function(area_k_exp_lambda_y);
//        P_=1-exp(-k_*area*exp(-lambda_*y_));

        area_=area;

};

void pvalues::get_P_error_using_splitting_method(
set_of_parameters &par_,
bool blast_,
double y_,
double m_,
double n_,

double &P_,
double &P_error_,

double &area_,

double a_normal_,
double b_normal_,
double h_normal_,
Int4 N_normal_,
double *p_normal_,

bool &area_is_1_flag_)
{
        Int4 dim=par_.m_LambdaSbs.size();
        if(dim==0)
        {
                throw error("Unexpected error in get_P_error_using_splitting_method\n",1);
        };

        P_=0;
        P_error_=0;


        vector<double> P_values(dim);

        Int4 i;
        for(i=0;i<dim;i++)
        {
                set_of_parameters par_tmp;

                par_tmp.a_I=par_.m_AISbs[i];
                par_tmp.a_I_error=0;

                par_tmp.a_J=par_.m_AJSbs[i];
                par_tmp.a_J_error=0;

                

                par_tmp.gapless_a=par_.gapless_a;
                par_tmp.gapless_a_error=par_.gapless_a_error;

                par_tmp.a=0.5*(par_tmp.a_I+par_tmp.a_J);
                par_tmp.a_error=0;


                par_tmp.sigma=par_.m_SigmaSbs[i];
                par_tmp.sigma_error=0;

                par_tmp.gapless_alpha=par_.gapless_alpha;
                par_tmp.gapless_alpha_error=par_.gapless_alpha_error;


                par_tmp.C=par_.m_CSbs[i];
                par_tmp.C_error=0;

                par_tmp.K=par_.m_KSbs[i];
                par_tmp.K_error=0;


                par_tmp.lambda=par_.m_LambdaSbs[i];
                par_tmp.lambda_error=0;

                par_tmp.alpha_I=par_.m_AlphaISbs[i];
                par_tmp.alpha_I_error=0;

                par_tmp.alpha_J=par_.m_AlphaJSbs[i];
                par_tmp.alpha_J_error=0;

                par_tmp.alpha=0.5*(par_tmp.alpha_I+par_tmp.alpha_J);
                par_tmp.alpha_error=0;

                par_tmp.G=par_.G;

                double P_tmp,P_tmp_error,area_tmp;

                get_appr_tail_prob_with_cov_without_errors(
                par_tmp,
                blast_,
                y_,
                m_,
                n_,

                P_tmp,
                P_tmp_error,

                area_tmp,

                a_normal_,
                b_normal_,
                h_normal_,
                N_normal_,
                p_normal_,

                area_is_1_flag_);

                P_values[i]=P_tmp;

                P_+=P_tmp;

        };

        if(dim<=1)
        {
                return;
        };


        if(P_<=0)
        {
                return;
        };

        P_/=(double)dim;

        for(i=0;i<dim;i++)
        {
                double tmp=P_values[i]/P_;
                P_error_+=tmp*tmp;
        };

        P_error_/=(double)dim;
        P_error_-=1;
        

        
        P_error_=P_*alp_reg::sqrt_for_errors(P_error_/(double)dim);

};


pvalues::pvalues()
{
        bool ee_error_flag=false;
        error ee_error("",0);
        
        p_normal=NULL;

        try
        {
        try
        {

                blast=false;
                eps=0.0001;
                a_normal=-10;
                b_normal=10;
                N_normal=NORMAL_DISTR_ARRAY_DIM;
                h_normal=(b_normal-a_normal)/(double)N_normal;
                p_normal=normal_distr_array_for_P_values_calculation;

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
                ee_error=error("Internal error in the program\n",1);
        };

        //memory release

        if(ee_error_flag)
        {
                this->~pvalues();
                throw error(ee_error.st,ee_error.error_code);
        };


};


pvalues::~pvalues()
{
        
};

void pvalues::calculate_P_values(
Int4 Score1,
Int4 Score2,
double Seq1Len,
double Seq2Len,
set_of_parameters &ParametersSet,
vector<double> &P_values,
vector<double> &P_values_errors)
{
        if(Score2<Score1)
        {
                throw error("Error - Score2<Score1\n",2);
        };

        if(Seq1Len<=0||Seq2Len<=0)
        {
                throw error("Error - Seq1Len<=0||Seq2Len<=0\n",2);
        };

        Int4 ymin_=Score1;//calculation of P-values in the range [ymin_,ymax_]
        Int4 ymax_=Score2;
        double m_=Seq1Len;//length of the first sequence
        double n_=Seq2Len;//length of the second sequence
        set_of_parameters &par_=ParametersSet;

        P_values.resize(ymax_-ymin_+1);
        P_values_errors.resize(ymax_-ymin_+1);



        Int4 y;
        for(y=ymin_;y<=ymax_;y++)
        {

                double P;
                double P_error;
                double area;
                bool area_is_1_flag=false;


                if(read_Sbs_par_flag)
                {
                        

                        get_appr_tail_prob_with_cov_without_errors(
                        par_,
                        blast,
                        (double)y,
                        m_,
                        n_,

                        P,
                        P_error,

                        area,
                        a_normal,
                        b_normal,
                        h_normal,
                        N_normal,
                        p_normal,
                        area_is_1_flag);


                        
                        double P_tmp,area_tmp;

                        if(par_.m_LambdaSbs.size()>0)
                        {
                                get_P_error_using_splitting_method(
                                par_,
                                blast,
                                (double)y,
                                m_,
                                n_,

                                P_tmp,
                                P_error,

                                area_tmp,
                                a_normal,
                                b_normal,
                                h_normal,
                                N_normal,
                                p_normal,
                                area_is_1_flag);


                                if(P_tmp>0)
                                {
                                        P_error=P_error/P_tmp*P;
                                };

                                P_values_errors[y-ymin_]=P_error;
                        }
                        else
                        {
                                P_values_errors[y-ymin_]=-DBL_MAX;
                        };

                        
                        
                }
                else
                {
                        get_appr_tail_prob_with_cov(
                        par_,
                        blast,
                        (double)y,
                        m_,
                        n_,

                        P,
                        P_error,

                        area,
                        a_normal,
                        b_normal,
                        h_normal,
                        N_normal,
                        p_normal,
                        area_is_1_flag);

                        P_values_errors[y-ymin_]=P_error;
                };
        

                // MCF: get E-values instead of p-values
                //P_values[y-ymin_]=P;
                P_values[y-ymin_]=area*par_.K*exp(-par_.lambda*y);
                


        };
        
};

