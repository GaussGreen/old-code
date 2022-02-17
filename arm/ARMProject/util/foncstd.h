/*
 * $Log: foncstd.h,v $
 * Revision 1.11  2003/10/21 12:19:04  jpriaudel
 * modif du define de Precision en ARM_Precision
 *
 * Revision 1.10  2003/01/21 15:57:02  mab
 * Formatting
 *
 * Revision 1.9  2002/11/25 16:23:33  mab
 * Formatting
 *
 */

/*-------------------------------------------------------------------------*/
#ifndef _FONCT_STD_H
#define _FONCT_STD_H



#include "armglob.h"
 



#define ARM_Precision                   1.0e-8
#define INFINIT                     1000000
#define alpha_start                       0.5
#define MAX_EXP                     700
#define ECART_RACINE                10000

#ifndef M_PI
#define M_PI                        3.141592654
#endif
 
                                                    
                       
static double Abs[49]=
{
    9.079405,
	8.417337,
	7.868063,
	7.376859,
	6.922853,
	6.495293,
	6.087727,
	5.695900,
	5.316811,
	4.948243,
	4.588493,
	4.236213,
	3.890307,
	3.549870,
	3.214133,
	2.882436,
	2.554201,
	2.228918,
	1.906125,
	1.585402,
	1.266362,
	0.948641,
	0.631893,
	0.315787,
	-0.000000,
	-0.315787,
	-0.631893,
	-0.948641,
	-1.266362,
	-1.585402,
	-1.906125,
	-2.228918,
	-2.554201,
	-2.882436,
	-3.214133,
	-3.549870,
	-3.890307,
	-4.236213,
	-4.588493,
	-4.948243,
	-5.316811,
	-5.695900,
	-6.087727,
	-6.495293,
	-6.922853,
	-7.376859,
	-7.868063,
	-8.417337,
	-9.079405
};
	
	
static double Pds[49]=
{
    0.000000,
    0.000000,
    0.000000,
    0.000000,
    0.000000,
    0.000000,
    0.000000,
    0.000000,
    0.000000,
    0.000000,
    0.000000,
    0.000000,
    0.000000,
    0.000001,
    0.000011,
    0.000081,
    0.000480,
    0.002254,
    0.008502,
    0.025900,
    0.064032,
    0.128966,
    0.212222,
    0.285911,
    0.315734,
    0.285911,
    0.212222,
    0.128966,
    0.064032,
    0.025900,
    0.008502,
    0.002254,
    0.000480,
    0.000081,
    0.000011,
    0.000001,
    0.000000,
    0.000000,
    0.000000,
    0.000000,
    0.000000,
    0.000000,
    0.000000,
    0.000000,
    0.000000,
    0.000000,
    0.000000,
    0.000000,
    0.000000
};
	
// déclaration de parametres generiques pour les S.O.

struct marginal_descr // work for Normal, LogNormal or other one parameter volatility marginal law
{
	int ModelType;
	double volatility;
};


struct sabr_descr:public marginal_descr
{
	// int ModelType;
	double alpha;
	double beta;
	double rho;
	double nu;
	int sabr_flag;
};


struct jointdistrib_descr
{
	int JointType;
};


struct one_parameter_bivariate:public jointdistrib_descr
{
	double correl;
};


struct bi_sabr:public jointdistrib_descr
{
	double rho_s1_s2;
	double rho_s1_vol2;
	double rho_s2_vol1;
	double rho_vol1_vol2;
};


struct bi_sabr_payIndex:public jointdistrib_descr
{
	bi_sabr* joint_s1_s2;
	double rho_s1_s3;
	double rho_s2_s3;
};


extern double generic_spreadoption_price(double S1, double S2,
                    marginal_descr *m1, marginal_descr *m2, jointdistrib_descr *j,
					double Rate, double K, double Mat, int pc, double w1, double w2,
					int ModelType, int ComputedFormula = 1, int NumSteps = 100);

extern double generic_DigitalSO_price(double S1, double S2,
                    marginal_descr *m1, marginal_descr *m2, jointdistrib_descr *j,
					double K, double Mat, int pc, double w1, double w2,
					int ModelType, int ComputedFormula = 1, int NumSteps = 100);
       
extern double generic_DigitalSO_PayingIndex_price(double S1, double S2, double S3,
                    marginal_descr *m1, marginal_descr *m2, marginal_descr *m3, jointdistrib_descr *j,
					double K, double Mat, int pc, double w1, double w2,
					int ModelType, int ComputedFormula = 1, int NumSteps = 100);
                
                 
// fonctions de calcul d'une double integrale (integration numerique)
extern double fonct_std_d0(double S1, double S2, double vol1, double vol2, 
                    double d1, double d2, double rho, double r, 
                    double K, double t,double x);

extern double fonct_std_d2(double S1, double S2, double vol1, 
                    double vol2, double d1, double d2, 
                    double rho, double r, double K, double t,double x);

extern double fonct_std_d1(double S1, double S2, double vol1,
                    double vol2, double d1, double d2, 
                    double rho, double r, double K,
                    double t,double x);

extern double Proba0(double S1, double S2, double vol1, double vol2,
              double d1, double d2, double rho, double r,
              double K, double t);                                                    
extern double Proba1(double S1, double S2, double vol1, double vol2,
              double d1, double d2, double rho, double r,
              double K, double t);                                                

extern double Proba2(double S1, double S2, double vol1, double vol2,
              double d1, double d2, double rho,
              double r, double K, double t);                                                 

// fonctions de calcul des racines           
extern double condition(double Var, double A1, double A2, double B1,
                 double B2, double  K);

extern double derivee_condition(double Var, double A1,double A2,
                         double B1, double B2);

extern double Newton(double Pt_depart, double epsilon, double A1, 
              double A2,double B1,double B2, double B);

extern double Dicho_croissante(double xmax, double epsilon, double A1,
                        double A2, double B1, double B2,double B);

extern double Dicho_decroissante(double xmax, double epsilon, double A1,
                          double A2, double B1, double B2,double B);

extern double maximum(double x, double y);

extern double depart(double var, double A1, double A2, double B1, double B2,
              double B);


// fonctions de calcul des probabilités d'exercice dans le cas ou
// rho vaut 1 en valeur absolue

extern double Probamin(double A1, double A2, double B1, double B2, double B); 

extern double Probamax(double A1, double A2, double B1, double B2, double B);

extern double Pay_Off_Spread(double S1, double S2, double K, double d1,
                      double d2, double Rate, double Mat,
                      double ProbaS1, double ProbaS2,
                      double ProbaK, int pc);

// fonctions de calcul du prix
extern double Spread_opt_norm(double S1, double S2, double vol1,
                       double vol2, double d1, double d2,double rho,
                       double r, double K, double t, int pc);                                                    
extern double Spread_opt_min(double S1, double S2, double vol1, double vol2,
                      double d1, double d2, double r, double K,
                      double t, int pc);                                                    
extern double Spread_opt_max(double S1, double S2, double vol1,
                      double vol2, double d1, double d2,
                      double r, double K, double t, int pc);                                                    


extern double MargrabeSpreadOption(double S1, double S2, 
                                   double vol1,double vol2,
							double d1, double d2, double Rho, 
                            double Rate, double Mat, int pc); 
                                               
extern double SpreadOption(double S1,double S2, double vol1, double vol2,
                    double d1, double d2,double rho, double r,
                    double K, double t, int pc, int ComputedFormula = 1 /* 0 formule d'Olivier, 1 Formule Classique*/);                              

#ifdef WIN32
extern double SpreadOption_ClosedForms(double S1, double S2,
									   double vol1, double vol2,
									   double Rho, double K,
									   double Mat, int pc);
#endif

extern double GaussianSpreadOption(double S1, double S2, double SpreadVol,
                            double d1, double d2, double Rho, 
							double Rate, double K, double Mat, int pc);



#undef MAXIT                                             



#endif
/*---------------------------------------------------------------------*/
/*---- End Of File ----*/
