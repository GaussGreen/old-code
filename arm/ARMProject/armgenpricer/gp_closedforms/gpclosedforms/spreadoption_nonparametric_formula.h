/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 03/15/2007
 *
 *  basic functions for the closed form framework 
 *
 *	\file spreadoption_nonparametric_formula.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date march 2007
 */
 
#ifndef _GP_CF_SPREADOPTION_NONPARAMETRIC_FORMULA_H
#define _GP_CF_SPREADOPTION_NONPARAMETRIC_FORMULA_H

#include "firsttoinc.h"
#include "gpbase/port.h"

#include "gpclosedforms/basic_distributions.h" /// for ArgumentList_Checking_Result
#include "gpclosedforms/powerspreadoption.h"		/// for instanciation of the templates
#include "gpclosedforms/smile_nonparametric.h"
#include "gpclosedforms/nonparametric_spline.h"
#include "gpclosedforms/bismile_templated_pricing.h"
#include "gpclosedforms/nonparametric_quantile.h"


CC_BEGIN_NAMESPACE(ARM)


//////////////////////////////////////////////////////////////////////
///  
///			class : ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula
///			Purpose : Evaluation of spread option cashflow=  (a1*S1+b1*S2-k1) if {a2*S1+b2*S2-k2>0},  and 0 otherwise
///			Assumptions: Extended Non parametric hypothesis fo the factors and a gaussian copula
///
///////////////////////////////////////////////////////////////////////


struct ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula
{
	typedef ARM_CF_Bismiled_PowerSpreadOption_Formula<NonParametric_LogSmile,NonParametric_LogSmile,GaussianCopula> structure;
	enum ArgumentType
	{
			FORWARD1,
			FORWARD2,
			INDEXBEGIN1,
			INDEXEND1,
			FLAGBEGIN1,
			FLAGEND1,
			INDEXBEGIN2,
			INDEXEND2,
			FLAGBEGIN2,
			FLAGEND2,
			CORRELATION,
			TIMETOEXPIRATION,
			A1,
			B1,
			K1,
			A2,
			B2,
			K2,
			NBSTEPS,
			ALGORITHM			/// internal
	};
	enum VectorArgumentType
	{
			STRIKEVECTOR1,
			VOLVECTOR1,
			STRIKEVECTOR2,
			VOLVECTOR2,
	};
	
	enum	TailExtrapolationMethod
	{	
		VOLATILITYTAILINDEX,
		VOLATILITYCONSTANT,
		VOLATILITYSLOGTAILINDEX,
		VOLATILITYLINEAR,
		DISTRIBUTIONSMOOTH
	};


	enum 
	{ 
		Nb_Parameters =20
	};
	enum 
	{ 
		Nb_Parameters2 =4
	};

	enum
	{ 
		Nb_Derivable_Parameters =18
	};
	static double value(const ArgumentList& a);
	static double value(int i,const ArgumentList& a, double s);
	static double value2(const ArgumentList& a);
	static double value(int i,int j,const ArgumentList& a, double s1,double s2);

	static double specific_shift(int i) ;
	static double value(int i,const ArgumentList& a);
	static ArgumentList_Checking_Result check_argument(const ArgumentList& a);
	static ArgumentList_Checking_Result check_dimension(int rank);
	///   Generic methods (relative to the smile distribution and the copula
	static double Generic_Pricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula,double t,int n,BivariateToDoubleFunc* option_payoff);
	static double PowerSpreadOption_Pricing(const ArgumentList& Underlying1,const ArgumentList& Underlying2,
				   const ArgumentList& copula,double t,int n,
				   double a10,double b10,double k10,double a20,double b20,double k20);
	static double Certitude(const ArgumentList& copula,double t,int n);
	static void Tail_1_Slope_Compute(std::vector<double>* x,std::vector<double>* y,double index_begin,double index_end,
							int beginflag,int endflag, double* beginderivative,double* endderivative );
	static void Tail_2_Slope_Compute(std::vector<double>* x,std::vector<double>* y,double index_begin,double index_end,
							int beginflag,int endflag, double* beginderivative,double* endderivative );


};

template<typename T>
class dummy_VC6
{
	public:
	static	double ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula_value(const ArgumentList& a);
};

template<typename A>
double dummy_VC6<A>::ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula_value(const ArgumentList& a)
{
	int argsize=a.size();
	int argsize2=a.size2();
	if (argsize!=ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::Nb_Parameters)
	{
		 throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::value  : bad argsize");
	 }
	 if (argsize2!=ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::Nb_Parameters2)
	 {
		 throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::value  : bad argsize2");
	 }
	 ///		Here the limitation of the dependence in the smile and the copula is that they 
	 ///	are described by alpha,beta,rho,nu   and the copula described by one number (called here correlation)
	 ///
	 ///		we use here only PowerSpreadOption_Pricing and not Basic_Pricing_With_Limits which is more efficient but which involve the computation of 
	 ///		DistributionZaInverseLimit<smile,copula> that is not simple and prone to error in the case of SABR 
	 ///
	 ArgumentList* CopulaArg=new ArgumentList(
		 a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::CORRELATION]	//CORRELATION
		 );
	GaussianCopula cop(CopulaArg);


	std::vector<double>* x1=a.V(ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::STRIKEVECTOR1);
	std::vector<double>* y1=a.V(ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLVECTOR1);
	int n1size=x1->size();
	std::vector<double> SecondDer1(n1size);
	double beginderivative1,endderivative1;

	A::Tail_1_Slope_Compute(x1,y1,
		a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::INDEXBEGIN1				],	//
		a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::INDEXEND1				],	//
		a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FLAGBEGIN1				],	//
		a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FLAGEND1					],	//
		&beginderivative1,&endderivative1 );

	cubicspline_precompute(
		*x1,
		*y1,
		beginderivative1,
		endderivative1,
		SecondDer1);
	
	std::vector<double>* x2=a.V(ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::STRIKEVECTOR2);
	std::vector<double>* y2=a.V(ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLVECTOR2);
	int n2size=x2->size();
	std::vector<double> SecondDer2(n2size);
	double beginderivative2,endderivative2;

	A::Tail_2_Slope_Compute(x2,y2,
		a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::INDEXBEGIN2				],	//
		a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::INDEXEND2				],	//
		a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FLAGBEGIN2				],	//
		a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FLAGEND2					],	//
		&beginderivative2,&endderivative2 );

	cubicspline_precompute(
		*x2,
		*y2,
		beginderivative2,
		endderivative2,
		SecondDer2);

	double val=A::structure::Pricing(

		 ArgumentList(
		 a.V(ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::STRIKEVECTOR1			),	//
		 a.V(ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLVECTOR1			),	//
		 &SecondDer1,	//
		 a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::INDEXBEGIN1				],	//
		 a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::INDEXEND1				],	//
		 a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FLAGBEGIN1				],	//
		 a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FLAGEND1				],	//
		 a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FORWARD1				]	//
		 ),

		 ArgumentList(
		 a.V(ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::STRIKEVECTOR2			),	//
		 a.V(ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLVECTOR2			),	//
		 &SecondDer2,	//
		 a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::INDEXBEGIN2				],	//
		 a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::INDEXEND2				],	//
		 a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FLAGBEGIN2				],	//
		 a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FLAGEND2				],	//
		 a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FORWARD2				]	//
		 ),

		 cop,				/// Copula

		 a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::TIMETOEXPIRATION],	//TIMETOEXPIRATION
		 a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::NBSTEPS			],
		 a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::A1				],	//A1
		 a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::B1				],	//B1
		 a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::K1				],	//K1
		 a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::A2				],	//A2
		 a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::B2				],	//B2 
		 a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::K2				],	//K2
		 a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::ALGORITHM		]	//Algorithm
		 );
		 
	 return val;
}


//////////////////////////////////////////////////////////////////////
///  
///			class : ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula
///			Purpose : Evaluation of spread option cashflow=  (a1*S1+b1*S2-k1) if {a2*S1+b2*S2-k2>0},  and 0 otherwise
///			Assumptions: Extended Non parametric hypothesis fo the factors and a gaussian copula
///
///////////////////////////////////////////////////////////////////////


struct ARM_CF_NonParametric_Gaussian_PowerSpreadOption2_Formula : public ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula
{
	typedef ARM_CF_Bismiled_PowerSpreadOption_Formula<NonParametric_NormalSmile,NonParametric_NormalSmile,GaussianCopula> structure;	
	static ArgumentList_Checking_Result check_argument(const ArgumentList& a);
	static double value(const ArgumentList& a);
	static void Tail_1_Slope_Compute(std::vector<double>* x,std::vector<double>* y,double index_begin,double index_end,
		int beginflag,int endflag, double* beginderivative,double* endderivative );
	static void Tail_2_Slope_Compute(std::vector<double>* x,std::vector<double>* y,double index_begin,double index_end,
		int beginflag,int endflag, double* beginderivative,double* endderivative );
	
	
};

struct ARM_CF_NonParametric_Gaussian_PowerSpreadOption3_Formula : public ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula
{
	typedef ARM_CF_Bismiled_PowerSpreadOption_Formula<NonParametric_LogSmile,NonParametric_NormalSmile,GaussianCopula> structure;	
	static ArgumentList_Checking_Result check_argument(const ArgumentList& a);
	static double value(const ArgumentList& a);
	
	static void Tail_1_Slope_Compute(std::vector<double>* x,std::vector<double>* y,double index_begin,double index_end,
		int beginflag,int endflag, double* beginderivative,double* endderivative );
	static void Tail_2_Slope_Compute(std::vector<double>* x,std::vector<double>* y,double index_begin,double index_end,
		int beginflag,int endflag, double* beginderivative,double* endderivative );
	
};

struct ARM_CF_NonParametric_Gaussian_PowerSpreadOption4_Formula : public ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula
{
	typedef ARM_CF_Bismiled_PowerSpreadOption_Formula<NonParametric_NormalSmile,NonParametric_LogSmile,GaussianCopula> structure;	
	static ArgumentList_Checking_Result check_argument(const ArgumentList& a);
	static double value(const ArgumentList& a);
	
	static void Tail_1_Slope_Compute(std::vector<double>* x,std::vector<double>* y,double index_begin,double index_end,
		int beginflag,int endflag, double* beginderivative,double* endderivative );
	static void Tail_2_Slope_Compute(std::vector<double>* x,std::vector<double>* y,double index_begin,double index_end,
		int beginflag,int endflag, double* beginderivative,double* endderivative );
	
};

CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

