/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 03/15/2007
 *
 *  basic functions for the closed form framework 
 *
 *	\file spreadoption_nonparametric_formula.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date march 2007
 */
#include "firsttoinc.h"
#include "gpbase/port.h"

#include "gpclosedforms/spreadoption_nonparametric_formula.h"
#include "gpclosedforms/smile_nonparametric.h"
#include "gpclosedforms/digital_templated_pricing.h"
#include "gpclosedforms/nonparametric_spline.h"
#include "gpclosedforms/nonparametric_quantile.h"
#include "gpclosedforms/bismile_templated_pricing.h"


#include "expt.h"




CC_BEGIN_NAMESPACE(ARM)




////////////////////////////////////////////////////////////////////////
///
///      Value of the formula (introduces the checking and the ability to automatically derive
///
////////////////////////////////////////////////////////////////////////


double ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::value(const ArgumentList& a)
{
	return dummy_VC6<ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula>::ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula_value(a);
}


void ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::Tail_1_Slope_Compute(std::vector<double>* x,std::vector<double>* y,double index_begin,double index_end,
																						  int beginflag,int endflag, double* beginderivative,double* endderivative )
{
	NonParametric_LogVolatility_TailSlope_Compute(x,y,index_begin,index_end,beginflag,endflag,beginderivative,endderivative);
}
void ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::Tail_2_Slope_Compute(std::vector<double>* x,std::vector<double>* y,double index_begin,double index_end,
																						  int beginflag,int endflag, double* beginderivative,double* endderivative )
{
	NonParametric_LogVolatility_TailSlope_Compute(x,y,index_begin,index_end,beginflag,endflag,beginderivative,endderivative);
}

double ARM_CF_NonParametric_Gaussian_PowerSpreadOption2_Formula::value(const ArgumentList& a)
{
	return dummy_VC6<ARM_CF_NonParametric_Gaussian_PowerSpreadOption2_Formula>::ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula_value(a);
}


void ARM_CF_NonParametric_Gaussian_PowerSpreadOption2_Formula::Tail_1_Slope_Compute(std::vector<double>* x,std::vector<double>* y,double index_begin,double index_end,
																					int beginflag,int endflag, double* beginderivative,double* endderivative )
{
	NonParametric_NormalVolatility_TailSlope_Compute(x,y,index_begin,index_end,beginflag,endflag,beginderivative,endderivative);
}
void ARM_CF_NonParametric_Gaussian_PowerSpreadOption2_Formula::Tail_2_Slope_Compute(std::vector<double>* x,std::vector<double>* y,double index_begin,double index_end,
																					int beginflag,int endflag, double* beginderivative,double* endderivative )
{
	NonParametric_NormalVolatility_TailSlope_Compute(x,y,index_begin,index_end,beginflag,endflag,beginderivative,endderivative);
}



double ARM_CF_NonParametric_Gaussian_PowerSpreadOption3_Formula::value(const ArgumentList& a)
{
	return dummy_VC6<ARM_CF_NonParametric_Gaussian_PowerSpreadOption3_Formula>::ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula_value(a);
}

void ARM_CF_NonParametric_Gaussian_PowerSpreadOption3_Formula::Tail_1_Slope_Compute(std::vector<double>* x,std::vector<double>* y,double index_begin,double index_end,
																					int beginflag,int endflag, double* beginderivative,double* endderivative )
{
	NonParametric_LogVolatility_TailSlope_Compute(x,y,index_begin,index_end,beginflag,endflag,beginderivative,endderivative);
}
void ARM_CF_NonParametric_Gaussian_PowerSpreadOption3_Formula::Tail_2_Slope_Compute(std::vector<double>* x,std::vector<double>* y,double index_begin,double index_end,
																					int beginflag,int endflag, double* beginderivative,double* endderivative )
{
	NonParametric_NormalVolatility_TailSlope_Compute(x,y,index_begin,index_end,beginflag,endflag,beginderivative,endderivative);
}


double ARM_CF_NonParametric_Gaussian_PowerSpreadOption4_Formula::value(const ArgumentList& a)
{
	return dummy_VC6<ARM_CF_NonParametric_Gaussian_PowerSpreadOption4_Formula>::ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula_value(a);
}

void ARM_CF_NonParametric_Gaussian_PowerSpreadOption4_Formula::Tail_1_Slope_Compute(std::vector<double>* x,std::vector<double>* y,double index_begin,double index_end,
																					int beginflag,int endflag, double* beginderivative,double* endderivative )
{
	NonParametric_NormalVolatility_TailSlope_Compute(x,y,index_begin,index_end,beginflag,endflag,beginderivative,endderivative);
}
void ARM_CF_NonParametric_Gaussian_PowerSpreadOption4_Formula::Tail_2_Slope_Compute(std::vector<double>* x,std::vector<double>* y,double index_begin,double index_end,
																					int beginflag,int endflag, double* beginderivative,double* endderivative )
{
	NonParametric_LogVolatility_TailSlope_Compute(x,y,index_begin,index_end,beginflag,endflag,beginderivative,endderivative);
}


////////////////////////////////////////////////////////////////////////
///
///      Value of the derivatives
///
////////////////////////////////////////////////////////////////////////


double ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::value(int i,const ArgumentList& a, double s)
///     i vaut:
///		
///		FORWARD1,
///		FORWARD2,
///		INDEXBEGIN1,
///		INDEXEND1,
///		FLAGBEGIN1,
///		FLAGEND1,
///		INDEXBEGIN2,
///		INDEXEND2,
///		FLAGBEGIN2,
///		FLAGEND2,
///		CORRELATION,
///		TIMETOEXPIRATION,
///		A1,
///		B1,
///		K1,
///		A2,
///		B2,
///		K2,


{
			return standard_first_derivative(ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::value,i,a,s);
}


////////////////////////////////////////////////////////////////////////
///
///      Second Derivative
///
////////////////////////////////////////////////////////////////////////


double ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::value(int i,int j,const ArgumentList& a, double s1,double s2)
{
	 		return standard_second_derivative(ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::value,i,j,a,s1,s2);
}



////////////////////////////////////////////////////////////////////////
///
///      Description of the shifting used for the derivationhttp://www.pixmania.com/fr/fr/tv-video/3/onglet.html?itag=9810
///
////////////////////////////////////////////////////////////////////////


double ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::specific_shift(int i) {
	switch(i)
	{
		///		Here the limitation of the dependence in the smile and the copula is that they 
		///	are described by alpha,beta,rho,nu   and the copula described by one number (called here correlation)
		///

	case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FORWARD1:
		return 0.00001;		// 
	case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FORWARD2:
		return 0.00001;		// 
	case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::INDEXBEGIN1:
		return 0.00001;		// 
	case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::INDEXEND1:
		return 0.00001;		// 
	case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FLAGBEGIN1:
		return 0.00001;		// 
	case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FLAGEND1:
		return 0.00001;		// 
	case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::INDEXBEGIN2:
		return 0.00001;		// 
	case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::INDEXEND2:
		return 0.00001;		// 
	case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FLAGBEGIN2:
		return 0.00001;		// 
	case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FLAGEND2:
		return 0.00001;		// 
	case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::CORRELATION:
		return 0.00001;		// 
	case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::TIMETOEXPIRATION:
		return 0.00001;		// 
	case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::A1:
		return 0.00001;		// 
	case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::B1:
		return 0.00001;		// 
	case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::K1:
		return 0.00001;		// 
	case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::A2:
		return 0.00001;		// 
	case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::B2:
		return 0.00001;		// 
	case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::K2:
		return 0.00001;		// 
	default :
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::specific_shift : incorrect input");
	}
}

ArgumentList_Checking_Result ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula_check_argument(
													const ArgumentList& a,int posforw1,int posforw2)
{
	///		Here the limitation of the dependence in the smile and the copula is that they 
	///	are described by alpha,beta,nu  all positive number and rho, included in {-1,1}
	///  and the copula described by one number (called here correlation)
	///
	
	if(a.size()!=ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::Nb_Parameters) 
	{
		return ArgumentList_Checking_Result(false,"Bad number of arguments ");
	}
///		INDEXBEGIN1,
///		INDEXEND1,
///		FLAGBEGIN1,
///		FLAGEND1,
///		INDEXBEGIN2,
///		INDEXEND2,
///		FLAGBEGIN2,
///		FLAGEND2,
///		CORRELATION,
///		TIMETOEXPIRATION,
///		A1,
///		B1,
///		K1,
///		A2,
///		B2,
///		K2,
if(posforw1>0)
{
	if(a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FORWARD1]<=0) return ArgumentList_Checking_Result(false," Forward1  not positive"); //			positivity of INDEX1
}
if(posforw2>0)
{
	if(a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FORWARD2]<=0) return ArgumentList_Checking_Result(false," Forward2  not positive"); //			positivity of INDEX2
}

	if(abs(a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::CORRELATION])>1.) return ArgumentList_Checking_Result(false," CORRELATION   not positive"); //	Abs(CORRELATION)<=1 
	if(a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::TIMETOEXPIRATION]<=0) return ArgumentList_Checking_Result(false," TIMETOEXPIRATION not positive"); //positivity of TIMETOEXPIRATION
  int a0=a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FLAGBEGIN1];
  int b0=ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYTAILINDEX;
	if (
		(a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FLAGBEGIN1]!=
		ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYTAILINDEX)&&
		(a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FLAGBEGIN1]!=
		ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYCONSTANT)&&
		(a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FLAGBEGIN1]!=
		ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYSLOGTAILINDEX)&&
		(a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FLAGBEGIN1]!=
		ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYLINEAR)&&
		(a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FLAGBEGIN1]!=
		ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::DISTRIBUTIONSMOOTH)
		)
	{
	return ArgumentList_Checking_Result(false," ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FLAGBEGIN1 should be 0,1,2,3"); //
	}
		if (
		(a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FLAGEND1]!=
		ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYTAILINDEX)&&
		(a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FLAGEND1]!=
		ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYCONSTANT)&&
		(a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FLAGEND1]!=
		ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYSLOGTAILINDEX)&&
		(a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FLAGEND1]!=
		ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYLINEAR)&&
		(a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FLAGEND1]!=
		ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::DISTRIBUTIONSMOOTH)
		)
	{
	return ArgumentList_Checking_Result(false," ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FLAGEND1 should be 0,1,2,3"); //
	}
			if (
		(a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FLAGBEGIN2]!=
		ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYTAILINDEX)&&
		(a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FLAGBEGIN2]!=
		ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYCONSTANT)&&
		(a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FLAGBEGIN2]!=
		ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYSLOGTAILINDEX)&&
		(a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FLAGBEGIN2]!=
		ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYLINEAR)&&
		(a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FLAGBEGIN2]!=
		ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::DISTRIBUTIONSMOOTH)
		)
	{
	return ArgumentList_Checking_Result(false," ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FLAGBEGIN2 should be 0,1,2,3"); //
	}
		if (
		(a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FLAGEND2]!=
		ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYTAILINDEX)&&
		(a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FLAGEND2]!=
		ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYCONSTANT)&&
		(a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FLAGEND2]!=
		ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYSLOGTAILINDEX)&&
		(a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FLAGEND2]!=
		ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYLINEAR)&&
		(a[ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FLAGEND2]!=
		ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::DISTRIBUTIONSMOOTH)
		)
	{
	return ArgumentList_Checking_Result(false," ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::FLAGEND2 should be 0,1,2,3"); //
	}
	return ArgumentList_Checking_Result(true,string(""));
}


////////////////////////////////////////////////////////////////////////
///
///      Checking of the arguments : description of the domain of validity
///
////////////////////////////////////////////////////////////////////////


ArgumentList_Checking_Result ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::check_argument(const ArgumentList& a)
{
	return ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula_check_argument(a,1,1);
}


////////////////////////////////////////////////////////////////////////
///
///      Checking of the dimension rank with respect to which one want ot derive
///
////////////////////////////////////////////////////////////////////////


ArgumentList_Checking_Result ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::check_dimension(int rank)
{
	///		Here the limitation of the dependence in the smile and the copula is that they 
	///	are described by alpha,beta,rho,nu   and the copula described by one number (called here correlation)
	///		so the number of parameters is 19
	if ((rank<0)||(rank>=ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::Nb_Derivable_Parameters)) return ArgumentList_Checking_Result(false," Invalide derivative dimension"); // Invalide derivative dimension	
	
	return ArgumentList_Checking_Result(true,string(""));
}





////////////////////////////////////////////////////////////////////////
///
///      Checking of the arguments : description of the domain of validity
///                  for the derived classes
///
////////////////////////////////////////////////////////////////////////

 
 ArgumentList_Checking_Result ARM_CF_NonParametric_Gaussian_PowerSpreadOption2_Formula::check_argument(const ArgumentList& a)
 {
	 return ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula_check_argument(a,0,0);
 }
 
 ArgumentList_Checking_Result ARM_CF_NonParametric_Gaussian_PowerSpreadOption3_Formula::check_argument(const ArgumentList& a)
 {
	 return ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula_check_argument(a,1,0);
 }
 
 ArgumentList_Checking_Result ARM_CF_NonParametric_Gaussian_PowerSpreadOption4_Formula::check_argument(const ArgumentList& a)
 {
	 return ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula_check_argument(a,0,1);
 }


#undef ARM_CF_SQRT_2_PI 
#undef ARM_CF_SQRT2 
#undef ARM_CF_INVSQRTPI 

CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/