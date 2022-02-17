/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file change_numeraire.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
#include "firsttoinc.h"
#include "gpbase/port.h"

#include "gpclosedforms/basic_distributions.h"

#include "gpclosedforms/tri_spreadoption_lognormal.h"
#include "gpclosedforms/tri_spreadoption_lognormal_formula.h"

#include "expt.h"

 

CC_BEGIN_NAMESPACE(ARM)


//////////////////////////////////////////////////////////////////////////////////////////
///
/// class : ARM_CF_TriSpreadDigitalOption_Formula
/// 	which pays  1  if A0+A1*S1+A2*S2+A3*S3>0
///
//////////////////////////////////////////////////////////////////////////////////////////


 double ARM_CF_TriSpreadDigitalOption_Formula::value(const ArgumentList& a)
	{
	 int argsize=a.size();
	 if (argsize!=ARM_CF_TriSpreadDigitalOption_Formula::Nb_Parameters) 
	 {
		 throw ("ARM_CF_TriSpreadDigitalOption_Formula::value  : bad argsize");
	 }
	 
	 double val= TriSpreadDigitalOption(
		 a[ARM_CF_TriSpreadDigitalOption_Formula::INDEX1],	//INDEX1
		 a[ARM_CF_TriSpreadDigitalOption_Formula::INDEX2],	//INDEX2
		 a[ARM_CF_TriSpreadDigitalOption_Formula::INDEX3],	//INDEX3
		 a[ARM_CF_TriSpreadDigitalOption_Formula::VOLATILITY1],	//VOLATILITY1
		 a[ARM_CF_TriSpreadDigitalOption_Formula::VOLATILITY2],	//VOLATILITY2
		 a[ARM_CF_TriSpreadDigitalOption_Formula::VOLATILITY3],	//VOLATILITY3
		 a[ARM_CF_TriSpreadDigitalOption_Formula::MU1],	//M1
		 a[ARM_CF_TriSpreadDigitalOption_Formula::MU2],	//M2
		 a[ARM_CF_TriSpreadDigitalOption_Formula::MU3],	//M3
		 a[ARM_CF_TriSpreadDigitalOption_Formula::CORRELATION12],	//CORRELATION12
		 a[ARM_CF_TriSpreadDigitalOption_Formula::CORRELATION13],	//CORRELATION13
		 a[ARM_CF_TriSpreadDigitalOption_Formula::CORRELATION23],	//CORRELATION23
		 a[ARM_CF_TriSpreadDigitalOption_Formula::A0],	//A0
		 a[ARM_CF_TriSpreadDigitalOption_Formula::A1],	//A1
		 a[ARM_CF_TriSpreadDigitalOption_Formula::A2],	//A2
		 a[ARM_CF_TriSpreadDigitalOption_Formula::A3],	//A3
		 a[ARM_CF_TriSpreadDigitalOption_Formula::TIMETORESET],	//t
		 a[ARM_CF_TriSpreadDigitalOption_Formula::CALLORPUT],	//CALL-PUT flag
		 a[ARM_CF_TriSpreadDigitalOption_Formula::NBSTEPS]	//NBSTEPS
		 );
	 return val;
	}




double ARM_CF_TriSpreadDigitalOption_Formula::value(int i,const ArgumentList& a, double s)
{
		return standard_first_derivative(ARM_CF_TriSpreadDigitalOption_Formula::value,i,a,s);
}


double ARM_CF_TriSpreadDigitalOption_Formula::value(int i,int j,const ArgumentList& a, double s1,double s2)
{
	
	return standard_second_derivative(ARM_CF_TriSpreadDigitalOption_Formula::value,i,j,a,s1,s2);

}



 double ARM_CF_TriSpreadDigitalOption_Formula::specific_shift(int i) {
		switch(i)
		{
		case ARM_CF_TriSpreadDigitalOption_Formula::TIMETORESET :
			return 0.00001;		// time to maturity
		case ARM_CF_TriSpreadDigitalOption_Formula::INDEX1 :
			return 0.00001;		// Forward index 1
		case ARM_CF_TriSpreadDigitalOption_Formula::INDEX2 :
			return 0.00001;		// Forward index 2
		case ARM_CF_TriSpreadDigitalOption_Formula::INDEX3 :
			return 0.00001;		// Forward index 3
		case ARM_CF_TriSpreadDigitalOption_Formula::VOLATILITY1 :
			return 0.00001;		// Total Volatility 1
		case ARM_CF_TriSpreadDigitalOption_Formula::VOLATILITY2 :
			return 0.00001;		// Total Volatility 2
		case ARM_CF_TriSpreadDigitalOption_Formula::VOLATILITY3 :
			return 0.00001;		// Total Volatility 3
		case ARM_CF_TriSpreadDigitalOption_Formula::CORRELATION12 :
			return 0.00001;		// Correlation between the two indexes
		case ARM_CF_TriSpreadDigitalOption_Formula::CORRELATION13 :
			return 0.00001;		// Correlation between the two indexes
		case ARM_CF_TriSpreadDigitalOption_Formula::CORRELATION23 :
			return 0.00001;		// Correlation between the two indexes
		case ARM_CF_TriSpreadDigitalOption_Formula::MU1 :
			return 0.00001;		// MU1
		case ARM_CF_TriSpreadDigitalOption_Formula::MU2 :
			return 0.00001;		// MU2
		case ARM_CF_TriSpreadDigitalOption_Formula::MU3 :
			return 0.00001;		// MU3
		default :
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_TriSpreadDigitalOption_Formula::specific_shift : incorrect input");
		}
	}



 ArgumentList_Checking_Result ARM_CF_TriSpreadDigitalOption_Formula::check_argument(const ArgumentList& a)
 {
	 if(a.size()!=ARM_CF_TriSpreadDigitalOption_Formula::Nb_Parameters) 
	 {
		 return ArgumentList_Checking_Result(false,"Bad number of arguments  ");
	 }
	 
	 if (a[ARM_CF_TriSpreadDigitalOption_Formula::TIMETORESET]<=0) return ArgumentList_Checking_Result(false," TIMETORESET not positive"); //		positivity of TIMETORESET

	 if (a[ARM_CF_TriSpreadDigitalOption_Formula::INDEX1]<=0) return ArgumentList_Checking_Result(false," INDEX1 not positive"); //			positivity of INDEX1
	 if (a[ARM_CF_TriSpreadDigitalOption_Formula::INDEX2]<=0) return ArgumentList_Checking_Result(false," INDEX2 not positive"); //			positivity of INDEX2
	 if (a[ARM_CF_TriSpreadDigitalOption_Formula::INDEX3]<=0) return ArgumentList_Checking_Result(false," INDEX3 not positive"); //			positivity of INDEX3
	
	 if (a[ARM_CF_TriSpreadDigitalOption_Formula::VOLATILITY1 ]<=0) return ArgumentList_Checking_Result(false," VOLATILITY1 not positive"); //		positivity of VOLATILITY1	
	 if (a[ARM_CF_TriSpreadDigitalOption_Formula::VOLATILITY2]<=0) return ArgumentList_Checking_Result(false," VOLATILITY2 not positive"); //		positivity of VOLATILITY2
	 if (a[ARM_CF_TriSpreadDigitalOption_Formula::VOLATILITY3]<=0) return ArgumentList_Checking_Result(false," VOLATILITY3 not positive"); //		positivity of VOLATILITY3
	 
	 if (abs(a[ARM_CF_TriSpreadDigitalOption_Formula::CORRELATION12])>1.) return ArgumentList_Checking_Result(false," Abs(CORRELATION12) not <=1"); //		Abs(CORRELATION12)<=1
	 if (abs(a[ARM_CF_TriSpreadDigitalOption_Formula::CORRELATION13])>1.) return ArgumentList_Checking_Result(false," Abs(CORRELATION13) not <=1"); //		Abs(CORRELATION13)<=1
	 if (abs(a[ARM_CF_TriSpreadDigitalOption_Formula::CORRELATION23])>1.) return ArgumentList_Checking_Result(false," Abs(CORRELATION23) not <=1"); //		Abs(CORRELATION23)<=1

	 if (abs(a[ARM_CF_TriSpreadDigitalOption_Formula::NBSTEPS])<=1) return ArgumentList_Checking_Result(false," NBSTEPS should be >= 1"); //		NBSTEPS should >=1
	
	if ((a[ARM_CF_TriSpreadDigitalOption_Formula::CALLORPUT]!=K_CALL) &&
		  (a[ARM_CF_TriSpreadDigitalOption_Formula::CALLORPUT]!=K_PUT))
		  return ArgumentList_Checking_Result(false," ARM_CF_TriSpreadDigitalOption_Formula:CALLORPUT should be CALL or PUT"); //			positivity of STRIKE
 	return ArgumentList_Checking_Result(true,string(""));
 }


 ArgumentList_Checking_Result ARM_CF_TriSpreadDigitalOption_Formula::check_dimension(int rank)
{
	if ((rank<0)||(rank>12)) return ArgumentList_Checking_Result(false," Invalide derivative dimension"); // Invalide derivative dimension	
	
	return ArgumentList_Checking_Result(true,string(""));
}

 //////////////////////////////////////////////////////////////////////////////////////////
///
/// class : ARM_CF_TriSpreadDigitalOption2_Formula
/// 	which pays  1  if A0+A2*S2+A3*S3>0 and B0+B1*S1>0
///
//////////////////////////////////////////////////////////////////////////////////////////


 double ARM_CF_TriSpreadDigitalOption2_Formula::value(const ArgumentList& a)
	{
	 int argsize=a.size();
	 if (argsize!=ARM_CF_TriSpreadDigitalOption2_Formula::Nb_Parameters) 
	 {
		 throw ("ARM_CF_TriSpreadDigitalOption2_Formula::value  : bad argsize");
	 }
	 
	 double val= TriSpreadDigitalCall2(
		 a[ARM_CF_TriSpreadDigitalOption2_Formula::INDEX1],	//INDEX1
		 a[ARM_CF_TriSpreadDigitalOption2_Formula::INDEX2],	//INDEX2
		 a[ARM_CF_TriSpreadDigitalOption2_Formula::INDEX3],	//INDEX3
		 a[ARM_CF_TriSpreadDigitalOption2_Formula::VOLATILITY1],	//VOLATILITY1
		 a[ARM_CF_TriSpreadDigitalOption2_Formula::VOLATILITY2],	//VOLATILITY2
		 a[ARM_CF_TriSpreadDigitalOption2_Formula::VOLATILITY3],	//VOLATILITY3
		 a[ARM_CF_TriSpreadDigitalOption2_Formula::MU1],	//M1
		 a[ARM_CF_TriSpreadDigitalOption2_Formula::MU2],	//M2
		 a[ARM_CF_TriSpreadDigitalOption2_Formula::MU3],	//M3
		 a[ARM_CF_TriSpreadDigitalOption2_Formula::CORRELATION12],	//CORRELATION12
		 a[ARM_CF_TriSpreadDigitalOption2_Formula::CORRELATION13],	//CORRELATION13
		 a[ARM_CF_TriSpreadDigitalOption2_Formula::CORRELATION23],	//CORRELATION23
		 a[ARM_CF_TriSpreadDigitalOption2_Formula::A0],	//A0
		 a[ARM_CF_TriSpreadDigitalOption2_Formula::A2],	//A2
		 a[ARM_CF_TriSpreadDigitalOption2_Formula::A3],	//A3
		 a[ARM_CF_TriSpreadDigitalOption2_Formula::B0],	//A3
		 a[ARM_CF_TriSpreadDigitalOption2_Formula::B1],	//A3
		 a[ARM_CF_TriSpreadDigitalOption2_Formula::TIMETORESET],	//t
		 a[ARM_CF_TriSpreadDigitalOption2_Formula::NBSTEPS]	//NBSTEPS
		 );
	 return val;
	}




double ARM_CF_TriSpreadDigitalOption2_Formula::value(int i,const ArgumentList& a, double s)
{
		return standard_first_derivative(ARM_CF_TriSpreadDigitalOption2_Formula::value,i,a,s);
}


double ARM_CF_TriSpreadDigitalOption2_Formula::value(int i,int j,const ArgumentList& a, double s1,double s2)
{
	
	return standard_second_derivative(ARM_CF_TriSpreadDigitalOption2_Formula::value,i,j,a,s1,s2);

}



 double ARM_CF_TriSpreadDigitalOption2_Formula::specific_shift(int i) {
		switch(i)
		{
		case ARM_CF_TriSpreadDigitalOption2_Formula::TIMETORESET :
			return 0.00001;		// time to maturity
		case ARM_CF_TriSpreadDigitalOption2_Formula::INDEX1 :
			return 0.00001;		// Forward index 1
		case ARM_CF_TriSpreadDigitalOption2_Formula::INDEX2 :
			return 0.00001;		// Forward index 2
		case ARM_CF_TriSpreadDigitalOption2_Formula::INDEX3 :
			return 0.00001;		// Forward index 3
		case ARM_CF_TriSpreadDigitalOption2_Formula::VOLATILITY1 :
			return 0.00001;		// Total Volatility 1
		case ARM_CF_TriSpreadDigitalOption2_Formula::VOLATILITY2 :
			return 0.00001;		// Total Volatility 2
		case ARM_CF_TriSpreadDigitalOption2_Formula::VOLATILITY3 :
			return 0.00001;		// Total Volatility 3
		case ARM_CF_TriSpreadDigitalOption2_Formula::CORRELATION12 :
			return 0.00001;		// Correlation between the two indexes
		case ARM_CF_TriSpreadDigitalOption2_Formula::CORRELATION13 :
			return 0.00001;		// Correlation between the two indexes
		case ARM_CF_TriSpreadDigitalOption2_Formula::CORRELATION23 :
			return 0.00001;		// Correlation between the two indexes
		case ARM_CF_TriSpreadDigitalOption2_Formula::MU1 :
			return 0.00001;		// MU1
		case ARM_CF_TriSpreadDigitalOption2_Formula::MU2 :
			return 0.00001;		// MU2
		case ARM_CF_TriSpreadDigitalOption2_Formula::MU3 :
			return 0.00001;		// MU3
		default :
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_TriSpreadDigitalOption2_Formula::specific_shift : incorrect input");
		}
	}



 ArgumentList_Checking_Result ARM_CF_TriSpreadDigitalOption2_Formula::check_argument(const ArgumentList& a)
 {
	 if(a.size()!=ARM_CF_TriSpreadDigitalOption2_Formula::Nb_Parameters) 
	 {
		 return ArgumentList_Checking_Result(false,"Bad number of arguments  ");
	 }
	 
	 if (a[ARM_CF_TriSpreadDigitalOption2_Formula::TIMETORESET]<=0) return ArgumentList_Checking_Result(false," TIMETORESET not positive"); //		positivity of TIMETORESET

	 if (a[ARM_CF_TriSpreadDigitalOption2_Formula::INDEX1]<=0) return ArgumentList_Checking_Result(false," INDEX1 not positive"); //			positivity of INDEX1
	 if (a[ARM_CF_TriSpreadDigitalOption2_Formula::INDEX2]<=0) return ArgumentList_Checking_Result(false," INDEX2 not positive"); //			positivity of INDEX2
	 if (a[ARM_CF_TriSpreadDigitalOption2_Formula::INDEX3]<=0) return ArgumentList_Checking_Result(false," INDEX3 not positive"); //			positivity of INDEX3
	
	 if (a[ARM_CF_TriSpreadDigitalOption2_Formula::VOLATILITY1 ]<=0) return ArgumentList_Checking_Result(false," VOLATILITY1 not positive"); //		positivity of VOLATILITY1	
	 if (a[ARM_CF_TriSpreadDigitalOption2_Formula::VOLATILITY2]<=0) return ArgumentList_Checking_Result(false," VOLATILITY2 not positive"); //		positivity of VOLATILITY2
	 if (a[ARM_CF_TriSpreadDigitalOption2_Formula::VOLATILITY3]<=0) return ArgumentList_Checking_Result(false," VOLATILITY3 not positive"); //		positivity of VOLATILITY3
	 
	 if (abs(a[ARM_CF_TriSpreadDigitalOption2_Formula::CORRELATION12])>1.) return ArgumentList_Checking_Result(false," Abs(CORRELATION12) not <=1"); //		Abs(CORRELATION12)<=1
	 if (abs(a[ARM_CF_TriSpreadDigitalOption2_Formula::CORRELATION13])>1.) return ArgumentList_Checking_Result(false," Abs(CORRELATION13) not <=1"); //		Abs(CORRELATION13)<=1
	 if (abs(a[ARM_CF_TriSpreadDigitalOption2_Formula::CORRELATION23])>1.) return ArgumentList_Checking_Result(false," Abs(CORRELATION23) not <=1"); //		Abs(CORRELATION23)<=1

	 if (abs(a[ARM_CF_TriSpreadDigitalOption2_Formula::NBSTEPS])<=1) return ArgumentList_Checking_Result(false," NBSTEPS should be >= 1"); //		NBSTEPS should >=1
 	return ArgumentList_Checking_Result(true,string(""));
 }


 ArgumentList_Checking_Result ARM_CF_TriSpreadDigitalOption2_Formula::check_dimension(int rank)
{
	if ((rank<0)||(rank>12)) return ArgumentList_Checking_Result(false," Invalide derivative dimension"); // Invalide derivative dimension	
	
	return ArgumentList_Checking_Result(true,string(""));
}




CC_END_NAMESPACE()
 
#undef ARM_CF_SQRT_2_PI 
#undef ARM_CF_SQRT2 
#undef ARM_CF_INVSQRTPI

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
