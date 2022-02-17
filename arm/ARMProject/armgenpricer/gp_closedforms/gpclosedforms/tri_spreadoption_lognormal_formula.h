/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file spreadoption_lognormal.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_TRISPREADOPTION_LOGNORMAL_FORMULA_H
#define _GP_CF_TRISPREADOPTION_LOGNORMAL_FORMULA_H

#include "firsttoinc.h"
#include "gpbase/port.h"

#include "basic_distributions.h" /// for ArgumentList_Checking_Result

#include "product.h"			/// for typedef manipulation !
#include "sum.h"
#include "difference.h"
#include "add_submodel.h"
#include "change_numeraire.h"


CC_BEGIN_NAMESPACE(ARM)



///////////////////////////////////////////////////////////////////////
///  
///			class : ARM_CF_TriSpreadDigitalOption_Formula
///			Purpose : Evaluation of spread option cashflow= {A0+A1*S1+A2*S2+A3*S3}^+
///			Assumptions: Lognormal hypothesis
///
///////////////////////////////////////////////////////////////////////

struct ARM_CF_TriSpreadDigitalOption_Formula
{
	enum ArgumentType
	{
		TIMETORESET,
		INDEX1,
		INDEX2,
		INDEX3,
		VOLATILITY1,
		VOLATILITY2,
		VOLATILITY3,
		CORRELATION12,
		CORRELATION13,
		CORRELATION23,
		MU1,
		MU2,
		MU3,
        A0,
		A1,
		A2,
		A3,
		CALLORPUT,
		NBSTEPS
	};
	enum 
	{ 
		Nb_Parameters =19
	};
		enum
	{ 
		Nb_Derivable_Parameters =17
	};
	
	static double value(const ArgumentList& a);
	static double value(int i,const ArgumentList& a, double s);
	static double value(int i,int j,const ArgumentList& a, double s1,double s2);

	static double specific_shift(int i) ;
	static double value(int i,const ArgumentList& a);
	static ArgumentList_Checking_Result check_argument(const ArgumentList& a);
	static ArgumentList_Checking_Result check_dimension(int rank);
};

///////////////////////////////////////////////////////////////////////
///  
///			class : ARM_CF_TriSpreadDigitalOption_Formula
///			Purpose : Evaluation of spread option cashflow= 1 if
///				 A0+A2*S2+A3*S3>=0 if B0+B1*S1>=0
///			Assumptions: Lognormal hypothesis
///
///////////////////////////////////////////////////////////////////////

struct ARM_CF_TriSpreadDigitalOption2_Formula
{
	enum ArgumentType
	{
		TIMETORESET,
		INDEX1,
		INDEX2,
		INDEX3,
		VOLATILITY1,
		VOLATILITY2,
		VOLATILITY3,
		CORRELATION12,
		CORRELATION13,
		CORRELATION23,
		MU1,
		MU2,
		MU3,
        A0,
		A2,
		A3,
		B0,
		B1,
		NBSTEPS
	};
	enum 
	{ 
		Nb_Parameters =19
	};
		enum
	{ 
		Nb_Derivable_Parameters =17
	};
	
	static double value(const ArgumentList& a);
	static double value(int i,const ArgumentList& a, double s);
	static double value(int i,int j,const ArgumentList& a, double s1,double s2);

	static double specific_shift(int i) ;
	static double value(int i,const ArgumentList& a);
	static ArgumentList_Checking_Result check_argument(const ArgumentList& a);
	static ArgumentList_Checking_Result check_dimension(int rank);
};



///////////////////////////////////////////////////////////////////////
///  
///			typedef : ARM_CF_Index1Paying_TriSpreadDigitalOption_Formula
///			Purpose : Evaluation of  cashflow call =  S1*1{a0+a1*S1+a2*S2+a3*S3>0}
///									 cashflow put  =  S1*1{a0+a1*S1+a2*S2+a3*S3<0}
///			Assumptions: Lognormal hypothesis
///			no smile taken into account
///
///////////////////////////////////////////////////////////////////////


typedef 	Product<Add_SubModel<Add_SubModel<Add_SubModel<ARM_CF_TriSpreadDigitalOption_Formula,
				Change_Numeraire1<		
					ARM_CF_TriSpreadDigitalOption_Formula::INDEX1,			// on modifie l'input S1		
					ARM_CF_TriSpreadDigitalOption_Formula::VOLATILITY1,	// multiplicativement par l'exp(sig1^2*T)
					ARM_CF_TriSpreadDigitalOption_Formula::TIMETORESET,
					1> 
				>,
				Change_Numeraire2<
					ARM_CF_TriSpreadDigitalOption_Formula::INDEX2,			// on modifie l'input S2
					ARM_CF_TriSpreadDigitalOption_Formula::CORRELATION12,	// multiplicativement par exp(sig1*sig2*corr12*T)
					ARM_CF_TriSpreadDigitalOption_Formula::VOLATILITY1,
					ARM_CF_TriSpreadDigitalOption_Formula::VOLATILITY2,
					ARM_CF_TriSpreadDigitalOption_Formula::TIMETORESET,
					1> 
				>,
				Change_Numeraire2<
					ARM_CF_TriSpreadDigitalOption_Formula::INDEX3,			// on modifie l'input S3
					ARM_CF_TriSpreadDigitalOption_Formula::CORRELATION13,	// multiplicativement par exp(sig1*sig3*corr13*T)
					ARM_CF_TriSpreadDigitalOption_Formula::VOLATILITY1,
					ARM_CF_TriSpreadDigitalOption_Formula::VOLATILITY3,
					ARM_CF_TriSpreadDigitalOption_Formula::TIMETORESET,
					1> 
				>,
				ARM_CF_TriSpreadDigitalOption_Formula::INDEX1>
								
	   ARM_CF_Index1Paying_TriSpreadDigitalOption_Formula;



///////////////////////////////////////////////////////////////////////
///  
///			typedef : ARM_CF_Index2Paying_TriSpreadDigitalOption_Formula
///			Purpose : Evaluation of  cashflow call =  S2*1{a0+a1*S1+a2*S2+a3*S3>0}
///									 cashflow put  =  S2*1{a0+a1*S1+a2*S2+a3*S3<0}
///			Assumptions: Lognormal hypothesis
///			no smile taken into account
///
///////////////////////////////////////////////////////////////////////



typedef 	Product<Add_SubModel<Add_SubModel<Add_SubModel<ARM_CF_TriSpreadDigitalOption_Formula,
				Change_Numeraire1<		
					ARM_CF_TriSpreadDigitalOption_Formula::INDEX2,			// on modifie l'input S2		
					ARM_CF_TriSpreadDigitalOption_Formula::VOLATILITY2,	// multiplicativement par l'exp(sig2^2*T)
					ARM_CF_TriSpreadDigitalOption_Formula::TIMETORESET,
					1> 
				>,
				Change_Numeraire2<
					ARM_CF_TriSpreadDigitalOption_Formula::INDEX1,			// on modifie l'input S1
					ARM_CF_TriSpreadDigitalOption_Formula::CORRELATION12,	// multiplicativement par exp(sig1*sig2*corr12*T)
					ARM_CF_TriSpreadDigitalOption_Formula::VOLATILITY1,
					ARM_CF_TriSpreadDigitalOption_Formula::VOLATILITY2,
					ARM_CF_TriSpreadDigitalOption_Formula::TIMETORESET,
					1> 
				>,
				Change_Numeraire2<
					ARM_CF_TriSpreadDigitalOption_Formula::INDEX3,			// on modifie l'input S3
					ARM_CF_TriSpreadDigitalOption_Formula::CORRELATION23,	// multiplicativement par exp(sig2*sig3*corr23*T)
					ARM_CF_TriSpreadDigitalOption_Formula::VOLATILITY2,
					ARM_CF_TriSpreadDigitalOption_Formula::VOLATILITY3,
					ARM_CF_TriSpreadDigitalOption_Formula::TIMETORESET,
					1> 
				>,
				ARM_CF_TriSpreadDigitalOption_Formula::INDEX1>
								
	   ARM_CF_Index2Paying_TriSpreadDigitalOption_Formula;

///////////////////////////////////////////////////////////////////////
///  
///			typedef : ARM_CF_Index3Paying_TriSpreadDigitalOption_Formula
///			Purpose : Evaluation of  cashflow call =  S3*1{a0+a1*S1+a2*S2+a3*S3>0}
///									 cashflow put  =  S3*1{a0+a1*S1+a2*S2+a3*S3<0}
///			Assumptions: Lognormal hypothesis
///			no smile taken into account
///
///////////////////////////////////////////////////////////////////////



typedef 	Product<Add_SubModel<Add_SubModel<Add_SubModel<ARM_CF_TriSpreadDigitalOption_Formula,
				Change_Numeraire1<		
					ARM_CF_TriSpreadDigitalOption_Formula::INDEX3,			// on modifie l'input S3		
					ARM_CF_TriSpreadDigitalOption_Formula::VOLATILITY3,	// multiplicativement par l'exp(sig3^2*T)
					ARM_CF_TriSpreadDigitalOption_Formula::TIMETORESET,
					1> 
				>,
				Change_Numeraire2<
					ARM_CF_TriSpreadDigitalOption_Formula::INDEX1,			// on modifie l'input S1
					ARM_CF_TriSpreadDigitalOption_Formula::CORRELATION13,	// multiplicativement par exp(sig1*sig3*corr13*T)
					ARM_CF_TriSpreadDigitalOption_Formula::VOLATILITY1,
					ARM_CF_TriSpreadDigitalOption_Formula::VOLATILITY3,
					ARM_CF_TriSpreadDigitalOption_Formula::TIMETORESET,
					1> 
				>,
				Change_Numeraire2<
					ARM_CF_TriSpreadDigitalOption_Formula::INDEX2,			// on modifie l'input S2
					ARM_CF_TriSpreadDigitalOption_Formula::CORRELATION23,	// multiplicativement par exp(sig2*sig3*corr23*T)
					ARM_CF_TriSpreadDigitalOption_Formula::VOLATILITY2,
					ARM_CF_TriSpreadDigitalOption_Formula::VOLATILITY3,
					ARM_CF_TriSpreadDigitalOption_Formula::TIMETORESET,
					1> 
				>,
				ARM_CF_TriSpreadDigitalOption_Formula::INDEX1>
								
	   ARM_CF_Index3Paying_TriSpreadDigitalOption_Formula;


///////////////////////////////////////////////////////////////////////
///  
///			typedef : ARM_CF_TriSpreadOption_Formula
///			Purpose : Evaluation of  cashflow call =  (a0+a1*S1+a2*S2+a3*S3)*1{a0+a1*S1+a2*S2+a3*S3>0}
///									 cashflow put  =  -(a0+a1*S1+a2*S2+a3*S3)*1{a0+a1*S1+a2*S2+a3*S3<0}
///			Assumptions: Lognormal hypothesis
///
///////////////////////////////////////////////////////////////////////

 
typedef		Product<																		/// The formula reflects the payoff
					Sum<																	/// it's a linear combinaison
						Sum<Product<
									ARM_CF_TriSpreadDigitalOption_Formula,
									ARM_CF_TriSpreadDigitalOption_Formula::A0>,
							Product<
									ARM_CF_Index1Paying_TriSpreadDigitalOption_Formula,
									ARM_CF_TriSpreadDigitalOption_Formula::A1>
							>,
						Sum<Product<
									ARM_CF_Index2Paying_TriSpreadDigitalOption_Formula,
									ARM_CF_TriSpreadDigitalOption_Formula::A2>,
							Product<
									ARM_CF_Index3Paying_TriSpreadDigitalOption_Formula,
									ARM_CF_TriSpreadDigitalOption_Formula::A3>
							>
						>,
					ARM_CF_TriSpreadDigitalOption_Formula::CALLORPUT>

			
			ARM_CF_TriSpreadOption_Formula;
		

///////////////////////////////////////////////////////////////////////
///  
///			typedef : ARM_CF_Index1Paying_TriSpreadDigitalOption2_Formula
///			Purpose : Evaluation of  cashflow call =  S1*1{a0+a2*S2+a3*S3>0}*1{b0+a1*S1>0}
///									
///			Assumptions: Lognormal hypothesis
///			no smile taken into account
///
///////////////////////////////////////////////////////////////////////


typedef 	Product<Add_SubModel<Add_SubModel<Add_SubModel<ARM_CF_TriSpreadDigitalOption2_Formula,
				Change_Numeraire1<		
					ARM_CF_TriSpreadDigitalOption2_Formula::INDEX1,			// on modifie l'input S1		
					ARM_CF_TriSpreadDigitalOption2_Formula::VOLATILITY1,	// multiplicativement par l'exp(sig1^2*T)
					ARM_CF_TriSpreadDigitalOption2_Formula::TIMETORESET,
					1> 
				>,
				Change_Numeraire2<
					ARM_CF_TriSpreadDigitalOption2_Formula::INDEX2,			// on modifie l'input S2
					ARM_CF_TriSpreadDigitalOption2_Formula::CORRELATION12,	// multiplicativement par exp(sig1*sig2*corr12*T)
					ARM_CF_TriSpreadDigitalOption2_Formula::VOLATILITY1,
					ARM_CF_TriSpreadDigitalOption2_Formula::VOLATILITY2,
					ARM_CF_TriSpreadDigitalOption2_Formula::TIMETORESET,
					1> 
				>,
				Change_Numeraire2<
					ARM_CF_TriSpreadDigitalOption2_Formula::INDEX3,			// on modifie l'input S3
					ARM_CF_TriSpreadDigitalOption2_Formula::CORRELATION13,	// multiplicativement par exp(sig1*sig3*corr13*T)
					ARM_CF_TriSpreadDigitalOption2_Formula::VOLATILITY1,
					ARM_CF_TriSpreadDigitalOption2_Formula::VOLATILITY3,
					ARM_CF_TriSpreadDigitalOption2_Formula::TIMETORESET,
					1> 
				>,
				ARM_CF_TriSpreadDigitalOption2_Formula::INDEX1>
								
	   ARM_CF_Index1Paying_TriSpreadDigitalOption2_Formula;



///////////////////////////////////////////////////////////////////////
///  
///			typedef : ARM_CF_Index2Paying_TriSpreadDigitalOption_Formula
///			Purpose : Evaluation of  cashflow call =  S2*1{a0+a2*S2+a3*S3>0}*1{b0+a1*S1>0}
///									 
///			Assumptions: Lognormal hypothesis
///			no smile taken into account
///
///////////////////////////////////////////////////////////////////////



typedef 	Product<Add_SubModel<Add_SubModel<Add_SubModel<ARM_CF_TriSpreadDigitalOption2_Formula,
				Change_Numeraire1<		
					ARM_CF_TriSpreadDigitalOption2_Formula::INDEX2,			// on modifie l'input S2		
					ARM_CF_TriSpreadDigitalOption2_Formula::VOLATILITY2,	// multiplicativement par l'exp(sig2^2*T)
					ARM_CF_TriSpreadDigitalOption2_Formula::TIMETORESET,
					1> 
				>,
				Change_Numeraire2<
					ARM_CF_TriSpreadDigitalOption2_Formula::INDEX1,			// on modifie l'input S1
					ARM_CF_TriSpreadDigitalOption2_Formula::CORRELATION12,	// multiplicativement par exp(sig1*sig2*corr12*T)
					ARM_CF_TriSpreadDigitalOption2_Formula::VOLATILITY1,
					ARM_CF_TriSpreadDigitalOption2_Formula::VOLATILITY2,
					ARM_CF_TriSpreadDigitalOption2_Formula::TIMETORESET,
					1> 
				>,
				Change_Numeraire2<
					ARM_CF_TriSpreadDigitalOption2_Formula::INDEX3,			// on modifie l'input S3
					ARM_CF_TriSpreadDigitalOption2_Formula::CORRELATION23,	// multiplicativement par exp(sig2*sig3*corr23*T)
					ARM_CF_TriSpreadDigitalOption2_Formula::VOLATILITY2,
					ARM_CF_TriSpreadDigitalOption2_Formula::VOLATILITY3,
					ARM_CF_TriSpreadDigitalOption2_Formula::TIMETORESET,
					1> 
				>,
				ARM_CF_TriSpreadDigitalOption2_Formula::INDEX1>
								
	   ARM_CF_Index2Paying_TriSpreadDigitalOption2_Formula;

///////////////////////////////////////////////////////////////////////
///  
///			typedef : ARM_CF_Index3Paying_TriSpreadDigitalOption2_Formula
///			Purpose : Evaluation of  cashflow call =  S3*1{a0+a2*S2+a3*S3>0}*1{b0+a1*S1>0}
///								
///			Assumptions: Lognormal hypothesis
///			no smile taken into account
///
///////////////////////////////////////////////////////////////////////



typedef 	Product<Add_SubModel<Add_SubModel<Add_SubModel<ARM_CF_TriSpreadDigitalOption2_Formula,
				Change_Numeraire1<		
					ARM_CF_TriSpreadDigitalOption2_Formula::INDEX3,			// on modifie l'input S3		
					ARM_CF_TriSpreadDigitalOption2_Formula::VOLATILITY3,	// multiplicativement par l'exp(sig3^2*T)
					ARM_CF_TriSpreadDigitalOption2_Formula::TIMETORESET,
					1> 
				>,
				Change_Numeraire2<
					ARM_CF_TriSpreadDigitalOption2_Formula::INDEX1,			// on modifie l'input S1
					ARM_CF_TriSpreadDigitalOption2_Formula::CORRELATION13,	// multiplicativement par exp(sig1*sig3*corr13*T)
					ARM_CF_TriSpreadDigitalOption2_Formula::VOLATILITY1,
					ARM_CF_TriSpreadDigitalOption2_Formula::VOLATILITY3,
					ARM_CF_TriSpreadDigitalOption2_Formula::TIMETORESET,
					1> 
				>,
				Change_Numeraire2<
					ARM_CF_TriSpreadDigitalOption2_Formula::INDEX2,			// on modifie l'input S2
					ARM_CF_TriSpreadDigitalOption2_Formula::CORRELATION23,	// multiplicativement par exp(sig2*sig3*corr23*T)
					ARM_CF_TriSpreadDigitalOption2_Formula::VOLATILITY2,
					ARM_CF_TriSpreadDigitalOption2_Formula::VOLATILITY3,
					ARM_CF_TriSpreadDigitalOption2_Formula::TIMETORESET,
					1> 
				>,
				ARM_CF_TriSpreadDigitalOption2_Formula::INDEX1>
								
	   ARM_CF_Index3Paying_TriSpreadDigitalOption2_Formula;

///////////////////////////////////////////////////////////////////////
///  
///			typedef : ARM_CF_TriSpreadOption2_Formula
///			Purpose : Evaluation of  cashflow  =  (a0+a2*S2+a3*S3)*1{a0+a2*S2+a3*S3>0}*1{b0+a1*S1>0}
///									 
///			Assumptions: Lognormal hypothesis
///
///////////////////////////////////////////////////////////////////////

 
typedef																				/// The formula reflects the payoff
					Sum3<																	/// it's a linear combinaison
						Product<
									ARM_CF_TriSpreadDigitalOption2_Formula,
									ARM_CF_TriSpreadDigitalOption2_Formula::A0>,
						Product<
									ARM_CF_Index2Paying_TriSpreadDigitalOption2_Formula,
									ARM_CF_TriSpreadDigitalOption2_Formula::A2>,
						Product<
									ARM_CF_Index3Paying_TriSpreadDigitalOption2_Formula,
									ARM_CF_TriSpreadDigitalOption2_Formula::A3>
						>
			
			ARM_CF_TriSpreadOption2_Formula;

///////////////////////////////////////////////////////////////////////
///  
///			typedef : ARM_CF_TriSpreadOption3_Formula
///			Purpose : Evaluation of  cashflow  =  (b0+b1*S1)*1{a0+a2*S2+a3*S3>0}*1{b0+a1*S1>0}
///									 
///			Assumptions: Lognormal hypothesis
///
///////////////////////////////////////////////////////////////////////

 
typedef																				/// The formula reflects the payoff
					Sum<																	/// it's a linear combinaison
						Product<
									ARM_CF_TriSpreadDigitalOption2_Formula,
									ARM_CF_TriSpreadDigitalOption2_Formula::B0>,
						Product<
									ARM_CF_Index1Paying_TriSpreadDigitalOption2_Formula,
									ARM_CF_TriSpreadDigitalOption2_Formula::B1>
						>
			
			ARM_CF_TriSpreadOption3_Formula;

CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/


