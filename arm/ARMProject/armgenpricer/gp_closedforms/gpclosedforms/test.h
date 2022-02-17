/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file product.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _ARM_CF_TEST_H
#define _ARM_CF_TEST_H

#include "gpbase/port.h"
#include <vector>
#include "product.h"
#include "sum.h"
#include "difference.h"
#include "add_submodel.h"
#include "change_numeraire.h"
#include "spreadoption_lognormal.h"



CC_BEGIN_NAMESPACE(ARM)

template<int RANK1,int RANK2,int RTOTAL>
class TestFormula 
{
public:		
	static double value (ArgumentList * a);
	static double value(int i, ArgumentList * a, double s);
	static double value(int i, int j, ArgumentList * a, double s1,double s2);
	static int nb_arguments();
	static double specific_shift(int i);
	static ArgumentList_Checking_Result check_argument(ArgumentList * a);
};

template<int RANK1,int RANK2,int RTOTAL>
ArgumentList_Checking_Result TestFormula<RANK1,RANK2,RTOTAL>::check_argument(ArgumentList * a)
{
	return ArgumentList_Checking_Result(true,"");
}

template<int RANK1,int RANK2,int RTOTAL>
double TestFormula<RANK1,RANK2,RTOTAL>::specific_shift(int i)
{
	return 0.000001;
}

template<int RANK1,int RANK2,int RTOTAL>
int TestFormula<RANK1,RANK2,RTOTAL>::nb_arguments()
{
	return RTOTAL;
}

template<int RANK1,int RANK2,int RTOTAL>
double TestFormula<RANK1,RANK2,RTOTAL>::value(ArgumentList * a)
{
	double val=(*a)[RANK1]*(*a)[RANK2] ;				;
	return val;
}


template<int RANK1,int RANK2,int RTOTAL>
double TestFormula<RANK1,RANK2,RTOTAL>::value(int i, ArgumentList * a,double s)
{
	double val;
	switch (i)
	{
	case RANK1 :
		{
			if (RANK1==RANK2) val=2.*(*a)[RANK2];
			else val=(*a)[RANK2];
		}
		break;
	default :
		{
			val=0.;
		}
		break;
	}
	return val;
}
	

template<int RANK1,int RANK2,int RTOTAL>
double TestFormula<RANK1,RANK2,RTOTAL>::value(int i,int j, ArgumentList * a,double s1,double s2)
{
	double val;
	switch (i)
	{
	case RANK1 :
		{
			switch (j)
			{
			case RANK2 :
				{
					if (RANK1==RANK2) val=2.;
					else val=1.;
				}
				break;
			default :
				{
					val=0.;
				}
				break;
			}
		}
		break;
	case RANK2 :
		{
			switch (j)
			{
			case RANK1 :
				{
					if (RANK1==RANK2) val=2.;
					else val=1.;
				}
				break;
			default :
				{
					val=0.;
				}
				break;
			}
		}
		break;
	default :
		{
			val=0.;
		}
		break;
	}
	return val;
}


typedef		Add_SubModel<TestFormula<
						ARM_CF_SpreadDigitalOption_Formula::INDEX1,
						ARM_CF_SpreadDigitalOption_Formula::VOLATILITY1,
						8>,
				Change_Numeraire1<
					ARM_CF_SpreadDigitalOption_Formula::INDEX1,
					ARM_CF_SpreadDigitalOption_Formula::VOLATILITY1,
					ARM_CF_SpreadDigitalOption_Formula::TIMETORESET,
					1>  >
			Test_Formula_1;

typedef		Add_SubModel<TestFormula<
						ARM_CF_SpreadDigitalOption_Formula::INDEX1,
						ARM_CF_SpreadDigitalOption_Formula::VOLATILITY1,
						8>,
				Change_Numeraire2<
					ARM_CF_SpreadDigitalOption_Formula::INDEX1,
					ARM_CF_SpreadDigitalOption_Formula::CORRELATION,
					ARM_CF_SpreadDigitalOption_Formula::VOLATILITY1,
					ARM_CF_SpreadDigitalOption_Formula::VOLATILITY2,
					ARM_CF_SpreadDigitalOption_Formula::TIMETORESET,
					1>   >
			Test_Formula_2;


CC_END_NAMESPACE()

#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

