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
 
#ifndef _GP_CF_PRODUCT_H
#define _GP_CF_PRODUCT_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"

#include "basic_distributions.h"
#include "input_extender.h"

CC_BEGIN_NAMESPACE(ARM)



///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Product : extender class
/// Formula : V<-F*V
///////////////////////////////////////////////////////////////////////////////////////
template<class M,int i_forw>
class Product 
{
public:		
	enum
	{
		Nb_Parameters=M::Nb_Parameters
	};
	static double value (const ArgumentList& a);
	static double value(int i,const  ArgumentList& a, double s);
	static double value(int i, int j,const  ArgumentList& a, double s1,double s2);
	static int affected_index();
	static double specific_shift(int i);
	static ArgumentList_Checking_Result check_argument(const ArgumentList& a);
	static ArgumentList_Checking_Result check_dimension(int i);
};

///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Product
///	Routine: value
/// Argument: ArgumentList * a 
///	Returns: ArgumentList *
///	Action :  model evaluation
///////////////////////////////////////////////////////////////////////////////////////

template<class M,int i_forw>
double Product<M,i_forw>::value(const ArgumentList& a)
{
	double val= a[i_forw]*M::value(a);				;
	return val;
}



template <class M ,int i_forw>
double Product<M,i_forw>::value(int i,const ArgumentList& a, double s)
{

	double val;
		switch (i)
	{
	case i_forw : 
		{
			val=a[i_forw]* M::value(i,a,s)+M::value(a);
			break;
		}
	default :
		{
			val=a[i_forw]* M::value(i,a,s);
			break;

		}
	}
	return val;

}

template <class M ,int i_forw>
double Product<M,i_forw>::value(int i,int j,const ArgumentList& a, double s1,double s2)
{
	double val;
	switch (i)
	{
	case i_forw : 
		{
			switch (j)
			{
			case i_forw : 
				{
					val=2*M::value(i,a,s1)+a[i_forw]* M::value(i,i,a,s1,s1);
					break;
				}
			default :
				{
					val=a[i_forw]* M::value(i,j,a,s1,s2)+M::value(j,a,s2);
					break;		
				}
			}
			break;
		}
	default :
		{
			switch (j)
			{
			case i_forw : 
				{
					val=a[i_forw]* M::value(i,j,a,s1,s2)+M::value(i,a,s1);
					break;
				}
			default :
				{
					val=a[i_forw]* M::value(i,j,a,s1,s2);
					break;
				}
			}
			break;
			
			
		}
	}
	return val;
}



template <class M ,int i_forw>
double Product<M,i_forw>::specific_shift(int i)
{
	if (i<M::Nb_Parameters)
	{
		return M::specific_shift(i);
	}
	else
	{
		return ARM_CF_DEFAULT_SHIFT;
	}
}

template<class M,int i_forw>
ArgumentList_Checking_Result Product<M,i_forw>::check_argument(const ArgumentList& a)
{
	return M::check_argument(a);
}

template<class M,int i_forw>
ArgumentList_Checking_Result Product<M,i_forw>::check_dimension(int rank)
{
	return M::check_dimension(rank);
}


CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

