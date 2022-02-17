/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file difference.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_DIFFERENCE_H
#define _GP_CF_DIFFERENCE_H

#include "firsttoinc.h"
#include "gpbase/port.h"



CC_BEGIN_NAMESPACE(ARM)



///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Difference : extender class
/// Formula : V<-V1+V2
///////////////////////////////////////////////////////////////////////////////////////
template<class V1,class V2>
class Difference 
{
public:	
	enum
	{
		Nb_Parameters=V1::Nb_Parameters
	};
	static double value (const ArgumentList& a);
	static double value(int i,const ArgumentList& a, double s);
	static double value(int i, int j,const ArgumentList& a, double s1,double s2);
	static int affected_index();
	static double specific_shift(int i);
	static ArgumentList_Checking_Result check_argument(const ArgumentList& a);
	static ArgumentList_Checking_Result check_dimension(int i);
};



///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Difference
///	Routine: value
/// Argument: ArgumentList * a 
///	Returns: ArgumentList *
///	Action :  model evaluation
///////////////////////////////////////////////////////////////////////////////////////

template<class V1,class V2>
double Difference<V1,V2>::value(const ArgumentList& a)
{
	double val= V1::value(a)-V2::value(a);				;
	return val;
}




///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Difference
///	Routine: value
/// Argument:  int i, ArgumentList * a, double s
///	Returns: double
///	Action :  first derivative with respect to the ith argument
///////////////////////////////////////////////////////////////////////////////////////

template<class V1,class V2>
double Difference<V1,V2>::value(int i,const ArgumentList& a, double s)
{
	return V1::value(i,a,s)-V2::value(i,a,s);
}


///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Difference
///	Routine: value
/// Argument:  int i, int j, ArgumentList * a, double s
///	Returns: double
///	Action :  second derivative with respect to the ith argument and the jth argument
///////////////////////////////////////////////////////////////////////////////////////

template<class V1,class V2>
double Difference<V1,V2>::value(int i,int j,const ArgumentList& a, double s1,double s2)
{
	return V1::value(i,j,a,s1,s2)-V2::value(i,j,a,s1,s2);
}



///////////////////////////////////////////////////////////////////////////////////////
///	Class  : SDifferenceum
///	Routine: specific_shift
/// Argument:  int i
///	Returns: double
///	Action :  specific shift to compute the derivatives with respect to the ith dimension
///////////////////////////////////////////////////////////////////////////////////////

template<class V1,class V2>
double Difference<V1,V2>::specific_shift(int i)
{
		return V1::specific_shift(i);
}

///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Difference
///	Routine: check_argument
/// Argument:  ArgumentList * a
///	Returns: ArgumentList_Checking_Result
///	Action :  check the range of the argument list
///////////////////////////////////////////////////////////////////////////////////////
template<class V1,class V2>
ArgumentList_Checking_Result Difference<V1,V2>::check_argument(const ArgumentList& a)
{
	return (V1::check_argument(a))&&(V2::check_argument(a));
}

///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Difference
///	Routine: check_dimension
/// Argument:  ArgumentList * a
///	Returns: ArgumentList_Checking_Result
///	Action :  check the range of the dimension
///////////////////////////////////////////////////////////////////////////////////////


 template<class V1,class V2>
ArgumentList_Checking_Result Difference<V1,V2>::check_dimension(int rank)
{
	return (V1::check_dimension(rank))&&(V2::check_dimension(rank));
}

CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

