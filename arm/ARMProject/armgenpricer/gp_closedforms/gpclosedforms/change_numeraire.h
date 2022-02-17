/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file change_numeraire.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_CHANGE_NUMERAIRE_H
#define _GP_CF_CHANGE_NUMERAIRE_H

#include "firsttoinc.h"
#include "gpbase/port.h"

#include <cmath>

CC_BEGIN_NAMESPACE(ARM)

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
///
///     Log normal Assumption for the underlyings 
///
///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Change_Numeraire1 : extender class
/// Formula : F1=exp(sig1*sig1/2)*F
///////////////////////////////////////////////////////////////////////////////////////
template<int i_forw, int i_vol1,int i_t,int mul>
class Change_Numeraire1 
{
public:
	static	ArgumentList * prepare_model_input(const ArgumentList& a,double submodel_result);		
	static double value (const ArgumentList& a);
	static double value(int i,const ArgumentList& a, double s);
	static double value(int i, int j,const ArgumentList& a, double s1,double s2);
	static int affected_index();
	static ArgumentList_Checking_Result check_argument(const ArgumentList& a);
	static ArgumentList_Checking_Result check_dimension(int i);

};

///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Change_Numeraire2 : extender class
/// Formula : F1=exp(rho*sig1*sig2)*F
///////////////////////////////////////////////////////////////////////////////////////
template<int i_forw, int i_vol1,int i_vol2, int i_corr,int i_t,int mul>
class Change_Numeraire2 
{
public:
	static	ArgumentList * prepare_model_input(const ArgumentList& a,double submodel_result);		
	static double value (const ArgumentList& a);
	static double value(int i,const ArgumentList& a, double s);
	static double value(int i, int j,const ArgumentList& a, double s1,double s2);
	static int affected_index();
	static ArgumentList_Checking_Result check_argument(const ArgumentList& a);
	static ArgumentList_Checking_Result check_dimension(int i);

};


///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Change_Numeraire1 
///	Routine: affected_index
/// Argument: 
///	Returns: int
///	Action : return the index in the arglist which carries the foward affected by 
/// the change of numeraire
///////////////////////////////////////////////////////////////////////////////////////
template<int i_forw, int i_vol1,int i_t,int mul>
inline int Change_Numeraire1<i_forw,i_vol1,i_t,mul>::affected_index()
{
	return i_forw;
}

///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Change_Numeraire2 
///	Routine: affected_index
/// Argument: 
///	Returns: int
///	Action : return the index in the arglist which carries the foward affected by 
/// the change of numeraire
///////////////////////////////////////////////////////////////////////////////////////
template<int i_forw, int i_vol1,int i_vol2, int i_corr,int i_t,int mul>
inline int Change_Numeraire2<i_forw,i_vol1,i_vol2,i_corr,i_t,mul>::affected_index()
{
	return i_forw;
}


///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Change_Numeraire1
///	Routine: value
/// Argument: ArgumentList * a ,double submodel_result
///	Returns: ArgumentList *
///	Action :  model evaluation
///////////////////////////////////////////////////////////////////////////////////////

template<int i_forw, int i_vol1,int i_t,int mul>
double Change_Numeraire1<i_forw,i_vol1,i_t,mul>::value(const ArgumentList& a)
{
	double val= (a[i_forw])*
		exp(a[i_vol1]*a[i_vol1]*a[i_t]*mul);
	return val;
}

///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Change_Numeraire2
///	Routine: value
/// Argument: ArgumentList * a ,double submodel_result
///	Returns: ArgumentList *
///	Action :  model evaluation
///////////////////////////////////////////////////////////////////////////////////////

template<int i_forw, int i_vol1,int i_vol2, int i_corr,int i_t,int mul>
double Change_Numeraire2<i_forw,i_vol1,i_vol2,i_corr,i_t,mul>::value(const ArgumentList& a)
{
	double val= (a[i_forw])*
		exp(a[i_vol1]*a[i_vol2]*a[i_corr]*a[i_t]*mul);
	return val;
}

///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Change_Numeraire1
///	Routine: value
/// Argument: ArgumentList * a ,double submodel_result
///	Returns: int i, ArgumentList * a, double s
///	Action : returns the derivative of the computed value with respect to the ith input
///////////////////////////////////////////////////////////////////////////////////////

template<int i_forw, int i_vol1, int i_t,int mul>
double Change_Numeraire1<i_forw,i_vol1,i_t,mul>::value(int i,const ArgumentList& a, double s)
{
	double val;
	double prodall=a[i_vol1]*a[i_vol1]*a[i_t]*mul;
	switch (i)
	{
	case i_t : 
		{
			val=mul*	a[i_vol1]*
						a[i_vol1]*
						a[i_forw];
			break;
		}
	case i_forw : 
		{
			val=1;
			break;
		}
	case i_vol1 : 
		{
			val=2.0*mul*a[i_vol1]*
						a[i_forw]*
						a[i_t];
			break;
		}
	default :
		{
			val=0.0;
			break;
			
		}
	}
	return val*exp(prodall);
}

///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Change_Numeraire2
///	Routine: value
/// Argument: ArgumentList * a ,double submodel_result
///	Returns: int i, ArgumentList * a, double s
///	Action : returns the derivative of the computed value with respect to the ith input
///////////////////////////////////////////////////////////////////////////////////////

template<int i_forw, int i_vol1,int i_vol2, int i_corr,int i_t,int mul>
double Change_Numeraire2<i_forw,i_vol1,i_vol2,i_corr,i_t,mul>::value(int i,const ArgumentList& a, double s)
{
	double val;
	double prodall=a[i_vol1]*a[i_vol2]*a[i_corr]*a[i_t]*mul;
	switch (i)
	{
	case i_t : 
	case i_vol1 : 
	case i_vol2 : 
	case i_corr : 
		{
			val=prodall/a[i]*
				a[i_forw];
			break;
		}
	case i_forw : 
		{
			val=1;
			break;
		}
	default :
		{
			val=0.0;
			break;
			
		}
	}
	return val*exp(prodall);
}


///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Change_Numeraire1 
///	Routine: value
/// Argument: ArgumentList * a ,double submodel_result
///	Returns: int i, int j,ArgumentList * a, double s
///	Action : returns the 2nd derivative of the computed value 
/// with respect to the ith input  and jth input
///////////////////////////////////////////////////////////////////////////////////////

template<int i_forw, int i_vol1, int i_t,int mul>
double Change_Numeraire1<i_forw,i_vol1,i_t,mul>::value(int i,int j,const ArgumentList& a, double s1,double s2)
{
	double val;
	double prodall=a[i_vol1]*a[i_vol1]*a[i_t]*mul;
	switch (i)	{
	case i_t :
		{
			switch (j)	{
			case i_t : 
				val=prodall*prodall/(a[i_t]*a[i_t])*a[i_forw];
				break;
			case i_forw : 
				val=prodall/a[i_t];
				break;
			case i_vol1 : 
				val=mul*2.0*(1+	prodall)*a[i_vol1]*a[i_forw];
				break;
			default :
				val=0.0;
				break;
			}
		}
		break;
	case i_forw : 
		{
			switch (j)	{
			case i_t : 
				val=prodall/a[i_t];
				break;
			case i_forw : 
				val=0.0;
				break;
			case i_vol1 : 
				val=2.0*mul*a[i_vol1]*a[i_t];
				break;
			default :
				val=0.0;
				break;
			}
		}
		break;
	case i_vol1 : 
		{
			switch (j)	{
			case i_t : 
				val=2.0*mul*a[i_forw]*(1+prodall)*a[i_vol1];
				break;
			case i_forw : 
				val=2.0*mul*a[i_vol1]*a[i_t];
				break;
			case i_vol1 : 
				val=2.0*mul*a[i_forw]*(1+2.0*prodall)*a[i_t];
				break;
			default :
				val=0.0;
				break;
			}
		}
		break;
	default :
		val=0.0;
		break;
	}
	return val*exp(prodall);
}
///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Change_Numeraire2
///	Routine: value
/// Argument: ArgumentList * a ,double submodel_result
///	Returns: int i, int j,ArgumentList * a, double s
///	Action : returns the 2nd derivative of the computed value 
/// with respect to the ith input  and jth input
///////////////////////////////////////////////////////////////////////////////////////

template<int i_forw, int i_vol1,int i_vol2, int i_corr,int i_t,int mul>
double Change_Numeraire2<i_forw,i_vol1,i_vol2,i_corr,i_t,mul>::value(int i,int j,const ArgumentList& a, double s1,double s2)
{
	double val;
	double prodall=a[i_vol1]*a[i_vol2]*a[i_corr]*a[i_t]*mul;
	switch (i)	{
	case i_t :
		{
			switch (j)	{
			case i_t : 
				val=prodall*prodall/(a[j]*a[j])*a[i_forw];
				break;
			case i_forw : 
				val=mul*prodall/a[i];
				break;
			case i_vol1 : 
			case i_vol2 : 
			case i_corr : 	
				val=mul*(1+	prodall)*prodall*a[i_forw]/(a[j]*a[i]);
				break;
			default :
				val=0.0;
				break;	
			}
			break;
		}
	case i_forw : 
		{
			switch (j)	{
			case i_forw : 
				val=0.0;
				break;
			case i_t : 
			case i_vol1 : 
			case i_vol2 :
			case i_corr : 
				val=mul*prodall/a[j];
				break;
			default :
				val=0.0;
				break;
			}
			break;
		}
	case i_vol1 : 
		{
			switch (j)	{
			case i_t : 
			case i_vol2 : 
			case i_corr : 
				val=mul*(1+	prodall)*prodall*a[i_forw]/(a[j]*a[i]);
				break;
			case i_forw : 
				val=mul*prodall/a[i];
				break;
			case i_vol1 : 
				val=prodall*prodall/(a[j]*a[j])*a[i_forw];
				break;
			default :
				val=0.0;
				break;
			}
			break;
		}
	case i_vol2 : 
		{
			switch (j)	{
			case i_t : 
			case i_vol1 :
			case i_corr : 
				val=mul*(1+	prodall)*prodall*a[i_forw]/(a[j]*a[i]);
				break;
			case i_forw : 
				val=mul*prodall/a[i];
				break;
			case i_vol2 : 
				val=prodall*prodall/(a[j]*a[j])*a[i_forw];
				break;
			default :
				val=0.0;
				break;
			}
			break;
		}
	case i_corr : 
		{
			switch (j)	{
			case i_t : 
			case i_vol1 : 
			case i_vol2 : 
				val=mul*(1+	prodall)*prodall*a[i_forw]/(a[j]*a[i]);
				break;
			case i_forw : 
				val=mul*prodall/a[i];
				break;
			case i_corr : 
				val=prodall*prodall/(a[j]*a[j])*a[i_forw];
				break;
			default :
				val=0.0;
				break;
			}
			break;
			default :
				val=0.0;
				break;
		}
		break;
	}
	return val*exp(prodall);
}

///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Change_Numeraire1 
///	Routine: prepare_model_input
/// Argument: ArgumentList * a ,double submodel_result
///	Returns: ArgumentList *
///	Action : prepare the input for a model evaluation
///////////////////////////////////////////////////////////////////////////////////////

template<int i_forw, int i_vol1, int i_t,int mul>
ArgumentList * Change_Numeraire1<i_forw,i_vol1,i_t,mul>::prepare_model_input(const ArgumentList& a,double submodel_result)
{
	ArgumentList * newarg=new ArgumentList(a);
	newarg->set_nth(i_forw,submodel_result);
	return newarg;
}
///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Change_Numeraire2
///	Routine: prepare_model_input
/// Argument: ArgumentList * a ,double submodel_result
///	Returns: ArgumentList *
///	Action : prepare the input for a model evaluation
///////////////////////////////////////////////////////////////////////////////////////

template<int i_forw, int i_vol1,int i_vol2, int i_corr,int i_t,int mul>
ArgumentList * Change_Numeraire2<i_forw,i_vol1,i_vol2,i_corr,i_t,mul>::prepare_model_input(const ArgumentList& a,double submodel_result)
{
	ArgumentList * newarg=new ArgumentList(a);
	newarg->set_nth(i_forw,submodel_result);
	return newarg;
}


///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Change_Numeraire1
///	Routine: check_argument
/// Argument: ArgumentList * a 
///	Returns: ArgumentList_Checking_Result 
///	Action : for now does no checking further on
///////////////////////////////////////////////////////////////////////////////////////

template<int i_forw, int i_vol1, int i_t,int mul>
ArgumentList_Checking_Result Change_Numeraire1<i_forw,i_vol1,i_t,mul>::check_argument(const ArgumentList& a)
{
	return ArgumentList_Checking_Result(true,"");
}

///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Change_Numeraire1
///	Routine: check_argument
/// Argument: ArgumentList * a 
///	Returns: ArgumentList_Checking_Result 
///	Action : for now does no checking further on
///////////////////////////////////////////////////////////////////////////////////////

template<int i_forw, int i_vol1,int i_vol2, int i_corr,int i_t,int mul>
ArgumentList_Checking_Result Change_Numeraire2<i_forw,i_vol1,i_vol2,i_corr,i_t,mul>::check_argument(const ArgumentList& a)
{
	return ArgumentList_Checking_Result(true,"");
}

///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Change_Numeraire1
///	Routine: check_dimension
/// Argument: int i 
///	Returns: ArgumentList_Checking_Result 
///	Action : for now does no checking further on
///////////////////////////////////////////////////////////////////////////////////////

template<int i_forw, int i_vol1, int i_t,int mul>
ArgumentList_Checking_Result Change_Numeraire1<i_forw,i_vol1,i_t,mul>::check_dimension(int i)
{
	return ArgumentList_Checking_Result(true,"");
}

///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Change_Numeraire2
///	Routine: check_dimension
/// Argument: int i
///	Returns: ArgumentList_Checking_Result
///	Action : for now does no checking further on
///////////////////////////////////////////////////////////////////////////////////////

template<int i_forw, int i_vol1,int i_vol2, int i_corr,int i_t,int mul>
ArgumentList_Checking_Result Change_Numeraire2<i_forw,i_vol1,i_vol2,i_corr,i_t,mul>::check_dimension(int i)
{
	return ArgumentList_Checking_Result(true,"");
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///     Normal Assumption for the underlyings 
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Change_Numeraire1_N : extender class
/// Formula : F1=exp(sig1*sig1/2)*F
///////////////////////////////////////////////////////////////////////////////////////
template<int i_forw, int i_vol1,int i_t,int mul>
class Change_Numeraire1_N 
{
public:
	static	ArgumentList * prepare_model_input(const ArgumentList& a,double submodel_result);		
	static double value (const ArgumentList& a);
	static double value(int i,const ArgumentList& a, double s);
	static double value(int i, int j,const ArgumentList& a, double s1,double s2);
	static int affected_index();
	static ArgumentList_Checking_Result check_argument(const ArgumentList& a);
	static ArgumentList_Checking_Result check_dimension(int i);

};

///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Change_Numeraire2_N : extender class
/// Formula : F1=exp(rho*sig1*sig2)*F
///////////////////////////////////////////////////////////////////////////////////////
template<int i_forw, int i_vol1,int i_vol2, int i_corr,int i_t,int mul>
class Change_Numeraire2_N 
{
public:
	static	ArgumentList * prepare_model_input(const ArgumentList& a,double submodel_result);		
	static double value (const ArgumentList& a);
	static double value(int i,const ArgumentList& a, double s);
	static double value(int i, int j,const ArgumentList& a, double s1,double s2);
	static int affected_index();
	static ArgumentList_Checking_Result check_argument(const ArgumentList& a);
	static ArgumentList_Checking_Result check_dimension(int i);

};


///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Change_Numeraire1_N 
///	Routine: affected_index
/// Argument: 
///	Returns: int
///	Action : return the index in the arglist which carries the foward affected by 
/// the change of numeraire
///////////////////////////////////////////////////////////////////////////////////////
template<int i_forw, int i_vol1,int i_t,int mul>
inline int Change_Numeraire1_N<i_forw,i_vol1,i_t,mul>::affected_index()
{
	return i_forw;
}

///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Change_Numeraire2_N 
///	Routine: affected_index
/// Argument: 
///	Returns: int
///	Action : return the index in the arglist which carries the foward affected by 
/// the change of numeraire
///////////////////////////////////////////////////////////////////////////////////////
template<int i_forw, int i_vol1,int i_vol2, int i_corr,int i_t,int mul>
inline int Change_Numeraire2_N<i_forw,i_vol1,i_vol2,i_corr,i_t,mul>::affected_index()
{
	return i_forw;
}


///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Change_Numeraire1_N
///	Routine: value
/// Argument: ArgumentList * a ,double submodel_result
///	Returns: ArgumentList *
///	Action :  model evaluation
///////////////////////////////////////////////////////////////////////////////////////

template<int i_forw, int i_vol1,int i_t,int mul>
double Change_Numeraire1_N<i_forw,i_vol1,i_t,mul>::value(const ArgumentList& a)
{
	double val= (a[i_forw])*
		(1.+a[i_vol1]*a[i_vol1]*a[i_t]*mul);
	return val;
}

///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Change_Numeraire2_N
///	Routine: value
/// Argument: ArgumentList * a ,double submodel_result
///	Returns: ArgumentList *
///	Action :  model evaluation
///////////////////////////////////////////////////////////////////////////////////////

template<int i_forw, int i_vol1,int i_vol2, int i_corr,int i_t,int mul>
double Change_Numeraire2_N<i_forw,i_vol1,i_vol2,i_corr,i_t,mul>::value(const ArgumentList& a)
{
	double val= (a[i_forw])*
		(1.+a[i_vol1]*a[i_vol2]*a[i_corr]*a[i_t]*mul);
	return val;
}

///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Change_Numeraire1_N
///	Routine: value
/// Argument: ArgumentList * a ,double submodel_result
///	Returns: int i, ArgumentList * a, double s
///	Action : returns the derivative of the computed value with respect to the ith input
///////////////////////////////////////////////////////////////////////////////////////

template<int i_forw, int i_vol1, int i_t,int mul>
double Change_Numeraire1_N<i_forw,i_vol1,i_t,mul>::value(int i,const ArgumentList& a, double s)
{
	double val;
	switch (i)
	{
	case i_t : 
		{
			val=mul*	a[i_vol1]*
						a[i_vol1]*
						a[i_forw];
			break;
		}
	case i_forw : 
		{
			val=1+mul*	a[i_vol1]*
						a[i_vol1]*
						a[i_forw]*
						a[i_t];
			break;
		}
	case i_vol1 : 
		{
			val=2.0*mul*a[i_vol1]*
						a[i_forw]*
						a[i_t];
			break;
		}
	default :
		{
			val=0.0;
			break;
			
		}
	}
	return val;
}

///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Change_Numeraire2_N
///	Routine: value
/// Argument: ArgumentList * a ,double submodel_result
///	Returns: int i, ArgumentList * a, double s
///	Action : returns the derivative of the computed value with respect to the ith input
///////////////////////////////////////////////////////////////////////////////////////

template<int i_forw, int i_vol1,int i_vol2, int i_corr,int i_t,int mul>
double Change_Numeraire2_N<i_forw,i_vol1,i_vol2,i_corr,i_t,mul>::value(int i,const ArgumentList& a, double s)
{
	double val;
	double prodall=a[i_vol1]*a[i_vol2]*a[i_corr]*a[i_t]*mul;
	switch (i)
	{
	case i_t : 
	case i_vol1 : 
	case i_vol2 : 
	case i_corr : 
		{
			val=prodall/a[i]*
				a[i_forw];
			break;
		}
	case i_forw : 
		{
			val=1+prodall;
			break;
		}
	default :
		{
			val=0.0;
			break;
			
		}
	}
	return val;
}


///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Change_Numeraire1_N 
///	Routine: value
/// Argument: ArgumentList * a ,double submodel_result
///	Returns: int i, int j,ArgumentList * a, double s
///	Action : returns the 2nd derivative of the computed value 
/// with respect to the ith input  and jth input
///////////////////////////////////////////////////////////////////////////////////////

template<int i_forw, int i_vol1, int i_t,int mul>
double Change_Numeraire1_N<i_forw,i_vol1,i_t,mul>::value(int i,int j,const ArgumentList& a, double s1,double s2)
{
	double val;
	double prodall=a[i_vol1]*a[i_vol1]*a[i_t]*mul;
	switch (i)	{
	case i_t :
		{
			switch (j)	{
			case i_t : 
				val=0.;
				break;
			case i_forw : 
				val=prodall/a[i_t];
				break;
			case i_vol1 : 
				val=mul*2.0*a[i_vol1]*a[i_forw];
				break;
			default :
				val=0.0;
				break;
			}
		}
		break;
	case i_forw : 
		{
			switch (j)	{
			case i_t : 
				val=prodall/a[i_t];
				break;
			case i_forw : 
				val=0.0;
				break;
			case i_vol1 : 
				val=2.0*mul*a[i_vol1]*a[i_t];
				break;
			default :
				val=0.0;
				break;
			}
		}
		break;
	case i_vol1 : 
		{
			switch (j)	{
			case i_t : 
				val=2.0*mul*a[i_forw]*a[i_vol1];
				break;
			case i_forw : 
				val=2.0*mul*a[i_vol1]*a[i_t];
				break;
			case i_vol1 : 
				val=2.0*mul*a[i_forw]*a[i_t];
				break;
			default :
				val=0.0;
				break;
			}
		}
		break;
	default :
		val=0.0;
		break;
	}
	return val;
}
///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Change_Numeraire2_N
///	Routine: value
/// Argument: ArgumentList * a ,double submodel_result
///	Returns: int i, int j,ArgumentList * a, double s
///	Action : returns the 2nd derivative of the computed value 
/// with respect to the ith input  and jth input
///////////////////////////////////////////////////////////////////////////////////////

template<int i_forw, int i_vol1,int i_vol2, int i_corr,int i_t,int mul>
double Change_Numeraire2_N<i_forw,i_vol1,i_vol2,i_corr,i_t,mul>::value(int i,int j,const ArgumentList& a, double s1,double s2)
{
	double val;
	double prodall=a[i_vol1]*a[i_vol2]*a[i_corr]*a[i_t]*mul;
	switch (i)	{
	case i_t :
		{
			switch (j)	{
			case i_t : 
				val=0.;
				break;
			case i_forw : 
				val=prodall/a[i];
				break;
			case i_vol1 : 
			case i_vol2 : 
			case i_corr : 	
				val=prodall*a[i_forw]/(a[j]*a[i]);
				break;
			default :
				val=0.0;
				break;	
			}
			break;
		}
	case i_forw : 
		{
			switch (j)	{
			case i_forw : 
				val=0.0;
				break;
			case i_t : 
			case i_vol1 : 
			case i_vol2 :
			case i_corr : 
				val=prodall/a[j];
				break;
			default :
				val=0.0;
				break;
			}
			break;
		}
	case i_vol1 : 
		{
			switch (j)	{
			case i_t : 
			case i_vol2 : 
			case i_corr : 
				val=prodall*a[i_forw]/(a[j]*a[i]);
				break;
			case i_forw : 
				val=prodall/a[i];
				break;
			case i_vol1 : 
				val=0.;
				break;
			default :
				val=0.0;
				break;
			}
			break;
		}
	case i_vol2 : 
		{
			switch (j)	{
			case i_t : 
			case i_vol1 :
			case i_corr : 
				val=prodall*a[i_forw]/(a[j]*a[i]);
				break;
			case i_forw : 
				val=prodall/a[i];
				break;
			case i_vol2 : 
				val=0.;
				break;
			default :
				val=0.0;
				break;
			}
			break;
		}
	case i_corr : 
		{
			switch (j)	{
			case i_t : 
			case i_vol1 : 
			case i_vol2 : 
				val=prodall*a[i_forw]/(a[j]*a[i]);
				break;
			case i_forw : 
				val=prodall/a[i];
				break;
			case i_corr : 
				val=0.;
				break;
			default :
				val=0.0;
				break;
			}
			break;
			default :
				val=0.0;
				break;
		}
		break;
	}
	return val*exp(prodall);
}

///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Change_Numeraire1_N
///	Routine: prepare_model_input
/// Argument: ArgumentList * a ,double submodel_result
///	Returns: ArgumentList *
///	Action : prepare the input for a model evaluation
///////////////////////////////////////////////////////////////////////////////////////

template<int i_forw, int i_vol1, int i_t,int mul>
ArgumentList * Change_Numeraire1_N<i_forw,i_vol1,i_t,mul>::prepare_model_input(const ArgumentList& a,double submodel_result)
{
	ArgumentList * newarg=new ArgumentList(a);
	newarg->set_nth(i_forw,submodel_result);
	return newarg;
}
///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Change_Numeraire2_N
///	Routine: prepare_model_input
/// Argument: ArgumentList * a ,double submodel_result
///	Returns: ArgumentList *
///	Action : prepare the input for a model evaluation
///////////////////////////////////////////////////////////////////////////////////////

template<int i_forw, int i_vol1,int i_vol2, int i_corr,int i_t,int mul>
ArgumentList * Change_Numeraire2_N<i_forw,i_vol1,i_vol2,i_corr,i_t,mul>::prepare_model_input(const ArgumentList& a,double submodel_result)
{
	ArgumentList * newarg=new ArgumentList(a);
	newarg->set_nth(i_forw,submodel_result);
	return newarg;
}


///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Change_Numeraire1_N
///	Routine: check_argument
/// Argument: ArgumentList * a 
///	Returns: ArgumentList_Checking_Result 
///	Action : for now does no checking further on
///////////////////////////////////////////////////////////////////////////////////////

template<int i_forw, int i_vol1, int i_t,int mul>
ArgumentList_Checking_Result Change_Numeraire1_N<i_forw,i_vol1,i_t,mul>::check_argument(const ArgumentList& a)
{
	return ArgumentList_Checking_Result(true,"");
}

///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Change_Numeraire2_N
///	Routine: check_argument
/// Argument: ArgumentList * a 
///	Returns: ArgumentList_Checking_Result 
///	Action : for now does no checking further on
///////////////////////////////////////////////////////////////////////////////////////

template<int i_forw, int i_vol1,int i_vol2, int i_corr,int i_t,int mul>
ArgumentList_Checking_Result Change_Numeraire2_N<i_forw,i_vol1,i_vol2,i_corr,i_t,mul>::check_argument(const ArgumentList& a)
{
	return ArgumentList_Checking_Result(true,"");
}

///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Change_Numeraire1_N
///	Routine: check_dimension
/// Argument: int i 
///	Returns: ArgumentList_Checking_Result 
///	Action : for now does no checking further on
///////////////////////////////////////////////////////////////////////////////////////

template<int i_forw, int i_vol1, int i_t,int mul>
ArgumentList_Checking_Result Change_Numeraire1_N<i_forw,i_vol1,i_t,mul>::check_dimension(int i)
{
	return ArgumentList_Checking_Result(true,"");
}

///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Change_Numeraire2_N
///	Routine: check_dimension
/// Argument: int i
///	Returns: ArgumentList_Checking_Result
///	Action : for now does no checking further on
///////////////////////////////////////////////////////////////////////////////////////

template<int i_forw, int i_vol1,int i_vol2, int i_corr,int i_t,int mul>
ArgumentList_Checking_Result Change_Numeraire2_N<i_forw,i_vol1,i_vol2,i_corr,i_t,mul>::check_dimension(int i)
{
	return ArgumentList_Checking_Result(true,"");
}


CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

