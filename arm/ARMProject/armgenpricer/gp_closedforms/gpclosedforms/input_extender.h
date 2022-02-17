/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  input extender for the closed form framework 
 *
 *	\file input_extender.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_INPUT_EXTENDER_H
#define _GP_CF_INPUT_EXTENDER_H


#define ARM_CF_DEFAULT_SHIFT  0.0001

#include "firsttoinc.h"
#include "gpbase/port.h"

#include "basic_distributions.h"

CC_BEGIN_NAMESPACE(ARM)


template <class model ,int n>
class Add_Input 
{
public:
	enum
	{
		Nb_Parameters=model::Nb_Parameters
	};
	static double value(const ArgumentList& a);
	static double value(int i,const ArgumentList& a, double s);
	static double value(int i,int j,const ArgumentList& a, double s1,double s2);
	static double specific_shift(int i) ;
	static ArgumentList_Checking_Result check_argument(const ArgumentList& a);
	static ArgumentList_Checking_Result check_dimension(int i);
};


template <class M ,int n>
double Add_Input<M,n>::value(const ArgumentList& a)
{
	return M::value(*a.sublist(a.size()-n));
}

template <class M ,int n>
double Add_Input<M,n>::value(int i,const ArgumentList& a, double s)
{
	return M::value(i,*a.sublist(a.size()-n),s);
}

template <class M ,int n>
double Add_Input<M,n>::value(int i,int j,const ArgumentList& a, double s1,double s2)
{
	return M::value(i,j,*a.sublist(a.size()-n),s1,s2);
}


template <class M ,int n>
double Add_Input<M,n>::specific_shift(int i)
{
	if (i>=M::Nb_Parameters) return ARM_CF_DEFAULT_SHIFT;
	else
	{
	
		return M::specific_shift(i);
	}
}

template <class M ,int n>
ArgumentList_Checking_Result Add_Input<M,n>::check_argument(const ArgumentList& a)
{
	return M::check_argument(*a.sublist(M::Nb_Parameters));
}

template <class M ,int n>
ArgumentList_Checking_Result Add_Input<M,n>::check_dimension(int rank)
{
	if ((rank<0)||(rank>(M::Nb_Parameters+n-1))) return ArgumentList_Checking_Result(false," Invalide derivative dimension"); // Invalide derivative dimension	
	
	return ArgumentList_Checking_Result(true,string(""));
}


CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

