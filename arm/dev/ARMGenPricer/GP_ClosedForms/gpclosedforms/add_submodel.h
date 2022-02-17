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
 
#ifndef _GP_CF_ADD_SUBMODEL_H
#define _GP_CF_ADD_SUBMODEL_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"

CC_BEGIN_NAMESPACE(ARM)


////////////////////////////////////////////////////////////////////
/// \class Add_SubModel
/// \brief
/// Add_SubModel class specifies that the input indexed by i is extended
/// by a function specified by the class extender
////////////////////////////////////////////////////////////////////


template <class model ,class submodel>
class Add_SubModel 
{
public:
	enum
	{
		Nb_Parameters=model::Nb_Parameters
	};
	static double value(const ArgumentList& a);
	static double value(int i, const ArgumentList& a, double s);
	static double value(int i,int j, const ArgumentList& a, double s1,double s2);
	static double specific_shift(int i) ;
	static double value(int i,const ArgumentList& a);
	static ArgumentList_Checking_Result check_argument(const ArgumentList& a);
	static ArgumentList_Checking_Result check_dimension(int i);

};


///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Add_SubModel 
///	Routine: value
/// Argument: ArgumentList * a 
///	Returns: double
///	Action : Computation of a value
///////////////////////////////////////////////////////////////////////////////////////

template<class M,class E>
double Add_SubModel<M,E>::value(const ArgumentList& a)
{
	double Modified_Input=E::value(a);
	ArgumentList * modified_arg=E::prepare_model_input (a,Modified_Input );
	double val= M::value (*modified_arg);
	delete modified_arg;
	return val;
}


///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Add_SubModel 
///	Routine: value
/// Argument: int i, ArgumentList * a, double s
///	Returns: double
///	Action : Computation of a first order derivative
///////////////////////////////////////////////////////////////////////////////////////


template<class M,class E>
double Add_SubModel<M,E>::value(int i, const ArgumentList& a, double s)
{
	int i_aff=E::affected_index();
	
	double Modified_Input=E::value(a);

	ArgumentList * modified_arg = E::prepare_model_input (a, Modified_Input );

	double deltabar_i_f=1.;
	if(i==i_aff) deltabar_i_f=0.;

	double result= deltabar_i_f*M::value(i,*modified_arg,s)+M::value(i_aff,*modified_arg,s)*E::value(i,a,s);

	delete modified_arg;

	return result;
	

}

///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Add_SubModel 
///	Routine: value
/// Argument: int i,int j, ArgumentList * a, double s1,double s2
///	Returns: double
///	Action : Computation of a second order derivative
///////////////////////////////////////////////////////////////////////////////////////


template <class M ,class E>
double Add_SubModel<M,E>::value(int i,int j,const  ArgumentList& a, double s1,double s2)
{

	
	int i_aff=E::affected_index();
	double s_aff=M::specific_shift(i_aff);
	
	double Modified_Input=E::value(a);

	ArgumentList * modified_arg = E::prepare_model_input (a, Modified_Input );

	double deltabar_i_f=1.;
	if(i==i_aff) deltabar_i_f=0.;

	double deltabar_j_f=1.;
	if(j==i_aff) deltabar_j_f=0.;

	double result= deltabar_i_f		*(deltabar_j_f*M::value(i,j,*modified_arg,s1,s2)+
		E::value(j,a,s2)*M::value(i,i_aff,*modified_arg,s1,s_aff))+
		E::value(i,a,s1)	*(deltabar_j_f*M::value(j,i_aff,*modified_arg,s2,s_aff)+	
		E::value(j,a,s1)*M::value(i_aff,i_aff,*modified_arg,s_aff,s_aff))+
		M::value(i_aff,*modified_arg,s_aff)*E::value(i,j,a,s1,s2);
	delete modified_arg;
	return result;
}

///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Add_SubModel 
///	Routine: specific_shift
/// Argument: int i
///	Returns: double
///	Action : return the shift in case of numarical derivation
///////////////////////////////////////////////////////////////////////////////////////

template <class M ,class E>
double Add_SubModel<M,E>::specific_shift(int i) 
{
		return M::specific_shift(i);
}

///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Add_SubModel 
///	Routine: check_argument
/// Argument: ArgumentList * a
///	Returns: double
///	Action : check the range of the arguments and the number of arguments
///////////////////////////////////////////////////////////////////////////////////////

template <class M ,class E>
 ArgumentList_Checking_Result Add_SubModel<M,E>::check_argument(const ArgumentList& a)
{
		return (M::check_argument(a))&&(E::check_argument(a));
}


///////////////////////////////////////////////////////////////////////////////////////
///	Class  : Add_SubModel 
///	Routine: check_dimension
/// Argument: int i
///	Returns: ArgumentList_Checking_Result
///	Action : check the dimension that we want to derive
///////////////////////////////////////////////////////////////////////////////////////
template <class M ,class E>
 ArgumentList_Checking_Result Add_SubModel<M,E>::check_dimension(int i)
{
		return (M::check_dimension(i))&&(E::check_dimension(i));
}





	
CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

