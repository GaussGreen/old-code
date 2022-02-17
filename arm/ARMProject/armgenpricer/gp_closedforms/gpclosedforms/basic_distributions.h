/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file basic_distribution.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_BASIC_DISTRIBUTIONS_H
#define _GP_CF_BASIC_DISTRIBUTIONS_H

#include "firsttoinc.h"
#include "gpbase/port.h"

#include "expt.h"
#include <vector>
#include <functional>
#include "gpbase/gpmatrix.h"
#include "gpbase/gpvector.h"
#include "gpbase/gplinalgtypedef.h"
//#include "gpbase/removenagwarning.h"

CC_USING_NS(std,vector)
CC_USING_NS(std,string)

CC_BEGIN_NAMESPACE(ARM)

////////////////////////////////////////////////////////////////////
///
/// by convention the time to expiration of  the options is the first parameter
/// of the argument list
/// it is necessary for the change of numeraire to work
///
///////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
/// \class ArgumentList
/// \brief
/// ArgumentList class specifies the argument of a pricing function
/// it can built from an int, or between 1 and 15 double
/// if it is an int, it build an argument list of the specified size 
////////////////////////////////////////////////////////////////////
class ArgumentList 
{

public:
	ArgumentList(int n);
	ArgumentList(int n,int p);
	ArgumentList(std::vector<double>*  x,std::vector<double>*  y,double x0);
	ArgumentList(std::vector<double>*  x,std::vector<double>*  y,double x0,double x1);
	ArgumentList(std::vector<double>*  x,std::vector<double>*  y,double x0,double x1,double x2);
	ArgumentList(std::vector<double>*  x,std::vector<double>*  y,double x0,double x1,double x2,double x3);
	ArgumentList(std::vector<double>*  x,std::vector<double>*  y,double x0,double x1,double x2,double x3,double x4);
	ArgumentList(std::vector<double>*  x,std::vector<double>*  y,double x0,double x1,double x2,double x3,double x4,double x5);
	ArgumentList(std::vector<double>*  x,std::vector<double>*  y,double x0,double x1,double x2,double x3,double x4,double x5,double x6);
	ArgumentList(std::vector<double>*  x,std::vector<double>*  y,std::vector<double>*  xv,std::vector<double>*  y2,double x0);
	ArgumentList(std::vector<double>*  x,std::vector<double>*  y,std::vector<double>*  xv,std::vector<double>*  y2,double x0,double x1);
	ArgumentList(std::vector<double>*  x,std::vector<double>*  y,std::vector<double>*  xv,std::vector<double>*  y2,double x0,double x1,double x2);
	ArgumentList(std::vector<double>*  x,std::vector<double>*  y,std::vector<double>*  xv,std::vector<double>*  y2,double x0,double x1,double x2,double x3);
	ArgumentList(std::vector<double>*  x,std::vector<double>*  y,std::vector<double>*  xv,std::vector<double>*  y2,double x0,double x1,double x2,double x3,double x4);
	ArgumentList(std::vector<double>*  x,std::vector<double>*  y,std::vector<double>*  xv,std::vector<double>*  y2,double x0,double x1,double x2,double x3,double x4,double x5);
	ArgumentList(std::vector<double>*  x,std::vector<double>*  y,std::vector<double>*  xv,std::vector<double>*  y2,double x0,double x1,double x2,double x3,double x4,double x5,double x6);
	ArgumentList(std::vector<double>*  x,std::vector<double>*  y,std::vector<double>*  xv,std::vector<double>*  y2,double x0,double x1,double x2,double x3,double x4,double x5,double x6,
		double x7,double x8,double x9,double x10,double x11,double x12,double x13,double x14,double x15,double x16,double x17);
	ArgumentList(std::vector<double>*  x,std::vector<double>*  y,std::vector<double>*  xv,std::vector<double>*  y2,double x0,double x1,double x2,double x3,double x4,double x5,double x6,
		double x7,double x8,double x9,double x10,double x11,double x12,double x13,double x14,double x15,double x16,double x17,double x18);
	ArgumentList(std::vector<double>*  x,std::vector<double>*  y,std::vector<double>*  xv,std::vector<double>*  y2,double x0,double x1,double x2,double x3,double x4,double x5,double x6,
		double x7,double x8,double x9,double x10,double x11,double x12,double x13,double x14,double x15,double x16,double x17,double x18,double x19);

	ArgumentList(std::vector<double>*  x,std::vector<double>*  y,std::vector<double>*  xv,double x0);
	ArgumentList(std::vector<double>*  x,std::vector<double>*  y,std::vector<double>*  xv,double x0,double x1);
	ArgumentList(std::vector<double>*  x,std::vector<double>*  y,std::vector<double>*  xv,double x0,double x1,double x2);
	ArgumentList(std::vector<double>*  x,std::vector<double>*  y,std::vector<double>*  xv,double x0,double x1,double x2,double x3);
	ArgumentList(std::vector<double>*  x,std::vector<double>*  y,std::vector<double>*  xv,double x0,double x1,double x2,double x3,double x4);
	ArgumentList(std::vector<double>*  x,std::vector<double>*  y,std::vector<double>*  xv,double x0,double x1,double x2,double x3,double x4,double x5);
	ArgumentList(std::vector<double>*  x,std::vector<double>*  y,std::vector<double>*  xv,double x0,double x1,double x2,double x3,double x4,double x5,double x6);
	ArgumentList(double x0);
	ArgumentList(double x0,double x1);
	ArgumentList(double x0,double x1,double x2);
	ArgumentList(double x0,double x1,double x2,double x3);
	ArgumentList(double x0,double x1,double x2,double x3,double x4);
	ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5);
	ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,double x6);
	ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,double x6,double x7);
	ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,double x6,double x7,double x8);
	ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,double x6,double x7,double x8,double x9);
	ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,double x6,double x7,double x8,double x9,double x10);
	ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,double x6,double x7,double x8,double x9,double x10,double x11);
	ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,double x6,double x7,double x8,double x9,double x10,double x11,double x12);
	ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,double x6,double x7,double x8,double x9,double x10,double x11,double x12,double x13);
	ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,double x6,double x7,double x8,double x9,double x10,double x11,double x12,double x13,double x14);
	ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,double x6,double x7,double x8,double x9,double x10,double x11,double x12,double x13,double x14,double x15);
	ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,double x6,double x7,double x8,double x9,double x10,double x11,double x12,double x13,double x14,double x15,double x16);
	ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,double x6,double x7,double x8,double x9,double x10,double x11,double x12,double x13,double x14,double x15,double x16,double x17);
	ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,double x6,double x7,double x8,double x9,double x10,double x11,double x12,double x13,double x14,double x15,double x16,double x17,double x18);
	ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,double x6,double x7,double x8,double x9,double x10,double x11,double x12,double x13,double x14,double x15,double x16,double x17,double x18,double x19);
	ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,double x6,double x7,double x8,double x9,double x10,double x11,double x12,double x13,double x14,double x15,double x16,double x17,double x18,double x19,double x20);
	ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,double x6,double x7,double x8,double x9,double x10,double x11,double x12,double x13,double x14,double x15,double x16,double x17,double x18,double x19,double x20,double x21);
	ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,double x6,double x7,double x8,double x9,double x10,double x11,double x12,double x13,double x14,double x15,double x16,double x17,double x18,double x19,double x20,
		double x21,double x22);
	ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,double x6,double x7,double x8,double x9,double x10,double x11,double x12,double x13,double x14,double x15,double x16,double x17,double x18,double x19,double x20,
		double x21,double x22,double x23);
	ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,double x6,double x7,double x8,double x9,double x10,double x11,double x12,double x13,double x14,double x15,double x16,double x17,double x18,double x19,double x20,
		double x21,double x22,double x23,double x24);
	ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,double x6,double x7,double x8,double x9,double x10,double x11,double x12,double x13,double x14,double x15,double x16,double x17,double x18,double x19,double x20,
		double x21,double x22,double x23,double x24,double x25);
	ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,double x6,double x7,double x8,double x9,double x10,double x11,double x12,double x13,double x14,double x15,double x16,double x17,double x18,double x19,double x20,
		double x21,double x22,double x23,double x24,double x25,double x26);
	ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,double x6,double x7,double x8,double x9,double x10,double x11,double x12,double x13,double x14,double x15,double x16,double x17,double x18,double x19,double x20,
		double x21,double x22,double x23,double x24,double x25,double x26,double x27);
	ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,double x6,double x7,double x8,double x9,double x10,double x11,double x12,double x13,double x14,double x15,double x16,double x17,double x18,double x19,double x20,
		double x21,double x22,double x23,double x24,double x25,double x26,double x27,double x28);
	ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,double x6,double x7,double x8,double x9,double x10,double x11,double x12,double x13,double x14,double x15,double x16,double x17,double x18,double x19,double x20,
		double x21,double x22,double x23,double x24,double x25,double x26,double x27,double x28,double x29);
	ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,double x6,double x7,double x8,double x9,double x10,double x11,double x12,double x13,double x14,double x15,double x16,double x17,double x18,double x19,double x20,
		double x21,double x22,double x23,double x24,double x25,double x26,double x27,double x28,double x29,double x30);
	ArgumentList(double x0,double x1,double x2,double x3,double x4,double x5,double x6,double x7,double x8,double x9,double x10,double x11,double x12,double x13,double x14,double x15,double x16,double x17,double x18,double x19,double x20,
		double x21,double x22,double x23,double x24,double x25,double x26,double x27,double x28,double x29,double x30,double x31);



	inline size_t size() const { return itsContainer.size();}
	inline size_t size2() const { return itsVectorContainer.size();}
	inline double& nth(int i) { return itsContainer[i];}
	inline double nth(int i) const { return itsContainer[i];}
	inline void set_nth(int i, double x) {itsContainer[i]=x;}
	const vector<double>& get_itsContainer() const { return itsContainer; }
	vector<double>& get_itsContainer() { return itsContainer; }
	const vector<std::vector<double>*>& get_itsVectorContainer() const { return itsVectorContainer; }
	vector<std::vector<double>*>& get_itsVectorContainer() { return itsVectorContainer; }
	inline double operator[](int i) const { return itsContainer[i];}
	inline double& operator[](int i) { return itsContainer[i];}
	std::vector<double>*  vectnth(int i) { return itsVectorContainer[i]; }
	std::vector<double>*  V(int i) const { return itsVectorContainer[i]; }
	std::vector<double>*& V(int i)  { return itsVectorContainer[i]; }
	inline void set_vectnth(int i,std::vector<double>*  x) {itsVectorContainer[i]=x;}
	ArgumentList* sublist(int n) const;
private:
	vector<double> itsContainer;
	vector<std::vector<double>*> itsVectorContainer;

};

////////////////////////////////////////////////////////////////////
/// \class Power_ArgumentList
/// \brief
/// Power_ArgumentList class creates an Argument list of the size
/// given by the type given in argument of the template
/// it allows to avoid to have to specify it.
////////////////////////////////////////////////////////////////////
template<typename A>
class Power_ArgumentList : public ArgumentList
{
public:
	Power_ArgumentList();
};



template<typename A>
Power_ArgumentList<A>::Power_ArgumentList():
	ArgumentList(A::nb_arguments())
{}

////////////////////////////////////////////////////////////////////
/// \class Expression
/// \brief
/// Expression class is an abstract class for the algorithms
/// to be computed
////////////////////////////////////////////////////////////////////
class Expression : public CC_NS(std,unary_function)<void,double> 
{
public:
	virtual double operator() (const ArgumentList& a) const=0;
private:
	int nb_arg;
};


////////////////////////////////////////////////////////////////////
/// \class ArgumentList_Checking_Result
/// \brief
/// ArgumentList_Checking_Result class is a container to return the 
/// result of range checking of the arguments
////////////////////////////////////////////////////////////////////
struct ArgumentList_Checking_Result
{
	bool	checkresult;
	string	trouble;
	ArgumentList_Checking_Result(bool res, string trb);
};

ArgumentList_Checking_Result operator&&(ArgumentList_Checking_Result arg1,ArgumentList_Checking_Result arg2);
ArgumentList_Checking_Result operator||(ArgumentList_Checking_Result arg1,ArgumentList_Checking_Result arg2);


////////////////////////////////////////////////////////////////////
/// \class Power_Expression
/// \brief
/// Power_Expression class is the class to use for creating closed 
/// forms. It creates an instance of a formula specified as the template
/// argument.
/// it has an operator interpretation that needs 1,2 or 3 arguments
/// it can compute the first and second derivative of the formula
////////////////////////////////////////////////////////////////////
template <typename A>
class Power_Expression: public Expression
{
public:
	double operator() (const ArgumentList& a) const
	{
		ArgumentList_Checking_Result argchk=A::check_argument(a);
		if (argchk.checkresult) return A::value(a);
		else
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,argchk.trouble );
		}
	}
	
	double operator() (int i, const ArgumentList& a) const
	{
		ArgumentList_Checking_Result argchk=A::check_argument(a);
		ArgumentList_Checking_Result argchk1=A::check_dimension(i);
		if ((argchk.checkresult) && (argchk1.checkresult))
		{

			double s=A::specific_shift(i);
			return A::value(i,a,s);
		}
		else
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,argchk.trouble+argchk1.trouble );
		}
		
	}

	double operator() (int i, int j,const ArgumentList& a) const
	{
		ArgumentList_Checking_Result argchk=A::check_argument(a);
		ArgumentList_Checking_Result argchk1=A::check_dimension(i);
		ArgumentList_Checking_Result argchk2=A::check_dimension(j);
		if ((argchk.checkresult) && (argchk1.checkresult) && (argchk2.checkresult))
		{
			double s1=A::specific_shift(i);
			double s2=A::specific_shift(j);
			return A::value(i,j,a,s1,s2);
		}
		else
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,argchk.trouble+argchk1.trouble+argchk2.trouble );
		}
		
	}		
};

double standard_first_derivative(double(*f)(const ArgumentList& a) ,int i,const ArgumentList& a,double s);
double standard_second_derivative(double(*f)(const ArgumentList& a) ,int i,int j,const ArgumentList& a,double s1,double s2);

double centred_standard_first_derivative(double(* f)(const ArgumentList& a),int i,const ArgumentList& a,double s);

CC_END_NAMESPACE()

#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

