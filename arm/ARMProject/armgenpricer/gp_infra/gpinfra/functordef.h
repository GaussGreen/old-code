/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: functordef.h,v $
 * Revision 1.1  2003/10/21 07:52:17  ebenhamou
 * Initial revision
 *
 *
 */

/*! \file functordef.h
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */

#ifndef _INGPINFRA_FUNCTORDEF_H
#define _INGPINFRA_FUNCTORDEF_H

#include "gpbase/port.h"
#include <functional>

#include <string>
CC_USING_NS( std, string )

///////////////////////////////////////////////////
/// these functions are largelly implemented in terms of the 
/// STL library <functional>
/// The difference is that
///		- the various predicate function returns double
///		instead of bool for consistency reason
///		- toString member function has been added to 
///		provides access to the name!
///     - access to template data is by value rather than
///     const reference (as implemented in the standard STL).
///     - modulus also casts into int...
///
/// EXCEPTIONALLY we do not use the namespace ARM
/// but rather ARM_GP to avoid confusion with the STL
/// library <functional> functions.!
///////////////////////////////////////////////////

CC_BEGIN_NAMESPACE( ARM_GP )

/// TEMPLATE STRUCT negate
template<typename T>
	struct negate : CC_NS( std, unary_function )<T,T> {
		T operator()(T _X) const
			{return (-_X); }
		string toString() const{ return "negate"; }
	};

/// TEMPLATE STRUCT logical_not
template<typename T>
	struct logical_not : CC_NS( std, unary_function )<T,double>{
		/// not bool but rather double
		/// for consistency reason
		double operator()(T _X) const
			{return (!_X); }
		string toString() const{ return "logical_not"; }
	};


/// TEMPLATE STRUCT plus
template<typename T>
	struct plus : CC_NS( std, binary_function )<T,T,T>{
		T operator()(T _X, T _Y) const
			{return (_X + _Y); }
		string toString() const{ return "plus"; }
	};

/// TEMPLATE STRUCT minus
template<typename T>
	struct minus : CC_NS( std, binary_function )<T,T,T>{
		T operator()(T _X, T _Y) const
			{return (_X - _Y); }
		string toString() const{ return "minus"; }
	};

/// TEMPLATE STRUCT multiplies
template<typename T>
	struct multiplies : CC_NS( std, binary_function )<T,T,T>{
		T operator()(T _X, T _Y) const
			{return (_X * _Y); }
		string toString() const{ return "multiplies"; }
	};

/// TEMPLATE STRUCT divides
template<typename T>
	struct divides : CC_NS( std, binary_function )<T,T,T>{
		T operator()(T _X, T _Y) const
			{return (_X / _Y); }
		string toString() const{ return "divides"; }
	};

/// TEMPLATE STRUCT modulus
template<typename T>
	struct modulus : CC_NS( std, binary_function )<T,T,double> {
	/// compared to the standard modulus function
	/// this casts the type to int and return a double correpsonding to it
	double operator()(T _X, T _Y) const
		{return ( int(_X) % int(_Y) ); }
	string toString() const{ return "modulus"; }
	};

/// TEMPLATE STRUCT equal_to
template<typename T>
	struct equal_to : CC_NS( std, binary_function )<T,T,double> {
	/// not bool but rather double
	/// for consistency reason
	double operator()(T _X, T _Y) const
		{return (_X == _Y); }
	string toString() const{ return "equal_to"; }
	};

/// TEMPLATE STRUCT not_equal_to
template<typename T>
	struct not_equal_to : CC_NS( std, binary_function )<T,T,double>  {
	/// not bool but rather double
	/// for consistency reason
	double operator()(T _X, T _Y) const
		{return (_X != _Y); }
	string toString() const{ return "not_equal_to"; }
	};

/// TEMPLATE STRUCT greater
template<typename T>
	struct greater : CC_NS( std, binary_function )<T,T,double> {
	/// not bool but rather double
	/// for consistency reason
	double operator()(T _X, T _Y) const
		{return (_X > _Y); }
	string toString() const{ return "greater"; }
	};

/// TEMPLATE STRUCT less
template<typename T>
	struct less : CC_NS( std, binary_function )<T,T,double> {
	/// not bool but rather double
	/// for consistency reason
	double operator()(T _X, T _Y) const
		{return (_X < _Y); }
	string toString() const{ return "less"; }
	};

/// TEMPLATE STRUCT greater_equal
template<typename T>
	struct greater_equal : CC_NS( std, binary_function )<T,T,double> {
	/// not bool but rather double
	/// for consistency reason
	double operator()(T _X, T _Y) const
		{return (_X >= _Y); }
	string toString() const{ return "greater_equal"; }
	};

/// TEMPLATE STRUCT less_equal
template<typename T>
	struct less_equal : CC_NS( std, binary_function )<T,T,double> {
	/// not bool but rather double
	/// for consistency reason
	double operator()(T _X, T _Y) const
		{return (_X <= _Y); }
	string toString() const{ return "less_equal"; }
	};

/// TEMPLATE STRUCT logical_and
template<typename T>
	struct logical_and : CC_NS( std, binary_function )<T,T,double> {
	/// not bool but rather double
	/// for consistency reason
	double operator()(T _X, T _Y) const
		{return (_X && _Y); }
	string toString() const{ return "logical_and"; }
	};

/// TEMPLATE STRUCT logical_or
template<typename T>
	struct logical_or : CC_NS( std, binary_function )<T,T,double> {
	/// not bool but rather double
	/// for consistency reason
	double operator()(T _X, T _Y) const
		{return (_X || _Y); }
	string toString() const{ return "logical_or"; }
	};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/




