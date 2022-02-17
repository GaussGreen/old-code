/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 * 
 *  \file assignop.h
 *  \brief Macro to generate automatically assignement operator
 * 
 *	\author  A. Chaix
 *	\version 1.0
 *	\date November 2005
 */

#ifndef _INGPBASE_ASSIGNOP_H
#define _INGPBASE_ASSIGNOP_H

///
///	--> Generates automatically the code for operator = 
/// -----------------------------------------------------
/// operator = will be based upon the copy constructor
/// 
/// Should be used as follows:
/// - - - - - -  - - - - - - -
///		class MyClass
///		{
///			private:
///				/// attributes ...
///
///			public:
///				MyClass();
///				MyClass(const MyClass& rhs);
///				virtual ~MyClass();
///				
///				ASSIGN_OPERATOR(MyClass)
///		};
///				
///

// FIXMEFRED: mig.vc8 (31/05/2007 10:26:37): problem in define
#define ASSIGN_OPERATOR(CLASS_NAME)					\
CLASS_NAME & operator = (const CLASS_NAME & rhs)	\
{													\
	if (&rhs != this)								\
	{												\
		this->~CLASS_NAME();						\
		new (this) CLASS_NAME (rhs);				\
	}												\
	return *this;									\
}															

	
#endif ///_INGPBASE_ASSIGNOP_H