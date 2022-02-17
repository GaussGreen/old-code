/*!
 *
 * Copyright (c) IXIS CIB December 2004 Paris
 *
 *	\file gensecurity.h
 *
 *  \brief 
 *
 *	\author  A. Schauly
 *	\version 1.0
 *	\date December 2004
 */


#ifndef _INGPINFRA_GENSECMANIP_H
#define _INGPINFRA_GENSECMANIP_H

CC_BEGIN_NAMESPACE( ARM )

// ARM_GenSecurity Forward declaration
class ARM_GenSecurity;

class ARM_GenSecManipulator
{
public: 
	// 
	void ChangeAmericanIntoTrigger( ARM_GenSecurity& gensec, ARM_DealDescriptionPtr& dealdesc );

	// default constructor
	ARM_GenSecManipulator() {}

};
CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/