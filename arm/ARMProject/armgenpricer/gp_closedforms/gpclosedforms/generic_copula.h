/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file gaussian_copula.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_GENERICCOPULADESQCRIPTION_H
#define _GP_CF_GENERICCOPULADESQCRIPTION_H

#include "firsttoinc.h"
#include "gpbase/port.h"

#include "basic_distributions.h"

CC_BEGIN_NAMESPACE(ARM)

///////////////////////////////////////////////////////////////////////////////////////
///
///	Class  : Generic Copula : 
/// 
///////////////////////////////////////////////////////////////////////////////////////

struct GenericCopula
{
		/// this the generic interface to copula 
	ArgumentList* structure;
	GenericCopula(ArgumentList* structure0):structure(structure0) {}
};



CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

