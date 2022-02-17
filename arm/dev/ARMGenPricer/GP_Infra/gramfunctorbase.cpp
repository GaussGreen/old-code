/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: gramfunctorbase.cpp,v $
 * Revision 1.1  2003/15/08 19:43:31  ebenhamou
 * Initial revision
 *
 */

/*! \file gramfunctorbase.cpp
 *
 *  \brief gramfunctorbase is the base class for all the grammar functors
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */

#include "gpinfra/gramfunctorbase.h"
#include "gpbase/ostringstream.h"
#include <glob/expt.h>

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_GramFctor
///	Member : itsFuncName is more for debugging output
////////////////////////////////////////////////////

string ARM_GramFctor::itsFuncName   = "ARM_GramFctor abstract Class"; 

////////////////////////////////////////////////////
///	Class  : ARM_GramFctor
///	Routine: GetUsedTimeLags
///	Returns: ARM_VectorPtr of time lags
///	Action : computes time lag given the input!
////////////////////////////////////////////////////

ARM_NodeInfo ARM_GramFctor::GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	double evalDate, vector<ARM_ExpNodePtr>& nodes )
{
	CC_Ostringstream os;
	os << "this function  " << itsFuncName << " does not use time lags";
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
}



CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

