/*Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *  \file meanrevfinder.cpp
 *
 *  \brief
 *
 *	\author  H. BAKHTRI 
 *	\version 1.0
 *	\date October 2005
 *
 * 
 * Brief file to obtain the mean reversion level corresponding to a given price. 
 * It is a basic newton raphson solver.
 */

#include "gpbase/removeidentifiedwarning.h"
#include "gpcalculators/meanrevfinder.h"
#include <iomanip>	
//GP_Base
#include "gpbase/curve.h"

/// gpinfra
#include "gpinfra/argconvdefault.h"
#include "gpinfra/gensecurity.h"
#include "gpinfra/dealdescription.h"
#include "gpinfra/gensecmanipulator.h"
#include "gpinfra/mktdatamanagerrep.h"
#include "gpinfra/pricingadviser.h"
#include "gpinfra/gramfunctorargdict.h"
#include "gpinfra/cstmanager.h"
#include "gpinfra/pricingmodel.h"

/// gpcalib
#include "gpcalib/calibmethod.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_MeanRevFinder
///	Routine: Constructor
///	Returns: 
///	Action : Build the MeanRevFinderObject
////////////////////////////////////////////////////
ARM_MeanRevFinder::ARM_MeanRevFinder(ARM_BermudaSwaptionCalculator* calculator, double mrs):
itsCalculator(calculator),
itsMrs(mrs)
{

}

////////////////////////////////////////////////////
///	Class  : ARM_MeanRevFinder
///	Routine: Operator =
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_MeanRevFinder& ARM_MeanRevFinder::operator=( const ARM_MeanRevFinder& rhs )
{
	if (this != &rhs)
	{
		itsCalculator		= rhs.itsCalculator;
		itsMrs				= rhs.itsMrs;
	}

	return *this;
}

////////////////////////////////////////////////////
///	Class  : ARM_MeanRevFinder
///	Routine: Desctructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_MeanRevFinder::~ARM_MeanRevFinder()
{

}

////////////////////////////////////////////////////
///	Class  : ARM_MeanRevFinder
///	Routine: Operator()
///	Returns: double
///	Action : 
////////////////////////////////////////////////////
double ARM_MeanRevFinder::operator() (double mrs) const  
{
	itsCalculator->SetInitMrs(mrs);
	itsCalculator->CreateAndSetModel();

	return itsCalculator->Price();
}


CC_END_NAMESPACE()

//-----------------------------------------------------------------------------
/*---- End of file ----*/
