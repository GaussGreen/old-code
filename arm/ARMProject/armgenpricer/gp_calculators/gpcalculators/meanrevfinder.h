/*Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *  \file meanrevfinder.h
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

#ifndef _INGPCALCULATORS_MEANREVFINDER_H
#define _INGPCALCULATORS_MEANREVFINDER_H

#include "bermudaswaptioncalculator.h"
#include "gpbase/port.h"
#include "gpcalib/typedef.h"
#include "gpbase/rootobject.h"

CC_BEGIN_NAMESPACE( ARM )


class ARM_MeanRevFinder : public ARM_GP::UnaryFunc<double,double>
{
public:

	ARM_MeanRevFinder(ARM_BermudaSwaptionCalculator* calculator, double mrs);
	ARM_MeanRevFinder& operator = ( const ARM_MeanRevFinder& rhs );
	virtual ~ARM_MeanRevFinder();

	virtual double operator() (double mrs) const; 

private:

	ARM_BermudaSwaptionCalculator*	itsCalculator;
	double							itsMrs;
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/


