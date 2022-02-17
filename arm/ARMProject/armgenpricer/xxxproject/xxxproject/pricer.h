/*!
 *
 * Copyright (c) IXIS-CIB March 2006
 *
 *	\file pricer.h
 *  \brief file for the abstract pricer
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date March 2006
 */

#ifndef INXXXPROJECT_PRICER
#define INXXXPROJECT_PRICER

#include <gpbase/port.h>
#include "xxxproject/mktdatas.h"

CC_BEGIN_NAMESPACE( ARM )

class ARM_MktData;

// This abstract class define a pricer : any kind of intrument which can be
// priced with the MktDatas
class ARM_Pricer
{
public:

	virtual double				Price ()= 0;

	virtual void				SetMkt(	ARM_MktData * )=0;

	virtual ARM_Observator*		GetObservator(	 ) const =0;
};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/