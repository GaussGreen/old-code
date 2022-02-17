/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file impsampleropt.h
 *
 *  \brief
 *
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date Dec 2005
 */


#ifndef _INGPNUMMETHODS_IMPSAMPLEROPT_H
#define _INGPNUMMETHODS_IMPSAMPLEROPT_H

/// gpbase
#include "gpbase/port.h"
#include "gpbase/functor.h"
#include "gpbase/curve.h"
#include "typedef.h"

CC_BEGIN_NAMESPACE( ARM )

// With this class you can optimize the alpa curve of the importance
// for a better a variance reduction.

class ARM_ImpSamplerOpt
{
public:
	static ARM_PricingModel* ComputeModel(
		ARM_GenSecurity* genSec,
		ARM_PricingModel* model,
		const ARM_MultiCurve& initGuess,
		const ARM_MultiCurve& lowerBound,
		const ARM_MultiCurve& upperBound,
		bool withMC,
		long nbSteps,
		bool bootstrap);
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

