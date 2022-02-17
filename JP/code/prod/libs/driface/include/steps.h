/*===================================================================
** PROGRAM:		STEPS.H
**
** AUTHOR:		Joris Klaentschi.	22 September 97.
**
** CONTAINS:	Tree building / evaluation function components.
**=================================================================*/

#ifndef STEPS_H
#define STEPS_H

#include "cgeneral.h"

GTO_EXPORT(int) DrTreesOptGen1DSingleStep(
		double*		vAssetEnd,			/* (I) underlying asset value at the end of the step */
		double*		vPayoffEnd,			/* (I) option payoff on underlying asset values at the end of the step */
		long		nEnd,				/* (I) number of nodes at end of step */
		double*		vAssetStart,		/* (I) underlying asset values at the start of the step */
		long		nStart,				/* (I) number of nodes at start of step */			
		double		stdev,				/* (I) standard deviation across step - NOT same as ChaseTools volatility! */
		double		hiLimit,			/* (I) if specified, value of high integration limit to use */	
		double		loLimit,			/* (I) if specified, value of low integration limit to use */	
		double*		vPayoffStart);		/* (O) calculated vector payoff values at the start of the step, not allocated here */

GTO_EXPORT(int) DrTreesSplit1DSingleStep(
		double*		vAssetEnd,			/* (I) underlying asset value at the end of the step */
		double*		vPayoffEnd,			/* (I) option payoff on underlying asset values at the end of the step */
		long		nEnd,				/* (I) number of nodes at end of step */
		long		nEndOne,			/* (I) number of nodes at end in middle slice. */
		long		nEndTwo,			/* (I) number of nodes at end in bottom slice. */
		double*		vAssetStart,		/* (I) underlying asset values at the start of the step */
		long		nStart,				/* (I) number of nodes at start of step */			
		double		stdev,				/* (I) standard deviation across step - NOT same as ChaseTools volatility! */
		double		hiLimit,			/* (I) if specified, value of high integration limit to use */	
		double		loLimit,			/* (I) if specified, value of low integration limit to use */	
		double      oneLimit,			/* (I) first middle limit to break integral at */
		double      twoLimit,			/* (I) second middle limit to break integral at */
		double*		vPayoffStart);		/* (O) calculated vector payoff values at the start of the step, not allocated here */

#endif /* STEPS_H */
