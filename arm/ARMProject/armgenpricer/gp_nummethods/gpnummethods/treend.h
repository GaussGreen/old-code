/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file tree3d.h
 *
 *  \brief
 *
 *	\author  JM Prie, E. Benhamou
 *	\version 1.0
 *	\date November 2004
 */


#ifndef _INGPNUMMETHODS_TREEND_H
#define _INGPNUMMETHODS_TREEND_H

#include "gpbase/port.h"
#include "treebase.h"

#include <cmath>


CC_BEGIN_NAMESPACE( ARM )

class ARM_TreeND : public ARM_TreeBase
{
	size_t itsDim;

protected:
    virtual void ExerciseSmoothing(const ARM_GP_VectorPtr& exerFct, ARM_GP_VectorPtr& smoothValues, ARM_GP_VectorPtr& exerStates) const;

public:
	ARM_TreeND( const ARM_SamplerBase* sampler, const ARM_TruncatorBase* truncator, const ARM_ReconnectorBase* reconnector, size_t dim, const ARM_SmootherBase* smoother, bool computeSpotProbas=false );
    ARM_TreeND( const ARM_TreeND& rhs);
    ARM_TreeND& operator=(const ARM_TreeND& rhs );
    virtual ~ARM_TreeND() {}

    /// Trinomial ND transition
// FIXMEFRED: mig.vc8 (25/05/2007 15:55:43):cast
	virtual size_t GetTransitionSize() const { return pow(3.,static_cast<int>(itsDim)); }

    virtual void ComputeStateVariables( int timeIdx, ARM_PricingStatesPtr& states ) const;

    virtual size_t dim() const { return itsDim; }

	/// Standard ARM support
	virtual ARM_Object* Clone() const { return new ARM_TreeND(*this); }
    virtual string toString(const string& indent="",const string& nextIndent="") const;
};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

