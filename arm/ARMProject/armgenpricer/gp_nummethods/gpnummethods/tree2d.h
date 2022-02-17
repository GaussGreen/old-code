/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file tree2d.h
 *
 *  \brief
 *
 *	\author  JM Prie, E. Benhamou
 *	\version 1.0
 *	\date November 2004
 */


#ifndef _INGPNUMMETHODS_TREE2D_H
#define _INGPNUMMETHODS_TREE2D_H

#include "gpbase/port.h"
#include "treebase.h"


CC_BEGIN_NAMESPACE( ARM )

class ARM_Tree2D : public ARM_TreeBase
{
protected:
    virtual void ExerciseSmoothing(const ARM_GP_VectorPtr& exerFct, ARM_GP_VectorPtr& smoothValues, ARM_GP_VectorPtr& exerStates) const;

public:
	ARM_Tree2D( const ARM_SamplerBase* sampler, const ARM_TruncatorBase* truncator, const ARM_ReconnectorBase* reconnector, const ARM_SmootherBase* smoother, bool computeSpotProbas=false);
    ARM_Tree2D(const ARM_Tree2D& rhs);
    ARM_Tree2D& operator=(const ARM_Tree2D& rhs );
    virtual ~ARM_Tree2D() {}

    /// Trinomial 2D transition
    virtual size_t GetTransitionSize() const { return 9; }

    virtual void ComputeStateVariables( int timeIdx, ARM_PricingStatesPtr& states ) const;

    virtual size_t dim() const { return 2; }

	/// Standard ARM support
	virtual ARM_Object* Clone() const { return new ARM_Tree2D(*this); }
    virtual string toString(const string& indent="",const string& nextIndent="") const;
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

