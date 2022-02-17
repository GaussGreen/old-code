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


#ifndef _INGPNUMMETHODS_TREE3D_H
#define _INGPNUMMETHODS_TREE3D_H

#include "gpbase/port.h"
#include "treebase.h"


CC_BEGIN_NAMESPACE( ARM )

class ARM_Tree3D : public ARM_TreeBase
{
    static int SmoothingDirection[3][3];

protected:
    virtual void ExerciseSmoothing(const ARM_GP_VectorPtr& exerFct, ARM_GP_VectorPtr& smoothValues, ARM_GP_VectorPtr& exerStates) const;

public:
	ARM_Tree3D( const ARM_SamplerBase* sampler, const ARM_TruncatorBase* truncator, const ARM_ReconnectorBase* reconnector, const ARM_SmootherBase* smoother, bool computeSpotProbas=false);
    ARM_Tree3D(const ARM_Tree3D& rhs);
    ARM_Tree3D& operator=(const ARM_Tree3D& rhs );
    virtual ~ARM_Tree3D() {}

    /// Trinomial 3D transition
    virtual size_t GetTransitionSize() const { return 27; }

    virtual void ComputeStateVariables( int timeIdx, ARM_PricingStatesPtr& states ) const;

    virtual size_t dim() const { return 3; }

	/// Standard ARM support
	virtual ARM_Object* Clone() const { return new ARM_Tree3D(*this); }
    virtual string toString(const string& indent="",const string& nextIndent="") const;
};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

