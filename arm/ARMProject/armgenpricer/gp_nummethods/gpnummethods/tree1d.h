/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file tree1d.h
 *
 *  \brief
 *
 *	\author  JM Prie, E. Benhamou
 *	\version 1.0
 *	\date November 2004
 */


#ifndef _INGPNUMMETHODS_TREE1D_H
#define _INGPNUMMETHODS_TREE1D_H

#include "gpbase/port.h"
#include "treebase.h"

CC_BEGIN_NAMESPACE( ARM )

class ARM_Tree1D : public ARM_TreeBase
{
public:
	ARM_Tree1D( const ARM_SamplerBase* sampler, const ARM_TruncatorBase* truncator, const ARM_ReconnectorBase* reconnector, const ARM_SmootherBase* smoother, bool computeSpotProbas=false);
    ARM_Tree1D(const ARM_Tree1D& rhs);
    ARM_Tree1D& operator=(const ARM_Tree1D& rhs );
    virtual ~ARM_Tree1D() {}

    /// Trinomial 1D transition
    virtual size_t GetTransitionSize() const { return 3; }

    virtual void ComputeStateVariables( int timeIdx, ARM_PricingStatesPtr& states ) const;

    virtual size_t dim() const { return 1; }

	/// Standard ARM support
	virtual ARM_Object* Clone() const { return new ARM_Tree1D(*this); }
    virtual string toString(const string& indent="",const string& nextIndent="") const;
};



CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

