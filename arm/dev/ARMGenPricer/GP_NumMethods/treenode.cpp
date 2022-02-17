/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file treenode.cpp
 *
 *  \brief
 *
 *	\author  JM Prie E Benhamou
 *	\version 1.0
 *	\date November 2004
 */


#include "gpnummethods/treenode.h"

/// gpbase
#include "gpbase/env.h"
#include "gpbase/cloneutilityfunc.h"

/// gpinfra
#include "gpinfra/pricingstates.h"

/// gpnummethods
#include "gpnummethods/slice.h"

/// ARM Kernel
#include <glob/expt.h>


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
////////////////////////////////////////////////////
///	ARM_NodeBase
////////////////////////////////////////////////////
////////////////////////////////////////////////////

ARM_Node1D* ARM_NodeBase::ToNode1D()
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": cannot cast to Node1D" );
}

ARM_NodeND* ARM_NodeBase::ToNodeND()
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": cannot cast to NodeND" );
}

//////////////////////////////////////
/// node N dimensional
//////////////////////////////////////


////////////////////////////////////////////////////
///	Class  : ARM_NodeND
///	Routine: copy constructor assignement operator, destructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_NodeND::ARM_NodeND( const ARM_NodeND& rhs )
:ARM_NodeBase(rhs), itsNodes( rhs.itsNodes.size() )
{
	DuplicateCloneablePointorAndNullVectorInPlace<ARM_Node1D>( rhs.itsNodes, itsNodes );
}

ARM_NodeND& ARM_NodeND::operator=( const ARM_NodeND& rhs )
{
	if( this != &rhs )
	{
		ARM_NodeBase::operator=(rhs);
		DeletePointorVector<ARM_Node1D>(itsNodes);
		DuplicateCloneablePointorAndNullVectorInPlace<ARM_Node1D>( rhs.itsNodes, itsNodes );
	}
	return *this;
}

ARM_NodeND::~ARM_NodeND()
{
	DeletePointorVector<ARM_Node1D>(itsNodes);
}

CC_END_NAMESPACE()



/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

