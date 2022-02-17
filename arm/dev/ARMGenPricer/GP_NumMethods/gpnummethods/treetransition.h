/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: treetransition.h,v $
 * Revision 1.1  2003/10/13 07:52:17  jmprie
 * Initial revision
 *
 *
 */



/*! \file treetransition.h
 *
 *  \brief
 *
 *	\author  JM Prie
 *	\version 1.0
 *	\date October 2003
 */


#ifndef _INGPINFRA_TREETRANSITION_H
#define _INGPINFRA_TREETRANSITION_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpinfra/typedef.h"
#include "gpbase/rootobject.h"
#include "gridindex.h"


CC_BEGIN_NAMESPACE( ARM )

//////////////////////////////////////////////
/// \class ARM_TreeTransition
/// \brief
/// ARM_TreeTransition class defines
/// the connection between a current state
/// to the next ones.
//////////////////////////////////////////////
class ARM_TreeTransition : public ARM_RootObject
{
private:
    /// The central connected state in the next
    /// slice defined by its tree index
    ARM_GridIndex itsIndex;

    /// 1D transition probabilities list
    /// for each dimension
    ARM_VectorVector itsProbas;

    void CopyNoCleanUp(const ARM_TreeTransition& rhs);
    void CleanUp();

public:
	ARM_TreeTransition();
	ARM_TreeTransition(const ARM_TreeTransition& rhs);
	ARM_TreeTransition(ARM_GridIndex& treeIndex,ARM_VectorVector& probas);
	virtual ~ARM_TreeTransition();

	ARM_TreeTransition& operator=(const ARM_TreeTransition& rhs);

    const ARM_GridIndex& GetIndex() const {return itsIndex;}
    const ARM_VectorVector& GetProbas() const {return itsProbas;}

	/// Standard ARM support
	virtual ARM_Object* Clone() const;
    virtual string toString(const string& indent="", const string& nextIndent="") const { return "ARM_TreeTransition"; }
};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

