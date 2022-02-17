/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *	\file treetransition.cpp
 *
 *  \brief
 *
 *	\author  JM Prie
 *	\version 1.0
 *	\date October 2003
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpnummethods/treetransition.h"
#include "gpbase/gpvector.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_TreeTransition
///	Routine: CopyNoCleanUp
///	Returns: 
///	Action : Arguments copy
////////////////////////////////////////////////////
void ARM_TreeTransition::CopyNoCleanUp(const ARM_TreeTransition& rhs)
{
	itsProbas.resize(rhs.itsProbas.size());
    for(size_t i=0;i<itsProbas.size();++i)
        itsProbas[i] = rhs.itsProbas[i] ? (ARM_GP_Vector*) rhs.itsProbas[i]->Clone() : NULL;
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeTransition
///	Routine: CleanUp
///	Returns: 
///	Action : Arguments destruction
////////////////////////////////////////////////////
void ARM_TreeTransition::CleanUp()
{
    for(size_t i=0;i<itsProbas.size();++i)
    {
        delete itsProbas[i];
        itsProbas[i]=NULL;
    }
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeTransition
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_TreeTransition::ARM_TreeTransition()
{}


////////////////////////////////////////////////////
///	Class  : ARM_TreeTransition
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_TreeTransition::ARM_TreeTransition(const ARM_TreeTransition& rhs)
: ARM_RootObject(rhs)
{
    CopyNoCleanUp(rhs);
}

////////////////////////////////////////////////////
///	Class  : ARM_TreeTransition
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising its mid index
///          and its transition probability matrix
////////////////////////////////////////////////////
ARM_TreeTransition::ARM_TreeTransition(ARM_GridIndex& treeIndex,ARM_VectorVector& probas)
: itsIndex(treeIndex)
{
	itsProbas.resize(probas.size());
    for(size_t i=0;i<itsProbas.size();++i)
        itsProbas[i] = probas[i] ? (ARM_GP_Vector*) probas[i]->Clone() : NULL;
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeTransition
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_TreeTransition::~ARM_TreeTransition()
{
    CleanUp();
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeTransition
///	Routine: operator =
///	Returns: 
///	Action : Affectation
////////////////////////////////////////////////////
ARM_TreeTransition& ARM_TreeTransition::operator = (const ARM_TreeTransition& rhs)
{
	if(this != &rhs)
	{
		ARM_RootObject::operator=(rhs);
		CleanUp();
		CopyNoCleanUp(rhs);
	}
	return *this;
}

////////////////////////////////////////////////////
///	Class  : ARM_TreeTransition
///	Routine: Copy, Clone, View
///	Returns: 
///	Action : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_TreeTransition::Clone() const
{
	return new ARM_TreeTransition(*this);
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

