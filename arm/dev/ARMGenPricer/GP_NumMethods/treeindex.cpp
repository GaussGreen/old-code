/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file treeindex.cpp
 *
 *  \brief
 *
 *	\author  JM Prie
 *	\version 1.0
 *	\date Decembre 2004
 */


#include "gpnummethods/treeindex.h"
#include "gpbase/ostringstream.h"
#include <glob/expt.h>

#include <iomanip> /// for setprecision()

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_TreeIndex
///	Routine: CopyNoCleanUp
///	Returns: 
///	Action : Arguments copy
////////////////////////////////////////////////////
void ARM_TreeIndex::CopyNoCleanUp(const ARM_TreeIndex& rhs)
{
	itsIndex = rhs.itsIndex;
	itsMinIndex = rhs.itsMinIndex;
	itsMaxIndex = rhs.itsMaxIndex;
	itsPosition = rhs.itsPosition;
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeIndex
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_TreeIndex::ARM_TreeIndex(const ARM_TreeIndex& rhs)
:	ARM_RootObject( rhs )	
{
    CopyNoCleanUp(rhs);
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeIndex
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_TreeIndex::~ARM_TreeIndex()
{}


////////////////////////////////////////////////////
///	Class  : ARM_TreeIndex
///	Routine: Operator =
///	Returns: 
///	Action : Affectation
////////////////////////////////////////////////////
ARM_TreeIndex& ARM_TreeIndex::operator = (const ARM_TreeIndex& rhs)
{
	if(this != &rhs)
	{
		ARM_RootObject::operator=(rhs);
		CopyNoCleanUp(rhs);
	}
	return *this;
}

////////////////////////////////////////////////////
///	Class  : ARM_TreeIndex
///	Routine: SetIndex
///	Returns: 
///	Action : Set values from a rhs int vector. It is
///          fast version with no update of
///          itsMaxIndex & itsPosition
////////////////////////////////////////////////////
void ARM_TreeIndex::SetIndex(const ARM_IntVector& index)
{
    int mySize=size();
    int newSize=index.size();
    if(mySize>0 && newSize != mySize)
    {
        ARM_THROW(ERR_INVALID_ARGUMENT,ARM_USERNAME + " : TreeIndex mismatch in index setting" );
    }
    else if(mySize==0)
    {
        itsIndex.resize(newSize);
        itsMinIndex.resize(newSize);
        itsMaxIndex.resize(newSize);
    }

    for(int i=0;i<newSize;++i)
    {
#if defined( __GP_STRICT_VALIDATION)
        if(index[i] < itsMinIndex[i] || index[i] > itsMaxIndex[i])
            ARM_THROW(ERR_INVALID_ARGUMENT,ARM_USERNAME + " : try to set out of [min,max] in TreeIndex " );
#endif
        itsIndex[i]=index[i];
    }

    /// Update position
    itsPosition = Position();
}


////////////////////////////////////////////////////
///	Class   : ARM_TreeIndex
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_TreeIndex::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;

    os << CC_NS(std,setw)(5);

    os << "[";
    for(int i=0;i<itsMaxIndex.size()-1;++i)
        os << itsMaxIndex[i] << ", ";
    os << itsMaxIndex[i] << "]";  /// no \n because very basic object and the calling function does it
    
    return os.str();
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

