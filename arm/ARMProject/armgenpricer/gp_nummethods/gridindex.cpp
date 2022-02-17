/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file gridindex.cpp
 *
 *  \brief
 *
 *	\author  JM Prie
 *	\version 1.0
 *	\date October 2003
 */


#include "gpnummethods/gridindex.h"
#include "gpbase/ostringstream.h"
#include <glob/expt.h>

#include <iomanip> /// for setprecision()

CC_BEGIN_NAMESPACE( ARM )


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// ARM_GridIndex
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
////////////////////////////////////////////////////
///	Class  : ARM_GridIndex
///	Routine: CopyNoCleanUp
///	Returns: 
///	Action : Arguments copy
////////////////////////////////////////////////////
void ARM_GridIndex::CopyNoCleanUp(const ARM_GridIndex& rhs)
{
	itsIndex = rhs.itsIndex;
	itsMaxIndex = rhs.itsMaxIndex;
	itsPosition = rhs.itsPosition;
}


////////////////////////////////////////////////////
///	Class  : ARM_GridIndex
///	Routine: default constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_GridIndex::ARM_GridIndex()
: itsPosition(0)
{}


////////////////////////////////////////////////////
///	Class  : ARM_GridIndex
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_GridIndex::ARM_GridIndex(const ARM_GridIndex& rhs)
:	ARM_RootObject( rhs )	
{
    CopyNoCleanUp(rhs);
}


////////////////////////////////////////////////////
///	Class  : ARM_GridIndex
///	Routine: Constructor
///	Returns: 
///	Action : Constructor intializing its size
////////////////////////////////////////////////////
ARM_GridIndex::ARM_GridIndex(int size)
: itsPosition(0)
{
    itsIndex.resize(size);
    itsMaxIndex.resize(size);
    for(int i=0;i<size;++i)
    {
        itsIndex[i]=0;
        itsMaxIndex[i]=0;
    }
}


////////////////////////////////////////////////////
///	Class  : ARM_GridIndex
///	Routine: Constructor
///	Returns: 
///	Action : Constructor intializing its index
////////////////////////////////////////////////////
ARM_GridIndex::ARM_GridIndex(const ARM_IntVector& index)
: itsIndex(index),itsMaxIndex(index)
{
    itsPosition=Position();
}


////////////////////////////////////////////////////
///	Class  : ARM_GridIndex
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_GridIndex::~ARM_GridIndex()
{}


////////////////////////////////////////////////////
///	Class  : ARM_GridIndex
///	Routine: Operator =
///	Returns: 
///	Action : Affectation
////////////////////////////////////////////////////
ARM_GridIndex& ARM_GridIndex::operator = (const ARM_GridIndex& rhs)
{
	if(this != &rhs)
	{
		ARM_RootObject::operator=(rhs);
		CopyNoCleanUp(rhs);
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_GridIndex
///	Routine: SetIndex
///	Returns: 
///	Action : Set values from a rhs index
////////////////////////////////////////////////////
void ARM_GridIndex::SetIndex(const ARM_GridIndex& rhs)
{
    int mySize=size();
    int newSize=rhs.size();
    if(mySize>0 && newSize != mySize)
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
       "GridIndex mismatch in index & max index settings");
    }
    else if(mySize==0)
    {
        itsIndex.resize(newSize);
        itsMaxIndex.resize(newSize);
    }

    /// Not through CopyNoCleanUp()
    /// just affect values... for efficiency
    for(int i=0;i<newSize;++i)
    {
        itsIndex[i]=rhs.itsIndex[i];
        itsMaxIndex[i]=rhs.itsMaxIndex[i];
    }
    itsPosition=rhs.itsPosition;
}

////////////////////////////////////////////////////
///	Class  : ARM_GridIndex
///	Routine: SetIndex
///	Returns: 
///	Action : Set values from a rhs int vector. It is
///          fast version with no update of
///          itsMaxIndex & itsPosition
////////////////////////////////////////////////////
void ARM_GridIndex::SetIndex(const ARM_IntVector& rhs)
{
    int mySize=size();
    int newSize=rhs.size();
    if(mySize>0 && newSize != mySize)
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
       "GridIndex mismatch in index setting");
    }
    else if(mySize==0)
    {
        itsIndex.resize(newSize);
        itsMaxIndex.resize(newSize);
    }

    for(int i=0;i<newSize;++i)
        itsIndex[i]=rhs[i];
}


////////////////////////////////////////////////////
///	Class  : ARM_GridIndex
///	Routine: SetMaxIndex
///	Returns: 
///	Action : Set max values from a rhs index
////////////////////////////////////////////////////
void ARM_GridIndex::SetMaxIndex(const ARM_GridIndex& maxIndex)
{
    int mySize=size();
    int newSize=maxIndex.size();
    if(mySize>0 && newSize != mySize)
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
       "GridIndex mismatch in max index setting");
    }
    else if(mySize==0)
    {
        itsIndex.resize(newSize);
        itsMaxIndex.resize(newSize);
    }

    for(int i=0;i<newSize;++i)
        itsMaxIndex[i]=maxIndex[i];

    Reset();
}


////////////////////////////////////////////////////
///	Class  : ARM_GridIndex
///	Routine: Max
///	Returns: 
///	Action : The index is set to the maximum index
///          in absolute value of the list in input
////////////////////////////////////////////////////
void ARM_GridIndex::SetToMax(const ARM_GridIndexVector& indexList)
{
    int i,iv;
    for(i=0;i<itsIndex.size();++i)
        itsIndex[i] = 0;
    
    for(int j=0;j<indexList.size();++j)
    {
        for(i=0;i<itsIndex.size();++i)
        {
            iv=(*(indexList[j])).itsIndex[i];
            iv = (iv < 0 ? -iv : iv);
            if(iv > itsIndex[i])
                itsIndex[i]=iv;
        }
    }

    itsPosition=1;
    for(i=0;i<itsIndex.size();++i)
    {
        itsMaxIndex[i]=itsIndex[i];
        itsPosition *= (itsMaxIndex[i]<<1) + 1;
    }
    --itsPosition;
}


////////////////////////////////////////////////////
///	Class  : ARM_GridIndex
///	Routine: RelativeIndex
///	Returns: 
///	Action : Modify the index to be relative to a
///          reference index and separeted by the
///          "shift" amount on each coordinate
///          Max value may limit the shift
////////////////////////////////////////////////////
void ARM_GridIndex::RelativeIndex(const ARM_GridIndex& RefIndex,int shift)
{
    size_t offset=1;
    int j;
    itsPosition=0;
    for(int i=0;i<itsIndex.size();++i)
    {
        j=RefIndex.itsIndex[i]+shift;
        itsIndex[i] = ( (j < -RefIndex.MaxIndex(i)) ? -RefIndex.MaxIndex(i) :
                       (j > RefIndex.MaxIndex(i) ? RefIndex.MaxIndex(i) : j) )
                      - RefIndex[i];

        /// Max indexes are also relative
        itsMaxIndex[i] = (itsIndex[i] > 0 ? itsIndex[i] : -itsIndex[i]);

        itsPosition += offset*(itsIndex[i]+itsMaxIndex[i]);

        offset *= (itsMaxIndex[i]<<1) + 1;
    }
}


////////////////////////////////////////////////////
///	Class  : ARM_GridIndex
///	Routine: RelativeIndex
///	Returns: 
///	Action : Modify the index to be relative to a
///          reference index and separeted by the
///          "shift" amount on each coordinate
///          Min and Max values may limit the shift
////////////////////////////////////////////////////
void ARM_GridIndex::RelativeIndex(const ARM_UnsymGridIndex& RefIndex,int shift)
{
    size_t offset=1;
    int j;
    itsPosition=0;
    for(int i=0;i<itsIndex.size();++i)
    {
        j=RefIndex[i]+shift;
        itsIndex[i] = ( (j < RefIndex.MinIndex(i)) ? RefIndex.MinIndex(i) :
                       (j > RefIndex.MaxIndex(i) ? RefIndex.MaxIndex(i) : j) )
                      - RefIndex[i];

        /// Max indexes are also relative
        itsMaxIndex[i] = (itsIndex[i] > 0 ? itsIndex[i] : -itsIndex[i]);

        itsPosition += offset*(itsIndex[i]+itsMaxIndex[i]);

        offset *= (itsMaxIndex[i]<<1) + 1;
    }
}


////////////////////////////////////////////////////
///	Class   : ARM_GridIndex
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_GridIndex::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;

    os << CC_NS(std,setw)(5);

    os << "[";
	int i;
    for(i=0;i<itsMaxIndex.size()-1;++i)
        os << itsMaxIndex[i] << ", ";
    os << itsMaxIndex[i] << "]";  /// no \n because very basic object and the calling function does it
    
    return os.str();
}

////////////////////////////////////////////////////
///	Class   : ARM_GridIndex
///	Routines: Clone,View
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_GridIndex::Clone() const
{
	return new ARM_GridIndex(*this);
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// ARM_UnsymGridIndex
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
////////////////////////////////////////////////////
///	Class  : ARM_UnsymGridIndex
///	Routine: CopyNoCleanUp
///	Returns: 
///	Action : Arguments copy
////////////////////////////////////////////////////
void ARM_UnsymGridIndex::CopyNoCleanUp(const ARM_UnsymGridIndex& rhs)
{
	itsMinIndex = rhs.itsMinIndex;
}


////////////////////////////////////////////////////
///	Class  : ARM_UnsymGridIndex
///	Routine: default constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_UnsymGridIndex::ARM_UnsymGridIndex()
{}


////////////////////////////////////////////////////
///	Class  : ARM_UnsymGridIndex
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_UnsymGridIndex::ARM_UnsymGridIndex(const ARM_UnsymGridIndex& rhs)
:	ARM_GridIndex( rhs )	
{
    CopyNoCleanUp(rhs);
}


////////////////////////////////////////////////////
///	Class  : ARM_UnsymGridIndex
///	Routine: Constructor
///	Returns: 
///	Action : Constructor intializing its size
////////////////////////////////////////////////////
ARM_UnsymGridIndex::ARM_UnsymGridIndex(int size)
: ARM_GridIndex(size)
{
    itsMinIndex.resize(size);
    for(int i=0;i<size;++i)
    {
        itsMinIndex[i]=0;
    }
}


////////////////////////////////////////////////////
///	Class  : ARM_UnsymGridIndex
///	Routine: Constructor
///	Returns: 
///	Action : Constructor intializing its indexes
////////////////////////////////////////////////////
ARM_UnsymGridIndex::ARM_UnsymGridIndex(const ARM_IntVector& minIndex,const ARM_IntVector& maxIndex)
: ARM_GridIndex(maxIndex),itsMinIndex(minIndex)
{
    itsPosition=Position();
}


////////////////////////////////////////////////////
///	Class  : ARM_UnsymGridIndex
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_UnsymGridIndex::~ARM_UnsymGridIndex()
{}


////////////////////////////////////////////////////
///	Class  : ARM_UnsymGridIndex
///	Routine: Operator =
///	Returns: 
///	Action : Affectation
////////////////////////////////////////////////////
ARM_UnsymGridIndex& ARM_UnsymGridIndex::operator = (const ARM_UnsymGridIndex& rhs)
{
	if(this != &rhs)
	{
		ARM_GridIndex::operator=(rhs);
		CopyNoCleanUp(rhs);
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_UnsymGridIndex
///	Routine: SetMinIndex
///	Returns: 
///	Action : Set min values from a rhs index
////////////////////////////////////////////////////
void ARM_UnsymGridIndex::SetMinIndex(const ARM_GridIndex& minIndex)
{
    int mySize=size();
    if(mySize>0 && minIndex.size() != mySize)
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
       "GridIndex mismatch in min index setting");
    }
    else if(minIndex.size() != mySize)
        SetMaxIndex(minIndex);
        
    for(int i=0;i<size();++i)
        itsMinIndex[i]=minIndex[i];

    Reset();
}


////////////////////////////////////////////////////
///	Class  : ARM_UnsymGridIndex
///	Routine: SetMinIndex
///	Returns: 
///	Action : Set min values from a rhs index
////////////////////////////////////////////////////
void ARM_UnsymGridIndex::SetRangeIndex(const ARM_GridIndex& minIndex,const ARM_GridIndex& maxIndex)
{
    for(int i=0;i<size();++i)
    {
        itsMinIndex[i]=minIndex[i];
        MaxIndex(i)=maxIndex[i];
    }
    Reset();
}


////////////////////////////////////////////////////
///	Class   : ARM_UnsymGridIndex
///	Routines: Clone,View
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_UnsymGridIndex::Clone() const
{
	return new ARM_UnsymGridIndex(*this);
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

