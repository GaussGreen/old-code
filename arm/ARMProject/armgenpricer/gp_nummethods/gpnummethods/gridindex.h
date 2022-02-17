/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 * \file gridindex.h
 *
 *  \brief
 *
 *	\author  JM Prie
 *	\version 1.0
 *	\date October 2003
 */


#ifndef _INGPNUMMETHODS_GRIDINDEX_H
#define _INGPNUMMETHODS_GRIDINDEX_H

#include "gpbase/port.h"
#include "gpbase/rootobject.h"

#include "typedef.h"

#include "gpbase/gpvector.h"

CC_BEGIN_NAMESPACE( ARM )

class ARM_UnsymGridIndex;

//////////////////////////////////////////////
/// \class ARM_GridIndex
/// \brief
/// ARM_GridIndex class is used to locate
/// a point in a multi-dimensional grid
/// in a symetrical domain [-max,+max]
//////////////////////////////////////////////
class ARM_GridIndex : public ARM_RootObject
{
private:
    ARM_IntVector   itsIndex;
    ARM_IntVector   itsMaxIndex;

    inline virtual size_t Position() const;

    void CopyNoCleanUp(const ARM_GridIndex& rhs);

protected :
    /// For efficiency purpose
    size_t          itsPosition;

public:
    ARM_GridIndex();
    ARM_GridIndex(const ARM_GridIndex& rhs);
    ARM_GridIndex(int size);
    ARM_GridIndex(const ARM_IntVector& index);
    virtual ~ARM_GridIndex();
    ARM_GridIndex& operator = (const ARM_GridIndex& rhs);

    size_t GetPosition() const {return itsPosition;};
    virtual void UpdatePosition() {itsPosition=Position();}

    int size() const {return itsIndex.size();};

    void SetIndex(const ARM_GridIndex& rhs);
    void SetIndex(const ARM_IntVector& rhs);

    void SetMaxIndex(const ARM_GridIndex& maxIndex);
    int& MaxIndex(int i) {return itsMaxIndex[i];}
    int MaxIndex(int i) const {return itsMaxIndex[i];}

    inline virtual size_t Range() const;
    inline virtual void Reset();
    void SetToMax(const ARM_GridIndexVector& indexList);

    void RelativeIndex(const ARM_GridIndex& RefIndex,int shift);
    void RelativeIndex(const ARM_UnsymGridIndex& RefIndex,int shift);

    inline size_t AbsolutePosition(const ARM_GridIndex& RefIndex) const;
    inline size_t AbsolutePosition(const ARM_UnsymGridIndex& RefIndex) const;

    inline virtual ARM_GridIndex& operator ++ ();
    inline ARM_GridIndex operator - ();
    inline bool operator <= (const ARM_GridIndex& rhs) const;
    inline bool operator < (int index) const;
    int& operator [] (int i) {return itsIndex[i];}
    int operator [] (int i) const {return itsIndex[i];}

	/// Standard ARM Support
	virtual ARM_Object* Clone() const;
    virtual string toString(const string& indent="",const string& nextIndent="") const;
};


//////////////////////////////////////////////
/// \class ARM_UnsymGridIndex
/// \brief
/// ARM_UnsymGridIndex class is used to locate
/// a point in a multi-dimensional grid
/// in a non symetrical domain [min,max]
//////////////////////////////////////////////
class ARM_UnsymGridIndex : public ARM_GridIndex
{
private:
    ARM_IntVector   itsMinIndex;

    inline virtual size_t Position() const;

    void CopyNoCleanUp(const ARM_UnsymGridIndex& rhs);

public:
    ARM_UnsymGridIndex();
    ARM_UnsymGridIndex(const ARM_UnsymGridIndex& rhs);
    ARM_UnsymGridIndex(int size);
    ARM_UnsymGridIndex(const ARM_IntVector& minIndex, const ARM_IntVector& maxIndex);
    virtual ~ARM_UnsymGridIndex();
    ARM_UnsymGridIndex& operator = (const ARM_UnsymGridIndex& rhs);

    void SetMinIndex(const ARM_GridIndex& minIndex);
    void SetRangeIndex(const ARM_GridIndex& minIndex,const ARM_GridIndex& maxIndex);

    int& MinIndex(int i) {return itsMinIndex[i];}
    int MinIndex(int i) const {return itsMinIndex[i];}

    inline virtual size_t Range() const;
    inline virtual void Reset();

    inline virtual ARM_GridIndex& operator ++ ();
    inline virtual ARM_GridIndex operator - ();

	/// Standard ARM Support
	virtual ARM_Object* Clone() const;
};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// ARM_GridIndex
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
////////////////////////////////////////////////////
///	Class  : ARM_GridIndex
///	Routine: Position
///	Returns: Its position
///	Action : Compute the position in a 1D vector
///          saving all the index of a grid
////////////////////////////////////////////////////
inline size_t ARM_GridIndex::Position() const
{
    size_t offset=1;
    int relPos;
    size_t nb=0;
    for(int i=0;i<itsIndex.size();++i)
    {
        if((relPos=itsIndex[i]+itsMaxIndex[i]) > 0)
            nb += offset*relPos;
        if(i<itsIndex.size()-1)
            offset *= (itsMaxIndex[i]<<1) + 1;
    }

    return nb;
}


////////////////////////////////////////////////////
///	Class  : ARM_GridIndex
///	Routine: Range
///	Returns: Its range
///	Action : Compute the maximum value of its
///          position
////////////////////////////////////////////////////
inline size_t ARM_GridIndex::Range() const
{
    size_t nb=1;
    for(int i=0;i<itsIndex.size();++i)
        nb *= (itsMaxIndex[i]<<1) + 1;

    return nb;
}


////////////////////////////////////////////////////
///	Class  : ARM_GridIndex
///	Routine: Reset
///	Returns: 
///	Action : Reset the index to its mimimum value
////////////////////////////////////////////////////
inline void ARM_GridIndex::Reset()
{
    for(int i=0;i<itsIndex.size();++i)
        itsIndex[i]=-itsMaxIndex[i];
    itsPosition=0;
}


////////////////////////////////////////////////////
///	Class  : ARM_GridIndex
///	Routine: AbsolutePosition
///	Returns: 
///	Action : Give the absolute position of the index
///          assumed to be relative to a reference one
////////////////////////////////////////////////////
inline size_t ARM_GridIndex::AbsolutePosition(const ARM_GridIndex& RefIndex) const
{
    size_t offset=1;
    int relPos,nb=0;
    for(int i=0;i<itsIndex.size();++i)
    {
        if((relPos = itsIndex[i]+RefIndex[i]+RefIndex.MaxIndex(i)) > 0)
            nb += offset*relPos;

        if(i<itsIndex.size()-1)
            offset *= (RefIndex.MaxIndex(i)<<1) + 1;
    }

    return nb;
}


////////////////////////////////////////////////////
///	Class  : ARM_GridIndex
///	Routine: AbsolutePosition
///	Returns: 
///	Action : Give the absolute position of the index
///          assumed to be relative to an unsymetrical
///          reference index
////////////////////////////////////////////////////
inline size_t ARM_GridIndex::AbsolutePosition(const ARM_UnsymGridIndex& RefIndex) const
{
    size_t offset=1;
    int relPos,nb=0;
    for(int i=0;i<itsIndex.size();++i)
    {
        if((relPos = itsIndex[i]+RefIndex[i]-RefIndex.MinIndex(i)) > 0)
            nb += offset*relPos;

        if(i<itsIndex.size()-1)
            offset *= RefIndex.MaxIndex(i) - RefIndex.MinIndex(i) + 1;
    }

    return nb;
}


////////////////////////////////////////////////////
///	Class  : ARM_GridIndex
///	Routine: operator ++
///	Returns: itself
///	Action : Increment by one the index
////////////////////////////////////////////////////
inline ARM_GridIndex& ARM_GridIndex::operator ++ ()
{
    for(int i=0;i<itsIndex.size();++i)
    {
        if(itsIndex[i]<itsMaxIndex[i] || i==itsIndex.size()-1)
        {
            ++(itsIndex[i]);
            break;
        }
        else
            itsIndex[i] = -itsMaxIndex[i];
    }

    ++itsPosition;

    return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_GridIndex
///	Routine: operator -
///	Returns: itself
///	Action : Return its opposite
////////////////////////////////////////////////////
inline ARM_GridIndex ARM_GridIndex::operator - ()
{
    ARM_GridIndex oppIndex(*this);
    for(int i=0;i<itsIndex.size();++i)
        oppIndex[i]=-itsIndex[i];

    oppIndex.UpdatePosition();

    return oppIndex;
}


////////////////////////////////////////////////////
///	Class  : ARM_GridIndex
///	Routine: operator <=
///	Returns: itself
///	Action : Test if the index is lower than another
///          Indexes are ordered from left to right
////////////////////////////////////////////////////
inline bool ARM_GridIndex::operator <= (const ARM_GridIndex& rhs) const
{
    if(itsIndex.size() == rhs.size())
    {
        for(int i=itsIndex.size()-1;i>=0;--i)
        {
            if(itsIndex[i]<rhs[i])
                return true;
            else if(itsIndex[i]>rhs[i])
                return false;
        }
        return true;
    }
    else if(itsIndex.size() < rhs.size())
        return true;
    else
        return false;
}


////////////////////////////////////////////////////
///	Class  : ARM_GridIndex
///	Routine: operator <
///	Returns: itself
///	Action : Test if the index is lower than a value
///          identical for all directions
////////////////////////////////////////////////////
inline bool ARM_GridIndex::operator < (int index) const
{
    for(int i=0;i<itsIndex.size();++i)
    {
        if(itsIndex[i] >= index)
            return false;
    }
    return true;
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
///	Routine: Position
///	Returns: Its position
///	Action : Compute the position in a 1D vector
///          saving all the index of a grid
////////////////////////////////////////////////////
inline size_t ARM_UnsymGridIndex::Position() const
{
    size_t offset=1;
    int relPos;
    size_t nb=0;
    for(int i=0;i<size();++i)
    {
        if((relPos=(*this)[i]-itsMinIndex[i]) > 0)
            nb += offset*relPos;
        if(i<size()-1)
            offset *= MaxIndex(i) - itsMinIndex[i] + 1;
    }

    return nb;
}


////////////////////////////////////////////////////
///	Class  : ARM_UnsymGridIndex
///	Routine: Range
///	Returns: Its range
///	Action : Compute the maximum value of its
///          position
////////////////////////////////////////////////////
inline size_t ARM_UnsymGridIndex::Range() const
{
    size_t nb=1;
    for(int i=0;i<size();++i)
        nb *= MaxIndex(i) - itsMinIndex[i] + 1;

    return nb;
}


////////////////////////////////////////////////////
///	Class  : ARM_UnsymGridIndex
///	Routine: Reset
///	Returns: 
///	Action : Reset the index to its mimimum value
////////////////////////////////////////////////////
inline void ARM_UnsymGridIndex::Reset()
{
    for(int i=0;i<size();++i)
        (*this)[i]=itsMinIndex[i];
    itsPosition=0;
}


////////////////////////////////////////////////////
///	Class  : ARM_UnsymGridIndex
///	Routine: operator ++
///	Returns: itself
///	Action : Increment by one the index
////////////////////////////////////////////////////
inline ARM_GridIndex& ARM_UnsymGridIndex::operator ++ ()
{
    for(int i=0;i<size();++i)
    {
        if((*this)[i]<MaxIndex(i) || i==size()-1)
        {
            ++((*this)[i]);
            break;
        }
        else
            (*this)[i] = itsMinIndex[i];
    }

    ++itsPosition;

    return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_UnsymGridIndex
///	Routine: operator -
///	Returns: itself
///	Action : Return its opposite with min/max
///          constraint
////////////////////////////////////////////////////
inline ARM_GridIndex ARM_UnsymGridIndex::operator - ()
{
    ARM_UnsymGridIndex oppIndex(*this);
    for(int i=0;i<size();++i)
    {
        oppIndex[i]=-(*this)[i];
        if(oppIndex[i] < itsMinIndex[i])
            oppIndex[i] = itsMinIndex[i];
        if(oppIndex[i] > MaxIndex(i))
            oppIndex[i] = MaxIndex(i);
    }

    oppIndex.UpdatePosition();

    return oppIndex;
}



CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

