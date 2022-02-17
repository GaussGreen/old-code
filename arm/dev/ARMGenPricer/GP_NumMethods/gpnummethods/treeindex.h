/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 * \file treeindex.h
 *
 *  \brief
 *
 *	\author  JM Prie
 *	\version 1.0
 *	\date December 2004
 */


#ifndef _INGPNUMMETHODS_TREEINDEX_H
#define _INGPNUMMETHODS_TREEINDEX_H

#include "gpbase/port.h"
#include "gpbase/rootobject.h"
#include "gpbase/typedef.h"
#include "gpbase/gpvector.h"

CC_BEGIN_NAMESPACE( ARM )

//////////////////////////////////////////////
/// \class ARM_TreeIndex
/// \brief
/// ARM_TreeIndex class is used to locate
/// a point in a multi-dimensional grid
/// in a symetrical domain [-max,+max]
//////////////////////////////////////////////
class ARM_TreeIndex : public ARM_RootObject
{
private:
    ARM_IntVector   itsIndex;
    ARM_IntVector   itsMinIndex;
    ARM_IntVector   itsMaxIndex;
    size_t          itsPosition;

    inline size_t Position() const;

    void CopyNoCleanUp(const ARM_TreeIndex& rhs);

public:
    ARM_TreeIndex() : itsPosition(0) {};
    ARM_TreeIndex(int size)
        : itsIndex(size,0),itsMinIndex(size,0),itsMaxIndex(size,0),itsPosition(0) {}
    ARM_TreeIndex(const ARM_IntVector& minIndex,const ARM_IntVector& maxIndex)
        : itsIndex(minIndex),itsMinIndex(minIndex),itsMaxIndex(maxIndex),itsPosition(0) {}
    ARM_TreeIndex(const ARM_TreeIndex& rhs);
    ARM_TreeIndex& operator = (const ARM_TreeIndex& rhs);
    virtual ~ARM_TreeIndex();

    size_t GetPosition() const {return itsPosition;};
    virtual void UpdatePosition() {itsPosition=Position();}

    int size() const {return itsIndex.size();};

    void SetIndex(const ARM_IntVector& index);
    const ARM_IntVector& GetIndex() const { return itsIndex; }

    inline size_t Range() const;
    inline void Reset();
    inline bool IsOnHedge();
    inline size_t  NthIncr(size_t n);

    inline ARM_TreeIndex& operator ++ ();
    inline ARM_TreeIndex operator - ();
    inline bool operator <= (const ARM_TreeIndex& rhs) const;
    inline bool operator < (int index) const;
    int& operator [] (int i) {return itsIndex[i];}
    int operator [] (int i) const {return itsIndex[i];}

	/// Standard ARM Support
    virtual ARM_Object* Clone() const { return new ARM_TreeIndex(*this); }
    virtual string toString(const string& indent="",const string& nextIndent="") const;
};

////////////////////////////////////////////////////
///	Class  : ARM_TreeIndex
///	Routine: Position
///	Returns: Its position
///	Action : Compute the position in a 1D vector
///          saving all the index of a grid
////////////////////////////////////////////////////
inline size_t ARM_TreeIndex::Position() const
{
    int relPos;
	size_t nb=0;
    for(size_t i=0,offset=1;i<itsIndex.size();++i)
    {
        if((relPos=itsIndex[i]-itsMinIndex[i]) > 0)
            nb += offset*relPos;
        if(i<itsIndex.size()-1)
            offset *= itsMaxIndex[i] - itsMinIndex[i] + 1;
    }

    return nb;
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeIndex
///	Routine: Range
///	Returns: Its range
///	Action : Compute the maximum value of its
///          position
////////////////////////////////////////////////////
inline size_t ARM_TreeIndex::Range() const
{
    for(size_t i=0,nb=1;i<itsIndex.size();++i)
        nb *= itsMaxIndex[i] - itsMinIndex[i] + 1;

    return nb;
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeIndex
///	Routine: Reset
///	Returns: 
///	Action : Reset the index to its mimimum value
////////////////////////////////////////////////////
inline void ARM_TreeIndex::Reset()
{
    for(int i=0;i<itsIndex.size();++i)
        itsIndex[i]=itsMinIndex[i];
    itsPosition=0;
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeIndex
///	Routine: IsOnHedge
///	Returns: bool
///	Action : Test if the index is on the ND hedge
////////////////////////////////////////////////////
inline bool ARM_TreeIndex::IsOnHedge()
{
    for(size_t i=0;i<itsIndex.size();++i)
    {
        if(itsIndex[i] == itsMinIndex[i] || itsIndex[i] == itsMaxIndex[i])
            return true;
    }
    return false;
}

////////////////////////////////////////////////////
///	Class  : ARM_TreeIndex
///	Routine: NthIncr
///	Returns: int
///	Action : Incr by 1 the nth index
////////////////////////////////////////////////////
inline size_t ARM_TreeIndex::NthIncr(size_t n)
{
    if(n>size()-1) n=size()-1;
    itsIndex[n] += 1;
    for(size_t i=0,nb=1;i<n;++i)
        nb *= itsMaxIndex[i] - itsMinIndex[i] + 1;

    itsPosition += nb;

    return itsPosition;
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeIndex
///	Routine: operator ++
///	Returns: itself
///	Action : Increment by one the index
////////////////////////////////////////////////////
inline ARM_TreeIndex& ARM_TreeIndex::operator ++ ()
{
    for(int i=0;i<itsIndex.size();++i)
    {
        if(itsIndex[i]<itsMaxIndex[i] || i==itsIndex.size()-1)
        {
            ++(itsIndex[i]);
            break;
        }
        else
            itsIndex[i] = itsMinIndex[i];
    }

    ++itsPosition;

    return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeIndex
///	Routine: operator -
///	Returns: itself
///	Action : Return its opposite
////////////////////////////////////////////////////
inline ARM_TreeIndex ARM_TreeIndex::operator - ()
{
    ARM_TreeIndex oppIndex(*this);
    for(int i=0;i<itsIndex.size();++i)
        oppIndex[i]=-itsIndex[i];

    oppIndex.UpdatePosition();

    return oppIndex;
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeIndex
///	Routine: operator <=
///	Returns: itself
///	Action : Test if the index is lower than another
///          Indexes are ordered from left to right
////////////////////////////////////////////////////
inline bool ARM_TreeIndex::operator <= (const ARM_TreeIndex& rhs) const
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
///	Class  : ARM_TreeIndex
///	Routine: operator <
///	Returns: itself
///	Action : Test if the index is lower than a value
///          identical for all directions
////////////////////////////////////////////////////
inline bool ARM_TreeIndex::operator < (int index) const
{
    for(int i=0;i<itsIndex.size();++i)
    {
        if(itsIndex[i] >= index)
            return false;
    }
    return true;
}

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

