/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file treenode.h
 *
 *  \brief
 *
 *	\author  JM Prie, E. Benhamou
 *	\version 1.0
 *	\date November 2004
 */


#ifndef _INGPNUMMETHODS_TREENODE_H
#define _INGPNUMMETHODS_TREENODE_H

#include "gpbase/port.h"
#include "gpbase/rootobject.h"
#include "expt.h"
#include <vector>
CC_USING_NS(std,vector)

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////
/// NODE PART
/// a slice of a tree.
////////////////////////////////////

class ARM_Node1D;
class ARM_NodeND;

struct ARM_NodeBase : public ARM_RootObject
{
    /// downcast to the appriate type
    virtual ARM_Node1D* ToNode1D();
    virtual ARM_NodeND* ToNodeND();
	virtual ~ARM_NodeBase() {};
};


//// node 1 dimensional
class ARM_Node1D : public ARM_NodeBase
{
private:
	double itsProbaDown;
	double itsProbaUp;
    int itsNextNodeIndex;
public:
	/// constructor, copy constructor assignement operator, destructor
	ARM_Node1D( double probaUp = 0.0, double probaDown=0.0, int nextNodeIndex=0 )
	: itsProbaUp(probaUp), itsProbaDown(probaDown), itsNextNodeIndex(nextNodeIndex) 
	{}

	ARM_Node1D( const ARM_Node1D& rhs )
	:	ARM_NodeBase(rhs),  itsProbaDown( rhs.itsProbaDown ), itsProbaUp( rhs.itsProbaUp ), 
		itsNextNodeIndex( rhs.itsNextNodeIndex )
	{}

	ARM_Node1D& operator=( const ARM_Node1D& rhs )
	{
		if( this != &rhs )
		{
			ARM_NodeBase::operator=(rhs);
			itsProbaDown	= rhs.itsProbaDown;
			itsProbaUp		= rhs.itsProbaUp;
			itsNextNodeIndex= rhs.itsNextNodeIndex;
		}
		return *this;

	}
	virtual ~ARM_Node1D() {};


    /// accessor
    double GetProbaDown() const { return itsProbaDown; }
    void SetProbaDown(double proba) { itsProbaDown = proba; }
    double GetProbaUp() const { return itsProbaUp; }
    void SetProbaUp(double proba) { itsProbaUp = proba; }
    int GetNextNodeIndex() const { return itsNextNodeIndex; }
    void SetNextNodeIndex(int nextNodeIndex) { itsNextNodeIndex=nextNodeIndex; }

	/// downcast
    virtual ARM_Node1D* ToNode1D() { return this; }

	/// Standard ARM support
	virtual ARM_Object* Clone() const { return new ARM_Node1D(*this); }
    virtual string toString(const string& indent="",const string& nextIndent="") const { return "ARM_Node1D"; }
};

//////////////////////////////////////
/// node N dimensional
//////////////////////////////////////

class ARM_NodeND : public ARM_NodeBase
{
private:
	vector<ARM_Node1D*> itsNodes;
public:
	/// constructor, copy constructor assignement operator, destructor
    ARM_NodeND( size_t nbDims ): itsNodes(nbDims) { for(size_t i=0;i<nbDims;++i) itsNodes[i] = new ARM_Node1D; }
	ARM_NodeND( const ARM_NodeND& rhs );
	ARM_NodeND& operator=( const ARM_NodeND& rhs );
	virtual ~ARM_NodeND();

    /// accessor
    double GetProbaDown(size_t i) const { return itsNodes[i]->GetProbaDown(); }
    void SetProbaDown(size_t i, double proba) { itsNodes[i]->SetProbaDown(proba); }
    double GetProbaUp(size_t i ) const { return itsNodes[i]->GetProbaUp(); }
    void SetProbaUp(size_t i, double proba) { itsNodes[i]->SetProbaUp(proba); }
    int GetNextNodeIndex(size_t i ) const { return itsNodes[i]->GetNextNodeIndex(); }
    void SetNextNodeIndex(size_t i, int nextNodeIndex) { itsNodes[i]->SetNextNodeIndex(nextNodeIndex); }
	size_t dim() const { return itsNodes.size(); }

	ARM_Node1D* GetNode1D(size_t i) const
	{ 
#if defined(__ARM_VECTOR_NO_RANGE_CHECK)
		if( i>=itsNodes.size() )
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " Out of Range! in GetNode1D!" );
#endif
		return itsNodes[i]; 
	}

	/// downcast
    virtual ARM_NodeND* ToNodeND() { return this; }

	/// Standard ARM support
	virtual ARM_Object* Clone() const { return new ARM_NodeND(*this); }
    virtual string toString(const string& indent="",const string& nextIndent="") const { return "ARM_NodeND"; }
};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

