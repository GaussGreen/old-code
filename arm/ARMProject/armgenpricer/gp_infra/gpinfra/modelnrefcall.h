/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file modelnrefcall.h
 *  \brief model and reference call classes
 *         useful for the pricing adviser
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */


#ifndef _INGPINFRA_MODELNREFCALL_H
#define _INGPINFRA_MODELNREFCALL_H

/// this headers has to come first
/// as it contains pragma for stupid warning
#include "gpbase/removeidentifiedwarning.h"

/// kernel
#include <glob/dates.h>

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix
#include "gpbase/env.h"
#include "gpbase/port.h"

/// unix does not support sstream
#include "gpbase/ostringstream.h"
#include "gpbase/pair.h"

#include "typedef.h"
#include "gramnode.h"	/// for accessors!

#include <vector>
CC_USING_NS( std, vector )
CC_USING_NS( std, pair )
//CC_USING_NS_BI_T( std, pair, string, ARM::ARM_ExpNodePtr )

#include <map>
CC_USING_NS( std, map )
CC_USING_NS( std, less )

CC_BEGIN_NAMESPACE( ARM )

///	Some helper classes
/// There isn't much in it but it gives less
/// work to the pricing Adviser which can 
/// focuss on more fundamental things


///////////////////////////////////////////////////////
/// \class ARM_CellRefCall
/// \brief  Very simple class for cell referencing
///		the design is to have a vector of parent nodes
///		that all share the child node..
///////////////////////////////////////////////////////
class ARM_CellRefCall
{
private:
	ARM_ExpNodePtr itsRefNode;
	int itsPayoffPosition;
	
public:
	/// inline for fast access
	inline ARM_CellRefCall( 
		const ARM_ExpNodePtr&	refNode, 
		int position			= -1 )
	:
		itsRefNode( refNode ), 
		itsPayoffPosition( position )
	{}

	ARM_ExpNodePtr AddRerenceNGetRefNode( const ARM_RowColCoords& pCoords );
	
	/// inline for fast access!
	inline void SetPayoffPosition( int nb ) {itsPayoffPosition=nb;}

	/// casting in the good type
	ARM_ExpNodeRef* GetCastedRefNode() const { return  (ARM_ExpNodeRef*) &*itsRefNode; }

	/// accessors
	inline ARM_Date GetChildDate() const { return GetCastedRefNode()->GetChildDate();}
	inline ARM_Date GetParentDate() const { return GetCastedRefNode()->GetParentDate();}
	inline int GetPayoffPosition() const { return itsPayoffPosition; }
	inline ARM_ExpNodePtr GetRefNode() const { return itsRefNode; }

	/// stringify the object!
	string toString( const string& indent = "" ) const;
};

///////////////////////////////////////////////////////
/// operator== to compare two cellRefCalls
///////////////////////////////////////////////////////
bool operator==( const ARM_CellRefCall& lhs, const ARM_CellRefCall& rhs );



///////////////////////////////////////////////////////
/// \class ARM_ModelCall
/// \brief
/// Very simple nested class for model call
///////////////////////////////////////////////////////
class ARM_ModelCall
{
public: 
	/// public to allow the operator< to use this type
	typedef pair< string, ARM_ExpNodePtr > strNodePair;

private:
	ARM_Date itsCallDate;
	typedef vector< strNodePair > strNodePairVec;
	typedef map< string, strNodePairVec, less< string > > strStrNodePairVecMap;
	strStrNodePairVecMap itsModelFunctionNNodeMap;

public:
	/// constructor
	ARM_ModelCall( const ARM_Date& cDate, const string& modelName, 
		const string& pricingFunc, const ARM_ExpNodePtr fNode );

	/// add a model checking for the modelName
	void AddModelCall( const ARM_Date& cDate, const string& modelName, 
		const string& pricingFunc, const ARM_ExpNodePtr fNode );

	/// function to sort function name and remove duplicates!
	void SortAndRemoveDuplicates();

	/// function that loops over the various model names and pricing function
	/// to get the various used time lags!
	/// the result is not sorted and with potential duplicates!
	pair<ARM_VectorPtr,ARM_AdditionalTimeInfoPtrVector> GetUsedTimeLags( ARM_PricingModel* model );

	// Reset the node after the GetUsedTimeLags
	void Reset( ARM_ExpNode::ResetLevel resetParam );

	/// function to stringify the object!
	string toString( const string& indent = "" ) const;

	inline ARM_Date GetCallDate() const { return itsCallDate; }
};
   

CC_END_NAMESPACE()


CC_BEGIN_NAMESPACE( std )

/// operator for the sorting in the std namespace!
inline bool operator<( const CC_NS(ARM,ARM_ModelCall)::strNodePair& lhs, const CC_NS(ARM,ARM_ModelCall)::strNodePair& rhs )
{
	return lhs.first < rhs.first;
}

CC_END_NAMESPACE()


#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

