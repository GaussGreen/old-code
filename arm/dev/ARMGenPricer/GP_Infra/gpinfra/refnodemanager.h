/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file refnodemanager.h
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date December 2003
 */


#ifndef _INGPINFRA_REFNODEMANAGER_H
#define _INGPINFRA_REFNODEMANAGER_H

/// this headers has to come first
/// as it contains pragma for stupid warning
#include "gpbase/removeidentifiedwarning.h"

#include "gpbase/port.h"
#include "typedef.h"

#include "gramnode.h"
#include <glob/dates.h>
#include "modelnrefcall.h"


CC_BEGIN_NAMESPACE( ARM )


//////////////////////////////////////////////////
/// \struct the reference node manager
/// object is designed to keep track of
///		1) all the reference nodes sorted by the child node coordinate of the date of the parent node caller
///			this data structure is a map to benefit from the balanced binary tree performance
///			for lookup.. this structure keeps growing when parsing the generic security and allows
///			to share nodes!
///		2) the potential circular references for the current row of a deal description
//////////////////////////////////////////////////
class ARM_RefNNode_Manager
{
public:
	/// a few simple typedefs and corresponding functor
	/// for the data type of Shared Nodes!

	/// sharedNodeInfo is a pair of the child node cell coordinate and the corresponding 
	///	date of the parent node caller
	/// typedef CC_NS( std, pair )< ARM_RowColCoords, ARM_Date > ARM_SharedNodeInfo;
	/// defined in typedef.h

	/// this functor is used to sort our set of shared node info
	/// it says that a ARM_SharedNodeInfo is smaller if
	///		- its child node coordinate is strictly smaller
	///		- its child node coordinate is the same but its date of its parent node caller is smaller
	struct SharedNodeLess : CC_NS(std,binary_function)<ARM_SharedNodeInfo, ARM_SharedNodeInfo, bool>
	{
		inline bool operator()(const ARM_SharedNodeInfo& lhs, const ARM_SharedNodeInfo& rhs) const
		{	
			return	lhs.first < rhs.first ||
				(lhs.first == rhs.first && lhs.second < rhs.second); 	
		}
	};

	/// map for all shared nodes
	typedef map<ARM_SharedNodeInfo,ARM_CellRefCallPtr,SharedNodeLess> ARM_SharedNodeMap;
	typedef CC_NS( std, pair )< const ARM_SharedNodeInfo, ARM_CellRefCallPtr > ARM_SharedNodeElem;

	/// set of potential circular references!
	typedef set<ARM_RowColCoords> ARM_RowColCoordsSet;
	
	/// constructor
	ARM_RefNNode_Manager( const ARM_SharedNodeMap& sharedNodeMap = ARM_SharedNodeMap(),
		const ARM_RowColCoordsSet& circularRefCoords = ARM_RowColCoordsSet() )
	:
		itsSharedNodeMap( sharedNodeMap ),
		itsCircularRefCoords( circularRefCoords )
	{}

	/// function to reset the set of potential circular references!
	inline void ResetCircularReferences() {	itsCircularRefCoords	= ARM_RowColCoordsSet(); }

	/// Get accessors
	inline ARM_SharedNodeMap::iterator GetSharedNodeMapBegin() { return itsSharedNodeMap.begin(); }
	inline ARM_SharedNodeMap::iterator GetSharedNodeMapEnd() { return itsSharedNodeMap.end(); }
	/// to avoid copy we use this reference return ... This is only for performance reason
	/// and should by no mean used as a reference of good programming style!
	inline ARM_RowColCoordsSet& GetCircularRefCoords() { return itsCircularRefCoords; }


	/// functions for sharedNodes
	inline ARM_SharedNodeMap::iterator findSharedNode( const ARM_SharedNodeInfo& sharedNodeInfo ) { return itsSharedNodeMap.find( sharedNodeInfo); }
	inline bool InsertSharedNode( const ARM_SharedNodeInfo& sharedNodeInfo, const ARM_CellRefCallPtr& cellRefCall )
		{ return itsSharedNodeMap.insert( ARM_SharedNodeElem( sharedNodeInfo,cellRefCall) ).second; }

	/// functions for the circular references!
	inline void SetCircularRefCoords( const set<ARM_RowColCoords>& crc ) { itsCircularRefCoords= crc; }
	inline bool InsertCircularRefCoords( const ARM_RowColCoords& coords ) { return itsCircularRefCoords.insert( coords ).second; }

private:
	/// the various member data
	ARM_SharedNodeMap		itsSharedNodeMap;
	ARM_RowColCoordsSet		itsCircularRefCoords;

	/// disallow copy construction nor assignment!
	ARM_RefNNode_Manager( const ARM_RefNNode_Manager& );
	ARM_RefNNode_Manager& operator=( const ARM_RefNNode_Manager& );
};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

