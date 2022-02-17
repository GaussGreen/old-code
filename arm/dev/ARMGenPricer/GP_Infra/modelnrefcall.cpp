/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file modelnrefcall.cpp
 *  \brief classes to handle model calls
 *		as well as other cell refence calls
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */

/// this headers has to come first
/// as it contains pragma for stupid warning
#include "gpbase/removeidentifiedwarning.h"
#include "gpinfra/modelnrefcall.h"
#include "gpinfra/gramnode.h"
#include "gpinfra/additionaltimeinfo.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
////////////////////////////////////////////////////
/////			ARM_CellRefCall				////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_CellRefCall
///	Routine: AddRerenceNGetRefNode
///	Returns: ARM_ExpNodePtr
///	Action : Says that the reference node has one more reference
///			and return it;
////////////////////////////////////////////////////
ARM_ExpNodePtr ARM_CellRefCall::AddRerenceNGetRefNode( const ARM_RowColCoords& pCoords )
{ 
	GetCastedRefNode()->IncCountor(pCoords);
	return itsRefNode;
}


////////////////////////////////////////////////////
///	Class  : ARM_CellRefCall
///	Routine: toString
///	Returns: string
///	Action : stringify the object
////////////////////////////////////////////////////
string ARM_CellRefCall::toString( const string& indent ) const
{
	ARM_ExpNodeRef* castedRefNode = GetCastedRefNode();
	CC_Ostringstream os;
	os << indent <<    "Location	: "	<< castedRefNode->GetChildLocationString() << "\n";
	os << indent << "Shared Nodes	: "	<<  castedRefNode->GetCountor() << "\n";
	os << indent << "Parent Date	: " << castedRefNode->GetParentDate().toString() << "\n";

#if defined(__GP_SHOW_SHARED_NODE_COORDINATES)
	os << indent << "Parent Coords    : total nb=" << castedRefNode->GetParentCoords().size() << " with coordinates=";
	for( size_t i=0; i<castedRefNode->GetParentCoords().size(); ++i )
		os << " (" << castedRefNode->GetParentCoords()[i].itsElem1 << "," << castedRefNode->GetParentCoords()[i].itsElem2 << ")";
	os << "\n";
#endif
	os << indent <<  "Child Date	: "	<< castedRefNode->GetChildDate().toString() << "\n";
	return os.str();
}


////////////////////////////////////////////////////
///	Routine: operator==
///	Returns: bool
///	Action : checks whether two cell reference calls are
///			the same!
////////////////////////////////////////////////////
bool operator==( const ARM_CellRefCall& lhs, const ARM_CellRefCall& rhs )
{
	return ( lhs.GetRefNode()  == rhs.GetRefNode() );
}



////////////////////////////////////////////////////
////////////////////////////////////////////////////
/////			ARM_ModelCall				////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////


///////////////////////////////////////////////////
///	Class  : ARM_ModelCall
///	Routine: ARM_ModelCall (constructor)
///	Returns: 
///	Action : build the object
////////////////////////////////////////////////////
ARM_ModelCall::ARM_ModelCall( const ARM_Date& cDate, const string& modelName, 
		const string& pricingFunc, const ARM_ExpNodePtr fNode )
	:	itsCallDate( cDate ), itsModelFunctionNNodeMap()
{
	/// this insert cannot fail because the map is empty
	/// hence no testing!
	itsModelFunctionNNodeMap.insert( pair< const string, strNodePairVec >
		( modelName, strNodePairVec( 1, strNodePair( pricingFunc, fNode ) ) ) );
}




///////////////////////////////////////////////////
///	Class  : ARM_ModelCall
///	Routine: AddModelCall
///	Returns: void
///	Action : Add a model checking for the modelName
////////////////////////////////////////////////////
void ARM_ModelCall::AddModelCall( const ARM_Date& cDate, const string& modelName, 
	const string& pricingFunc, const ARM_ExpNodePtr fNode )
{
#ifdef __GP_STRICT_VALIDATION
	if( cDate != itsCallDate )
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			"Trying to add a model call with not the same event date!" );
#endif

	strStrNodePairVecMap::iterator iter =itsModelFunctionNNodeMap.find( modelName ),
		end = itsModelFunctionNNodeMap.end();

	if( iter != end )
		((*iter).second).push_back( strNodePair( pricingFunc, fNode ) );
	else
		itsModelFunctionNNodeMap.insert( pair< const string, strNodePairVec >
			( modelName, strNodePairVec( 1, strNodePair( pricingFunc, fNode ) ) ) );
}




///////////////////////////////////////////////////
///	Class  : ARM_ModelCall
///	Routine: SortAndRemoveDuplicates
///	Returns: void
///	Action : function to sort function name and remove duplicates!
////////////////////////////////////////////////////
void ARM_ModelCall::SortAndRemoveDuplicates()
{
	strStrNodePairVecMap::iterator 
		iter = itsModelFunctionNNodeMap.begin(),
		end = itsModelFunctionNNodeMap.end();

	/// loop over the model call
	while( iter!= end)
	{
		CC_NS(std,sort)( ((*iter).second).begin(),((*iter).second).end() );
		strNodePairVec::iterator pos = CC_NS(std,unique)( ((*iter).second).begin(), ((*iter).second).end() );
		((*iter).second).erase( pos, ((*iter).second).end() );
		++iter;
	}
}

///////////////////////////////////////////////////
///	Class  : ARM_ModelCall
///	Routine: GetUsedTimeLags
///	Returns: ARM_VectorPtr
///	Action : function that loops over the various model 
///				names and pricing function
///				to get the various used time lags!
///				the result is not sorted and with 
///				potential duplicates!
///////////////////////////////////////////////////
pair<ARM_VectorPtr,ARM_AdditionalTimeInfoPtrVector> ARM_ModelCall::GetUsedTimeLags( ARM_PricingModel* model )
{
	vector<double> result;
	ARM_AdditionalTimeInfoPtrVector AdditionalTimeInfos;
	
	strStrNodePairVecMap::iterator
		iter = itsModelFunctionNNodeMap.begin(),
		end = itsModelFunctionNNodeMap.end();

	/// loop over the various model names!
	while( iter!= end)
	{
		size_t i;
		/// loop over the various model call
		for(i=0; i<(*iter).second.size(); ++i)
		{
			ARM_NodeInfo nodeInfo = ((ARM_ExpNodeFunc*)	&*((*iter).second)[i].second)->GetUsedTimeLags( model );
			ARM_VectorPtr tmpUsedTimeLags = nodeInfo.first;
			ARM_AdditionalTimeInfoPtr additionalTimeInfo = nodeInfo.second;
            result.insert( result.end(), tmpUsedTimeLags->begin(), tmpUsedTimeLags->end() );

			if( additionalTimeInfo != ARM_AdditionalTimeInfoPtr(NULL) ) 
				AdditionalTimeInfos.push_back( additionalTimeInfo );
		}
		++iter;
	}
	return pair<ARM_VectorPtr,ARM_AdditionalTimeInfoPtrVector>( ARM_VectorPtr( new ARM_GP_Vector( result ) ), AdditionalTimeInfos );
}


///////////////////////////////////////////////////
///	Class  : ARM_ModelCall
///	Routine: Reset
///	Returns: void
///	Action : Reset the node after the GetUsedTimeLags
///////////////////////////////////////////////////
void ARM_ModelCall::Reset( ARM_ExpNode::ResetLevel resetParam )
{
	strStrNodePairVecMap::iterator
		iter = itsModelFunctionNNodeMap.begin(),
		end = itsModelFunctionNNodeMap.end();

	while( iter!= end)
	{
		size_t i;
		/// loop over the various model call
		for(i=0; i<(*iter).second.size(); ++i)
		{
			/// do not store the Node Eval result as the caching as been done with dumStates!
			/// compared to Reset... do not reset as well the GramFunction!
			((ARM_ExpNodeFunc*)	&*((*iter).second)[i].second)->Reset( resetParam );
		}
		++iter;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelCall
///	Routine: toString
///	Returns: string
///	Action : stringify the object
////////////////////////////////////////////////////
string ARM_ModelCall::toString( const string& indent ) const
{
	CC_Ostringstream os;
	os << indent <<   "Call Date  : " << itsCallDate.toString() << "\n";
	strStrNodePairVecMap::const_iterator 
		iter = itsModelFunctionNNodeMap.begin(),
		end = itsModelFunctionNNodeMap.end();

	while(iter!=end)
	{
		os << indent << "Model Name : " << (*iter).first << "\n";
		os << indent << "Function(s): \n";
		size_t i;
		for(i=0; i<(*iter).second.size(); ++i)
			os << indent+"\t" << ((*iter).second)[i].first << "\n";
		++iter;
	}
	return os.str();
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
