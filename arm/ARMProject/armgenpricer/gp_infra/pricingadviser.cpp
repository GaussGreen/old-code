/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file pricingadviser.cpp
 *
 *  \brief object to handle the pricing information 
 *		in the generic security
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */

/// gpbase
/// this headers has to come first
/// as it contains pragma for stupid warning
#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/env.h"
#include "gpinfra/typedef.h"

/// gpinfra
#include "gpinfra/pricingadviser.h"

#include "gpinfra/gensecurity.h"
#include "gpinfra/timeinfo.h"
#include "gpinfra/gramnode.h"
#include "gpinfra/pricingstatesprocess.h"
#include "gpinfra/modelnrefcall.h"
#include "gpinfra/dealdescription.h"
#include "gpinfra/typedef.h"

/// STL
#include <algorithm>
#include <functional>
#include <ostream>

CC_USING_NS(std,ostream)


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_PricingAdviser
///	Const  : NoValidationPb
/// Meaning: the validation problem in this case, no problem
///		as it says
////////////////////////////////////////////////////

const string ARM_PricingAdviser::NoValidationPb = "No Problem";


////////////////////////////////////////////////////
///	Class  : ARM_PricingAdviser
///	Const  : itsSeperator
/// Meaning: separator when stringifying
////////////////////////////////////////////////////

const string ARM_PricingAdviser::itsSeperator = "\t";

/// Dates sorting routines!
///		- Inc for increasing order
///		- Dec for decreasing order
///		- In for date to put in
///		- Out for dates to get out

bool DateIn_Inc_Sort( const ARM_CellRefCallPtr& lhs, const ARM_CellRefCallPtr& rhs )
{
	return lhs->GetChildDate() < rhs->GetChildDate(); 
}

bool DateIn_Dec_Sort( const ARM_CellRefCallPtr& lhs, const ARM_CellRefCallPtr& rhs )
{
	return lhs->GetChildDate()  > rhs->GetChildDate() ;
}

bool DateOut_Inc_Sort( const ARM_CellRefCallPtr& lhs, const ARM_CellRefCallPtr& rhs )
{
	return lhs->GetParentDate() < rhs->GetParentDate();
}

bool DateOut_Dec_Sort( const ARM_CellRefCallPtr& lhs, const ARM_CellRefCallPtr& rhs )
{
	return lhs->GetParentDate() > rhs->GetParentDate();
}


////////////////////////////////////////////////////
///	Class  : ARM_PricingAdviser
///	Routine: Constructor 
///	Returns: built object
///	Action : build the object
////////////////////////////////////////////////////
ARM_PricingAdviser::ARM_PricingAdviser()
:	
itsFwdCoords(), 
itsBckwdCoords(), 
itsModelCalls(),
itsPParseTree(NULL), 
itsRowDates(NULL), 
itsPVNodes(NULL),  
itsPricingDir(),
itsFromTimes(),
itsFirstFromTime(),
itsFirstToTime(), 
itsFirstNodeIdx(), 
itsLastNodeIdx(), 
itsIncrement()
{ 
	itsPVNodes = new CC_STL_VECTOR( ARM_ExpNodePtr )(0); 
}
	


////////////////////////////////////////////////////
///	Class  : ARM_PricingAdviser
///	Routine: Desctructor
///	Returns: built object
///	Action : desctruct the object
////////////////////////////////////////////////////
ARM_PricingAdviser::~ARM_PricingAdviser()
{
	delete itsPVNodes;
	itsPVNodes = NULL;
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingAdviser
///	Routine: Copy Constructor 
///	Returns: built object
///	Action : constructor
////////////////////////////////////////////////////
ARM_PricingAdviser::ARM_PricingAdviser( const ARM_PricingAdviser& rhs)
:
	itsFwdCoords(rhs.itsFwdCoords),
	itsBckwdCoords(rhs.itsBckwdCoords),
	itsModelCalls(rhs.itsModelCalls),
	itsPricingDir(rhs.itsPricingDir),
	itsFromTimes(rhs.itsFromTimes),
	itsFirstFromTime(rhs.itsFirstFromTime),
	itsLastFromTime(rhs.itsLastFromTime),
	itsToTimes(rhs.itsToTimes),
	itsFirstToTime(rhs.itsFirstToTime),
	itsLastToTime(rhs.itsLastToTime),
	itsFirstNodeIdx(rhs.itsFirstNodeIdx),
    itsLastNodeIdx(rhs.itsLastNodeIdx),
	itsIncrement(rhs.itsIncrement),
	itsFWDInRefCalls(rhs.itsFWDInRefCalls),
	itsFWDOutRefCalls(rhs.itsFWDOutRefCalls),
	itsBCKWDInRefCalls(rhs.itsBCKWDInRefCalls),
	itsBCKWDOutRefCalls(rhs.itsBCKWDOutRefCalls)
{
	itsPParseTree	= ARM_GP_NodeMatrixPtr ((ARM_GP_NodeMatrix*) (rhs.itsPParseTree->Clone()));
	itsRowDates		= ARM_GP_VectorPtr((ARM_GP_Vector*) rhs.itsRowDates->Clone());
	itsPVNodes      = new CC_STL_VECTOR( ARM_ExpNodePtr ) (*rhs.itsPVNodes);
}


////////////////////////////////////////////////////
///	Class  : ARM_PricingAdviser
///	Routine: Operator = 
///	Returns: built object
///	Action : constructor
////////////////////////////////////////////////////
ARM_PricingAdviser& ARM_PricingAdviser::operator = ( const ARM_PricingAdviser& rhs )
{
	if( this != &rhs )
	{
		itsFwdCoords		= (rhs.itsFwdCoords);
		itsBckwdCoords		= (rhs.itsBckwdCoords);
		itsModelCalls		= (rhs.itsModelCalls);
		itsPricingDir		= (rhs.itsPricingDir);
		itsFromTimes		= (rhs.itsFromTimes);
		itsFirstFromTime	= (rhs.itsFirstFromTime);
		itsLastFromTime		= (rhs.itsLastFromTime);
		itsToTimes			= (rhs.itsToTimes);
		itsFirstToTime		= (rhs.itsFirstToTime);
		itsLastToTime		= (rhs.itsLastToTime);
		itsFirstNodeIdx		= (rhs.itsFirstNodeIdx);
		itsLastNodeIdx		= (rhs.itsLastNodeIdx);
		itsIncrement		= (rhs.itsIncrement);
		itsFWDInRefCalls	= (rhs.itsFWDInRefCalls);
		itsFWDOutRefCalls	= (rhs.itsFWDOutRefCalls);
		itsBCKWDInRefCalls	= (rhs.itsBCKWDInRefCalls);
		itsBCKWDOutRefCalls = (rhs.itsBCKWDOutRefCalls);

		//Pointors
		itsPParseTree	= ARM_GP_NodeMatrixPtr ((ARM_GP_NodeMatrix*) (rhs.itsPParseTree->Clone()));
		itsRowDates		= ARM_GP_VectorPtr ((ARM_GP_Vector*) rhs.itsRowDates->Clone());
		itsPVNodes      = new CC_STL_VECTOR( ARM_ExpNodePtr ) (*rhs.itsPVNodes);
	}
	return *this;
}



////////////////////////////////////////////////////
///	Class  : ARM_PricingAdviser::CellRefCalls
///	Routine: Constructor of cellRefCalls
///	Returns: built object
///	Action : constructor
////////////////////////////////////////////////////

ARM_PricingAdviser::CellRefCalls::CellRefCalls(
	const vector< ARM_CellRefCallPtr >& refCalls, DateType Type )
:	itsCoords( refCalls ), itsIdx(0), itsDateType(Type)
{
	/// sorting of the refCall according to the Date type
	switch(itsDateType)
	{
	case DATEIN_INC:
		{
			CC_NS( std, sort )( itsCoords.begin(), itsCoords.end(), DateIn_Inc_Sort );
			break;
		}
	case DATEIN_DEC:
		{
			CC_NS( std, sort )( itsCoords.begin(), itsCoords.end(), DateIn_Dec_Sort );
			break;
		}
	case DATEOUT_INC:
		{
			CC_NS( std, sort )( itsCoords.begin(), itsCoords.end(), DateOut_Inc_Sort );
			break;
		}
	case DATEOUT_DEC:
		{
			CC_NS( std, sort )( itsCoords.begin(), itsCoords.end(), DateOut_Dec_Sort );
			break;
		}
	default:
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Unknow date type!" );
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingAdviser::CellRefCalls
///	Routine: GetExpNodeAndIncIdx
///	Returns: a vector<ARM_ExpNodePtr> 
///	Action : basically all the cells corresponding 
///				to the date with an accuracy of 0.5 day!
///				and increase the index
////////////////////////////////////////////////////
vector<ARM_CellRefCallPtr> ARM_PricingAdviser::CellRefCalls::GetCellRefCallAndIncIdx( const ARM_Date& date )
{
	vector<ARM_CellRefCallPtr> result;
	const double ACCURACY = 0.5; /// half a day difference!

	switch(itsDateType)
	{
	case DATEIN_DEC:
		{
			while( itsIdx < itsCoords.size() && 
				( date.GetJulian()- itsCoords[itsIdx]->GetChildDate().GetJulian() ) < ACCURACY )
			{
				result.push_back( itsCoords[itsIdx] );
				++itsIdx;
			}
			break;
		}
	case DATEIN_INC:
		{
			while( itsIdx < itsCoords.size() && 
				( itsCoords[itsIdx]->GetChildDate().GetJulian() - date.GetJulian() ) < ACCURACY )
			{
				result.push_back( itsCoords[itsIdx] );
				++itsIdx;
			}
			break;
		}
		
	case DATEOUT_DEC:
		{
			while( itsIdx < itsCoords.size() && 
				( date.GetJulian()- itsCoords[itsIdx]->GetParentDate().GetJulian() ) < ACCURACY )
			{
				result.push_back( itsCoords[itsIdx] );
				++itsIdx;
			}
			break;
		}
	case DATEOUT_INC:	
		{
			while( itsIdx < itsCoords.size() && 
				( itsCoords[itsIdx]->GetParentDate().GetJulian() - date.GetJulian() ) < ACCURACY )
			{
				result.push_back( itsCoords[itsIdx] );
				++itsIdx;
			}
			break;
		}
	default:
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Unknow date type!" );
	}

	return result;
}


////////////////////////////////////////////////////
///	Class  : ARM_PricingAdviser
///	Routine: AddCellRefCall
///	Returns: void
///	Action : add Cell reference call to 
///		-itsBckwdCoords: if the call 
///	 is from a cell with a higher row (= higher date) 
///		-itsFwdCoords: if the call 
///	is from a cell with a lower row (= lower date) 
////////////////////////////////////////////////////

void ARM_PricingAdviser::AddCellRefCall( const ARM_CellRefCallPtr& cellRefCall )
{
    /// in the future so backwdlooking
	if( cellRefCall->GetParentDate() < cellRefCall->GetChildDate() )
		itsBckwdCoords.push_back( cellRefCall );

	/// case ParentDate > ChildDate
	/// in the past so fwdlooking
	else if( cellRefCall->GetParentDate() > cellRefCall->GetChildDate() )
			itsFwdCoords.push_back( cellRefCall );

	/// for nodes that are on the same rows in the deal description
	/// the pricing adviser do not need to store them!
}


////////////////////////////////////////////////////
///	Class  : ARM_PricingAdviser
///	Routine: AddModelCall
///	Returns: void
///	Action : adds model call 
////////////////////////////////////////////////////

void ARM_PricingAdviser::AddModelCall( const ARM_ModelCallPtr& modelCall )
{
	itsModelCalls.push_back( modelCall );
}


///////////////////////////////////////////////////
///	Class  : ARM_PricingAdviser
///	Routine: toString,detailledString,toStringCommon
///	Returns: stringify the object
///	Action : 
////////////////////////////////////////////////////
string ARM_PricingAdviser::toString(const string& indent, const string& nextIndent ) const
{ 
#ifdef _DEBUG
	return toStringCommon(indent,nextIndent,true); 
#else
	return toStringCommon(indent,nextIndent,false); 
#endif
}


string ARM_PricingAdviser::detailledString(const string& indent, const string& nextIndent ) const
{ return toStringCommon(indent,nextIndent,true); }


string ARM_PricingAdviser::toStringCommon(const string& indent, const string& nextIndent, bool detailledMode) const
{
	CC_Ostringstream os;
	os << "ARM_PricingAdviser Object\n";
	int PricingDirEnum = GetPricingDir();
	os << "Pricing Direction	: " 
		<< ARM_NumMethod::GP_PricingDirectionTxt[ PricingDirEnum  ] << "\n";
	
	if( detailledMode )
	{
		if( !itsFwdCoords.empty() )
		{
			os << "\nForward Coordinates	: " << itsFwdCoords.size() << "\n";
			for( int i=0; i<itsFwdCoords.size(); ++i )
			{
				os << i << "\n";
				os << itsFwdCoords[i]->toString( itsSeperator );
			}
			
		}
		
		if( !itsBckwdCoords.empty() )
		{
			os << "\nBackward Coordinates	: " << itsBckwdCoords.size() << "\n";
			for( int i=0; i<itsBckwdCoords.size(); ++i )
			{
				os << i << "\n";
				os << itsBckwdCoords[i]->toString( itsSeperator );
			}
		}
		

		if( !itsModelCalls.empty() )
		{
			os << "\nModel Calls		: " << itsModelCalls.size() << "\n";
			for( int i=0; i<itsModelCalls.size(); ++i )
			{
				os << i << "\n";
				os << itsModelCalls[i]->toString( itsSeperator );
			}
		}
	}
	
	return os.str();
}


////////////////////////////////////////////////////
///	Class  : ARM_PricingAdviser
///	Routine: ARM_PricingAdviser
///	Returns: ARM_NumMethod::GP_PricingDirection
///	Action : says whether the deal is 
///		- forward looking
///		- backward looking
///		- forward backward looking
///		- ambiguous
////////////////////////////////////////////////////

ARM_NumMethod::GP_PricingDirection ARM_PricingAdviser::GetPricingDir( ) const
{
	if( !itsFwdCoords.empty() && !itsBckwdCoords.empty() )
		return ARM_NumMethod::GP_FWDBCKWDLOOKING;
	if( !itsFwdCoords.empty() )
		return ARM_NumMethod::GP_FWDLOOKING;
	if( !itsBckwdCoords.empty() )
		return ARM_NumMethod::GP_BCKWDLOOKING;
	//// other cases are ambiguous
	return ARM_NumMethod::GP_AMBIGUOUS;
}


////////////////////////////////////////////////////
///	simple predictate to sort the model calls
////////////////////////////////////////////////////
struct ModelCallPtrPred : CC_NS( std, binary_function )< ARM_ModelCallPtr, ARM_ModelCallPtr, bool>
{  
	bool operator()( const ARM_ModelCallPtr& lhs, const ARM_ModelCallPtr& rhs) const
    {
        return lhs->GetCallDate() < rhs->GetCallDate();
    }
};


////////////////////////////////////////////////////
///	Class  : ARM_PricingAdviser
///	Routine: GetModelAndEventTimeInfos
///	Returns: ARM_TimeInfoPtrVector
///	Action : Compute the timeInfo from the model call information
///		using the TimeInfoEventTimes as the basic event times!
///		this is not const because it sorts model calls! (in sort)
////////////////////////////////////////////////////
ARM_TimeInfoPtrVector ARM_PricingAdviser::GetModelAndEventTimeInfos( const ARM_VectorPtr& TimeInfoEventTimes,
	ARM_PricingModel* model )
{
	size_t TimeInfoSize = TimeInfoEventTimes->size();
	ARM_TimeInfoPtrVector result, tmpResult;

   	const double ACCURACY = 0.5; /// half a day difference!
	tmpResult.reserve(TimeInfoSize);
    
    /// FIRST get all the model calls
	/// sort the model calls for easy lookup of the corresponding dates!
    CC_NS( std, sort )( itsModelCalls.begin(), itsModelCalls.end(), ModelCallPtrPred() );

    if( !itsModelCalls.empty() )
    {
        double eventTime;
        size_t pos;
        size_t endPos = itsModelCalls.size();

		for (pos = 0; pos < endPos; ++pos)
		{
			itsModelCalls[pos]->Reset(ARM_ExpNode::ResetCountor);
		}

        for( pos = 0; pos < endPos; ++pos )
		{
			eventTime = model->GetTimeFromDate( itsModelCalls[pos]->GetCallDate() );
            vector<double> payTimes;
			ARM_AdditionalTimeInfoPtrVector AdditionalTimeInfos;

			/// find all the model call with the same date
            while( pos < endPos &&
                fabs( model->GetTimeFromDate( itsModelCalls[pos]->GetCallDate() ) 
					- eventTime ) < ACCURACY )
			{ 
				pair<ARM_VectorPtr, ARM_AdditionalTimeInfoPtrVector> PayTimesAndAdditionalTimeInfos = itsModelCalls[pos]->GetUsedTimeLags(model);
                ARM_VectorPtr currentPayTimes = PayTimesAndAdditionalTimeInfos.first;
				ARM_AdditionalTimeInfoPtrVector currentAdditionalTimeInfos = PayTimesAndAdditionalTimeInfos.second;
				
                payTimes.insert( payTimes.end(), currentPayTimes->begin(), currentPayTimes->end() );
				AdditionalTimeInfos.insert( AdditionalTimeInfos.end(), currentAdditionalTimeInfos.begin(), currentAdditionalTimeInfos.end() );
                ++pos;
            }

			/// pos is one step after the last model call with the same eventtime
			/// hence the need to remove it!
			if(pos>=1)
				--pos;

			/// sort the payTimes and remove duplicates
			CC_NS( std, sort )( payTimes.begin(), payTimes.end() );
			CC_STL_VECTOR( double )::iterator lastUniquePos = CC_NS( std, unique )( payTimes.begin(), payTimes.end() );
			payTimes.resize( CC_NS( std, CC_DISTANCE )( payTimes.begin(), lastUniquePos ) );

            /// create the corresponding ARM_InfoPtr
            tmpResult.push_back( ARM_TimeInfoPtr( new ARM_TimeInfo( eventTime, ARM_VectorPtr( new std::vector<double>( payTimes ) ),AdditionalTimeInfos ) ) );	
        }

		for (pos = 0; pos < endPos; ++pos)
		{
			itsModelCalls[pos]->Reset(ARM_ExpNode::ResetBucket);
		}
    }

	/// SECOND add all the induct times. We need to merge
	/// time info from tmpResult and from FromTimes
	/// 1) clones the from time to sort it without any problem!
	ARM_VectorPtr tmpTimes( new std::vector<double>( *TimeInfoEventTimes) );
	CC_NS( std, sort )( tmpTimes->begin(), tmpTimes->end() );
	std::vector<double>::iterator times		= tmpTimes->begin();
	std::vector<double>::iterator endTimes	= tmpTimes->end();

	/// 2) Get the the model Time info
	ARM_TimeInfoPtrVector::iterator modTimeInfos	= tmpResult.begin();
	ARM_TimeInfoPtrVector::iterator endmodTimeInfos	= tmpResult.end();
	result.reserve(TimeInfoSize);

	/// THIRD times contains all modelTimes so do a loop over times
	/// for each times test whether we are in a modelTimes event time
	/// if this is the case, result is filled with the corresponding modelTimes
	/// otherwise creates a time info with only an event time
	for( ; times != endTimes; ++times )
	{
		if( modTimeInfos != endmodTimeInfos 
			&& (*modTimeInfos)->GetEventTime() == (*times) )
		{
			result.push_back( (*modTimeInfos) );
			++modTimeInfos;
		}
		else
			result.push_back( ARM_TimeInfoPtr( new ARM_TimeInfo( *times, ARM_VectorPtr( new std::vector<double>(1,*times) ) ) ) );
	}

	return result;
}


////////////////////////////////////////////////////
///	Class  : ARM_PricingAdviser
///	Routine: GetTimeInfos
///	Returns: ARM_TimeInfoPtrVector
///	Action : Compute the timeInfo from the information
///		the pricing adviser has!
///		this is not const because it sorts model calls! (in GetModelAndEventTimeInfos)
////////////////////////////////////////////////////
ARM_TimeInfoPtrVector ARM_PricingAdviser::GetTimeInfos( ARM_PricingModel* model )
{
	ARM_TimeInfoPtrVector result;
    switch( itsPricingDir )
	{
	case ARM_NumMethod::GP_BCKWDLOOKING: case ARM_NumMethod::GP_FWDLOOKING: case ARM_NumMethod::GP_FWDBCKWDLOOKING:
		{
			result = GetModelAndEventTimeInfos( itsFromTimes, model );
			break;
		}

	//// Unknown is for closed forms
	case ARM_NumMethod::GP_UNKNOWN:
		{
			if(itsFromTimes->size()==1 && (*itsFromTimes)[0] < K_NEW_DOUBLE_TOL)
            {
                result.reserve(1);
			    result.push_back( ARM_TimeInfoPtr( new ARM_TimeInfo(0.0, ARM_VectorPtr( NULL ) ) ) );
            }
            else
            {
			    throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                "Not set or unknown numerical method only allowed to price at as of date");
            }
			break;
        }
	default:
		{
			CC_Ostringstream os;
			os << "Invalid direction " << ARM_USERNAME << ": please advise! ";
			throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
		}

	}

	return result;
}


////////////////////////////////////////////////////
///	Class  : ARM_PricingAdviser
///	Routine: CompatibilityPb
///	Returns: string
///	Action : return a string for compatibility problem
////////////////////////////////////////////////////
string ARM_PricingAdviser::CompatibilityPb( ARM_PricingModel* model ) const
{
	/// checks that we have a numerical method
    if( model->GetNumMethod() == ARM_NumMethodPtr(NULL) )
    {
		CC_Ostringstream os;
		os << "no numerical method set with the model "
			<< ARM_USERNAME << " please advise!";
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
    }

	string result( ARM_PricingAdviser::NoValidationPb );
	ARM_NumMethod::GP_PricingDirection mPricingDir =
		model->GetNumMethod()->GetPricingDirection();
	
	/// tests the direction of the numerical model
	switch(mPricingDir)
	{
	/// fine, it can do everything!
	case ARM_NumMethod::GP_FWDBCKWDLOOKING:
		break;
		
	/// compatible if only same method or GP_AMBIGUOUS
	case ARM_NumMethod::GP_FWDLOOKING: case ARM_NumMethod::GP_BCKWDLOOKING:
		{
			if(	mPricingDir!= GetPricingDir()
				&&  ARM_NumMethod::GP_AMBIGUOUS != GetPricingDir() )
			{
				CC_Ostringstream os;
				os << "Incompatible pricing direction: "
					<< "Security of type " << ARM_NumMethod::GP_PricingDirectionTxt[ GetPricingDir() ]
					<< " while numerical method of type: " << ARM_NumMethod::GP_PricingDirectionTxt[ mPricingDir ]
					<< "  " << ARM_USERNAME << ": please advise!  ";
				result = os.str();
			}
			break;
		}
	
	/// this is nonsense so throw exception!
	case ARM_NumMethod::GP_AMBIGUOUS: default:
		{
			CC_Ostringstream os;
			os << "A model with its associated numerical method has to know in which direction it can price\n"
				<< "Please advise!\n";
			result = os.str();
		}
	}

	/// Should we also validate model calls?
	return result;
}


////////////////////////////////////////////////////
///	Class  : ARM_PricingAdviser
///	Routine: IsCompatible
///	Returns: bool
///	Action : the bool tells whether 
////////////////////////////////////////////////////

bool ARM_PricingAdviser::IsCompatible(ARM_PricingModel* model) const
{
	return CompatibilityPb(model) == ARM_PricingAdviser::NoValidationPb;
}


////////////////////////////////////////////////////
///	Class  : ARM_PricingAdviser
///	Routine: SetComputedNodesCommon
///	Returns: size_t
///	Action : This routine gets the precomputed nodes
///			sets it to the corresponding parent nodes
///			and then remove the precomputed payoff from 
///			the full payoff matrix.. it also updates
///			the payoff position of other nodes
////////////////////////////////////////////////////

size_t ARM_PricingAdviser::SetComputedNodesCommon( 
	ARM_PricingAdviser::CellRefCalls& inCellRefCalls,
	ARM_PricingAdviser::CellRefCalls& outCellRefCalls,
	const ARM_Date& date, ARM_PricingModel* model,
	ARM_PricingStatesPtr& states )
{
	vector<ARM_CellRefCallPtr> cellCalls = outCellRefCalls.GetCellRefCallAndIncIdx(date);
	size_t i;
	size_t nbPayoffsToRemove= cellCalls.size();

    if( nbPayoffsToRemove )
    {
	    /// 1) set the parent node with the computed payoff
	    for( i=0; i<nbPayoffsToRemove; ++i )
	    {
		    ARM_VectorPtr computedPayoff = states->GetPayoffVec( cellCalls[i]->GetPayoffPosition() );
		    ARM_GramFctorArg val( computedPayoff );
		    cellCalls[i]->GetCastedRefNode()->SetValue(val, model, states );
	    }

	    /// 2) shift the payoff from one and build the indexToRemoveVec
	    vector<size_t> indexToRemoveVec;
	    indexToRemoveVec.reserve(nbPayoffsToRemove);
	    size_t k, allInRefIndexSize = inCellRefCalls.size();
	    const int UNASSIGNED = -1;

	    /// get the list of all the indexes to remove
	    for (size_t i=0; i<nbPayoffsToRemove; ++i )
	    {
		    int oldpPosition = cellCalls[i]->GetPayoffPosition();
    #ifdef __GP_STRICT_VALIDATION
		    if(oldpPosition == UNASSIGNED)
			    throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				    "trying to remove unassigned payoff!" );
    #endif
		    indexToRemoveVec.push_back(oldpPosition);
	    }

	    /// sort them to avoid problem of unsorted indexes!
	    CC_NS( std,sort )( indexToRemoveVec.begin(), indexToRemoveVec.end() );

	    for (size_t i=0; i<nbPayoffsToRemove; ++i )
	    {
		    /// needs to do minus i because it is shifted by minus 1 each time!
		    int oldpPosition = indexToRemoveVec[i]-i;
		    
		    /// remove payoff 
		    /// 1st: for any cellRef Decrease by 1 if higher than oldpPosition!
		    double currentPayoffPosition;
		    for( k=0; k<allInRefIndexSize; ++k )
		    {
			    currentPayoffPosition = inCellRefCalls.itsCoords[k]->GetPayoffPosition();
			    if( currentPayoffPosition > oldpPosition )
				    inCellRefCalls.itsCoords[k]->SetPayoffPosition(currentPayoffPosition-1);
			    else
				    if( inCellRefCalls.itsCoords[k]->GetPayoffPosition() == oldpPosition )
					    inCellRefCalls.itsCoords[k]->SetPayoffPosition(UNASSIGNED);
		    }
	    }

	    /// 3) remove the index from the Payoff matrix!
	    states->RemovePayoffVec( indexToRemoveVec );
    }
	return nbPayoffsToRemove;
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingAdviser
///	Routine: SetComputedNodes
///	Returns: size_t
///	Action : This method apply the SetComputeNodes method
/// to FWD and BCKWD refcalls
////////////////////////////////////////////////////

size_t ARM_PricingAdviser::SetComputedNodes( 
	const ARM_Date& date, ARM_PricingModel* model,
	ARM_PricingStatesPtr& states )
{
	size_t nbPayoffsToRemove = SetComputedNodesCommon(
		itsFWDInRefCalls,
		itsFWDOutRefCalls,
		date, model,
		states);
	nbPayoffsToRemove += SetComputedNodesCommon(
		itsBCKWDInRefCalls,
		itsBCKWDOutRefCalls,
		date, model,
		states);

	return nbPayoffsToRemove;
}


////////////////////////////////////////////////////
///	Class  : ARM_PricingAdviser
///	Routine: GetAuxiliaryNodesCommon
///	Returns: void
///	Action : This method evaluates the auxiliary nodes 
/// based on the in refcalls.
////////////////////////////////////////////////////

size_t ARM_PricingAdviser::GetAuxiliaryNodesCommon( 
	ARM_PricingAdviser::CellRefCalls& cellRefCall,
	const ARM_Date& date, ARM_PricingModel* model,
	ARM_PricingStatesPtr& states )
{
	vector<ARM_CellRefCallPtr> cellCalls = cellRefCall.GetCellRefCallAndIncIdx( date );
	for( size_t i=0; i<cellCalls.size(); ++i )
	{
		cellCalls[i]->SetPayoffPosition( states->GetPayoffsSize() );
		
		/// then evaluate and set the corresponding child node
		EvalChildNodeAndProcess( cellCalls[i]->GetRefNode(), model, states, -1,
            PushBackAndSetFtor(), model->GetTimeFromDate( cellCalls[i]->GetChildDate() ), true );
	}
	return cellCalls.size();
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingAdviser
///	Routine: GetAuxiliaryNodesCommon
///	Returns: void
///	Action : This method evaluates the auxliary node
/// of both FWD and BKWD ref calls.
////////////////////////////////////////////////////

size_t ARM_PricingAdviser::GetAuxiliaryNodes( 
	const ARM_Date& date, ARM_PricingModel* model,
	ARM_PricingStatesPtr& states )
{
	size_t nbAuxliaryNodes = GetAuxiliaryNodesCommon(
		itsFWDInRefCalls,
		date, model,
		states);

	nbAuxliaryNodes += GetAuxiliaryNodesCommon(
		itsBCKWDInRefCalls,
		date, model,
		states);

	return nbAuxliaryNodes;
}



////////////////////////////////////////////////////
///	Class  : ARM_PricingAdviser
///	Routine: ResetPrecomputedNodes
///	Returns: void
///	Action : browse the parse tree and remove the precomputed nodes!
////////////////////////////////////////////////////
void ARM_PricingAdviser::ResetPrecomputedNodes(ARM_ExpNode::ResetLevel partialReset)
{
	ExpNodeIter 
		pnode	= itsPParseTree->begin(),
		nodeEnd	= itsPParseTree->end();
	
	while( pnode != nodeEnd )
	{
		(*pnode)->Reset(partialReset);
		++pnode;
	}
}



////////////////////////////////////////////////////
///	Class  : ARM_PricingAdviser
///	Routine: DecidePricing
///	Returns: void
///	Action : Check if the pricing will be possible with numerical
/// method and initialize the fromTimes and the toTimes
////////////////////////////////////////////////////

void ARM_PricingAdviser::DecidePricing(ARM_PricingModel* model)
{
	ARM_NumMethodPtr numMethod=model->GetNumMethod();
    PricingDir mPricingDir;

	if (numMethod != ARM_NumMethodPtr(NULL) )
	{
		/// test if the model is compatible with the pricing adviser!
		if( !IsCompatible( model ) )
			throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, CompatibilityPb( model) );

		/// NumMethod modifies the ParseTree. For Example AMC will change MAX(a,PV(b)) into Exercise(a,b). 
		numMethod->ModifyParseTree( this );
	    mPricingDir = numMethod->GetPricingDirection();
	}
    else
        mPricingDir = ARM_NumMethod::GP_UNKNOWN;

	/// get the times in the correct order
	size_t rowDateSize = itsRowDates->size();

	int i;

	itsPricingDir	= mPricingDir;

	switch(mPricingDir)
	{
	case ARM_NumMethod::GP_FWDLOOKING: case ARM_NumMethod::GP_FWDBCKWDLOOKING: 
		{
			int realSize = 1;

			if (mPricingDir == ARM_NumMethod::GP_FWDLOOKING)
			{
				// With this filtering we get the shorter numeraire time !!!!

				/// all of this is to filter remaining zero in the payoff
				/// we use dummed pricing states and model since no pricing is required
				
				int rowNb = itsPParseTree->GetRowsNb();
				int colNb = itsPParseTree->GetColsNb();

				/// at least one!
				ARM_ExpNodeRef* currentNodeRef=NULL;
				ARM_ExpNodeDateOrDouble* currentNode=NULL;
				ARM_PricingStatesPtr dumPricingStates(NULL);
				ARM_PricingModel* dumModel= NULL;

				for(size_t row=0;row<rowNb;++row )
				{
					bool nonNullRow = false;

					
					for (size_t col=0;col<colNb; ++col)
					{
						currentNodeRef = dynamic_cast<ARM_ExpNodeRef*>(&*(*itsPParseTree)(row,col));
						if (currentNodeRef)
						{
							/// a little ugly but the only way to get above the smart pointor stuff
							currentNode = dynamic_cast<ARM_ExpNodeDateOrDouble*>(&*currentNodeRef->GetChildNode());
							if( !currentNode || currentNode->Eval( dumModel, dumPricingStates ).GetDouble() != 0.0 )
								nonNullRow = true;
						}
					}

					if (nonNullRow)
						realSize=row+1;
				}	
			}
			else
			{
				realSize = rowDateSize;
			}

			std::vector<double> tmpToTimes(realSize);
			std::vector<double> tmpFromTimes(realSize);
			
			tmpFromTimes[0]	= model->GetTimeFromDate( (ARM_Date) (*itsRowDates)[0] );
			for(i=1; i<realSize; ++i)
			{
				tmpFromTimes[i]		= model->GetTimeFromDate( (ARM_Date) (*itsRowDates)[i] );
				tmpToTimes[i-1]		= tmpFromTimes[i];
			}
			tmpToTimes[realSize-1]	= model->GetTimeFromDate( (ARM_Date) (*itsRowDates)[realSize-1] );

			itsToTimes		= ARM_VectorPtr( new std::vector<double>( tmpToTimes ) );
			itsFromTimes	= ARM_VectorPtr( new std::vector<double>( tmpFromTimes ) );

			break;
		}
	case ARM_NumMethod::GP_BCKWDLOOKING:
		{
			std::vector<double> tmpToTimes(rowDateSize);
			std::vector<double> tmpFromTimes(rowDateSize);
			
			tmpFromTimes[rowDateSize-1]	= model->GetTimeFromDate( (ARM_Date) (*itsRowDates)[rowDateSize-1] );
			for(i=rowDateSize-2; i>=0; --i)
			{
				tmpFromTimes[i]		= model->GetTimeFromDate( (ARM_Date) (*itsRowDates)[i] );
				tmpToTimes[i+1]		= tmpFromTimes[i];
			}
			tmpToTimes[0]	= 0.0;

			itsToTimes		= ARM_VectorPtr( new std::vector<double>( tmpToTimes ) );
			itsFromTimes	= ARM_VectorPtr( new std::vector<double>( tmpFromTimes ) );

			break;
		}
	case ARM_NumMethod::GP_UNKNOWN:
		{
            if(rowDateSize == 1 && model->GetTimeFromDate( (ARM_Date) (*itsRowDates)[0]) < K_NEW_DOUBLE_TOL)
            {
				std::vector<double> tmpToTimes(rowDateSize);
				std::vector<double> tmpFromTimes(rowDateSize);

                tmpToTimes[0]   = 0.0;
                tmpFromTimes[0] = 0.0;

			    itsToTimes		= ARM_VectorPtr( new std::vector<double>( tmpToTimes ) );
			    itsFromTimes	= ARM_VectorPtr( new std::vector<double>( tmpFromTimes ) );
            }
            else
            {
			    throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                "Can't decide a direction to price because the numerical method is not set or unknown");
            }
            break;
        }
	default:
		{
			CC_Ostringstream os;
			os << "A model with its associated numerical method has to know in which direction to price";
			throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
		}
	}

	ResetPrecomputedNodes( ARM_ExpNode::ResetTotal );

}

////////////////////////////////////////////////////
///	Class  : ARM_PricingAdviser
///	Routine: InitPricing
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

void ARM_PricingAdviser::InitPricing(ARM_PricingModel* model)
{
	ARM_NumMethodPtr numMethod=model->GetNumMethod();
	PricingDir pricingDir;
    PricingDir pricingDirCurrLoop;
    if(numMethod != ARM_NumMethodPtr(NULL))
    {
		pricingDir = numMethod->GetPricingDirection();
	    pricingDirCurrLoop = numMethod->GetPricingDirCurrLoop();
    }
    else
	{
        pricingDir = pricingDirCurrLoop = ARM_NumMethod::GP_UNKNOWN;
	}

	/// get the times in the correct order
	size_t rowDateSize = itsRowDates->size();

	/// tests the direction of the numerical model
	switch(pricingDirCurrLoop)
	{
	case ARM_NumMethod::GP_FWDLOOKING:
		{
			itsPricingDir	= ARM_NumMethod::GP_FWDLOOKING;
			int realSize = itsFromTimes->size();
			itsFirstNodeIdx	= 0;						/// nodes for the computation of the payoff.
			itsLastNodeIdx	= itsPParseTree->GetRowsNb()-(rowDateSize-realSize);	// last node minus the removed nodes
			itsFirstFromTime	= itsFromTimes->begin();
			itsLastFromTime	= itsFromTimes->end();
			itsFirstToTime	= itsToTimes->begin();
			itsLastToTime	= itsToTimes->end();
			itsIncrement	= 1;

			itsFWDInRefCalls	= CellRefCalls( itsFwdCoords, CellRefCalls::DATEIN_INC  );
			itsFWDOutRefCalls	= CellRefCalls( itsFwdCoords, CellRefCalls::DATEOUT_INC );

			vector< ARM_CellRefCallPtr > dumCoords;
			
			if (pricingDir == ARM_NumMethod::GP_FWDBCKWDLOOKING)
			{
				itsBCKWDInRefCalls	= CellRefCalls( itsBckwdCoords, CellRefCalls::DATEIN_INC  );
				itsBCKWDOutRefCalls	= CellRefCalls( dumCoords, CellRefCalls::DATEOUT_INC  );
			}
			else
			{
				itsBCKWDInRefCalls	= CellRefCalls( dumCoords, CellRefCalls::DATEIN_INC  );
				itsBCKWDOutRefCalls	= CellRefCalls( dumCoords, CellRefCalls::DATEOUT_INC  );
			}

			break;
		}
	case ARM_NumMethod::GP_BCKWDLOOKING:
		{ 
			itsPricingDir	= ARM_NumMethod::GP_BCKWDLOOKING;
			itsFirstNodeIdx	= itsPParseTree->GetRowsNb()-1;
			itsLastNodeIdx		= -1;
			itsFirstFromTime	= itsFromTimes->end()-1;
			itsLastFromTime	= itsFromTimes->begin()-1;
			itsFirstToTime	= itsToTimes->end()-1;
			itsLastToTime	= itsToTimes->begin()-1;
			itsIncrement	= -1;

			vector< ARM_CellRefCallPtr > dumCoords;

			itsFWDInRefCalls	= CellRefCalls( dumCoords, CellRefCalls::DATEIN_INC  );
			itsFWDOutRefCalls	= CellRefCalls( dumCoords, CellRefCalls::DATEOUT_INC  );

			itsBCKWDInRefCalls	= CellRefCalls( itsBckwdCoords, CellRefCalls::DATEIN_DEC  );
			itsBCKWDOutRefCalls	= CellRefCalls( itsBckwdCoords, CellRefCalls::DATEOUT_DEC );
			
			break;
		}
	case ARM_NumMethod::GP_FWDBCKWDLOOKING:
		{
			throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                "BCKWDFWDLOOKING invalid loop direction.");
		}
	case ARM_NumMethod::GP_UNKNOWN:
		{
			if(rowDateSize == 1 && model->GetTimeFromDate((ARM_Date) (*itsRowDates)[0]) < K_NEW_DOUBLE_TOL)
            {
				/// The only one case allowed : a single event date at as of
			    itsPricingDir	= pricingDir;
			    itsFirstNodeIdx	= 0;
			    itsLastNodeIdx	= itsPParseTree->GetRowsNb();
				itsFirstFromTime	= itsFromTimes->end()-1;
				itsLastFromTime	= itsFromTimes->begin()-1;
				itsFirstToTime	= itsToTimes->end()-1;
				itsLastToTime	= itsToTimes->begin()-1;
			    itsIncrement	= 1;
            }
            else
            {
			    throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                "Can't decide a direction to price because the numerical method is not set or unknown");
            }
            break;
        }
	default:
		{
			CC_Ostringstream os;
			os << "A model with its associated numerical method has to know in which direction to price";
			throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
		}
	}
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_PricingAdviser 
///	Routine: Clone
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_PricingAdviser* ARM_PricingAdviser::Clone() const
{
	return new ARM_PricingAdviser(*this);
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
