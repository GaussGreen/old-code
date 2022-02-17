/*!
 *
 *	Copyright (c) CDC IXIS CM July 2003 Paris
 *	\file pricingadviser.h
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */

#ifndef _INGPINFRA_PRICINGADVISER_H
#define _INGPINFRA_PRICINGADVISER_H

/// this headers has to come first
/// as it contains pragma for stupid warning
#include "gpbase/removeidentifiedwarning.h"

#include <glob/dates.h>

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/gpvector.h"

#include "typedef.h"
#include "nummethod.h"
#include "modelnrefcall.h"

/// unix does not support sstream
#include "gpbase/ostringstream.h"

/// STL
#include <vector>
CC_USING_NS( std, vector )

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration 
class ARM_GenSecurity;


///	Some helper classes
/// There isn't much in it but it gives less
/// work to the pricing Adviser which can 
/// focuss on more fundamental things

///////////////////////////////////////////////////////
/// \class ARM_PricingAdviser
/// \brief
/// The pricing adviser is the middle
///		man between the generic security and the pricing
///		model. It gets pricing information
///		when parsing the generic security. It advises
///		the generic security how to be priced, hence
///		the name!
///////////////////////////////////////////////////////

class ARM_PricingAdviser
{
public:
	///////////////////////////////////////////////////////
	/// \class CellRefCalls
	/// \brief
	/// Very simple nested class to manage ref calls with
	/// meaning of the DateType enum is the following:
	///		-DATEIN_INC:	Use ChildDate in Increasing order to sort
	///		-DATEIN_DEC:	Use ChildDate in Decreasing order to sort
	///		-DATEOUT_INC:	Use ParentDate in Increasing order to sort
	///		-DATEOUT_DEC:	Use ParentDate in Decreasing order to sort
	///////////////////////////////////////////////////////

	//typedef ARM_GP_T_Matrix<ARM_ExpNodePtr> ARM_GP_NodeMatrix;
	//typedef	ARM_CountedPtr< ARM_GP_NodeMatrix >	ARM_GP_NodeMatrixPtr;

	struct CellRefCalls
	{
		enum DateType{ DATEIN_INC, DATEIN_DEC, DATEOUT_INC, DATEOUT_DEC	};
		CellRefCalls(){};
		CellRefCalls( const vector< ARM_CellRefCallPtr >& refCalls, DateType Type );
		
		vector<ARM_CellRefCallPtr> GetCellRefCallAndIncIdx( const ARM_Date& date );
		size_t size() const { return itsCoords.size(); }
		
		/// data members
		vector<ARM_CellRefCallPtr> itsCoords;
		size_t itsIdx;
		DateType itsDateType;
	};

	ARM_PricingAdviser( const ARM_PricingAdviser& rhs);
	ARM_PricingAdviser& operator = ( const ARM_PricingAdviser& rhs );

public:

	ARM_PricingAdviser();
	~ARM_PricingAdviser();

	/// typedef for shorter typing
	typedef ARM_NumMethod::GP_PricingDirection PricingDir;

	PricingDir GetPricingDir( ) const;
	
	/// function for the genSecurity
	void AddCellRefCall( const ARM_CellRefCallPtr& cellRefCall );
	void AddModelCall( const ARM_ModelCallPtr& modelCall );

	/// pricing information 
	void DecidePricing(ARM_PricingModel* model);

	void InitPricing(ARM_PricingModel* model);

	/// Loop processing
	typedef std::vector< ARM_ExpNodePtr >::iterator ExpNodeIter;
	inline size_t GetFirstNodeIdx() const { return itsFirstNodeIdx; }
	inline size_t GetLastNodeIdx() const { return itsLastNodeIdx; }
	inline void Increment( size_t& nodeIdx, ARM_GP_Vector::iterator& fromTime, ARM_GP_Vector::iterator& toTime ) const { nodeIdx+=itsIncrement; fromTime+=itsIncrement; toTime+=itsIncrement; }

	/// GetTimeInfos and toTime vector
	inline ARM_VectorPtr GetFromTimes() const { return itsFromTimes; }
	ARM_GP_Vector::iterator GetFirstFromTime() const { return itsFirstFromTime; }
	ARM_GP_Vector::iterator GetLastFromTime() const { return itsLastFromTime; }
	inline ARM_VectorPtr GetToTimes() const { return itsToTimes; }
	ARM_GP_Vector::iterator GetFirstToTime() const { return itsFirstToTime; }
	ARM_GP_Vector::iterator GetLastToTime() const { return itsLastToTime; }
	ARM_TimeInfoPtrVector GetTimeInfos( ARM_PricingModel* model );
	ARM_TimeInfoPtrVector GetModelAndEventTimeInfos( const ARM_VectorPtr& TimeInfoEventTimes, ARM_PricingModel* model );

	ARM_ExpNodePtr GetNode(size_t rowIdx, size_t colIdx) { return ((ARM_ExpNodePtr) ((*itsPParseTree)(rowIdx, colIdx))); };

	/// compatibility with a given numerical method
	bool IsCompatible( ARM_PricingModel*) const;
	string CompatibilityPb( ARM_PricingModel*) const;

	/// for the part on reference calls
	size_t SetComputedNodesCommon(
		CellRefCalls& inCellRefCalls,
		CellRefCalls& outCellRefCalls,
		const ARM_Date& date, ARM_PricingModel* model,
		ARM_PricingStatesPtr& states );
	size_t GetAuxiliaryNodesCommon( 
		CellRefCalls& cellRefCalls,	
		const ARM_Date& date, ARM_PricingModel* model,
		ARM_PricingStatesPtr& states );
	size_t SetComputedNodes( 
		const ARM_Date& date, ARM_PricingModel* model,
		ARM_PricingStatesPtr& states );
	size_t GetAuxiliaryNodes( 
		const ARM_Date& date, ARM_PricingModel* model,
		ARM_PricingStatesPtr& states );

	inline CC_STL_VECTOR( ARM_ExpNodePtr )* getPVNodes() { return itsPVNodes; }
	inline void AddPVNode( const ARM_ExpNodePtr& PVNode ) { itsPVNodes->push_back( PVNode ); }
	inline void ResetPVNodesList( ) { delete itsPVNodes; itsPVNodes = new CC_STL_VECTOR( ARM_ExpNodePtr )(0);}

	/// for easy debugging
    string toString(const string& indent="", const string& nextIndent="") const;
	string detailledString(const string& indent, const string& nextIndent ) const;

	/// because we use pre-computed technology, we need to remove precomputed node from
	/// previous computation!
	// By default we do a full reset
	void ResetPrecomputedNodes(ARM_ExpNode::ResetLevel partialReset );

private:	
	vector< ARM_CellRefCallPtr >	itsFwdCoords;
	vector< ARM_CellRefCallPtr >	itsBckwdCoords;
	vector< ARM_ModelCallPtr >		itsModelCalls;

	/// pointor to the parse tree!
	ARM_GP_NodeMatrixPtr			itsPParseTree;
	ARM_GP_VectorPtr				itsRowDates;

	CC_STL_VECTOR( ARM_ExpNodePtr )* itsPVNodes;
	
	/// this data are set after the DecidePricing method step
	PricingDir itsPricingDir;
	ARM_VectorPtr itsFromTimes;
	ARM_GP_Vector::iterator itsFirstFromTime;
	ARM_GP_Vector::iterator itsLastFromTime;
	ARM_VectorPtr itsToTimes;
	ARM_GP_Vector::iterator itsFirstToTime;
	ARM_GP_Vector::iterator itsLastToTime;
	size_t itsFirstNodeIdx;
	size_t itsLastNodeIdx;
	int itsIncrement;

	/// for cell call sorted once and for all
	CellRefCalls itsFWDInRefCalls;
	CellRefCalls itsFWDOutRefCalls;
	CellRefCalls itsBCKWDInRefCalls;
	CellRefCalls itsBCKWDOutRefCalls;

	/// some static variables
	static const string itsSeperator;
	static const string NoValidationPb;

	/// to allow ARM_GenSecurity to access these
	/// make it friend
	friend class ARM_GenSecurity;
	/// set methods to populate the pricing adviser
	/// function to set a pointor to the parse tree  shared... hence the vanilla pointor
	void SetPParseTree( ARM_GP_NodeMatrixPtr pParseTree ) { itsPParseTree=pParseTree;}
	/// similar for the dates!
	void SetPDates( ARM_GP_VectorPtr RowDates ) { itsRowDates= RowDates;}

    string toStringCommon(const string& indent="", const string& nextIndent="",bool detailledMode=false) const;

	ARM_PricingAdviser* Clone() const;
};




CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

