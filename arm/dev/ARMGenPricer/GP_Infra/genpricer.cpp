/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file genpricer.cpp
 *
 *  \brief generic pricer object does all the pricing loop
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */

/// this headers has to come first
/// as it contains pragma for stupid warning
#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/env.h"
#include "gpbase/stringmanip.h"
#include "gpbase/autocleaner.h"

/// this header should be second as
/// it includes some preprocessor constant for
#include "gpinfra/genpricer.h"
#include "gpinfra/gensecurity.h"
#include "gpinfra/pricingmodel.h"
#include "gpinfra/numeraire.h"
#include "gpinfra/pricingadviser.h"
#include "gpinfra/lexerdec.h"
#include "gpinfra/pricingstatesprocess.h"
#include "gpinfra/modelnamemap.h"
#include "gpinfra/lexerdec.h"
#include "gpinfra/timeinfo.h"
#include "gpinfra/modelnrefcall.h"
#include "gpinfra/pricerinfo.h"
#include "gpinfra/dealdescription.h"

CC_BEGIN_NAMESPACE( ARM )

const string ARM_GenPricer::CVTag = "Control Variate";

////////////////////////////////////////////////////
///	Class  : ARM_GenPricer
///	Routine: size
///	Returns: size of the deal description
///	Action : 
////////////////////////////////////////////////////

string ARM_GenPricer::toString(const string& indent, const string& nextIndent) const
{ 
	CC_Ostringstream os;
	
	os << string( "Generic Pricer Object\n\n" );

	os << "======================> Pricing Info\n\n";
	os << itsPricerInfo->toString();
	os << "======================> End Pricing Info\n\n";

	os << " Detail Mode : " << ( itsDetailMode? " On": "Off" );
	os << "\n\n\n";

	if(itsDetailMode)
	{
		if( itsGenSecurity )
			os <<  "Security		:" << itsGenSecurity->toString();

		if( itsPricingMod)
			os <<  "Pricing Model	:" << itsPricingMod->toString();
	}

	return os.str();
}


////////////////////////////////////////////////////
///	Class  : ARM_GenPricer
///	Routine: Clone method
///	Returns: a new copy of this
///	Action : Copy this
////////////////////////////////////////////////////
ARM_Object* ARM_GenPricer::Clone() const
{
	return new ARM_GenPricer( *this );
}


////////////////////////////////////////////////////
///	Class  : ARM_GenPricer
///	Routine: CheckColumnNameExists
///	Returns: checks that the name of a column is included in the generic security
///	Action : throw an exception if not present
////////////////////////////////////////////////////

void ARM_GenPricer::CheckColumnNameExists( const string& columnName, const ARM_StringVector& secColumnNames ) const
{
	bool found=false;
	for( size_t j=0; j<secColumnNames.size(); ++j )
	{
		if( stringGetUpper(secColumnNames[j]) == stringGetUpper(columnName) )
		{ found =true;break; }
	}

	if( !found )
		ARM_THROW( ERR_INVALID_ARGUMENT, " could not find the CVColumnName " + columnName + " in the generic security!" );
}


////////////////////////////////////////////////////
///	Class  : ARM_GenPricer
///	Routine: ARM_GenPricer
///	Returns: the built object
///	Action : Constructor
////////////////////////////////////////////////////
ARM_GenPricer::ARM_GenPricer( ARM_GenSecurity* sec, ARM_PricingModel* mod, 
	const ARM_StringVector& controlVariatesColumns,	const ARM_GP_Vector& controlVariatesPrices,	const string& refPriceColumn, 
	const ARM_GP_Vector& betas )
:	
	ARM_RootObject(),
	itsGenSecurity(sec),
	itsPricingMod(mod),
	itsPricerInfo( NULL ), 
	itsCVInfo(NULL),
	itsDetailMode(false)
{
	if( controlVariatesColumns.size() )
	{
		/// 1) check compatibility with pricerinfo
		ARM_PricerInfo* dummyPricerInfo = itsPricingMod->CreatePricerInfo();
		ARM_AutoCleaner<ARM_PricerInfo> holdPricerInfo(dummyPricerInfo);
		if( !dynamic_cast<ARM_PInfo_MultipleLoop*>(dummyPricerInfo) )
			ARM_THROW( ERR_INVALID_ARGUMENT, " the model with its numerical method is not compatible with the generic control variate!" );

		/// 2) check size
		if( betas.size() && betas.size() != controlVariatesColumns.size() )
			ARM_THROW( ERR_INVALID_ARGUMENT, " betas.size() && betas.size() != controlVariatesColumns.size() !" );

		if( controlVariatesPrices.size() != controlVariatesColumns.size() )
			ARM_THROW( ERR_INVALID_ARGUMENT, " controlVariatesPrices.size() != controlVariatesColumns.size()  !" );

		/// 3) check columns names
		ARM_StringVector secColumnNames;
		if( sec )
			secColumnNames = sec->GetDealDescription().GetPricedColumnNames();
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, " sec == NULL!" );

		for( size_t i=0; i<controlVariatesColumns.size(); ++i )
			CheckColumnNameExists( controlVariatesColumns[i], secColumnNames );

		/// if we do not give the ref column takes the last column
		string refPriceColumnTmp = refPriceColumn;
		if( refPriceColumnTmp == "" )
			refPriceColumnTmp = sec->GetDealDescription().GetElem( 0, sec->GetDealDescription().GetColsNb()-1);
		else 
			CheckColumnNameExists( refPriceColumnTmp, secColumnNames );

		itsCVInfo = new ARM_CVInfo(controlVariatesColumns, refPriceColumnTmp, controlVariatesPrices, betas );
	}

	CC_ARM_SETNAME( ARM_GENPRICER );
}


////////////////////////////////////////////////////
///	Class  : ARM_GenPricer
///	Routine: Copy Constructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_GenPricer::ARM_GenPricer( const ARM_GenPricer& rhs)
:
    ARM_RootObject( rhs ), 
	itsGenSecurity( rhs.itsGenSecurity ),
	itsPricingMod( rhs.itsPricingMod ),
    itsPricerInfo( rhs.itsPricerInfo? (ARM_PricerInfo*) rhs.itsPricerInfo->Clone() : NULL ),
	itsCVInfo( rhs.itsCVInfo? (ARM_CVInfo*) rhs.itsCVInfo->Clone() : NULL ),
    itsDetailMode(rhs.itsDetailMode)
{}


////////////////////////////////////////////////////
///	Class  : ARM_GenPricer
///	Routine: Destructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_GenPricer::~ARM_GenPricer()
{
	delete itsPricerInfo;
	delete itsCVInfo;
}

////////////////////////////////////////////////////
///	Class  : ARM_GenPricer
///	Routine: Assignment operator
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_GenPricer& ARM_GenPricer::operator=( const ARM_GenPricer& rhs )
{
	if( this !=	 &rhs )
	{
	    ARM_RootObject::operator=( rhs );
		itsGenSecurity				= rhs.itsGenSecurity;
		itsPricingMod				= rhs.itsPricingMod;
		delete itsPricerInfo;
        itsPricerInfo				= rhs.itsPricerInfo? (ARM_PricerInfo*) rhs.itsPricerInfo->Clone() : NULL;
		delete itsCVInfo;
		itsCVInfo		= rhs.itsCVInfo? (ARM_CVInfo*) rhs.itsCVInfo->Clone() : NULL;
		itsDetailMode				= rhs.itsDetailMode;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_GenPricer
///	Routine: Price
///	Returns: a double corresponding to the price of the GenericSecurity
///				according to the PricingModel
///	Action : Compute price!
////////////////////////////////////////////////////

double ARM_GenPricer::Price()
{
	/// not to be removed as the pricer info is created later on!
	/// cannot use timer!
	clock_t startTime = clock();

    const string payModelName = itsGenSecurity->GetPayModelName();
	itsPricingMod->SetPayModelName( payModelName );

	/// Get pricing Adviser and decides pricing!
	ARM_PricingAdviserPtr pricingAdviser = itsGenSecurity->GetPricingAdviser();
	pricingAdviser->DecidePricing( itsPricingMod );

	/// Init is a pre-treatment by the model of the security!
	ARM_TimeInfoPtrVector TimeInfos	= pricingAdviser->GetTimeInfos( itsPricingMod );

	ARM_PricingStatesPtr states		= itsPricingMod->Init(payModelName, TimeInfos );

	size_t NumPricedColumns = itsGenSecurity->GetDealDescription().GetNbPricedColumns(), Col;
	vector<string> PricedColumnNames = itsGenSecurity->GetDealDescription().GetPricedColumnNames();
	
	if (!states.IsNull())
	{
		states->resizePayoffStatesVector(NumPricedColumns);
		states->SetOtherPayoffsFlag(itsGenSecurity->GetOtherPayoffsFlag() && itsPricingMod->GetOtherPayoffsFlag());
		states->SetIVFlag(itsGenSecurity->GetOtherPayoffsFlag() && itsPricingMod->GetOtherPayoffsFlag() && itsGenSecurity->GetIVFlag());
	}
	pricingAdviser->InitPricing( itsPricingMod );

	/// Get the nodes to process as iterator!
	size_t firstNodeIdx = pricingAdviser->GetFirstNodeIdx();
	size_t lastNodeIdx = pricingAdviser->GetLastNodeIdx();
	size_t currentNodeIdx = firstNodeIdx;
	ARM_ExpNodePtr pNode;

	/// Get the induct time and the corresponding iterator
	ARM_VectorPtr toTimes				= pricingAdviser->GetToTimes();
	ARM_GP_Vector::iterator toTime		= pricingAdviser->GetFirstToTime();
	ARM_GP_Vector::iterator endToTime	= pricingAdviser->GetLastToTime();

	// from times and to times defines the time interval on which we induct
	ARM_VectorPtr fromTimes				= pricingAdviser->GetFromTimes();
	ARM_GP_Vector::iterator fromTime	= pricingAdviser->GetFirstFromTime();
	ARM_GP_Vector::iterator endFromTime	= pricingAdviser->GetLastFromTime();

	/// check for closed form solution ... in this case no need to loop
	/// and induct since everything is computed from day one!
	if(fabs( *fromTime - *toTime) < K_NEW_DOUBLE_TOL && *fromTime < K_NEW_DOUBLE_TOL)
	{

	    /// a closed form pricing is always with a single loop pricer info
        delete itsPricerInfo;

		if (NumPricedColumns > 1)
		{

			// Take care the pricer info copy constructor is not deep !
			vector<ARM_PricerInfo*> pricerInfos;
			for(size_t i = 0; i < NumPricedColumns; ++i)
				pricerInfos.push_back(itsPricingMod->CreatePricerInfo());

			itsPricerInfo = new ARM_PInfo_Map(PricedColumnNames, pricerInfos);
		}
		else
			itsPricerInfo = new ARM_PInfo_SingleLoop;

        /// start timer
	    itsPricerInfo->SetStartTime(startTime);

		for (Col = 0; Col < NumPricedColumns; ++Col)
		{
			pNode = pricingAdviser->GetNode(currentNodeIdx, Col);
			EvalChildNodeAndProcess( pNode, itsPricingMod, states, Col, CreateFirstPayoffClosedForms(), *fromTime );
		}
		itsPricerInfo->StoreCurrentLoopPrice( states, 1, 0, pricingAdviser->GetPricingDir() );

	}
    /// an only analytical model cannot induct!
	else if( itsPricingMod->IsCurrentlyOnlyAnalyticalModel() )
    {
        throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL, "Analytical pricing is only supported for payoff with eventdate = today!" );
    }
    /// other case for a numerical method
    else
    {
	    /// a closed form pricing is always with a single loop pricer info
        delete itsPricerInfo;

		if (NumPricedColumns > 1)
		{

			// Take care the pricer info copy constructor is not deep !
			vector<ARM_PricerInfo*> pricerInfos;
			for(size_t i = 0; i < NumPricedColumns; ++i)
				pricerInfos.push_back(itsPricingMod->CreatePricerInfo());

			itsPricerInfo = new ARM_PInfo_Map(PricedColumnNames, pricerInfos);
		}
		else
		{
			itsPricerInfo = itsPricingMod->CreatePricerInfo();
		}


        /// start timer
	    itsPricerInfo->SetStartTime(startTime);

	    /// reset in particular the numeraire! should be here after the pricing adviser initialisation!
	    itsPricerInfo->reset();

	    /// get number of iteration!
	    size_t bucketsNb = itsPricingMod->GetNumMethod()->GetBuckets()->size();

		// we first iterate on the buckets
	    for(size_t i=0; i<bucketsNb; ++i )
	    {
			/// do not loop in the case of empty bucket!
			///if( (*itsPricingMod->GetNumMethod()->GetBuckets())[i] )
			{
				size_t LoopNb = itsPricingMod->GetNumMethod()->GetLoopNb();

				// we then iterate on the "loops". This is useful for backward forward nummethods only. 
				for( size_t j=0 ; j<LoopNb ; ++j )
				{
					// We then iterate on nodes. 
				
					/// first evaluates the last payoff
					for (Col = 0; Col < NumPricedColumns; ++Col)
					{
						pNode = pricingAdviser->GetNode(currentNodeIdx, Col);
						EvalChildNodeAndProcess( pNode, itsPricingMod, states, Col, CreateFirstPayoff(payModelName), *fromTime );
					}
				
					/// get the correspnding auxiliary nodes corresponding to cell reference calls!
					size_t newNodes = pricingAdviser->GetAuxiliaryNodes( itsPricingMod->GetDateFromTime( *fromTime ), itsPricingMod, states );
				
					/// induct
					states = itsPricingMod->Induct( states, *toTime );
				
					/// add the other computed nodes
					size_t removedNodes = pricingAdviser->SetComputedNodes( itsPricingMod->GetDateFromTime( *toTime ), itsPricingMod, states );

					pricingAdviser->Increment(currentNodeIdx, fromTime, toTime);

					/// the pricingAdviser is responsible for pricing in the correct direction the deal!
					while( currentNodeIdx!=lastNodeIdx )
					{
			
					/// test for out of bound from and to times!
					/// only in strict validation mode
		#if defined(__GP_STRICT_VALIDATION)
					if( fromTime == endFromTime )
						throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL, "from time out of bound!" );
					if( toTime == endToTime )
						throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL, "to time out of bound!" );
		#endif
					
						/// evaluates the current payoff and add to current payoff
						for (Col = 0; Col < NumPricedColumns; ++Col)
						{
							pNode = pricingAdviser->GetNode(currentNodeIdx, Col);
							EvalChildNodeAndProcess( pNode, itsPricingMod, states, Col, AddToState0(payModelName), *fromTime );
						}
					
						/// get the correspnding auxiliary nodes corresponding to cell reference calls!
						size_t newNodes = pricingAdviser->GetAuxiliaryNodes( itsPricingMod->GetDateFromTime( *fromTime ), itsPricingMod, states );

						/// induct
						states = itsPricingMod->Induct( states, *toTime );
					
						/// add the other computed nodes
						size_t removedNodes = pricingAdviser->SetComputedNodes( itsPricingMod->GetDateFromTime( *toTime ), itsPricingMod, states );

						pricingAdviser->Increment(currentNodeIdx, fromTime, toTime);
					}

					if ( LoopNb > 1 && j < LoopNb-1 )
					{
						pricingAdviser->ResetPrecomputedNodes( ARM_ExpNode::ResetLoop );
						states	 = itsPricingMod->ReInitLoop();
						states->resizePayoffStatesVector(NumPricedColumns);
						states->SetOtherPayoffsFlag(itsGenSecurity->GetOtherPayoffsFlag() && itsPricingMod->GetOtherPayoffsFlag());
						states->SetIVFlag(itsGenSecurity->GetOtherPayoffsFlag() && itsPricingMod->GetOtherPayoffsFlag() && itsGenSecurity->GetIVFlag());
						pricingAdviser->InitPricing( itsPricingMod );
						firstNodeIdx = pricingAdviser->GetFirstNodeIdx();
						currentNodeIdx = firstNodeIdx;
						lastNodeIdx	  = pricingAdviser->GetLastNodeIdx();

						/// Get the induct time and the corresponding iterator
						toTimes		= pricingAdviser->GetToTimes();
						toTime		= pricingAdviser->GetFirstToTime();
						endToTime	= pricingAdviser->GetLastToTime();

						// from times and to times defines the time interval on which we induct
						fromTimes				= pricingAdviser->GetFromTimes();
						fromTime		= pricingAdviser->GetFirstFromTime();
						endFromTime		= pricingAdviser->GetLastFromTime();
					}

				}
				/// storing the current loop price
				itsPricerInfo->StoreCurrentLoopPrice( states, itsPricingMod->GetNumMethod()->GetBuckets()->Elt( itsPricingMod->GetNumMethod()->GetBucketIndex() ), i, itsPricingMod->GetNumMethod()->GetPricingDirection() );

				/// Reset the initial state!
				if(bucketsNb>1 && i<bucketsNb-1 )
				{
					pricingAdviser->ResetPrecomputedNodes( ARM_ExpNode::ResetBucket );
					states	 = itsPricingMod->ReInit();
					states->resizePayoffStatesVector(NumPricedColumns);
					states->SetOtherPayoffsFlag(itsGenSecurity->GetOtherPayoffsFlag() && itsPricingMod->GetOtherPayoffsFlag());
					states->SetIVFlag(itsGenSecurity->GetOtherPayoffsFlag() && itsPricingMod->GetOtherPayoffsFlag() && itsGenSecurity->GetIVFlag());
					pricingAdviser->InitPricing( itsPricingMod );
					firstNodeIdx = pricingAdviser->GetFirstNodeIdx();
					currentNodeIdx = firstNodeIdx;
					lastNodeIdx	  = pricingAdviser->GetLastNodeIdx();

					/// Get the induct time and the corresponding iterator
					toTimes		= pricingAdviser->GetToTimes();
					toTime		= pricingAdviser->GetFirstToTime();
					endToTime	= pricingAdviser->GetLastToTime();

					// from times and to times defines the time interval on which we induct
					fromTimes				= pricingAdviser->GetFromTimes();
					fromTime		= pricingAdviser->GetFirstFromTime();
					endFromTime		= pricingAdviser->GetLastFromTime();
				}
			}
		}
	}
	
	if( itsGenSecurity->GetExercBoundaryResetFlag() )
		pricingAdviser->ResetPrecomputedNodes( ARM_ExpNode::ResetTotal );
	else
		pricingAdviser->ResetPrecomputedNodes( ARM_ExpNode::ResetBucket );

    /// remove memory pool allocation
	states->FreeMemoryPool();

	itsPricerInfo->PostProcess(TimeInfos,itsPricingMod->GetNumMethod());

	/// Finalize pricingmodel (frees useless memory used by pricingmodel).
	itsPricingMod->Finalize();

	/// Take care of control variates
	if( itsCVInfo )
		ComputeCV();

	/// end timer
	itsPricerInfo->ClockEndTime();

	return itsPricerInfo->GetPrice();
}


////////////////////////////////////////////////////
///	Class  : ARM_GenPricer
///	Routine: Price
///	Returns: a double corresponding to the price of the GenericSecurity
///				according to the PricingModel
///	Action : Compute price!
////////////////////////////////////////////////////

void ARM_GenPricer::ComputeCV()
{
	ARM_PInfo_Map* pInfoMap = dynamic_cast<ARM_PInfo_Map*>( itsPricerInfo );
	if( !pInfoMap )
		ARM_THROW( ERR_INVALID_ARGUMENT, " the pricer info is not a price info map, dynamic_cast failes!" );

	ARM_PInfo_Cov* infoCov = new ARM_PInfo_Cov( pInfoMap, itsCVInfo->GetRefPriceColumn(), itsCVInfo->GetCVColumns(), itsCVInfo->GetCVPrices(), itsCVInfo->GetBetas() ) ;
	itsCVInfo->SetBetas (*(infoCov->GetBetas()));
	ARM_PricerInfoPtr pinfoCov = ARM_PricerInfoPtr( infoCov );

	/// swap first element and control variate (since the main result is the first element)
	pair<string,ARM_PricerInfoPtr> firstElem = (*pInfoMap)[0];
	(*pInfoMap)[0] = pair<string,ARM_PricerInfoPtr>( ARM_GenPricer::CVTag, pinfoCov );
	pInfoMap->push_back( firstElem.first, firstElem.second );
}

CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

