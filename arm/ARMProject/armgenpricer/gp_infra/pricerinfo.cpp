/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file pricerinfo.cpp
 *
 *  \brief pricer info is a dictionary to contain all the result
 *		of a pricing request
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date February 2004
 */

/// this headers has to come first
/// as it contains pragma for stupid warning
#include "gpbase/removeidentifiedwarning.h"

#include "gpbase/env.h"
#include "gpinfra/pricerinfo.h"
#include "gpbase/gpmatrixlinalg.h"
#include "gpbase/stringmanip.h"
#include "gpbase/badalloc.h"

/// gpbase
#include "gpbase/gpvector.h"
#include "gpbase/ostringstream.h"

/// gpinfra
#include "gpinfra/typedef.h"
#include "gpinfra/pricingmodel.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/nummethod.h"
#include "gpinfra/gramfunctorargdict.h"
#include "gpinfra/timeinfo.h"

// gpnumlib
#include "gpnumlib/stddevfunc.h"

#include <algorithm>
#include <utility>
#include <climits>
#include <new>
CC_USING_NS( std, make_pair )
CC_USING_NS( std, bad_alloc )


CC_BEGIN_NAMESPACE( ARM )

#define PRICER_INFO_MULTILOOP_RUNNING_STEP_CST 10

////////////////////////////////////////////////////
///	Class  : ARM_PricerInfo
///	Routine: constructor
////////////////////////////////////////////////////

/// timer part
ARM_PricerInfo::ARM_PricerInfo()
:	ARM_RootObject(),
	ARM_Timer()
{
	CC_ARM_SETNAME(ARM_PRICERINFO);
}


////////////////////////////////////////////////////
///	Class  : ARM_PricerInfo
///	Routine: copy constructor
////////////////////////////////////////////////////

/// timer part
ARM_PricerInfo::ARM_PricerInfo( const ARM_PricerInfo& rhs )
:	ARM_RootObject(rhs),
	/// timer part
	ARM_Timer(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_PricerInfo
///	Routine: assignment operator
////////////////////////////////////////////////////
ARM_PricerInfo& ARM_PricerInfo::operator=(const ARM_PricerInfo& rhs )
{
	if( this != &rhs )
	{
		ARM_RootObject::operator=(rhs);
		ARM_Timer::operator=(rhs);
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_PricerInfo
///	Routine: destructor
////////////////////////////////////////////////////
ARM_PricerInfo::~ARM_PricerInfo()
{}




////////////////////////////////////////////////////
///	Class  : ARM_PInfo_SingleLoop
///	Routine: toString
///	Returns: 
///	Action : returns the price and does a toString
////////////////////////////////////////////////////
string ARM_PricerInfo::toString(const string& indent, const string& nextIndent ) const
{
	CC_Ostringstream os;
	os << GetName() << " pricer\n";
	os << GetContents().toString();
	return os.str();
}




////////////////////////////////////////////////////
////////////////////////////////////////////////////
////		ARM_PInfo_SingleLoop				////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_PInfo_SingleLoop
///	Routine: Constructor,Copy constructor, assignment operator
///				Destructor and Clone
///	Returns: 
///	Action : builds the object,....
////////////////////////////////////////////////////
ARM_PInfo_SingleLoop::ARM_PInfo_SingleLoop()
: ARM_PricerInfo()
{}

ARM_PInfo_SingleLoop::ARM_PInfo_SingleLoop( const ARM_PInfo_SingleLoop& rhs )
:	ARM_PricerInfo(rhs)
{}

ARM_PInfo_SingleLoop& ARM_PInfo_SingleLoop::operator=( const ARM_PInfo_SingleLoop& rhs )
{
	if( this != &rhs )
	{
		ARM_PricerInfo::operator=(rhs);
	}
	return *this;
}

ARM_PInfo_SingleLoop::~ARM_PInfo_SingleLoop() 
{}


ARM_Object* ARM_PInfo_SingleLoop::Clone() const
{
	return new ARM_PInfo_SingleLoop( *this );
}


////////////////////////////////////////////////////
///	Class  : ARM_PInfo_SingleLoop
///	Routine: StoreCurrentLoopPrice,PostProcess,GetPrice
///	Returns: 
///	Action : store information about current loop price,
///		postprocess, get price and show the result as a string (toString)
////////////////////////////////////////////////////
void ARM_PInfo_SingleLoop::StoreCurrentLoopPrice( const ARM_PricingStatesPtr& states, size_t bucketSize, size_t currentBucketIndex, ARM_PricingAdviser::PricingDir pricingDir, int colIdx )
{
#if defined(__GP_STRICT_VALIDATION)
	if( bucketSize != 1 )
		
		Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			ARM_USERNAME + ": single loop method with bucket of size != 1 not supported" );
#endif
	itsPrice = states->GetPayoffStates(colIdx).GetPayoff(0);

	size_t nbIntermediatePayoffs = states->GetPayoffStates(colIdx).GetIntermediatePayoffsSize();

	itsIntermediatePrices = ARM_GP_VectorPtr(new ARM_GP_Vector(nbIntermediatePayoffs));

	bool isFwd = (pricingDir == ARM_NumMethod::GP_FWDLOOKING);

    int i;
	for (size_t i = 0;i < nbIntermediatePayoffs;++i)
		(*itsIntermediatePrices)[i] = states->GetPayoffStates(colIdx).GetIntermediatePayoff(0,(isFwd?i:nbIntermediatePayoffs-i-1));

    /// Payoff snapshots
	size_t nbSnapshots = states->GetPayoffStates(colIdx).GetPayoffSnapshotsSize();
    size_t nbStates,maxStates = 0;
	for (size_t i=0;i<nbSnapshots;++i)
    {
        nbStates = states->GetPayoffStates(colIdx).GetPayoffSnapshot(i)->size();
        maxStates = (maxStates < nbStates ? nbStates : maxStates);
    }
    itsPriceSnapshots = ARM_GP_MatrixPtr(new ARM_GP_Matrix(maxStates,nbSnapshots));

    ARM_GP_VectorPtr prices;
    int j;
	for (size_t i=0;i<nbSnapshots;++i)
    {
        prices = states->GetPayoffStates(colIdx).GetPayoffSnapshot(isFwd ? i : nbSnapshots-i-1);
        nbStates = prices->size();
	    for(j=0;j<nbStates;++j)
		    (*itsPriceSnapshots)(j,i) = (*prices)[j];
	    for(;j<maxStates;++j)
		    (*itsPriceSnapshots)(j,i) = 0.0;
    }
}

double ARM_PInfo_SingleLoop::GetPrice() const
{	
	return itsPrice;
}



void ARM_PInfo_SingleLoop::PostProcess(const ARM_TimeInfoPtrVector& timeInfo,const const ARM_NumMethodPtr& numMethod)
{
    /// Get probabilities to access to future states of numerical method
	std::vector<double> eventTimes(timeInfo.size());
    for(size_t i=0;i<timeInfo.size();++i)
        eventTimes[i]=timeInfo[i]->GetEventTime();

    if(numMethod != ARM_NumMethodPtr(NULL) )
        itsSpotProbabilities=numMethod->GetSpotProbabilities(eventTimes);
}


////////////////////////////////////////////////////
///	Class  : ARM_PInfo_SingleLoop
///	Routine: GetContents
///	Returns: 
///	Action : Creates the corresponding dicitionary
////////////////////////////////////////////////////
ARM_GramFctorArgDict ARM_PInfo_SingleLoop::GetContents(const string& colName) const
{
	ARM_GramFctorArgDict contents;
	contents[ "Price"		]			= itsPrice;
	contents[ "IntermediatePrices" ]	= itsIntermediatePrices;
	contents[ "PriceSnapshots" ]	    = itsPriceSnapshots;
	contents[ "SpotProbabilities" ]	    = itsSpotProbabilities;
	contents[ "Duration"    ]			= GetDuration();
	return contents;	
}

////////////////////////////////////////////////////
////////////////////////////////////////////////////
////		ARM_PInfo_SingleLoopPDE				////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_PInfo_SingleLoop
///	Routine: Constructor,Copy constructor, assignment operator
///				Destructor and Clone
///	Returns: 
///	Action : builds the object,....
////////////////////////////////////////////////////
ARM_PInfo_SingleLoopPDE::ARM_PInfo_SingleLoopPDE()
: ARM_PInfo_SingleLoop()
{}

ARM_PInfo_SingleLoopPDE::ARM_PInfo_SingleLoopPDE( const ARM_PInfo_SingleLoopPDE& rhs )
:	ARM_PInfo_SingleLoop(rhs), itsPriceIndex( rhs.itsPriceIndex )
{}

ARM_PInfo_SingleLoopPDE& ARM_PInfo_SingleLoopPDE::operator=( const ARM_PInfo_SingleLoopPDE& rhs )
{
	if( this != &rhs )
	{
		ARM_PricerInfo::operator=(rhs);
	}
	return *this;
}

ARM_PInfo_SingleLoopPDE::~ARM_PInfo_SingleLoopPDE() 
{}


ARM_Object* ARM_PInfo_SingleLoopPDE::Clone() const
{
	return new ARM_PInfo_SingleLoopPDE( *this );
}


////////////////////////////////////////////////////
///	Class  : ARM_PInfo_SingleLoop
///	Routine: StoreCurrentLoopPrice,PostProcess,GetPrice
///	Returns: 
///	Action : store information about current loop price,
///		postprocess, get price and show the result as a string (toString)
////////////////////////////////////////////////////
void ARM_PInfo_SingleLoopPDE::StoreCurrentLoopPrice( const ARM_PricingStatesPtr& states, size_t bucketSize, size_t currentBucketIndex, ARM_PricingAdviser::PricingDir pricingDir, int colIdx )
{
#if defined(__GP_STRICT_VALIDATION)
	if( bucketSize != 1 )
		
		Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			ARM_USERNAME + ": single loop method with bucket of size != 1 not supported" );
#endif
	itsPrice = states->GetPayoffStates(colIdx).GetPayoff(itsPriceIndex);

	size_t nbIntermediatePayoffs = states->GetPayoffStates(colIdx).GetIntermediatePayoffsSize();

	itsIntermediatePrices = ARM_GP_VectorPtr(new ARM_GP_Vector(nbIntermediatePayoffs));

	bool isFwd = (pricingDir == ARM_NumMethod::GP_FWDLOOKING);

 
	for (size_t i = 0;i < nbIntermediatePayoffs;++i)
		(*itsIntermediatePrices)[i] = states->GetPayoffStates(colIdx).GetIntermediatePayoff(itsPriceIndex,(isFwd?i:nbIntermediatePayoffs-i-1));

    /// Payoff snapshots
	size_t nbSnapshots = states->GetPayoffStates(colIdx).GetPayoffSnapshotsSize();
    size_t nbStates,maxStates = 0;
	for (size_t i=0;i<nbSnapshots;++i)
    {
        nbStates = states->GetPayoffStates(colIdx).GetPayoffSnapshot(i)->size();
        maxStates = (maxStates < nbStates ? nbStates : maxStates);
    }
    itsPriceSnapshots = ARM_GP_MatrixPtr(new ARM_GP_Matrix(maxStates,nbSnapshots));

    ARM_GP_VectorPtr prices;
    int j;
	for (size_t i=0;i<nbSnapshots;++i)
    {
        prices = states->GetPayoffStates(colIdx).GetPayoffSnapshot(isFwd ? i : nbSnapshots-i-1);
        nbStates = prices->size();
	    for(j=0;j<nbStates;++j)
		    (*itsPriceSnapshots)(j,i) = (*prices)[j];
	    for(;j<maxStates;++j)
		    (*itsPriceSnapshots)(j,i) = 0.0;
    }
}

double ARM_PInfo_SingleLoopPDE::GetPrice() const
{	
	return itsPrice;
}



void ARM_PInfo_SingleLoopPDE::PostProcess(const ARM_TimeInfoPtrVector& timeInfo,const const ARM_NumMethodPtr& numMethod)
{
    /// Get probabilities to access to future states of numerical method
	std::vector<double> eventTimes(timeInfo.size());
    for(size_t i=0;i<timeInfo.size();++i)
        eventTimes[i]=timeInfo[i]->GetEventTime();

    if(numMethod != ARM_NumMethodPtr(NULL) )
        itsSpotProbabilities=numMethod->GetSpotProbabilities(eventTimes);
}


////////////////////////////////////////////////////
///	Class  : ARM_PInfo_SingleLoop
///	Routine: GetContents
///	Returns: 
///	Action : Creates the corresponding dicitionary
////////////////////////////////////////////////////
ARM_GramFctorArgDict ARM_PInfo_SingleLoopPDE::GetContents(const string& colName) const
{
	ARM_GramFctorArgDict contents;
	contents[ "Price"		]			= itsPrice;
	contents[ "IntermediatePrices" ]	= itsIntermediatePrices;
	contents[ "PriceSnapshots" ]	    = itsPriceSnapshots;
	contents[ "SpotProbabilities" ]	    = itsSpotProbabilities;
	contents[ "Duration"    ]			= GetDuration();
	return contents;	
}

////////////////////////////////////////////////////
////////////////////////////////////////////////////
////		ARM_PInfo_MultipleLoop				////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_PInfo_MultipleLoop
///	Routine: Constructor,Copy constructor, assignment operator
///				Destructor and Clone
///	Returns: 
///	Action : builds the object
////////////////////////////////////////////////////
ARM_PInfo_MultipleLoop::ARM_PInfo_MultipleLoop(const ARM_PricingModel& mod, const ARM_MomentFuncPtr& momentFunc )
:
	ARM_PricerInfo(), 
	itsPrice(0),
	itsIntermediatePrices(NULL),
	itsStdDev(0),
	itsIntermediateStdDevs(NULL),
	itsValues(NULL),
	itsIntermediateValues(NULL),
	itsRunningAverage(NULL),
	itsIntermediateRunningAverages(NULL),
	itsEmpiricalStdDev(0),
	itsIntermediateEmpiricalStdDevs(NULL),
	itsBuckets(0),
	itsPricesPerBucket(0),
	itsStdDevPerBucket(0),
	itsValuesPerBucket(0),
	itsValuesIndex(0),
	itsMomentFunc(momentFunc),
	itsIntermediateArrayInitialized(false)
{
#if defined(__GP_STRICT_VALIDATION)
	/// not strictly necessary but for robustness
	if(mod.GetNumMethod() == ARM_NumMethodPtr(NULL) )	
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			ARM_USERNAME + ": no numerical method set for this pricing model" );
#endif
	itsBuckets		= mod.GetNumMethod()->GetBuckets();
}


ARM_PInfo_MultipleLoop::ARM_PInfo_MultipleLoop( const ARM_PInfo_MultipleLoop& rhs )
:
	ARM_PricerInfo(	rhs														), 
	itsPrice(rhs.itsPrice													),
	itsIntermediatePrices( rhs.itsIntermediatePrices						),
	itsStdDev( rhs.itsStdDev												),
	itsIntermediateStdDevs( rhs.itsIntermediateStdDevs						),
	itsValues( rhs.itsValues												),
	itsIntermediateValues( rhs.itsIntermediateValues						),
	itsRunningAverage( rhs.itsRunningAverage								),
	itsIntermediateRunningAverages(rhs.itsIntermediateRunningAverages		),
	itsEmpiricalStdDev( rhs.itsEmpiricalStdDev								),
	itsIntermediateEmpiricalStdDevs( rhs.itsIntermediateEmpiricalStdDevs	),
	itsRunningAvgStep( rhs.itsRunningAvgStep								),
	itsBuckets( rhs.itsBuckets												),
	itsValuesPerBucket( rhs.itsValuesPerBucket                              ),
	itsPricesPerBucket( rhs.itsPricesPerBucket                              ),
	itsStdDevPerBucket( rhs.itsStdDevPerBucket                              ),
	itsValuesIndex( rhs.itsValuesIndex										),
	itsMomentFunc( rhs.itsMomentFunc							),
	itsInterMomentFunc( rhs.itsInterMomentFunc				),
	itsIntermediateArrayInitialized( rhs.itsIntermediateArrayInitialized	)
{}


ARM_PInfo_MultipleLoop& ARM_PInfo_MultipleLoop::operator=( const ARM_PInfo_MultipleLoop& rhs )
{
	if( this != &rhs )
	{
		ARM_PricerInfo::operator=(rhs);
		itsPrice						= rhs.itsPrice;
		itsIntermediatePrices			= rhs.itsIntermediatePrices;
		itsStdDev						= rhs.itsStdDev;
		itsIntermediateStdDevs			= rhs.itsIntermediateStdDevs;
		itsValues						= rhs.itsValues;
		itsIntermediateValues			= rhs.itsIntermediateValues;
		itsRunningAverage				= rhs.itsRunningAverage;
		itsIntermediateRunningAverages	= rhs.itsIntermediateRunningAverages;
		itsEmpiricalStdDev				= rhs.itsEmpiricalStdDev;
		itsIntermediateEmpiricalStdDevs	= rhs.itsIntermediateEmpiricalStdDevs;
		itsRunningAvgStep				= rhs.itsRunningAvgStep;
		itsBuckets						= rhs.itsBuckets;
		itsPricesPerBucket              = rhs.itsPricesPerBucket;
		itsStdDevPerBucket              = rhs.itsStdDevPerBucket;
		itsValuesPerBucket              = rhs.itsValuesPerBucket;
		itsValuesIndex					= rhs.itsValuesIndex;
		itsMomentFunc					= rhs.itsMomentFunc;
		itsInterMomentFunc				= rhs.itsInterMomentFunc;
		itsIntermediateArrayInitialized	= rhs.itsIntermediateArrayInitialized;
	}
	return *this;
}


ARM_PInfo_MultipleLoop::~ARM_PInfo_MultipleLoop() 
{}


ARM_Object* ARM_PInfo_MultipleLoop::Clone() const
{
	return new ARM_PInfo_MultipleLoop( *this );
}


////////////////////////////////////////////////////
///	Class  : ARM_PInfo_MultipleLoop
///	Routine: Constructor,Destructor and Clone
///	Returns: 
///	Action : builds the object
////////////////////////////////////////////////////
void ARM_PInfo_MultipleLoop::reset()
{
	itsValuesIndex  = 0;
}


////////////////////////////////////////////////////
///	Class  : ARM_PInfo_MultipleLoop
///	Routine: StoreCurrentLoopPrice,PostProcess,GetPrice
///	Returns: 
///	Action : store information about current loop price,
///		postprocess, get price and show the result as a string (toString)
////////////////////////////////////////////////////
void ARM_PInfo_MultipleLoop::StoreCurrentLoopPrice( const ARM_PricingStatesPtr& states, size_t bucketSize, size_t currentBucketIndex, ARM_PricingAdviser::PricingDir pricingDir, int colIdx)
{
	size_t nbIntermediatePayoffs = states->GetPayoffStates(colIdx).GetIntermediatePayoffsSize();
	bool otherPayoffFlag = states->GetPayoffStates(colIdx).GetOtherPayoffsFlag();
	bool ivFlag = states->GetPayoffStates(colIdx).GetIVFlag();

	// This initialization should be done just once.
	if (!itsIntermediateArrayInitialized)
	{	size_t totalSize = 0;
		for(size_t i=0; i<itsBuckets->size(); ++i )
			totalSize += (*itsBuckets)[i];
		itsIntermediatePrices			= ARM_VectorPtr(new std::vector<double>(nbIntermediatePayoffs));
		itsIntermediateStdDevs			= ARM_VectorPtr(new std::vector<double>(nbIntermediatePayoffs));
		itsIntermediateEmpiricalStdDevs	= ARM_VectorPtr(new std::vector<double>(nbIntermediatePayoffs));
		itsIntermediateRunningAverages	= ARM_GP_MatrixPtr(new ARM_GP_Matrix);

		if (otherPayoffFlag)
		{
			itsValues				= ARM_VectorPtr( new std::vector<double>(totalSize,0.0) );
			if (ivFlag)
			{
				itsIntermediateValues	= ARM_GP_MatrixPtr(new ARM_GP_Matrix(nbIntermediatePayoffs,totalSize));
			}
		}

		// We can have 1000 points maximum in the running average vector
		itsRunningAvgStep = ceil((double)totalSize/1000);
		itsMomentFunc->init(totalSize,itsRunningAvgStep);
		
		for (size_t i = 0; i < nbIntermediatePayoffs; ++i)
		{
			itsInterMomentFunc.push_back(ARM_MomentFuncPtr(itsMomentFunc->Clone()));
			itsInterMomentFunc[i]->init(totalSize,itsRunningAvgStep);
		}

		itsIntermediateArrayInitialized = true;
	}

	std::vector<double> curBucketValues(bucketSize);
	ARM_GP_Matrix interCurBucketValues(nbIntermediatePayoffs,bucketSize);

	bool isFwd = (pricingDir == ARM_NumMethod::GP_FWDLOOKING);

	double currentPrice;

	size_t i,j;

	for (size_t i=0;i<bucketSize;++i)
	{
		currentPrice	= states->GetPayoffStates(colIdx).GetPayoff(i);

#if defined(__GP_STRICT_VALIDATION)
		/// to catch up out of  bound numbers
		if ( currentPrice < -1e+308 || currentPrice >1e+308 )
		if( currentPrice < -1.0e+300 || currentPrice >1.0e+300 )
			ARM_THROW( ERR_INVALID_ARGUMENT, "currentPrice is out of bound!" );
#endif

		curBucketValues[i] = currentPrice;
		if (otherPayoffFlag)
			(*itsValues)[itsValuesIndex] = currentPrice;
		for (j = 0; j < nbIntermediatePayoffs; ++j)
		{
			currentPrice = states->GetPayoffStates(colIdx).GetIntermediatePayoff(i,(isFwd?j:nbIntermediatePayoffs-j-1));
			interCurBucketValues(j,i) = currentPrice;
			if (otherPayoffFlag && ivFlag)
			{
				(*itsIntermediateValues)(j,itsValuesIndex) = currentPrice;
			}
		}
		

		itsValuesIndex++;
	}
	
	itsMomentFunc->addBucket(curBucketValues);
	for (size_t i=0;i<nbIntermediatePayoffs;++i)
	{
		ARM_GP_Vector* tmpRow = interCurBucketValues.GetRow(i);
		itsInterMomentFunc[i]->addBucket(tmpRow->GetValues());
		delete tmpRow;
	}
}


void ARM_PInfo_MultipleLoop::PostProcess(const ARM_TimeInfoPtrVector& timeInfo,const ARM_NumMethodPtr& numMethod)
{
	itsMomentFunc->finalize();
	
	size_t i;

	for (size_t i=0;i<itsInterMomentFunc.size();++i)
		itsInterMomentFunc[i]->finalize();

	size_t nbIntermediatePayoffs = itsIntermediatePrices->size();

	itsPrice = itsMomentFunc->GetMean();
	itsStdDev = itsMomentFunc->GetStdDev();
	itsRunningAverage = itsMomentFunc->GetRunningAverage();
	itsEmpiricalStdDev = itsMomentFunc->GetEmpiricalStdDev();
	itsPricesPerBucket = itsMomentFunc->GetMeanPerBucket();
	itsStdDevPerBucket = itsMomentFunc->GetStdDevPerBucket();

	for (size_t i = 0; i < itsInterMomentFunc.size(); ++i)
	{
		(*itsIntermediatePrices)[i] = itsInterMomentFunc[i]->GetMean();
		(*itsIntermediateStdDevs)[i] = itsInterMomentFunc[i]->GetStdDev();
		(*itsIntermediateEmpiricalStdDevs)[i] = itsInterMomentFunc[i]->GetEmpiricalStdDev();

		ARM_GP_VectorPtr runningAverage = itsInterMomentFunc[i]->GetRunningAverage();
		itsIntermediateRunningAverages->push_backRow(*runningAverage);
	}

#if defined(__GP_STRICT_VALIDATION)
	if( itsStdDev < -K_NEW_DOUBLE_TOL )
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			ARM_USERNAME + ": a variance cannot be negative!" );

	if( itsEmpiricalStdDev < -K_NEW_DOUBLE_TOL )
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			ARM_USERNAME + ": a variance cannot be negative!" );

	for (size_t j = 0; j < nbIntermediatePayoffs; ++j)
	{
		if( (*itsIntermediateStdDevs)[j] < -K_NEW_DOUBLE_TOL )
			throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				ARM_USERNAME + ": a variance cannot be negative!" );

		if( (*itsIntermediateEmpiricalStdDevs)[j] < -K_NEW_DOUBLE_TOL )
			throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				ARM_USERNAME + ": a variance cannot be negative!" );
	}
#endif
}


double ARM_PInfo_MultipleLoop::GetPrice() const
{	
	return itsPrice; 
}

////////////////////////////////////////////////////
///	Class  : ARM_PInfo_MultipleLoop
///	Routine: GetContents
///	Returns: 
///	Action : Creates the corresponding dicitionary
////////////////////////////////////////////////////
ARM_GramFctorArgDict ARM_PInfo_MultipleLoop::GetContents(const string& colName) const
{
	ARM_GramFctorArgDict contents;
	contents[ "EmpiricalStdDev"					] = itsEmpiricalStdDev;
	contents[ "IntermediateEmpiricalStdDevs"	] = itsIntermediateEmpiricalStdDevs;
	contents[ "Price"							] = itsPrice;
	contents[ "Prices Per Bucket"				] = itsPricesPerBucket;
	contents[ "StdDev Per Bucket"				] = itsStdDevPerBucket;
	contents[ "IntermediatePrices"				] = itsIntermediatePrices;
	contents[ "StdDev"							] = itsStdDev;
	contents[ "IntermediateStdDevs"				] = itsIntermediateStdDevs;
	contents[ "Values"							] = itsValues;
	contents[ "IntermediateValues"				] = itsIntermediateValues;
	contents[ "RunningAvg"						] = itsRunningAverage;
	contents[ "IntermediateRunningAvgs"			] = itsIntermediateRunningAverages;
	contents[ "RunningAvgStep"					] = itsRunningAvgStep;
	contents[ "Buckets"							] = itsBuckets;
	contents[ "Duration"						] = GetDuration();
	return contents;	
}


////////////////////////////////////////////////////
////////////////////////////////////////////////////
////		ARM_PInfo_Cov						////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_PInfo_Cov
///	Routine: ARM_PInfo_Cov
///	Returns: constructor
///	Action : Creates the object
////////////////////////////////////////////////////
ARM_PInfo_Cov::ARM_PInfo_Cov( ARM_PInfo_Map* pinfoMap, const string& itsRefColumn, const ARM_StringVector& itsCVColumns, const std::vector<double>& CVPrices, const std::vector<double>& userBetas )
:
	ARM_PricerInfo(),
	itsCVColumns( ARM_StringVectorPtr( new ARM_StringVector(itsCVColumns) ) ),
	itsRefColumn( itsRefColumn),
	itsCVPrices( ARM_GP_VectorPtr( new ARM_GP_Vector(CVPrices) ) )
{
	ARM_PInfo_MultipleLoop* pinfoMultiLoop = dynamic_cast<ARM_PInfo_MultipleLoop*>( &*((*pinfoMap)[0].second ) );
	if( !pinfoMultiLoop )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": dynamic_cast<ARM_PInfo_MultipleLoop*> failes!"  );

	ComputeBetas(pinfoMap, pinfoMultiLoop, userBetas );
	PostProcessCV(pinfoMap, pinfoMultiLoop, userBetas );
}



////////////////////////////////////////////////////
///	Class  : ARM_PInfo_Cov
///	Routine: ComputeBetas
///	Returns: void 
///	Action : computes if necessary the optimal betas
////////////////////////////////////////////////////
void ARM_PInfo_Cov::ComputeBetas( ARM_PInfo_Map* pinfoMap, ARM_PInfo_MultipleLoop* pinfoMultiLoop, const std::vector<double>& userBetas )
{
	if( userBetas.size() )
	{
#if defined(__GP_STRICT_VALIDATION)
		if( itsCVColumns->size() != userBetas.size() )
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " itsCVColumns->size() != userBetas.size()! " );
#endif
		itsBetas = ARM_GP_VectorPtr( new ARM_GP_Vector( userBetas ) );
	}
	else
	{
		/// 1) compute the varCV cov between control variates and the cov between the ref column and the control variate column
		itsVarCovCV = ARM_GP_MatrixPtr( new ARM_GP_Matrix( itsCVColumns->size(), itsCVColumns->size(), 0.0 ) );
		itsCovRefCV = ARM_GP_VectorPtr( new ARM_GP_Vector( itsCVColumns->size(), 0.0 ) );
		itsCorrRefCV= ARM_GP_VectorPtr( new ARM_GP_Vector( itsCVColumns->size(), 0.0 ) );

		ARM_GP_VectorPtr vecRefValues = pinfoMap->GetContents( itsRefColumn ).GetData( "Values" ).GetVector();
		double stdDevRef =  pinfoMap->GetContents( itsRefColumn ).GetData( "StdDev" ).GetDouble();
		double stdDevI;

		size_t totalSize = vecRefValues->size();
		size_t i,j;
		ARM_GP_VectorPtr vecValuesI, vecValuesJ;

		for( i=0; i<itsCVColumns->size(); ++i )
		{
			stdDevI = pinfoMap->GetContents( (*itsCVColumns)[i] ).GetData( "StdDev" ).GetDouble();
			(*itsVarCovCV)(i,i) = stdDevI*stdDevI;

			/// compute the varCV Covar matrix
			vecValuesI = pinfoMap->GetContents( (*itsCVColumns)[i] ).GetData( "Values" ).GetVector();

			for( j=0; j<i; ++j )
			{
				vecValuesJ = pinfoMap->GetContents( (*itsCVColumns)[j] ).GetData( "Values" ).GetVector();

	#ifdef __GP_STRICT_VALIDATION
				if( vecValuesI->size() != vecValuesJ->size() )
					ARM_THROW( ERR_INVALID_ARGUMENT, " vecValuesI.size() != vecValuesJ.size() !" );
	#endif
				(*itsVarCovCV)(i,j) = (*itsVarCovCV)(j,i) = pinfoMultiLoop->GetMomentFunc()->ComputeCov( vecValuesJ, vecValuesI );
			}

			/// part on the covariance between the ref column and the control variate
			(*itsCovRefCV)[i] = pinfoMultiLoop->GetMomentFunc()->ComputeCov( vecValuesI, vecRefValues );
			(*itsCorrRefCV)[i]= (*itsCovRefCV)[i]/(stdDevI * stdDevRef);
		}

		/// 2) determines the itsBetas of the control variate
		itsBetas = ARM_GP_VectorPtr(static_cast<ARM_GP_Vector*>( itsCovRefCV->Clone() ) );
		LinSolve( &*itsVarCovCV, &*itsBetas );
	}	
}


////////////////////////////////////////////////////
///	Class  : ARM_PInfo_Cov
///	Routine: PostProcessCV
///	Returns: void 
///	Action : given the beta, computes the control variate
////////////////////////////////////////////////////
void ARM_PInfo_Cov::PostProcessCV( ARM_PInfo_Map* pinfoMap, ARM_PInfo_MultipleLoop* pinfoMultiLoop, const std::vector<double>& userBetas ) 
{
	/*itsRunningAvgStep = pinfoMultiLoop->GetRunningAvgStep();
	ARM_GP_VectorPtr vecRefValues = ARM_GP_VectorPtr( static_cast<ARM_GP_Vector*>( pinfoMap->GetContents( itsRefColumn ).GetData( "Values" ).GetVector()->Clone() ) );

	size_t totalSize = vecRefValues->size();
	size_t i,j;

	/// substract the beta to the initial values
	for( j=0; j<itsCVColumns->size(); ++j )
	{
		ARM_GP_VectorPtr vecValuesJ= pinfoMap->GetContents( (*itsCVColumns)[j] ).GetData( "Values" ).GetVector();
		for( i=0; i<totalSize; ++i )
			(*vecRefValues)[i] += (*itsBetas)[j]*((*itsCVPrices)[j]-(*vecValuesJ)[i]);
	}

	itsRunningAverage	= ARM_GP_VectorPtr( new std::vector<double>( int( totalSize/itsRunningAvgStep ), 0.0 ) );
	itsValues			= vecRefValues;
	itsPrice			= 0.0;
	
	pinfoMultiLoop->GetMomentFunc()->init(totalSize,itsRunningAvgStep);
	pinfoMultiLoop->GetMomentFunc()->addBucket(*itsValues);
	pinfoMultiLoop->GetMomentFunc()->finalize();

	/// standard std dev
	itsPrice = pinfoMultiLoop->GetMomentFunc()->GetMean();
	itsStdDev = pinfoMultiLoop->GetMomentFunc()->GetStdDev();
	itsRunningAverage = pinfoMultiLoop->GetMomentFunc()->GetRunningAverage();

	if( userBetas.empty() )
	{
		itsTheoreticalStdDev= pinfoMap->GetContents( itsRefColumn ).GetData( "StdDev" ).GetDouble();
		itsTheoreticalStdDev*= itsTheoreticalStdDev;

		/// stdDev = Var( X - Sum Beta * Y ) = Var(X) + Sum Beta^2*Var(Y) - Sum 2*Beta*Cov(X,Y)
		for( j=0; j<itsCVColumns->size(); ++j )
			itsTheoreticalStdDev+= (*itsBetas)[j]*(*itsBetas)[j]* (*itsVarCovCV)(j,j) - 2.0*(*itsBetas)[j]*(*itsCovRefCV)[j];
		if( itsTheoreticalStdDev< 0.0 )
			itsTheoreticalStdDev = 0.0;
		else
			itsTheoreticalStdDev= sqrt(itsTheoreticalStdDev);
	}
	else itsTheoreticalStdDev=itsStdDev;*/

}
////////////////////////////////////////////////////
///	Class  : ARM_PInfo_Cov
///	Routine: ARM_PInfo_Cov
///	Returns: copy constructor
///	Action : Creates the object
////////////////////////////////////////////////////

ARM_PInfo_Cov::ARM_PInfo_Cov( const ARM_PInfo_Cov& rhs )
:	
	ARM_PricerInfo(rhs),
	itsPrice(rhs.itsPrice),
	itsStdDev(rhs.itsStdDev),
	itsValues(rhs.itsValues),
	itsBetas(rhs.itsBetas),
	itsRunningAverage(rhs.itsRunningAverage),
	itsRunningAvgStep(rhs.itsRunningAvgStep),
	itsVarCovCV(rhs.itsVarCovCV),
	itsCovRefCV(rhs.itsCovRefCV),
	itsCorrRefCV(rhs.itsCorrRefCV),
	itsCVColumns(rhs.itsCVColumns),
	itsRefColumn(rhs.itsRefColumn),
	itsTheoreticalStdDev(rhs.itsTheoreticalStdDev),
	itsCVPrices(rhs.itsCVPrices)
{}


////////////////////////////////////////////////////
///	Class  : ARM_PInfo_Cov
///	Routine: operator=
///	Returns: assignment operator
////////////////////////////////////////////////////

ARM_PInfo_Cov& ARM_PInfo_Cov::operator=(const ARM_PInfo_Cov & rhs )
{
	if( this != & rhs )
	{
		ARM_PricerInfo::operator=(rhs),
		itsPrice			= rhs.itsPrice;
		itsStdDev			= rhs.itsStdDev;
		itsValues			= rhs.itsValues;
		itsBetas				= rhs.itsBetas;
		itsRunningAverage	= rhs.itsRunningAverage;
		itsRunningAvgStep	= rhs.itsRunningAvgStep;
		itsVarCovCV			= rhs.itsVarCovCV;
		itsCovRefCV			= rhs.itsCovRefCV;
		itsCorrRefCV		= rhs.itsCorrRefCV;
		itsCVColumns			= rhs.itsCVColumns;
		itsRefColumn		= rhs.itsRefColumn;
		itsTheoreticalStdDev= rhs.itsTheoreticalStdDev;
		itsCVPrices		= rhs.itsCVPrices;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_PInfo_Cov
///	Routine: StoreCurrentLoopPrice,PostProcess
///	Returns: throw exception
////////////////////////////////////////////////////

void ARM_PInfo_Cov::StoreCurrentLoopPrice( const ARM_PricingStatesPtr& states, 
		size_t bucketSize, size_t currentBucketIndex, ARM_PricingAdviser::PricingDir pricingDir, int colIdx )
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": make no sense to call StoreCurrentLoopPrice!"  );
}

void ARM_PInfo_Cov::PostProcess(const ARM_TimeInfoPtrVector& timeInfo,const ARM_NumMethodPtr& numMethod)
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": make no sense to call PostProcess!"  );
}



////////////////////////////////////////////////////
///	Class  : ARM_PInfo_Cov
///	Routine: GetContents
///	Returns: 
///	Action : Creates the corresponding dicitionary
////////////////////////////////////////////////////

ARM_GramFctorArgDict ARM_PInfo_Cov::GetContents(const string& colName ) const
{
	ARM_GramFctorArgDict contents;
	contents[ "Price"			] = itsPrice;
	contents[ "StdDev"			] = itsStdDev;
	contents[ "TheoStdDev"		] = itsTheoreticalStdDev;
	contents[ "Values"			] = itsValues;
	contents[ "Beta"			] = itsBetas;
	contents[ "RunningAvg"		] = itsRunningAverage;
	contents[ "RunningAvgStep"	] = itsRunningAvgStep;
	contents[ "CovRefCV"		] = itsCovRefCV;
	contents[ "CorrRefCV"		] = itsCorrRefCV;
	contents[ "VarCovCV"		] = itsVarCovCV;
	contents[ "itsCVColumns"		] = itsCVColumns;
	contents[ "RefColumn"		] = itsRefColumn;
	contents[ "Duration"		] = GetDuration();
	contents[ "RefPrices"		] = itsCVPrices;
	return contents;	
}

////////////////////////////////////////////////////
////////////////////////////////////////////////////
////		ARM_PInfo_Map						////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_PInfo_Map
///	Routine: Constructor,Copy constructor, assignment operator
///				Destructor and Clone
///	Returns: 
///	Action : builds the object
////////////////////////////////////////////////////
ARM_PInfo_Map::ARM_PInfo_Map(const vector<string>& columnNames, const vector<ARM_PricerInfo*>& pricerInfos)
: ARM_PricerInfo()
{
	for (size_t i = 0; i < columnNames.size(); ++i)
		itsPricerInfoMap.push_back(make_pair(columnNames[i], ARM_PricerInfoPtr(pricerInfos[i])));
}

ARM_PInfo_Map::ARM_PInfo_Map(const ARM_PInfo_Map& rhs)
: ARM_PricerInfo(),
itsPricerInfoMap(rhs.itsPricerInfoMap)
{
}

ARM_PInfo_Map& ARM_PInfo_Map::operator=(const ARM_PInfo_Map& rhs)
{
	if (&rhs != this)
	{
		ARM_PInfo_Map::operator =(rhs);
		itsPricerInfoMap = rhs.itsPricerInfoMap;
	}

	return *this;
}

ARM_Object* ARM_PInfo_Map::Clone() const
{
	return new ARM_PInfo_Map(*this);
}

////////////////////////////////////////////////////
///	Class  : ARM_PInfo_Map
///	Routine: reset
///	Returns: 
///	Action : reset the pricer info
////////////////////////////////////////////////////
void ARM_PInfo_Map::reset()
{
	for (PricerInfoMap::iterator it = itsPricerInfoMap.begin(); it != itsPricerInfoMap.end(); ++it)
		(*it).second->reset();
}


////////////////////////////////////////////////////
///	Class  : ARM_PInfo_Map
///	Routine: StoreCurrentLoopPrice,PostProcess,GetPrice,toString
///	Returns: 
///	Action : store information about current loop price,
///		postprocess, get price and show the result as a string (toString)
////////////////////////////////////////////////////
void ARM_PInfo_Map::StoreCurrentLoopPrice( const ARM_PricingStatesPtr& states, size_t bucketSize, size_t currentBucketIndex, ARM_PricingAdviser::PricingDir pricingDir, int colIdx )
{
	for (size_t i = 0; i < itsPricerInfoMap.size(); i++)
		itsPricerInfoMap[i].second->StoreCurrentLoopPrice(states,bucketSize,currentBucketIndex,pricingDir, i);
}


void ARM_PInfo_Map::PostProcess(const ARM_TimeInfoPtrVector& timeInfo,const ARM_NumMethodPtr& numMethod)
{
	for (size_t i = 0; i < itsPricerInfoMap.size(); i++)
		itsPricerInfoMap[i].second->PostProcess(timeInfo, numMethod);

}


double ARM_PInfo_Map::GetPrice() const
{	
	for (size_t i = 0; i < itsPricerInfoMap.size(); i++)
		itsPricerInfoMap[i].second->SetDuration( GetDuration() );
	return (*itsPricerInfoMap.begin()).second->GetPrice(); 
}

////////////////////////////////////////////////////
///	Class  : ARM_PInfo_SingleLoop
///	Routine: GetContents
///	Returns: 
///	Action : Creates the corresponding dicitionary
////////////////////////////////////////////////////
ARM_GramFctorArgDict ARM_PInfo_Map::GetContents(const string& colName) const
{
	/// default behavior is to return the first pricer if no colName given
	if( "" == colName )
		return itsPricerInfoMap[0].second->GetContents();

	for (size_t i = 0; i < itsPricerInfoMap.size(); i++)
	{
		if ( stringGetUpper(itsPricerInfoMap[i].first) == stringGetUpper(colName) ) 
			return itsPricerInfoMap[i].second->GetContents();
	}

	throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			ARM_USERNAME + ": The column" + colName + "has not been priced!"  );
}


////////////////////////////////////////////////////
///	Class  : ARM_PInfo_Map
///	Routine: toString
///	Returns: 
///	Action : returns the price and does a toString
////////////////////////////////////////////////////
string ARM_PInfo_Map::toString(const string& indent, const string& nextIndent ) const
{
	CC_Ostringstream os;

	for (size_t i = 0; i < itsPricerInfoMap.size(); i++)
	{
		os << itsPricerInfoMap[i].first << " pricer info\n";
		os << itsPricerInfoMap[i].second->GetContents().toString();
	}
	return os.str();
}

CC_END_NAMESPACE()

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/


