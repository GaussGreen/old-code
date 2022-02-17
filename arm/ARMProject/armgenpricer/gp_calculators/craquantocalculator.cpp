/*!
 /** Copyright (c) CDC IXIS CM July 2003 Paris
 *	\file CraQuantoCalculator.cpp
 *  \brief file for the CRA quanto
 *	\author  J.M
 *	\version 1.0
 *	\date Apr 2007
 */

#include "gpcalculators/craquantocalculator.h"

/// gpbase
#include "gpbase/env.h"
#include "gpbase/typedef.h"
#include "gpbase/argconvdefault.h"
#include "gpbase/datestripcombiner.h"
#include "gpbase/datestrip.h"
#include "gpbase/stringconvert.h"
#include "gpbase/curveconvert.h"

/// gpinfra
#include "gpinfra/argconvdefault.h"
#include "gpinfra/cstmanager.h"
#include "gpinfra/curvemodelparam.h"

/// gpmodels
#include "gpmodels/HWHWQtoModel.h"

/// gpcalib
#include "gpcalib/vanillaswaption.h"
#include "gpcalib/calibmethod.h"

/// kernel
//#include <inst/swaption.h>





CC_BEGIN_NAMESPACE( ARM )

const double SIGMA_LOWER_BOUND		= 0.0005;
const double SIGMA_UPPER_BOUND		= 0.05;
const double SIGMA_DEFAULT_VALUE	= 0.005;
const double NON_CALL_FEE			= 1e15;

const double DEFAULT_PRICE			= 1.0e+100;
const double DEFAULT_PRECISION		= 1.0;
const double DEFAULT_WEIGHT			= 1.0;


const string ARM_CraQuantoCalculator::CraQuantoColNamesTable [] =
{  
	"Funding",
	"Corridor",
	"Swap"
};

const int ARM_CraQuantoCalculator::CraQuantoProductToPriceColumns [] =
{  
	Funding,
	Corridor,
	Swap
};

//////////////////////////////////////////////////////////////////////////////
/// Class  : ARM_CraQuantoCalculator
/// Routines: 
/// Returns :
/// Action  : standard constructor
//////////////////////////////////////////////////////////////////////////////

ARM_CraQuantoCalculator::ARM_CraQuantoCalculator(	const ARM_Currency& ccyDom,
								const ARM_Currency& ccyFor,
								ARM_Date&	startDate,
								ARM_Date&	endDate,
								int payRec,
								const ARM_ReferenceValue&  notional,
								int callFreq,
								int callNotice,
								const string& callCal,
								const ARM_ReferenceValue&  callFees,
								int fundFreq,
								int fundDayCount,
								const ARM_ReferenceValue& fundSpread,
								int cpnPayFreq,
								int cpnDayCount,
								const string& cpnResetCal,
								const string& cpnPayCal,
								int payIndex,
								int payResetTiming,
								const ARM_ReferenceValue& fixRate,
								const ARM_ReferenceValue& barrierDown,
								const ARM_ReferenceValue& barrierUp,
								int refResetFreq,
								int refResetTiming,
								int refIndex,
								const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice)
:	ARM_CallableQuantoCalculator(	ccyDom,
									ccyFor,
									startDate,
									endDate,
									payRec,
									notional,
									callFreq,
									callNotice,
									callCal,
									callFees,
									fundFreq,
									fundDayCount,
									fundSpread,
									cpnPayFreq,
									cpnDayCount,
									cpnResetCal,
									cpnPayCal,
									productsToPrice),
	itsPayIndex(payIndex),
	itsPayResetTiming(payResetTiming),
	itsFixRate(fixRate),
	itsBarrierDown(barrierDown),
	itsBarrierUp(barrierUp),
	itsRefResetFreq(refResetFreq),
	itsRefResetTiming(refResetTiming),
	itsRefIndex(refIndex)
{

	ARM_StringVector pricedColumns = PricedColumnNames();

	ARM_CstManagerPtr cstManagerPtr = CreateCstManager();

	CreateAndSetDealDescriptionAndTimeIt(GetKeys()[YcDomKey], pricedColumns, cstManagerPtr);
}

//////////////////////////////////////////////////////////////////////////////
/// Class  : ARM_CraQuantoCalculator
/// Routines: 
/// Returns :
/// Action  : copy constructor
//////////////////////////////////////////////////////////////////////////////
ARM_CraQuantoCalculator::ARM_CraQuantoCalculator(const ARM_CraQuantoCalculator& rhs)
: ARM_CallableQuantoCalculator(rhs)
{
	CopyNoCleanUp(rhs);
}
	
//////////////////////////////////////////////////////////////////////////////
/// Class  : ARM_CraQuantoCalculator
/// Routine: CopyNoCleanUp
/// Returns: void
/// Action : Arguments copy
//////////////////////////////////////////////////////////////////////////////
void ARM_CraQuantoCalculator::CopyNoCleanUp(const ARM_CraQuantoCalculator& rhs)
{
	itsPayIndex = rhs.itsPayIndex;
	itsPayResetTiming = rhs.itsPayResetTiming;
	itsFixRate = rhs.itsFixRate;
	itsBarrierDown = rhs.itsBarrierDown;
	itsBarrierUp = rhs.itsBarrierUp;
	itsRefResetFreq = rhs.itsRefResetFreq;
	itsRefResetTiming = rhs.itsRefResetTiming;
	itsRefIndex = rhs.itsRefIndex;
}

//////////////////////////////////////////////////////////////////////////////
/// Class  : ARM_CraQuantoCalculator
/// Routines: Clone
/// Returns :
/// Action  : Standard ARM object support
///           Call copy constructor
//////////////////////////////////////////////////////////////////////////////
ARM_Object*	ARM_CraQuantoCalculator::Clone() const
{
	return(new ARM_CraQuantoCalculator(*this));
}
	
//////////////////////////////////////////////////////////////////////////////
/// Class  : ARM_CraQuantoCalculator
/// Routines: operator=
/// Returns :
/// Action  : Call copy constructor
//////////////////////////////////////////////////////////////////////////////
ARM_CraQuantoCalculator& ARM_CraQuantoCalculator::operator=(const ARM_CraQuantoCalculator& rhs)
{
	if( this != & rhs )
	{
		ARM_CallableQuantoCalculator::operator=(rhs);
        CleanUp();
        CopyNoCleanUp(rhs);
	}
	return *this;
}

//////////////////////////////////////////////////////////////////////////////
/// Class  : ARM_CraQuantoCalculator
/// Routine: CleanUp
/// Returns: void
/// Action : delete pointors
//////////////////////////////////////////////////////////////////////////////
void ARM_CraQuantoCalculator::CleanUp()
{
}

//////////////////////////////////////////////////////////////////////////////
/// Class  : ARM_CraQuantoCalculator
/// Routine: Destructor
/// Returns: void
/// Action : destroys the object
//////////////////////////////////////////////////////////////////////////////
ARM_CraQuantoCalculator::~ARM_CraQuantoCalculator()
{
	CleanUp();
}

//////////////////////////////////////////////////////////////////////////////
/// Class  : ARM_CraQuantoCalculator
/// Routine: ExtraDescSize
/// Returns: void
/// Action : returns extra desc size
//////////////////////////////////////////////////////////////////////////////

size_t ARM_CraQuantoCalculator::ExtraDescSize() const
{
	return sizeof(CraQuantoColNamesTable)/sizeof(CraQuantoColNamesTable[0]);
};

/////////////////////////////////////////////////////////////////////////////
/// Class  : ARM_CraQuantoCalculator
/// Routine: ExtraNbProductsToPrice
/// Returns: void
/// Action : returns extra desc size
//////////////////////////////////////////////////////////////////////////////

size_t ARM_CraQuantoCalculator::ExtraNbProductsToPrice() const
{
	return NbProductsToPrice;
};


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CraQuantoCalculator
///	Routine: CreateCstManager
///	Returns: void
///	Action : create the const manager (static data).
/////////////////////////////////////////////////////////////////
ARM_CstManagerPtr ARM_CraQuantoCalculator::CreateCstManager()
{
	try
	{
		double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

		vector<string>	cstNames; 
		cstNames.push_back("FundSpreadCst");
		cstNames.push_back("NotioCst");
		cstNames.push_back("FixRateCst");
		cstNames.push_back("BarrierDownCst");
		cstNames.push_back("BarrierUpCst");
		

		vector<ARM_GramFctorArg> cstVector;
		
		//FundSpreadCst
		ARM_Curve* tmpFundSpread = RefValueToCurve(GetFundSpread(), asOfDate, ARM_Constants::rateBase);
		cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpFundSpread))));
		
		//NotionalCst
		ARM_Curve* tmpNotional = RefValueToCurve(GetNotional(), asOfDate);
		cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpNotional))));

		//FixRateCst
		ARM_Curve* tmpFix;
		tmpFix = RefValueToCurve(GetFixRate(), asOfDate, ARM_Constants::rateBase);
		cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpFix))));
	
		//BarrierDownCst
		ARM_Curve* tmpBarrierDown = RefValueToCurve(GetBarrierDown(), asOfDate, ARM_Constants::rateBase);
		cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpBarrierDown))));

		//BarrierUpCst
		ARM_Curve* tmpBarrierUp = RefValueToCurve(GetBarrierUp(), asOfDate, ARM_Constants::rateBase);
		cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpBarrierUp))));

		ARM_CstManagerPtr cstManagerPtr = ARM_CstManagerPtr(new ARM_CstManager(cstNames, cstVector));
		return cstManagerPtr;
	}
	catch (Exception& x)
	{
	    x.DebugPrint();

		throw x;
	}
	catch (...)
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in CraCalculator::CreateCstManager" );
	}
}

//////////////////////////////////////////////////////////////////////////////
/// Class  : ARM_CraQuantoCalculator
/// Routine: PricedColumnNames
/// Returns: ARM_StringVector
/// Action : cf name
//////////////////////////////////////////////////////////////////////////////
ARM_StringVector ARM_CraQuantoCalculator::PricedColumnNames() const
{
	ARM_StringVector pricedColumns;

	int colSize = MIN(NbProductsToPrice+1, GetProductsToPrice().size());

	pricedColumns.push_back(
		ARM_CallableQuantoCalculator::CallableQuantoColNamesTable[
			ARM_CallableQuantoCalculator::CallableQuantoProductToPriceColumns[0] ]);

	for (int i=1; i<colSize; i++)
	{
		if (GetProductsToPrice()[i])
			pricedColumns.push_back(CraQuantoColNamesTable[CraQuantoProductToPriceColumns[i-1]]);
	}
	return pricedColumns;
}

//////////////////////////////////////////////////////////////////////////////
/// Class  : ARM_CraQuantoCalculator
/// Routine: UpdateExtraDesc
/// Returns: void
/// Action : end of middlerows
//////////////////////////////////////////////////////////////////////////////
void ARM_CraQuantoCalculator::UpdateExtraColNames(vector<string>& colNamesVec) const
{
	size_t offsetCol = ARM_CallableQuantoCalculator::CallableQuantoColAlias::NbCols;

	for (size_t i = 0; i < ExtraDescSize(); i++)
	{
		colNamesVec[offsetCol + i] = CraQuantoColNamesTable[i];
	}
}

//////////////////////////////////////////////////////////////////////////////
/// Class  : ARM_CraQuantoCalculator
/// Routine: UpdateExtraDesc
/// Returns: void
/// Action : end of middlerows
//////////////////////////////////////////////////////////////////////////////
string ARM_CraQuantoCalculator::ColNames(size_t i) const
{
	size_t offsetCol = CallableQuantoProductToPriceAlias::NbProductsToPrice;
	if (i<offsetCol)
		return ARM_CallableQuantoCalculator::CallableQuantoColNamesTable[ ARM_CallableQuantoCalculator::CallableQuantoProductToPriceColumns[i] ] ;
	else
		return CraQuantoColNamesTable[CraQuantoProductToPriceColumns[i-offsetCol]];
}


//////////////////////////////////////////////////////////////////////////////
/// Class  : ARM_CraQuantoCalculator
/// Routine: UpdateExtraDesc
/// Returns: void
/// Action : end of middlerows
//////////////////////////////////////////////////////////////////////////////
void ARM_CraQuantoCalculator::UpdateExtraDesc(size_t eventIdx, const ARM_DateStripCombiner& datesStructure, vector<string>& rowDescVec, vector<ARM_GP_VALUE_TYPE>& rowTypeVec) const
{
try
	{
		double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
		
		size_t EXER_SCHED = ARM_CallableQuantoCalculator::ScheduleCount::ExerSched;

		size_t eventSize = datesStructure.GetDateStrip(EXER_SCHED)->GetResetDates()->size();
		std::vector<double>& exerciseDates = datesStructure.GetDateStrip(EXER_SCHED)->GetResetDates();
		size_t nbPastNoCall=0;
		while(nbPastNoCall < eventSize && (*exerciseDates)[nbPastNoCall] < asOfDate)
			++nbPastNoCall;

		bool isFirstEvent = (eventIdx==nbPastNoCall);
		bool isLastEvent  = (eventIdx==(eventSize-1));

		size_t offsetCol = ARM_CallableQuantoCalculator::CallableQuantoColAlias::NbCols;

		string payRec = ( GetPayRec() == -1 )? "P" : "R";

		string fundFreq		= ARM_ArgConvReverse_MatFrequency.GetString(GetFundFreq());
		string fundDayCount = ARM_ArgConvReverse_DayCount.GetString(GetFundDayCount());
		
		//FUNDING
		CC_Ostringstream fundingDesc;
		fundingDesc << "SWAP(" << GetKeys()[YcDomKey] << ", StartDate[i], EndDate[i], 0," << payRec << ", ";
		fundingDesc << fundFreq << ", " << fundDayCount << ", " << fundFreq << ", " << fundDayCount << ", ";
		fundingDesc << "FundSpreadCst, NotioCst" << ")";
		rowDescVec[offsetCol + Funding] = fundingDesc.str();
		rowTypeVec[offsetCol + Funding] = ARM_STRING;

		string cpnPayFreq	= ARM_ArgConvReverse_StdFrequency.GetString(GetCpnPayFreq());
		string cpnDayCount	= ARM_ArgConvReverse_DayCount.GetString(GetCpnDayCount());
		string refResetFreq	= ARM_ArgConvReverse_StdFrequency.GetString(GetRefResetFreq());

		int refIndexTypeInt;
		string refIndexTerm = FromIndexTypeToTermAndType(GetRefIndex(),refIndexTypeInt);
		string refIndexType = ARM_ArgConvReverse_IndexClass.GetString(refIndexTypeInt);
		int payIndexTypeInt;
		string payIndexTerm = FromIndexTypeToTermAndType(GetPayIndex(),payIndexTypeInt);
		string payIndexType = ARM_ArgConvReverse_IndexClass.GetString(payIndexTypeInt);

		string payResetTiming = ARM_ArgConvReverse_Timing.GetString(GetPayResetTiming());
		string refResetTiming = ARM_ArgConvReverse_Timing.GetString(GetRefResetTiming());
		
		//CORRIDOR
		CC_Ostringstream corridorDesc;
		corridorDesc << "CORRIDOR(" << LOCAL_MODEL_NAME << ",";		// Model
		corridorDesc << "StartDate[i]" << ",";				// StartDate
		corridorDesc << "EndDate[i]" << ",";				// EndDate
		corridorDesc << payRec << ",";						// PayRec
		corridorDesc << payIndexType << ",";				// PayIndex
		corridorDesc << "FixRateCst" << ",";				// Fix Value
		corridorDesc << 1 << ",";							// PayIndexMult Value
		corridorDesc << cpnPayFreq << ",";					// Payment Freq
		corridorDesc << payIndexTerm << ",";				// Pay Index Term
		corridorDesc << payResetTiming << ",";				// Pay Index Timing
		corridorDesc << cpnDayCount << ",";					// Pay Day Count
		corridorDesc << ",";								// Pay Gap
		corridorDesc << ",";								// Int Rule
		corridorDesc << 0 << ",";							// Spread
		corridorDesc << refResetFreq << ",";				// Fixing Freq
		corridorDesc << refResetTiming << ",";				// Fixing Timing
		corridorDesc << refIndexType << ",";				// CMS, LIBOR
		corridorDesc << refIndexTerm << ",";				// 5y, 10y, 3m...
		corridorDesc << 1 << ",";							// Coeff 1
		corridorDesc << "BarrierDownCst" << ",";			// Barrier Down
		corridorDesc << "BarrierUpCst" << ",";				// Barrier Up
		corridorDesc << "NotioCst" << ",";					// Notional
		corridorDesc << "CMS" << ",";									// Fake : simple corridor not ready yet!!!
		corridorDesc << "2Y" << ",";									// Fake
		corridorDesc << "0)";											// Fake
			
		rowDescVec[offsetCol + Corridor] = corridorDesc.str();
		rowTypeVec[offsetCol + Corridor] = ARM_STRING;

		//SWAP
		CC_Ostringstream swapDesc;
		swapDesc << "Funding[i] + Corridor[i]";
		rowDescVec[offsetCol + Swap] = swapDesc.str();
		rowTypeVec[offsetCol + Swap] = ARM_STRING;


	}

	catch (Exception& x)
	{
	    x.DebugPrint();

		throw x;
	}
	catch (...)
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in CraQuantoCalculator::UpdateExtraDesc" );
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CraQuantoCalculator
///	Routine: Calibrate
///	Returns: void
///	Action : calibrate the model 
/////////////////////////////////////////////////////////////////
void ARM_CraQuantoCalculator::Calibrate()
{
	GetCalibMethod()->Calibrate(&*GetPricingModel());
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CraQuantoCalculator
///	Routine: CreateAndSetCalibration
///	Returns: void
///	Action : create the calibration methods
/////////////////////////////////////////////////////////////////
void ARM_CraQuantoCalculator::CreateAndSetCalibration()
{
	
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CraQuantoCalculator
///	Routine: UpdateCalibration
///	Returns: void
///	Action : update the calibration datas w.r.t. MarketDataManager
/////////////////////////////////////////////////////////////////
void ARM_CraQuantoCalculator::UpdateCalibration(bool isUpdateStrike)
{
	CreateEmptyCalibration();

	UpdateSwaptionPrices();
}
	

//////////////////////////////////////////////////////////////////////////////
/// Class  : ARM_CraQuantoCalculator
/// Routine: CreateEmptyCalibration
/// Returns: void
/// Action : 
//////////////////////////////////////////////////////////////////////////////
void ARM_CraQuantoCalculator::CreateEmptyCalibration()
{
	std::vector<double> calibTimes(1, 0.0);
	std::vector<double> volatility(1, SIGMA_DEFAULT_VALUE);
	std::vector<double> lowerVolatility(1, SIGMA_LOWER_BOUND);
	std::vector<double> upperVolatility(1, SIGMA_UPPER_BOUND);
	ARM_CurveModelParam* volParam = new ARM_CurveModelParam(ARM_ModelParamType::Volatility, &volatility, &calibTimes ,"","STEPUPRIGHT",&lowerVolatility,&upperVolatility);
	ARM_ModelParamVector volParamVec(1,volParam);

	ARM_StdPortfolioPtr pfDomSwopt = CreateSwaptionPortfolio(GetDomesticCcy()); 
	ARM_StdPortfolioPtr pfForSwopt = CreateSwaptionPortfolio(GetForeignCcy()); 
	
	ARM_CalibMethod calibVolMethodDom(
						pfDomSwopt,
						volParamVec,
						ARM_CalibMethodType::Bootstrap1D,
						ARM_NumericConstants::ARM_GP_MAX_ITER,
						ARM_CalibrationTarget::PriceTarget,
						NULL,
						NULL,
						false,
						ARM_HWHWQtoModel::DomModel,
						1,
						true);

	ARM_CalibMethodPtr calibVolMethodFor = ARM_CalibMethodPtr(
		new ARM_CalibMethod(
						pfDomSwopt,
						volParamVec,
						ARM_CalibMethodType::Bootstrap1D,
						ARM_NumericConstants::ARM_GP_MAX_ITER,
						ARM_CalibrationTarget::PriceTarget,
						NULL,
						&calibVolMethodDom,
						false,
						ARM_HWHWQtoModel::ForModel,
						1,
						false)	// to prevent validate
						);

	SetCalibMethod(calibVolMethodFor);
	delete volParam;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CraQuantoCalculator
///	Routine: CreateSwaptionPortfolio
///	Returns: a portfolio
///	Action : create the list of diagonal swaptions
/////////////////////////////////////////////////////////////////
ARM_StdPortfolioPtr ARM_CraQuantoCalculator::CreateSwaptionPortfolio(ARM_Currency& ccy)
{
	ARM_DateStripCombiner datesStructure	= DatesStructure();
	size_t EXER_SCHED = ARM_CallableQuantoCalculator::ScheduleCount::ExerSched;

	std::vector<double>& resetDates = datesStructure.GetDateStrip(EXER_SCHED)->GetResetDates();
	std::vector<double>& startDates = datesStructure.GetDateStrip(EXER_SCHED)->GetFlowStartDates();

	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

	int nbFlows = startDates->size();

	list< ARM_Security* > swaptionList;

	ARM_Date endDate = GetEndDate();

	ARM_INDEX_TYPE indexType = ccy.GetVanillaIndexType();
	
	for (int i = 0; i < nbFlows; ++i)
	{
		ARM_Date startDate ((*startDates)[i]);

		double fees = itsCallFees.Interpolate((*resetDates)[i]);

		if (((*resetDates)[i] > asOfDate) && (fees < NON_CALL_FEE))
		{
			ARM_Swap stdSwap(
							startDate,
							endDate,
							indexType,
							0.0,
							K_MARKET_RATE,
							K_RCV,
							K_DEF_FREQ,
							K_DEF_FREQ,
							&ccy);

			ARM_Date expiryDate((*(stdSwap.GetFloatLeg()->GetResetDates()))[0]);
			ARM_Swaption* swaption = new ARM_Swaption(&stdSwap,K_RCV,K_EUROPEAN,K_MARKET_RATE,expiryDate);

			swaptionList.push_back(static_cast< ARM_Security* >(swaption));
		}
	}

	ARM_StdPortfolioPtr swaptionPF(new ARM_StdPortfolio(swaptionList));
    	
	for(i=0;i<swaptionPF->size();++i)
	{
		swaptionPF->SetPrecision(DEFAULT_PRECISION,i);
        swaptionPF->SetWeight(DEFAULT_WEIGHT,i);
		swaptionPF->SetPrice(DEFAULT_PRICE,i);
	}
	return swaptionPF;
}

////////////////////////////////////////////////////////////////
///	Class  : ARM_CraQuantoCalculator
///	Routine: UpdateSwaptionPrices
///	Returns: void
///	Action : 
/////////////////////////////////////////////////////////////////
void ARM_CraQuantoCalculator::UpdateSwaptionPrices()
{
    double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	
	ARM_StdPortfolioPtr pfFor = GetCalibMethod()->GetPortfolio();
	ARM_StdPortfolioPtr pfDom = GetCalibMethod()->GetPreviousMethod()->GetPortfolio();

	ARM_BSModel* bsModelFor = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswForModelKey]) );
	ARM_BSModel* bsModelDom = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswDomModelKey]) );
	
	ComputeSwaptionPrices(pfFor,bsModelFor);
	ComputeSwaptionPrices(pfDom,bsModelDom);
}

////////////////////////////////////////////////////////////////
///	Class  : ARM_CraQuantoCalculator
///	Routine: ComputeSwaptionPrices
///	Returns: void
///	Action : 
/////////////////////////////////////////////////////////////////
void ARM_CraQuantoCalculator::ComputeSwaptionPrices(ARM_StdPortfolioPtr& pf, ARM_BSModel* bsModel)
{
	for(int i=0;i<pf->GetSize();++i)
	{
		ARM_Swaption* swaption = static_cast< ARM_Swaption* >(pf->GetAsset(i));
		swaption->SetModel(bsModel);
		double price = swaption->ComputePrice();
		pf->SetPrice (price, i);
	}
}


CC_END_NAMESPACE()