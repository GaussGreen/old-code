/*!
 *
 * Copyright (c) IXIS CIB July 2007 Paris
 *
 *	\file fxsnowballkocalculator.cpp
 *
 *  \brief file for the Fx Snow Ball Ko
 *	\author  N. Belgrade
 *	\version 1.0
 *	\date July 2007
 */


/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalculators/fxsnowballkocalculator.h"
#include "gpcalculators/basisconverter.h"

/// gpbase
#include "gpbase/env.h"
#include "gpbase/argconvdefault.h"
#include "gpbase/autocleaner.h"
#include "gpbase/ostringstream.h"
#include "gpbase/datestripcombiner.h"
#include "gpbase/datestrip.h"
#include "gpbase/singleton.h"
#include "gpbase/utilityport.h"  
#include "gpbase/curveconvert.h"  
#include "gpbase/gplinalgconvert.h"
#include "gpbase/datestripconvert.h"
#include "gpbase/globalconstant.h"

/// gpnumlib
#include "gpnumlib/randomgenfactory.h"
#include "gpnumlib/transposer.h"

/// kernel
#include <crv/volflat.h>
#include <inst/forex.h>
#include <inst/swaption.h>
#include <inst/option.h>
#include <inst/portfolio.h>
#include <inst/swap.h>
#include <util/fromto.h>

/// STL
#include <iomanip> /// for setprecision()
#include <list>
CC_USING_NS(std,list)


CC_BEGIN_NAMESPACE( ARM )

/// Default MDM key names

const string YC_KEY_NAME		= "YC_";
const string YC_BASIS_KEY_NAME	= "YC_BASIS_";
const string FOREX_KEY_NAME		= "FOREX_";
const string OSWMODEL_KEY_NAME	= "OSWMOD_";
const string FXMODEL_KEY_NAME	= "FXMOD_";
const string CORREL_KEY_NAME	= "CORREL_";
const string MRS_KEY_NAME		= "MRS_";
const string QFX_KEY_NAME		= "Q_";


const string ARM_FxSnowBallKoCalculator::fxSbKoColNamesTable[] =
{
		"ResetDate",
		"NextResetDate",
		"StartDate",
		"EndDate",
		"FundingStartDate",
		"FundingEndDate",
		"IT",
		"Notional",
		"Strike",
		"Barrier",
		"Leverage",
		"Coeff",
		"FixedFunding",
		"Offset",
		"Width",
		"FxSpot",
		"IsAlive",
		"Option",
		"Coupon",
		"RealCoupon",
		"Funding",
		"DiscountFunding",
		"RealFunding",
		"DF",
		"Proba",
		"FXSBKO"
};

const int ARM_FxSnowBallKoCalculator::ProductToPriceColumns[] = 
{
	ARM_FxSnowBallKoCalculator::SpotFwdPrice,
	ARM_FxSnowBallKoCalculator::OptionPrice,
	ARM_FxSnowBallKoCalculator::IsAliveValue,
	ARM_FxSnowBallKoCalculator::ProbaValue,
	ARM_FxSnowBallKoCalculator::CouponPrice,
	ARM_FxSnowBallKoCalculator::RealCouponPrice,
	ARM_FxSnowBallKoCalculator::FundingPrice,
	ARM_FxSnowBallKoCalculator::RealFundingPrice,
	ARM_FxSnowBallKoCalculator::FXSBKOPrice,
	ARM_FxSnowBallKoCalculator::NbProductsToPrice
};

/////////////////////////////////////////////////////////////////
///	Class  : ARM_FxSnowBallKoCalculator
///	Routine: Copy constructor
///	Returns: a built object
///	Action : builds the object
/////////////////////////////////////////////////////////////////

ARM_FxSnowBallKoCalculator::ARM_FxSnowBallKoCalculator( const ARM_FxSnowBallKoCalculator& rhs )
: ARM_HybridIRFXCalculator(rhs)
{
	itsStrikeCv		= rhs.itsStrikeCv;
	itsBarrierCv	= rhs.itsBarrierCv;
	itsLeverageCv	= rhs.itsLeverageCv;
	itsCoeffCv		= rhs.itsCoeffCv;
}


ARM_FxSnowBallKoCalculator::ARM_FxSnowBallKoCalculator(const ARM_Date& asOfDate,
														 const ARM_Date& startDate,
														 const ARM_Date& endDate,
														 const ARM_Currency& DomCcy,
														 const vector<ARM_Currency>& ForCcy,
														 const ARM_Currency& CpnCcy,
														 const ARM_Currency* FundCcy,
														 int payRec,
														 int cpnDayCount,
														 int cpnFreq,
														 int cpnResetGap,
														 const string& cpnResetCal,
														 const string& cpnPayCal,
														 int stubRule,
														 int cpnTiming,
														 int cpnIntRule,
														 const ARM_Curve& cpnNominal,
														 int fundFreq,
														 int fundDayCount,
														 const ARM_Curve& fundNominal,
														 const ARM_Curve& fundSpread,
														 const ARM_Curve& strike,
														 const ARM_Curve& barrier,
														 const ARM_Curve& leverage,
														 const ARM_Curve& coeff,
														 int intermediatePrices,
														 const ARM_StringVector& columnsToPrice,
														 int	optionType,
														 bool isPerform,
														 ARM_FixingSched* pastFixings)
: ARM_HybridIRFXCalculator( asOfDate,
							startDate,
							endDate,
							DomCcy,
							ForCcy,
							CpnCcy,
							FundCcy,
							payRec,
							cpnDayCount,
							cpnFreq,
							cpnResetGap,
							cpnResetCal,
							cpnPayCal,
							stubRule,
							cpnTiming,
							cpnIntRule,
							&cpnNominal,
							NULL,
							NULL,
							NULL,
							NULL,
							NULL,
							fundFreq,
							fundDayCount,
							&fundNominal,
							&fundSpread,
							0,
							0,
							0,
							NULL,
							NULL,
							"",
							intermediatePrices,
							columnsToPrice,
							ARM_TARNFXPayoffType::SWITCHER,
							pastFixings ),
itsStrikeCv(strike), itsBarrierCv(barrier), itsLeverageCv(leverage), itsCoeffCv(coeff), 
itsOptionType(optionType), itsIsPerform(isPerform)
{
	/// Create the Generic Security
	CreateAndSetDealDescription(itsPayModelName, itsColumnsToPrice, ARM_CstManagerPtr(), false);
}
/*
/////////////////////////////////////////////////////////////////
///	Class  : ARM_FxSnowBallKoCalculator
///	Routine: CheckData
///	Returns: void
///	Action : check if FX SB KO data are consistent
/////////////////////////////////////////////////////////////////

void ARM_FxSnowBallKoCalculator::CheckData()
{
	if (itsStartDate >= itsEndDate )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The start date frequency should be before the end date.");
} 

/////////////////////////////////////////////////////////////////
///	Class  : ARM_FxSnowBallKoCalculator
///	Routine: CheckMktData
///	Returns: void
///	Action : check if FX SB KO data are consistent
/////////////////////////////////////////////////////////////////
void ARM_FxSnowBallKoCalculator::CheckMktData()
{
	// Dom Zc and Basis curves
	string YcDomKey = YC_KEY_NAME + itsDomesticCcy.GetCcyName();
	ARM_ZeroCurve* domCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcDomKey));
	if(!domCurve)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Yield Curve for key=" + YcDomKey + " is expected in the Market Data Manager");
	string domCcy(domCurve->GetCurrencyUnit()->GetCcyName());

	string YcBasisDomKey = YC_BASIS_KEY_NAME + itsDomesticCcy.GetCcyName();
	ARM_ZeroCurve* basisDomCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcBasisDomKey));
	if(!basisDomCurve)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Yield Curve for key=" + YcBasisDomKey + " is expected in the Market Data Manager");
	string basisDomCcy(basisDomCurve->GetCurrencyUnit()->GetCcyName());

	if(domCcy != basisDomCcy)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Domestic(=coupon) Basis Curve currency should consistent with reference curve");

	// Dom Mrs and Models
	string MrsDomKey = MRS_KEY_NAME + itsDomesticCcy.GetCcyName();
	ARM_ModelParam* mrsDomParam = dynamic_cast< ARM_ModelParam* >( GetMktDataManager()->GetData(MrsDomKey) );
	if(!mrsDomParam || mrsDomParam->GetType() != ARM_ModelParamType::MeanReversion)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : MRS Param for key=" + MrsDomKey + " is expected in the Market Data Manager");

	string OswDomModelKey = OSWMODEL_KEY_NAME + itsDomesticCcy.GetCcyName();
	ARM_BSModel* oswDomBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(OswDomModelKey) );
	if(!oswDomBSModel)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : domestic swaption B&S model for key=" + OswDomModelKey + " is expected in the Market Data Manager");

	int i(0);
	// Fgn Zc and Basis curves
	string YcForKey = YC_KEY_NAME + itsForeignCcy[i].GetCcyName();
	ARM_ZeroCurve* forCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcForKey));
	if(!forCurve)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Yield Curve for key=" + YcForKey + " is expected in the Market Data Manager");
	string forCcy(forCurve->GetCurrencyUnit()->GetCcyName());

	string YcBasisForKey = YC_BASIS_KEY_NAME + itsForeignCcy[i].GetCcyName();
	ARM_ZeroCurve* basisForCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcBasisForKey));
	if(!basisForCurve)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Yield Curve for key=" + YcBasisForKey + " is expected in the Market Data Manager");
	string basisForCcy(basisForCurve->GetCurrencyUnit()->GetCcyName());

	if(forCcy != basisForCcy)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Foreign Basis Curve currency should consistent with reference curve");

	// Dom Mrs and Models
	string MrsForKey = MRS_KEY_NAME + itsForeignCcy[i].GetCcyName();
	ARM_ModelParam* mrsForParam = dynamic_cast< ARM_ModelParam* >( GetMktDataManager()->GetData(MrsForKey) );
	if(!mrsForParam || mrsForParam->GetType() != ARM_ModelParamType::MeanReversion)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Foreign MRS Param for key=" + MrsForKey + " is expected in the Market Data Manager");

	string OswForModelKey = OSWMODEL_KEY_NAME + itsForeignCcy[i].GetCcyName();
	ARM_BSModel* oswForBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(OswForModelKey) );
	if(!oswForBSModel)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : foreign swaption B&S model for key=" + OswForModelKey + " is expected in the Market Data Manager");

	// Fx models and 
	string QFxKey = QFX_KEY_NAME + itsForeignCcy[i].GetCcyName() + "/" + itsDomesticCcy.GetCcyName();
	ARM_ModelParam* qFxParam = dynamic_cast< ARM_ModelParam* >( GetMktDataManager()->GetData(QFxKey) );
	if(!qFxParam || qFxParam->GetType() != ARM_ModelParamType::QParameter)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Foreign Q Fx Param for key=" + QFxKey + " is expected in the Market Data Manager");

	string ForexKey = FOREX_KEY_NAME + itsForeignCcy[i].GetCcyName() + "/" + itsDomesticCcy.GetCcyName();
	ARM_Forex* forex = dynamic_cast<ARM_Forex*>(GetMktDataManager()->GetData(ForexKey));
	if(!forex)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : forex=" + ForexKey + " is expected in the Market Data Manager");

	string FxModelKey = FXMODEL_KEY_NAME + itsForeignCcy[i].GetCcyName() + "/" + itsDomesticCcy.GetCcyName();
	ARM_MixtureModel_Fx* fxMixModel = dynamic_cast< ARM_MixtureModel_Fx* >( GetMktDataManager()->GetData(FxModelKey) );
	ARM_BSModel* fxBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(FxModelKey) );

	if (!fxMixModel && (itsModelType==Model1IRFX))
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : 1IRFX model needs fx Mixture model for key.");

	if(!fxMixModel && !fxBSModel)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : fx Mixture or BS model for key=" + FxModelKey + " is expected in the Market Data Manager");
	
	// Crossed currency correlations
	string CorrelMatrixKey = CORREL_KEY_NAME +"(" + itsDomesticCcy.GetCcyName() + "," + itsForeignCcy[i].GetCcyName() + "," + itsForeignCcy[i].GetCcyName() + "/" + itsDomesticCcy.GetCcyName() + ")";
} 


/////////////////////////////////////////////////////////////////
///	Class  : ARM_FxSnowBallKoCalculator
///	Routine: ColumnNames
///	Returns: ARM_RowInfo
///	Action : create the column names of the deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_FxSnowBallKoCalculator::ColumnNames() const
{
	size_t colNamesSize = sizeof(fxSbKoColNamesTable)/sizeof(fxSbKoColNamesTable[0]);
	vector<string>				colNamesVec(colNamesSize);
	vector<ARM_GP_VALUE_TYPE>	colTypeVec(colNamesSize, ARM_STRING); 
	
	for (size_t i=0; i<colNamesSize; i++)
	{
		colNamesVec[i] = fxSbKoColNamesTable[i];
	}

	ARM_RowInfo	rowInfo(colNamesVec, colTypeVec);

	return rowInfo;
}	

/////////////////////////////////////////////////////////////////
///	Class  : ARM_FxSnowBallKoCalculator
///	Routine: InitPriceableColumns
///	Returns: void
///	Action : initialise to 0 column to be able to sum
/////////////////////////////////////////////////////////////////
void ARM_FxSnowBallKoCalculator::InitPriceableColumns(vector<string>& rowDescVec, vector<ARM_GP_VALUE_TYPE>& rowTypeVec) const
{
	string zeroValue("0");

	for(int i=0; i<( sizeof(ProductToPriceColumns)/sizeof(ProductToPriceColumns[0]) ); i++)
	{
		rowDescVec[ProductToPriceColumns[i]] = zeroValue;
		rowTypeVec[ProductToPriceColumns[i]] = ARM_DOUBLE;
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_FxSnowBallKoCalculator
///	Routine: MiddleRows
///	Returns: ARM_RowInfo
///	Action : create a row of a deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_FxSnowBallKoCalculator::MiddleRows(size_t eventIdx, const ARM_DateStripCombiner& datesStructure ) const
{
	ARM_DateStripPtr vDateStrip = datesStructure.GetDateStrip(0);
	double asOfDate				= GetMktDataManager()->GetAsOfDate().GetJulian();
	size_t eventSize			= vDateStrip->GetResetDates()->size();
	size_t descSize				= sizeof(fxSbKoColNamesTable)/sizeof(fxSbKoColNamesTable[0]);

	vector<string> rowDescVec(descSize);
	vector<ARM_GP_VALUE_TYPE> rowTypeVec(descSize, ARM_MISSING_TYPE);

	// Set default 0 value for each column to be able to sum it
	InitPriceableColumns(rowDescVec, rowTypeVec);
	int i(0);

	/// WARNING : can have only one past fixing !
	if ( ( eventIdx < vDateStrip->GetResetDates()->size()-1 )
		&&
		 ( (*( vDateStrip->GetResetDates() ) )[eventIdx+1] < asOfDate) )
		return ARM_RowInfo();

	double resetDate = ( *( vDateStrip->GetResetDates() ) )[eventIdx];
	double payDate	 = ( *( vDateStrip->GetPaymentDates() ) )[eventIdx];

	// Timing management
	bool isArrearsCpn		= (itsCpnTiming==K_ARREARS);
	bool isRealFirstEvent	= (eventIdx==0);
	bool isFirstEvent		= (eventIdx==itsFirstEventIdx);
	bool isFirstCpnEvent	= ( isArrearsCpn ? ( eventIdx == (itsFirstEventIdx+1) ) : (eventIdx == itsFirstEventIdx) );
	bool isLastEvent		= (eventIdx==(eventSize-1));

	// Reset Date
	CC_Ostringstream resetDateDesc;
	double newResetDate(resetDate);

	resetDateDesc << CC_NS(std, fixed) << newResetDate;

	rowDescVec[ColAlias::ResetDate] = resetDateDesc.str();
	rowTypeVec[ColAlias::ResetDate] = ARM_DATE_TYPE;

	//__________________________________________IN ADVANCE CASE__________________________________________

	// Before last event
	if ( !(isArrearsCpn && isLastEvent) )
	{
		// Next Reset Date
		CC_Ostringstream nextResetDateDesc;
		if(!isArrearsCpn)
			nextResetDateDesc << CC_NS(std, fixed) << ( *( vDateStrip->GetResetDates() ) )[eventIdx];
		else if (!isLastEvent)
			nextResetDateDesc << CC_NS(std, fixed) << ( *( vDateStrip->GetResetDates() ) )[eventIdx+1];

		rowDescVec[ColAlias::NextResetDate] = nextResetDateDesc.str();
		rowTypeVec[ColAlias::NextResetDate] = ARM_DATE_TYPE;

		// Funding Start Date
		CC_Ostringstream fundingStartDateDesc;
		double fundingStart;
		if(!isArrearsCpn)
			fundingStart = ( *( vDateStrip->GetFlowStartDates() ) )[eventIdx];
		else if (!isLastEvent)
			fundingStart = ( *( vDateStrip->GetFlowStartDates() ) )[eventIdx+1];

		fundingStartDateDesc << CC_NS(std, fixed) << fundingStart;

		rowDescVec[ColAlias::FundingStartDate] = fundingStartDateDesc.str();
		rowTypeVec[ColAlias::FundingStartDate] = ARM_DATE_TYPE;

		// Funding End Date
		CC_Ostringstream fundingEndDateDesc;
		double fundingEnd;
		if(!isArrearsCpn)
			fundingEnd = ( *( vDateStrip->GetFlowEndDates() ) )[eventIdx];
		else if (!isLastEvent)
			fundingEnd = ( *( vDateStrip->GetFlowEndDates()  ) )[eventIdx+1];

		fundingEndDateDesc << CC_NS(std, fixed) << fundingEnd;

		rowDescVec[ColAlias::FundingEndDate] = fundingEndDateDesc.str();
		rowTypeVec[ColAlias::FundingEndDate] = ARM_DATE_TYPE;

		// Funding Fixed
		double IT(0.);
		double pastValue( const_cast<ARM_FxSnowBallKoCalculator*>(this)->PastLiborValue(vDateStrip, IT, eventIdx, rowDescVec, rowTypeVec) );

		CC_Ostringstream fundingFixedDesc;
		if (pastValue>0)
		{
			// Don't fill IT if no past value
			CC_Ostringstream ITDesc;
			ITDesc << CC_NS(std, fixed) << IT;
			rowDescVec[ColAlias::IT] = ITDesc.str();
			rowTypeVec[ColAlias::IT] = ARM_DOUBLE;

			fundingFixedDesc << pastValue;
		} else
			fundingFixedDesc << "0";

		rowDescVec[ColAlias::FixedFunding] = fundingFixedDesc.str();
		rowTypeVec[ColAlias::FixedFunding] = ARM_STRING;

		// Funding
		CC_Ostringstream fundingDesc;
		if ( (asOfDate==payDate) && ( const_cast<ARM_FxSnowBallKoCalculator*>(this)->GetDiscPricingMode() == ARM_DISC_ACCOUNTING_METH ) )
		{
			fundingDesc << "0";
		} else {
			if ( isRealFirstEvent && (pastValue>0) )
				fundingDesc << "FixedFunding[i]";
			else
				fundingDesc << "FixedFunding[i]+SWAP(" << YC_BASIS_KEY_NAME << itsDomesticCcy.GetCcyName() << ","
							<< ColNamesTable[ColAlias::FundingStartDate] << "[i]" << ","
							<< ColNamesTable[ColAlias::FundingEndDate] << "[i]" << ","
							<< "0" << ","
							<< ARM_ArgConvReverse_RcvOrPay.GetString(itsPayRec) << ","
							<< ARM_ArgConvReverse_LgNameFrequency.GetString(itsCpnFreq) << ","
							<< ARM_ArgConvReverse_LgNameDayCount.GetString(itsCpnDayCount) << ","
							<< ARM_ArgConvReverse_LgNameFrequency.GetString(itsFundFreq) << ","
							<< ARM_ArgConvReverse_LgNameDayCount.GetString(itsFundDayCount) << ","
							<< "FundMargin,FundNominal,,,FundDateStrip,FundDateStrip," 
							<< ColNamesTable[ColAlias::Offset] << "[i]," 
							<< ColNamesTable[ColAlias::Offset] << "[i],"
							<< ColNamesTable[ColAlias::Width] << "[i],"
							<< ColNamesTable[ColAlias::Width] << "[i])";
		}

		rowDescVec[ColAlias::Funding] = fundingDesc.str();
		rowTypeVec[ColAlias::Funding] = ARM_STRING;

		// Discount Funding
		CC_Ostringstream discountFundingDesc;

		discountFundingDesc << ColNamesTable[ColAlias::Funding] << "[i]";
		discountFundingDesc << "/" << "DF(" << YC_BASIS_KEY_NAME << itsDomesticCcy.GetCcyName() << ",";
		discountFundingDesc << ColNamesTable[ColAlias::NextResetDate] << "[i])";
		
		rowDescVec[ColAlias::DiscountFunding] = discountFundingDesc.str();
		rowTypeVec[ColAlias::DiscountFunding] = ARM_STRING;
		
		// OffSet
		CC_Ostringstream offsetDesc;

		offsetDesc << itsvFundIndex[eventIdx] - itsNbPastFixings << endl;
		rowDescVec[ColAlias::Offset] = offsetDesc .str();
		rowTypeVec[ColAlias::Offset] = ARM_DOUBLE;

		// Width
		CC_Ostringstream widthDesc;
		size_t width = itsvFundIndex[eventIdx+1] -itsvFundIndex[eventIdx];
		
		widthDesc << CC_NS(std, fixed) << width;

		rowDescVec[ColAlias::Width] = widthDesc.str();
		rowTypeVec[ColAlias::Width] = ARM_DOUBLE;
	}

	// After first event
	if ( !(isArrearsCpn && isFirstEvent) )
	{
		// Start Date
		CC_Ostringstream startDateDesc;
		double flowStartDate = ( *( vDateStrip->GetFlowStartDates() ) )[eventIdx];;

		startDateDesc << CC_NS(std, fixed) << flowStartDate;
		rowDescVec[ColAlias::StartDate] = startDateDesc.str();
		rowTypeVec[ColAlias::StartDate] = ARM_DATE_TYPE;

		// End Date
		CC_Ostringstream endDateDesc;
		double flowEndDate = ( *( vDateStrip->GetFlowEndDates() ) )[eventIdx];;

		endDateDesc << CC_NS(std, fixed) << flowEndDate;
		rowDescVec[ColAlias::EndDate] = endDateDesc.str();
		rowTypeVec[ColAlias::EndDate] = ARM_DATE_TYPE;

		// Interst term
		CC_Ostringstream ITDesc;
		double it = ( *( vDateStrip->GetInterestTerms() ) )[eventIdx];

		ITDesc << CC_NS(std, fixed) << it;
		rowDescVec[ColAlias::IT] = ITDesc.str();
		rowTypeVec[ColAlias::IT] = ARM_DOUBLE;

		// Notional
		CC_Ostringstream notionalDesc;
		double notional = itsCpnNominalCv.Interpolate(flowStartDate-asOfDate);
		
		notionalDesc << CC_NS(std, fixed) << notional;
		rowDescVec[ColAlias::Notional] = notionalDesc.str();
		rowTypeVec[ColAlias::Notional] = ARM_DOUBLE;

		// Strike
		CC_Ostringstream strikeDesc;
		double strike = itsStrikeCv.Interpolate(flowStartDate-asOfDate);
		
		strikeDesc << CC_NS(std, fixed) << strike;
		rowDescVec[ColAlias::Strike] = strikeDesc.str();
		rowTypeVec[ColAlias::Strike] = ARM_DOUBLE;

		// Barrier
		CC_Ostringstream barrierDesc;
		double barrier = itsBarrierCv.Interpolate(flowStartDate-asOfDate);
		
		barrierDesc << CC_NS(std, fixed) << barrier;
		rowDescVec[ColAlias::Barrier] = barrierDesc.str();
		rowTypeVec[ColAlias::Barrier] = ARM_DOUBLE;

		// Leverage
		CC_Ostringstream leverageDesc;
		double leverage = itsLeverageCv.Interpolate(flowStartDate-asOfDate);
		
		leverageDesc << CC_NS(std, fixed) << leverage;
		rowDescVec[ColAlias::Leverage] = leverageDesc.str();
		rowTypeVec[ColAlias::Leverage] = ARM_DOUBLE;

		// Coeff
		CC_Ostringstream coeffDesc;
		double coeff = itsCoeffCv.Interpolate(flowStartDate-asOfDate);
		
		coeffDesc << CC_NS(std, fixed) << coeff;
		rowDescVec[ColAlias::Coeff] = coeffDesc.str();
		rowTypeVec[ColAlias::Coeff] = ARM_DOUBLE;

		// Fx spot
		CC_Ostringstream fxSpotDesc;

		fxSpotDesc << const_cast<ARM_FxSnowBallKoCalculator*>(this)->FXValue(vDateStrip, eventIdx, i);
		rowDescVec[ColAlias::FxSpot] = fxSpotDesc.str();
		rowTypeVec[ColAlias::FxSpot] = ARM_STRING;
		

		// ?? I don't know what it does?
		const_cast<ARM_FxSnowBallKoCalculator*>(this)->DoPastReset(vDateStrip, eventIdx, rowDescVec, rowTypeVec);

		// IsAlive
		CC_Ostringstream isAliveDesc;

		isAliveDesc << "IF(" << ColNamesTable[ColAlias::FxSpot] << "[i]<"
					<< ColNamesTable[ColAlias::Barrier] << "[i],1,0)";
	
		if (!isFirstCpnEvent)
			isAliveDesc << "*" << ColNamesTable[ColAlias::IsAlive] << "[i-1]";
		
		rowDescVec[ColAlias::IsAlive] = isAliveDesc.str();
		rowTypeVec[ColAlias::IsAlive] = ARM_STRING;

		// Option
		CC_Ostringstream optionDesc;
		if (!itsIsPerform)
		{
			if (itsOptionType==K_CALL)
			{
				optionDesc  << "MAX(" << ColNamesTable[ColAlias::FxSpot] << "[i]" << "-"
							<< ColNamesTable[ColAlias::Strike] << ",0)";
			} else {
				optionDesc  << "MAX(" << ColNamesTable[ColAlias::Strike] << "[i]" << "-"
							<< ColNamesTable[ColAlias::FxSpot] << ",0)";
			}
		} else {
			if (itsOptionType==K_CALL)
			{
				optionDesc  << "MAX(1-" << ColNamesTable[ColAlias::FxSpot] << "[i]" << "*"
							<< ColNamesTable[ColAlias::Strike] << ",0)";
			} else {
				optionDesc  << "MAX(" << ColNamesTable[ColAlias::Strike] << "[i]" << "*"
							<< ColNamesTable[ColAlias::FxSpot] << "-1,0)";
			}
		}
		
		rowDescVec[ColAlias::Coupon] = optionDesc.str();
		rowTypeVec[ColAlias::Coupon] = ARM_STRING;

		// Snow Ball Coupon
		CC_Ostringstream couponDesc;
		couponDesc = optionDesc;

		if (!isFirstCpnEvent)
			couponDesc << "+" << ColNamesTable[ColAlias::Coupon] << "[i-1]";
		
		rowDescVec[ColAlias::Coupon] = couponDesc.str();
		rowTypeVec[ColAlias::Coupon] = ARM_STRING;

		// Real Coupon
		CC_Ostringstream realCouponDesc;
		couponDesc << ColNamesTable[ColAlias::Coupon] << "[i]";

		if (!isFirstCpnEvent)
			couponDesc << "*" << ColNamesTable[ColAlias::IsAlive] << "[i-1]";
		
		rowDescVec[ColAlias::RealCoupon] = realCouponDesc.str();
		rowTypeVec[ColAlias::RealCoupon] = ARM_STRING;

		// Real Funding
		CC_Ostringstream realFundingDesc;
		realFundingDesc << ColNamesTable[ColAlias::Funding] << "[i]";

		if (!isFirstCpnEvent)
			realFundingDesc << "*" << ColNamesTable[ColAlias::IsAlive] << "[i-1]";

		rowDescVec[ColAlias::RealFunding] = realFundingDesc.str();
		rowTypeVec[ColAlias::RealFunding] = ARM_STRING;

		// DF
		CC_Ostringstream dfDesc;
		dfDesc  << "DF(" << YC_KEY_NAME << itsDomesticCcy.GetCcyName() << ","
			<< ColNamesTable[ColAlias::EndDate] << ")" << "[i]";

		rowDescVec[ColAlias::DF] = dfDesc.str();
		rowTypeVec[ColAlias::DF] = ARM_STRING;

		// Proba
		CC_Ostringstream probaDesc;

		if (isFirstCpnEvent)
		{
			probaDesc	<< "UNPAY(1-" << ColNamesTable[ColAlias::IsAlive] << "[i]" << ")";
		} else 
		if (isLastEvent) {
			probaDesc	<< "UNPAY(" << ColNamesTable[ColAlias::IsAlive] << "[i-1]" << "*(1-"
						<< ColNamesTable[ColAlias::IsAlive] << "[i]" << ")"
						<< "+" << ColNamesTable[ColAlias::IsAlive] << "[i])";
		} else {
			probaDesc	<< "UNPAY(" << ColNamesTable[ColAlias::IsAlive] << "[i-1]" << "*(1-"
							<< ColNamesTable[ColAlias::IsAlive] << "[i]" << "))";
		}
		
		rowDescVec[ColAlias::Proba] = probaDesc.str();
		rowTypeVec[ColAlias::Proba] = ARM_STRING;

		// FX Snow Ball Knock Out
		CC_Ostringstream fxSbKoDesc;
		fxSbKoDesc	<< ColNamesTable[ColAlias::RealCoupon] << "[i]" << "-"
					<< ColNamesTable[ColAlias::RealFunding] << "[i]";

		rowDescVec[ColAlias::FXSBKO] = fxSbKoDesc.str();
		rowTypeVec[ColAlias::FXSBKO] = ARM_STRING;
	}

	return ARM_RowInfo(rowDescVec, rowTypeVec);
}

////////////////////////////////////////////////////
///	Class   : ARM_FxSnowBallKoCalculator
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////

string ARM_FxSnowBallKoCalculator::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os("Not filled yet!");

	return os.str();
}


//ARM_FxSnowBallKoCalculator

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/