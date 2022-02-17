/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file tarncalculator.cpp
 *
 *  \brief file for the calculator for TARN Snow Ball
 *
 *	\author  JP Riaudel
 *	\version 1.0
 *	\date Janvier 2005
 */

/// put this header first to remove the ugly warning from the STL name mangling generation
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalculators/tarncalculatorsnowball.h"

/// gpbase
#include "gpbase/gpmatrix.h"

/// gpinfra
#include "gpinfra/gensecurity.h"
#include "gpinfra/dealdescription.h"

#include "gpinfra/pricingmodel.h"




#include "gpinfra/pricingadviser.h"

/// gpcalib
#include "gpcalib/calibmethod.h"

CC_BEGIN_NAMESPACE( ARM )

// Minimum Strike Value for Cap Floor calibration
const double MIN_CF_STRIKE			= 0.0025;

ARM_TARNCalculatorSnowBall::ARM_TARNCalculatorSnowBall(const ARM_Date& startDate,
    const ARM_Date& endDate,
    const ARM_Curve& strike,
	double coupon0,
    int payRec,
    int cpnDayCount,
    int cpnFreq,
    int cpnTiming,
    const string& cpnIndexTerm,
    int cpnIndexDayCount,
    const string& cpnResetCal,
    const string& cpnPayCal,
    int cpnResetGap,
	int intRule,
	int stubRule,
	long reverse,
    const ARM_Curve& leverage,
	const ARM_Curve& cpnMin,
	const ARM_Curve& cpnMax,
	const ARM_Curve& levPrev,
	double lifeTimeCapTarget,
	bool globalCapFlag,
    double lifeTimeFloorTarget,
    const ARM_Curve& fundSpread,
    int fundFreq,
    int fundDayCount,
    const ARM_Curve& nominal,
	const ARM_Curve& fees,
	const ARM_Curve& fundNominal,
	const std::vector<double>& nbIterations,
	ARM_ModelType modelType,
	CalibrationMode capFloorCalibMode,
	CalibrationMode digitalCalibMode,
	bool oswCalibFlag,
	bool controlVariableFlag,
	bool digitalSmoothingFlag,
	bool smiledFRMRescallingFlag,
	const string& genType1,
	const string& genType2,
	int firstNbTimes,
	int firstNbDims,
	const string& pathScheme,
	const string& pathOrder,
	bool antithetic,
	std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
	const ARM_Currency& cpnCcy,
    const ARM_Currency& fundingCcy,
    const ARM_Currency& basisCcy,
    const ARM_Currency& domesticCcy,
    const ARM_Currency& foreignCcy,
    const ARM_MarketData_ManagerRep& mktDataManager,
    const ARM_StringVector& mdmKeys)
	:  	ARM_TARNCalculator(startDate,
						   endDate,
						   strike,
						   payRec,
						   cpnDayCount,
						   cpnFreq,
						   cpnTiming,
						   cpnIndexTerm,
						   cpnIndexDayCount,
						   cpnResetCal,
						   cpnPayCal,
						   cpnResetGap,
						   intRule,
						   stubRule,
						   leverage,
						   cpnMin,
						   cpnMax,
						   levPrev,
						   lifeTimeCapTarget,
						   globalCapFlag,
						   lifeTimeFloorTarget,
						   fundSpread,
						   fundFreq,
						   fundDayCount,
						   nominal,
						   fees,
						   fundNominal,
						   nbIterations,
						   modelType,
						   capFloorCalibMode,
						   digitalCalibMode,
						   oswCalibFlag,
						   controlVariableFlag,
						   digitalSmoothingFlag,
						   smiledFRMRescallingFlag,
						   genType1,
						   genType2,
						   firstNbTimes,
						   firstNbDims,
						   pathScheme,
						   pathOrder,
						   antithetic,
						   productsToPrice,
						   mktDataManager),

	itsCoupon0( coupon0)
{
	SetName(ARM_TARN_SNOWBALL);	
	
	itsReverse = reverse == K_YES ? true: false;
	/// Set the coupon/payment currency (inherited from ARM_Security)
    SetCurrencyUnit( &const_cast<ARM_Currency&>(cpnCcy) );

    /// Set funding & basis currencies
    SetFundingCcy(const_cast<ARM_Currency&> (fundingCcy));
    SetBasisCcy(const_cast<ARM_Currency&> (basisCcy));

    /// Set keys for MDM datas access
	if ( mdmKeys.size() == 0 )
	{
		string cpnCcyName(cpnCcy.GetCcyName());
	    string fundingCcyName(fundingCcy.GetCcyName());

		// create internal keys
		ARM_StringVector keys(NbKeys);
		keys[YcKey]				= YC_KEY_NAME			+ cpnCcyName;
		keys[OswModelKey]		= OSWMODEL_KEY_NAME		+ cpnCcyName;
		keys[CfModelKey]		= CFMODEL_KEY_NAME		+ cpnCcyName;
		keys[MrsKey]			= MRS_KEY_NAME			+ cpnCcyName;
		keys[BetaKey]			= BETA_KEY_NAME			+ cpnCcyName;
		keys[CorrelKey]			= CORREL_KEY_NAME		+ cpnCcyName;
		keys[HumpKey]			= HUMP_KEY_NAME			+ cpnCcyName;
		keys[BetaCorrelKey]		= BETACORREL_KEY_NAME   + cpnCcyName;
 
		if( fundingCcyName == cpnCcyName )
		{
			keys.resize(NbKeys-1);
			// no basis, no forex but keep compatibility !
			keys[FundingKey] = keys[YcKey];
			keys[BasisKey]   = keys[YcKey];
			keys[FundBasisKey]= keys[YcKey];
		}
		else
		{
			keys[FundingKey]	= YC_KEY_NAME       + fundingCcyName;
			keys[BasisKey]		= YC_BASIS_KEY_NAME + string(basisCcy.GetCcyName());
			keys[FundBasisKey]   = YC_BASIS_KEY_NAME + string(fundingCcy.GetCcyName());
			keys[ForexKey]		= FOREX_KEY_NAME + string(foreignCcy.GetCcyName()) + "_" + string(domesticCcy.GetCcyName());
		}

        SetKeys(keys);
	}
	else
	{
		if(mdmKeys.size() < (NbKeys-5))
		{
			/// To be compatible with SBGM version, should have hump, betacorrel, fund and bs keys
			ARM_StringVector newMdMKeys(mdmKeys);
			newMdMKeys.resize(NbKeys-1);
			newMdMKeys[HumpKey]			= newMdMKeys[YcKey];
			newMdMKeys[BetaCorrelKey]   = newMdMKeys[YcKey];
			newMdMKeys[FundingKey]		= newMdMKeys[YcKey];
			newMdMKeys[BasisKey]		= newMdMKeys[YcKey];
			newMdMKeys[FundBasisKey]	= newMdMKeys[YcKey];

			SetKeys(newMdMKeys);
		}
		else if(mdmKeys.size() < (NbKeys-1))
		{
			/// To be compatible with non basis version, should at least have fund and bs keys
			ARM_StringVector newMdMKeys(mdmKeys);
			newMdMKeys.resize(NbKeys-1);
			newMdMKeys[FundingKey]	= newMdMKeys[YcKey];
			newMdMKeys[BasisKey]	= newMdMKeys[YcKey];
			newMdMKeys[FundBasisKey]= newMdMKeys[YcKey];
			SetKeys(newMdMKeys);
		}
		else
			SetKeys(mdmKeys);
	}

	// 12M is 
	if (GetCpnIndexTerm() == "12M")
		SetCpnIndexTerm("1Y");
    
    /// Check input datas
    CheckData();

    
    /// Set the coupon & funding model name alias i.e. keys to access to
    /// the right model for each GP keyword (1st argument)
    SetModelKeys();

   itsNeedOtherPayoff = itsProductsToPrice[ExerciseStrikesPrice] || itsProductsToPrice[ExerciseProbasPrice]|| controlVariableFlag;

    /// Create the Generic Security
	CreateAndSetDealDescription((GetKeys()[itsCpnModelKey]), CreateColumnNames(), ARM_CstManagerPtr(NULL),false,itsNeedOtherPayoff);
	CreateAndSetModelAndTimeIt();
	CreateAndSetCalibrationAndTimeIt();

}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNSnowBallCalculator
///	Routine: Copy constructor
///	Returns: a built object
///	Action : builds the object
/////////////////////////////////////////////////////////////////
ARM_TARNCalculatorSnowBall::ARM_TARNCalculatorSnowBall(const ARM_TARNCalculatorSnowBall& rhs)
                           : ARM_TARNCalculator(rhs)
{
	itsReverse = rhs.itsReverse;
	itsCoupon0 = rhs.itsCoupon0;
	itsNeedOtherPayoff = rhs.itsNeedOtherPayoff;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNSnowBallCalculator
///	Routine: assignment operator
///	Returns: 
///	Action : 
/////////////////////////////////////////////////////////////////
ARM_TARNCalculatorSnowBall& ARM_TARNCalculatorSnowBall::operator=(const ARM_TARNCalculatorSnowBall& rhs)
{
	if (this != &rhs)
	{
		ARM_TARNCalculator::operator=(rhs);

		itsReverse			= rhs.itsReverse;
		itsCoupon0			= rhs.itsCoupon0;
		itsNeedOtherPayoff  = rhs.itsNeedOtherPayoff;
	}

	return *this;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNSnowBallCalculator
///	Routine: CreateCpnDescription
///	Returns: the coupon description of a deal
///	Action : 
/////////////////////////////////////////////////////////////////
string ARM_TARNCalculatorSnowBall::CreateCpnDescription(bool isFirstExer) const
{
	string couponDesc;
	if (!isFirstExer)
	{
		couponDesc = "MIN(MAX(" + TARNColNamesTable[LevPrev] + "[i]*" + TARNColNamesTable[CouponWithoutIT] + "[i-1]";
	}
	else
	{
		char sCoupon0[15];
		sprintf(sCoupon0,"%lf",itsCoupon0);
		couponDesc = "MIN(MAX(" + TARNColNamesTable[LevPrev] + "[i]" + "*" + (string)sCoupon0;
	}

	if (itsReverse)
		couponDesc += "+";
	else
		couponDesc += "-";

	couponDesc += "(" + TARNColNamesTable[Strike] + "[i]-";
	couponDesc += TARNColNamesTable[Leverage] + "[i]*" + TARNColNamesTable[CpnIndex] + "[i])," + TARNColNamesTable[CouponMin] + "[i])," + TARNColNamesTable[CouponMax] + "[i])";

	return couponDesc;
}



string ARM_TARNCalculatorSnowBall::CreateCFCpnDdescription( bool isFirstExer, const string& cpnModelName) const
{
	CC_Ostringstream CouponCFDesc;
	CouponCFDesc << TARNColNamesTable[Leverage]  << "[i]*Caplet(" << cpnModelName << "," << TARNColNamesTable[IndexStartDate] << "[i],";
	CouponCFDesc << GetCpnIndexTerm()  << ",(" << TARNColNamesTable[Strike] << "[i]";
	
	if( isFirstExer)
	{
		char sCoupon0[15];
		sprintf(sCoupon0,"%lf",itsCoupon0);
		if( itsReverse )
			CouponCFDesc << "+";
		else
			CouponCFDesc << "-";
		CouponCFDesc << sCoupon0;
	}
	CouponCFDesc << ")/" << TARNColNamesTable[Leverage] << "[i],";
	if( itsReverse )
		CouponCFDesc << "F";
	else
		CouponCFDesc << "C";
	CouponCFDesc << ")*" << TARNColNamesTable[PaidNb] <<  "[i]*" << TARNColNamesTable[Nominal] << "[i]";;
	return CouponCFDesc.str();
}


void ARM_TARNCalculatorSnowBall::SetCVSwapMiddleText( ARM_DealDescriptionPtr& dealDesc, size_t knockOutIndex ) const
{
	size_t i;
	for( i=0; i<knockOutIndex; ++i )
	{
		CC_Ostringstream couponCFDesc;
		couponCFDesc << knockOutIndex-i;
		dealDesc->SetElem(i+1,PaidNb,couponCFDesc.str(),ARM_STRING);
	}
	for( i=knockOutIndex; i<dealDesc->GetRowsNb()-1; ++i )
	{
		dealDesc->SetElem(i+1,ControlVariable2,"0",ARM_DOUBLE);
	}

}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNSnowBallCalculator
///	Routine: CreateExStrikeDescription
///	Returns: the exercise strike description of a deal
///	Action : 
/////////////////////////////////////////////////////////////////
string ARM_TARNCalculatorSnowBall::CreateExStrikeDescription(bool isFirstExer, const string& cpnModelName) const
{
	string exerciseStrikeDesc;
	string basisModelName = GetKeys()[BasisKey];

	exerciseStrikeDesc = "MAX((" + TARNColNamesTable[Strike] + "[i]";

	if (itsReverse)
		exerciseStrikeDesc += "+(";
	else
		exerciseStrikeDesc += "-(";

	if (!isFirstExer)
		exerciseStrikeDesc += TARNColNamesTable[CouponWithoutIT] + "[i-1]-";
	else
	{
		char sCoupon0[15];
		sprintf(sCoupon0,"%lf",itsCoupon0);
		exerciseStrikeDesc += (string)sCoupon0 + "-";
	}

	exerciseStrikeDesc += TARNColNamesTable[TargetCap] + "[i]/" + TARNColNamesTable[IT] + "[i]))";

	char minCF[15];
	sprintf(minCF,"%lf",MIN_CF_STRIKE);

	exerciseStrikeDesc += "/" + TARNColNamesTable[Leverage] + "[i]," + minCF + ")*";
	exerciseStrikeDesc += "DF(" + basisModelName + "," + TARNColNamesTable[IndexEndDate] + "[i])";

	return exerciseStrikeDesc;
}
/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculatorSnowBall
///	Routine: UpdateCalibration
///	Returns: void
///	Action : update the calibration datas w.r.t. MarketDataManager
///          The context is an hedge ratios computation
/////////////////////////////////////////////////////////////////
void ARM_TARNCalculatorSnowBall::UpdateCalibration(bool isUpdateStrike)
{
	// first to re-compute Basis
	ARM_TARNCalculator::UpdateCalibration(isUpdateStrike);

	/// re-update the funding spread if necessairy
	if(IsBasis()){	
		/// Fake FIX FIX FIX
		CreateAndSetDealDescription((GetKeys()[itsCpnModelKey]), CreateColumnNames(), ARM_CstManagerPtr(NULL),false,itsNeedOtherPayoff);
	}
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
