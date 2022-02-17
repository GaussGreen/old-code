/*!
 *
 * Copyright (c) IXIS CIB August 2006 Paris
 *
 *	\file tarnfxcalculator.cpp
 *
 *  \brief file for the TARN FX Calculator
 *	\author  A. Lekrafi && P. Lam
 *	\version 1.0
 *	\date August 2006
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalculators/tarnfxcalculator.h"
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
//#include "gpbase/datestripconvert.h"
#include "gpbase/globalconstant.h"
#include "gpbase/globalconstant.h"

/// gpinfra
#include "gpinfra/mktdatamanagerrep.h"
#include "gpinfra/argconvdefault.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/gensecurity.h"
#include "gpinfra/dealdescription.h"
#include "gpinfra/pricingmodel.h"
#include "gpinfra/pricingadviser.h"
#include "gpinfra/genpricer.h"
#include "gpinfra/numerairefactory.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/pricerinfo.h"
#include "gpinfra/cstmanager.h"
#include "gpinfra/discretisationscheme.h"

/// gpcalib
#include "gpcalib/calibmethod.h"
#include "gpcalib/modelfitter.h"
#include "gpcalib/vanillapricer.h"


/// gpmodels
#include "gpmodels/2irfxModel.h"
#include "gpmodels/1irfxModel.h"
#include "gpmodels/NP1IRNFX.h"
#include "gpmodels/q1f.h"
#include "gpmodels/q1f_fx.h"
#include "gpmodels/Smiled_fx.h"
#include "gpmodels/mixture_fx.h"
#include "gpmodels/modelparamsq1f.h"
#include "gpmodels/eqfx_modelfactory.h"
#include "gpmodels/ModelParams_EqFxBase.h"
#include "gpmodels/forwardmarginbasis.h"

/// gpnummethods
#include "gpnummethods/scheduler.h"
#include "gpnummethods/meanrevertingsampler.h"
#include "gpnummethods/pathscheme.h"
#include "gpnummethods/pathschemefactory.h"
#include "gpnummethods/mcmethod.h"
#include "gpnummethods/sampler.h"
#include "gpnummethods/schedulerfactory.h"
#include "gpnummethods/samplerfactory.h"

/// gpnumlib
#include "gpnumlib/randomgenfactory.h"
#include "gpnumlib/transposer.h"

/// kernel
////#include <crv/volflat.h>
//#include <inst/forex.h>
//#include <inst/swaption.h>
//#include <inst/option.h>
//#include <inst/portfolio.h>
//#include <inst/swap.h>
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


/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNFXCalculator
///	Routine: Copy constructor
///	Returns: a built object
///	Action : builds the object
/////////////////////////////////////////////////////////////////
ARM_TARNFXCalculator::ARM_TARNFXCalculator( const ARM_TARNFXCalculator& rhs )
	:
	ARM_HybridIRFXCalculator(rhs)
{
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNFXCalculator
///	Routine: Destructor
///	Returns: void
///	Action : destroys the object
/////////////////////////////////////////////////////////////////
ARM_TARNFXCalculator::~ARM_TARNFXCalculator()
{
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNFXCalculator
///	Routine: Constructor
///	Returns: void
///	Action : builds the object (stand alone version)
/////////////////////////////////////////////////////////////////
ARM_TARNFXCalculator::ARM_TARNFXCalculator(
							 const ARM_Date& asOfDate,
							 const ARM_Date& startDate,
							 const ARM_Date& endDate,
							 const ARM_Currency& DomCcy,
							 const vector<ARM_Currency>& ForCcy,
							 const ARM_Currency& FundCcy,
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
							 const vector<ARM_Curve>& domesticCpn,
							 const vector<ARM_Curve>& foreignCpn,
							 const ARM_Curve& MinCpn,
							 const ARM_Curve& MaxCpn,
							 const vector<ARM_Curve>& InitialFX,
							 int fundFreq,
							 int fundDayCount,
							 const ARM_Curve& fundNominal,
							 const ARM_Curve& fundSpread,
							 const ARM_Curve& target,
							 const string& FXChoice,
							 int redemptionType,
							 int redemptionGap,
							 std::vector<double>& redemptionStrike,
							 const ARM_Curve& fees,
							 int intermediatePrices,
							 const ARM_StringVector& columnsToPrice,
							 const ARM_FXTARNPayoffType& payOffName,
							 ARM_FixingSched* pastFixings,
							 char* refDate,
							 ARM_Date& effDate)
:	
	ARM_HybridIRFXCalculator(asOfDate,
							 startDate,
							 endDate,
							 DomCcy,
							 ForCcy,
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
							 cpnNominal,
							 domesticCpn,
							 foreignCpn,
							 MinCpn,
							 MaxCpn,
							 InitialFX,
							 fundFreq,
							 fundDayCount,
							 fundNominal,
							 fundSpread,
							 redemptionType,
							 redemptionGap,
							 redemptionStrike,
							 target,
							 fees,
							 FXChoice,
							 intermediatePrices,
							 columnsToPrice,
							 payOffName,
							 pastFixings,
							 refDate,
							 effDate)
{
	itsNbFX = ForCcy.size();
	itsRedemptionResetDate = itsEndDate;
	itsRedemptionResetDate.PreviousBusinessDay( abs(itsRedemptionGap), const_cast<char*>(itsCpnResetCal.c_str()) );

	/// Check input datas
    CheckDataAndTimeIt();

	bool	needOtherPayoff = (itsIntermediatePrices?true:false);

    /// Create the Generic Security
	CreateAndSetDealDescription(itsPayModelName, itsColumnsToPrice, CreateCstManager(), false, needOtherPayoff);
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNFXCalculator
///	Routine: CreateCstManager
///	Returns: void
///	Action : create the const manager (static data).
/////////////////////////////////////////////////////////////////
ARM_CstManagerPtr ARM_TARNFXCalculator::CreateCstManager()
{
	ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();

	vector<string>	cstNames(3,"");
	cstNames[0] = "FundDateStrip";
	cstNames[1] = "FundNominal";
	cstNames[2] = "FundMargin";

	size_t fundSize = itsFundDateStrip->size();

	vector<ARM_GramFctorArg> cstVector;
	
	ARM_DateStripPtr cstDateStrip(static_cast<ARM_DateStrip*>(itsFundDateStrip->Clone()));
	if (GetFixings())
		cstDateStrip->ResizeAndBuilt(1, fundSize); // keep fwd flows only

	cstVector.push_back(ARM_GramFctorArg(cstDateStrip));

	ARM_GP_VectorPtr fundNominal(new std::vector<double>(fundSize));
	ARM_GP_VectorPtr fundMargin(new std::vector<double>(fundSize));

	std::vector<double>& resetDates = itsFundDateStrip->GetResetDates();

	size_t i;
	for (i = 0; i < fundSize; ++i)
	{
		(*fundNominal)[i] = itsCpnNominalCv.Interpolate((*resetDates)[i]-asOfDate.GetJulian());
		(*fundMargin)[i] = itsFundSpreadCv.Interpolate((*resetDates)[i]-asOfDate.GetJulian());
	}

	cstVector.push_back(ARM_GramFctorArg(fundNominal));
	cstVector.push_back(ARM_GramFctorArg(fundMargin));

	ARM_CstManagerPtr cstManagerPtr = ARM_CstManagerPtr(new ARM_CstManager(cstNames, cstVector));

	return cstManagerPtr;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNFXCalculator
///	Routine: CheckData
///	Returns: void
///	Action : check if TARN FX data are consistent
/////////////////////////////////////////////////////////////////
void ARM_TARNFXCalculator::CheckData()
{
	if (itsStartDate >= itsEndDate)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The start date frequency should be before the end date.");

	if (itsCpnFreq > itsFundFreq)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The funding frequency should be greater or equal than the coupon frequency.");

	if (itsRedemptionType == ARM_PRCSRedemptionType::dualOptionRedemption)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Dual Option is not supported.");
	
	itsvFundIndex.clear();
/// last: compute indexes to relate exer - funding - cpn		
	std::vector<double>& exerciseDates = itsCpnDateStrip->GetResetDates() ;
	std::vector<double>& fundStartDates = itsFundDateStrip->GetFlowStartDates();
	std::vector<double>& fundResetDates = itsFundDateStrip->GetResetDates() ;
	
	size_t exSize = itsCpnDateStrip->size()-1;
	size_t i = 0;
	for (size_t k(0); k<exSize; k++)
	{
		while ((*exerciseDates)[k]>(*fundStartDates)[i]) i++;
		if((*exerciseDates)[k]>(*fundResetDates)[i])
			ARM_THROW( ERR_INVALID_ARGUMENT, "TARN calculator : no call date allowed between fixing dates and start dates");
		itsvFundIndex.push_back(i);
	}
	itsvFundIndex.push_back(fundStartDates->size());
	/* 2 call dates cannot belong to the same period */
	/* this may happen when freq (call) > freq (exotic) or freq (funding) */
	for (k=0; k<exSize; k++)
	{
		if (itsvFundIndex[k]>=itsvFundIndex[k+1])
			ARM_THROW( ERR_INVALID_ARGUMENT, "TARN calculator : 2 call dates found within the same funding leg period. Freq (call) > Freq (Funding) is forbidden");
	}

	/* after each call call date, we must have same start dates for exo and funding legs */
	std::vector<double>& cpnStartDates  = itsCpnDateStrip->GetFlowStartDates() ;
	for (k=0; k<exSize; k++) 
	{		
		if (fabs((*fundStartDates)[itsvFundIndex[k]]-(*cpnStartDates)[k+1])> ARM_GlobalConstant::ARM_SEVENDAYS_LAG)
			ARM_THROW( ERR_INVALID_ARGUMENT, "TARN calculator : call date  defines a swap with mismatching exo and floating start dates (no broken period allowed)");
	}
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNFXCalculator
///	Routine: CheckMktData
///	Returns: void
///	Action : check if TARN FX market data are consistent
/////////////////////////////////////////////////////////////////
void ARM_TARNFXCalculator::CheckMktData()
{
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

	string YcFundKey = YC_KEY_NAME + itsFundingCcy.GetCcyName();
	string YcBasisFundKey = YC_BASIS_KEY_NAME + itsFundingCcy.GetCcyName();

	ARM_ZeroCurve* fundCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcFundKey));
    if(!fundCurve)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Yield Curve for key=" + YcFundKey + " is expected in the Market Data Manager");
    string fundCcy(fundCurve->GetCurrencyUnit()->GetCcyName());

	ARM_ZeroCurve* basisFundCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcBasisFundKey));
    if(!basisFundCurve)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Yield Curve for key=" + YcBasisFundKey + " is expected in the Market Data Manager");
    string basisFundCcy(basisFundCurve->GetCurrencyUnit()->GetCcyName());

    if(fundCcy != basisFundCcy)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Funding Basis Curve currency should consistent with reference curve");

	string FundForexKey = FOREX_KEY_NAME + itsFundingCcy.GetCcyName() + "/" + itsDomesticCcy.GetCcyName();
	ARM_Forex* fundforex = dynamic_cast<ARM_Forex*>(GetMktDataManager()->GetData(FundForexKey));
	if(!fundforex)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : forex=" + FundForexKey + " is expected in the Market Data Manager");

	string MrsDomKey = MRS_KEY_NAME + itsDomesticCcy.GetCcyName();
	ARM_ModelParam* mrsDomParam = dynamic_cast< ARM_ModelParam* >( GetMktDataManager()->GetData(MrsDomKey) );
    if(!mrsDomParam || mrsDomParam->GetType() != ARM_ModelParamType::MeanReversion)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : MRS Param for key=" + MrsDomKey + " is expected in the Market Data Manager");

	string OswDomModelKey = OSWMODEL_KEY_NAME + itsDomesticCcy.GetCcyName();
    ARM_BSModel* oswDomBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(OswDomModelKey) );
    if(!oswDomBSModel)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : domestic swaption B&S model for key=" + OswDomModelKey + " is expected in the Market Data Manager");

	string YcForKey;
	string YcBasisForKey;
	string ForexKey;
	string MrsForKey;
	string QFxKey;
	string CorrelMatrixKey;
	string OswForModelKey;
	string FxModelKey;
	int i = 0;

	for (i=0; i<itsNbFX; i++)
	{
		YcForKey = YC_KEY_NAME + itsForeignCcy[i].GetCcyName();
		ARM_ZeroCurve* forCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcForKey));
		if(!forCurve)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Yield Curve for key=" + YcForKey + " is expected in the Market Data Manager");
		string forCcy(forCurve->GetCurrencyUnit()->GetCcyName());

		YcBasisForKey = YC_BASIS_KEY_NAME + itsForeignCcy[i].GetCcyName();
		ARM_ZeroCurve* basisForCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcBasisForKey));
		if(!basisForCurve)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Yield Curve for key=" + YcBasisForKey + " is expected in the Market Data Manager");
		string basisForCcy(basisForCurve->GetCurrencyUnit()->GetCcyName());

		if(forCcy != basisForCcy)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Foreign Basis Curve currency should consistent with reference curve");

		ForexKey = FOREX_KEY_NAME + itsForeignCcy[i].GetCcyName() + "/" + itsDomesticCcy.GetCcyName();
		ARM_Forex* forex = dynamic_cast<ARM_Forex*>(GetMktDataManager()->GetData(ForexKey));
		if(!forex)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : forex=" + ForexKey + " is expected in the Market Data Manager");

		MrsForKey = MRS_KEY_NAME + itsForeignCcy[i].GetCcyName();
		ARM_ModelParam* mrsForParam = dynamic_cast< ARM_ModelParam* >( GetMktDataManager()->GetData(MrsForKey) );
		if(!mrsForParam || mrsForParam->GetType() != ARM_ModelParamType::MeanReversion)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Foreign MRS Param for key=" + MrsForKey + " is expected in the Market Data Manager");

		QFxKey = QFX_KEY_NAME + itsForeignCcy[i].GetCcyName() + "/" + itsDomesticCcy.GetCcyName();
		ARM_ModelParam* qFxParam = dynamic_cast< ARM_ModelParam* >( GetMktDataManager()->GetData(QFxKey) );
		if(!qFxParam || qFxParam->GetType() != ARM_ModelParamType::QParameter)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Foreign Q Fx Param for key=" + QFxKey + " is expected in the Market Data Manager");

		CorrelMatrixKey = CORREL_KEY_NAME +"(" + itsDomesticCcy.GetCcyName() + "," + itsForeignCcy[i].GetCcyName() + "," + itsForeignCcy[i].GetCcyName() + "/" + itsDomesticCcy.GetCcyName() + ")";
		
		OswForModelKey = OSWMODEL_KEY_NAME + itsForeignCcy[i].GetCcyName();
		ARM_BSModel* oswForBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(OswForModelKey) );
		if(!oswForBSModel)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : foreign swaption B&S model for key=" + OswForModelKey + " is expected in the Market Data Manager");

		FxModelKey = FXMODEL_KEY_NAME + itsForeignCcy[i].GetCcyName() + "/" + itsDomesticCcy.GetCcyName();
		ARM_MixtureModel_Fx* fxMixModel = dynamic_cast< ARM_MixtureModel_Fx* >( GetMktDataManager()->GetData(FxModelKey) );
		ARM_BSModel* fxBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(FxModelKey) );

		if (!fxMixModel && (itsModelType==Model1IRFX))
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : 1IRFX model needs fx Mixture model for key.");

		if(!fxMixModel && !fxBSModel)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : fx Mixture or BS model for key=" + FxModelKey + " is expected in the Market Data Manager");
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNFXCalculator
///	Routine: InitPriceableColumns
///	Returns: void
///	Action : initialise to 0 column to be able to sum
/////////////////////////////////////////////////////////////////
void ARM_TARNFXCalculator::InitPriceableColumns(vector<string>& rowDescVec, vector<ARM_GP_VALUE_TYPE>& rowTypeVec) const
{
    string	zeroValue("0");

	rowDescVec[ARM_TARNFXCalculator::Funding] = zeroValue;
    rowTypeVec[ARM_TARNFXCalculator::Funding] = ARM_DOUBLE;

	rowDescVec[ARM_TARNFXCalculator::Coupon] = zeroValue;
    rowTypeVec[ARM_TARNFXCalculator::Coupon] = ARM_DOUBLE;

	rowDescVec[ARM_TARNFXCalculator::DF] = zeroValue;
    rowTypeVec[ARM_TARNFXCalculator::DF] = ARM_DOUBLE;

	rowDescVec[ARM_TARNFXCalculator::PaidCoupon] = zeroValue;
    rowTypeVec[ARM_TARNFXCalculator::PaidCoupon] = ARM_DOUBLE;

	rowDescVec[ARM_TARNFXCalculator::RealCoupon] = zeroValue;
    rowTypeVec[ARM_TARNFXCalculator::RealCoupon] = ARM_DOUBLE;

	rowDescVec[ARM_TARNFXCalculator::IsAlive] = zeroValue;
    rowTypeVec[ARM_TARNFXCalculator::IsAlive] = ARM_DOUBLE;

	rowDescVec[ARM_TARNFXCalculator::Cap] = zeroValue;
    rowTypeVec[ARM_TARNFXCalculator::Cap] = ARM_DOUBLE;

	rowDescVec[ARM_TARNFXCalculator::Floor] = zeroValue;
    rowTypeVec[ARM_TARNFXCalculator::Floor] = ARM_DOUBLE;

	rowDescVec[ARM_TARNFXCalculator::Redemption] = zeroValue;
    rowTypeVec[ARM_TARNFXCalculator::Redemption] = ARM_DOUBLE;

	rowDescVec[ARM_TARNFXCalculator::RealRedemption] = zeroValue;
    rowTypeVec[ARM_TARNFXCalculator::RealRedemption] = ARM_DOUBLE;

	rowDescVec[ARM_TARNFXCalculator::RealFunding] = zeroValue;
    rowTypeVec[ARM_TARNFXCalculator::RealFunding] = ARM_DOUBLE;

	rowDescVec[ARM_TARNFXCalculator::Duration] = zeroValue;
    rowTypeVec[ARM_TARNFXCalculator::Duration] = ARM_DOUBLE;

	rowDescVec[ARM_TARNFXCalculator::Proba] = zeroValue;
    rowTypeVec[ARM_TARNFXCalculator::Proba] = ARM_DOUBLE;

	rowDescVec[ARM_TARNFXCalculator::TARN] = zeroValue;
    rowTypeVec[ARM_TARNFXCalculator::TARN] = ARM_DOUBLE;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNFXCalculator
///	Routine: PastLiborValue
///	Returns: CC_Ostringstream
///	Action : fill the rate value to pass to FUNDING
///			 and update resetDate
/////////////////////////////////////////////////////////////////
double ARM_TARNFXCalculator::PastLiborValue(ARM_DateStripPtr dateStrip, double& IT, size_t eventIdx, vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec)
{
	double fixing = 0.0;
	double resetDate = (*(dateStrip->GetResetDates()))[eventIdx];
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	
	if (itsNbPastFixings)
	{
		fixing = GetFixings()->GetLiborFixing(ARM_Date(resetDate), 
											  itsFundingCcy.GetCcyName(), 
											  GetSummitCurveIndexFromCurrency(itsFundingCcy.GetCcyName()), 
											  ConvertYearTermToStringMatu(1.0/itsFundFreq) );

		if ((resetDate <= asOfDate) && (fixing == -1))
			ARM_THROW( ERR_INVALID_ARGUMENT, "TARN FX : past libor fixing expected at " + ARM_Date(resetDate).toString());
	}

	if (fixing > 0)
	{
		CC_Ostringstream	resetDateDesc;
		double newResetDate = resetDate;
		if (resetDate <= asOfDate)
		{
			newResetDate = asOfDate + 1;
			if (itsCpnIntRule == K_ADJUSTED)
				newResetDate = ARM_Date(newResetDate).GoodBusinessDay(K_FOLLOWING, (char*)itsCpnResetCal.c_str()).GetJulian();

			double fundingEnd;
			if (itsCpnTiming == K_ADVANCE)
				fundingEnd = (*(dateStrip->GetFlowEndDates()))[eventIdx];
			else
				fundingEnd = (*(dateStrip->GetFlowEndDates()))[eventIdx+1];

			IT = CountYears(itsCpnDayCount, ARM_Date(newResetDate), ARM_Date(fundingEnd));
		}

		resetDateDesc << CC_NS(std, fixed) << newResetDate;
		rowDescVec[ColAlias::ResetDate] = resetDateDesc.str();
		rowTypeVec[ColAlias::ResetDate] = ARM_DATE_TYPE;
	}

	return fixing;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNFXCalculator
///	Routine: DoPastReset
///	Returns: 
///	Action : define the past fixings to take into account
///          update the target
/////////////////////////////////////////////////////////////////
void ARM_TARNFXCalculator::DoPastReset(ARM_DateStripPtr dateStrip, size_t eventIdx, vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec)
{
	ARM_Curve target(itsTarget);

	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	double startDate = (*(dateStrip->GetFlowStartDates()))[eventIdx];
	if (itsNbPastFixings)
	{
		double fixedCoupon = 0.0;
		ARM_Vector fixing(itsNbFX, 0.0);

		int i=0;
		double prevReset = (*(dateStrip->GetResetDates()))[eventIdx-1];
		for (i; i<itsNbFX; i++)
		{
			fixing.Elt(i) = GetFixings()->GetFxFixing( ARM_Date(prevReset), 
													   itsForeignCcy[i].GetCcyName(), 
													   itsDomesticCcy.GetCcyName() );
		}

		if (fixing.Elt(0) > 0)
		{
			double endDate = (*(dateStrip->GetFlowEndDates()))[eventIdx];

			double domCpn = itsDomesticCpnCv[0].Interpolate(0);
			double fgnCpn = itsForeignCpnCv[0].Interpolate(0);
			double initFX = itsInitialFXCv[0].Interpolate(0);
			double minCpn = itsMinCpnCv.Interpolate(0);
			double maxCpn = itsMaxCpnCv.Interpolate(0);

			if ((itsNbFX == 1) || (itsFXChoice == "First") || itsPayOffName == ARM_TARNFXPayoffType::TARNFX)
			{
				fixedCoupon = MIN(MAX(fgnCpn * fixing.Elt(0) / initFX - domCpn, minCpn), maxCpn);
			}
			else 
			{
				double domCpn2 = itsDomesticCpnCv[1].Interpolate(0);
				double fgnCpn2 = itsForeignCpnCv[1].Interpolate(0);
				double initFX2 = itsInitialFXCv[1].Interpolate(0);

				if (itsFXChoice == "Second")
				{
					fixedCoupon = MIN(MAX(fgnCpn2 * fixing.Elt(1) / initFX2 - domCpn2, minCpn), maxCpn);
				}
				else if (itsFXChoice == "Worst")
				{
					fixedCoupon = MIN(MAX(MIN(fgnCpn * fixing.Elt(0) / initFX - domCpn, fgnCpn2 * fixing.Elt(1) / initFX2 - domCpn2), minCpn), maxCpn);
				}
				else if (itsFXChoice == "Best")
				{
					fixedCoupon = MIN(MAX(MAX(fgnCpn * fixing.Elt(0) / initFX - domCpn, fgnCpn2 * fixing.Elt(1) / initFX2 - domCpn2), minCpn), maxCpn);
				}
			}
			
			double IT = 0.0;
			if (prevReset <= asOfDate)
			{
				double newResetDate = asOfDate + 1;
				if (itsCpnIntRule == K_ADJUSTED)
					newResetDate = ARM_Date(newResetDate).GoodBusinessDay(K_FOLLOWING, (char*)itsCpnResetCal.c_str()).GetJulian();

				IT = CountYears(itsCpnDayCount, ARM_Date(newResetDate), ARM_Date(endDate));
			}
			else
				IT = CountYears(itsCpnDayCount, ARM_Date(startDate), ARM_Date(endDate));

			target -= fixedCoupon*IT; 
		}
	}

	double targetdble = target.Interpolate(startDate - asOfDate);
	CC_Ostringstream	targetRedemptionDesc; 
	targetRedemptionDesc << targetdble;
	rowDescVec[ColAlias::Target] = targetRedemptionDesc.str();
	rowTypeVec[ColAlias::Target] = ARM_STRING;

}


////////////////////////////////////////////////////
///	Class   : ARM_TARNFXCalculator
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_TARNFXCalculator::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream	os;
	
	if (!itsDomSwaptionPF.IsNull())
	{
		os << indent << "Domestic Swaption Portfolio" << endl;
		os << indent << itsDomSwaptionPF->toString() << endl;
	}

	size_t i;
	for (i=0; i<itsNbFX; i++)
	{
		char num[3];
		sprintf(num, "%d", i);

		if (!itsForSwaptionPF[i].IsNull())
		{
			os << indent  << "Foreign Swaption Portfolio " << string(num) << endl;
			os << indent  << itsForSwaptionPF[i]->toString() << endl;
		}

		if (!itsFxOptionPF[i].IsNull())
		{
			os << indent  << "Fx Option Portfolio" << string(num) << endl;
			os << indent  << itsFxOptionPF[i]->toString() << endl;
		}
	}

	// A hook to prevent crash in the view
	if (itsModelType == Model2IRFX)
	{
		ARM_MultiAssetsModel* hybridModel = dynamic_cast<ARM_MultiAssetsModel*>(&*GetPricingModel());
		ARM_QModel1F_Fx* fxModel = dynamic_cast<ARM_QModel1F_Fx*>(&*(*hybridModel->GetModelMap())[ARM_2IRFXModel::FxModel]->Model() );
		fxModel->SetIntegratedVersion(false);
	}

    os <<  ARM_GenCalculator::toString(indent,nextIndent);

	// A hook to prevent crash in the view
	if (itsModelType == Model2IRFX)
	{
		ARM_MultiAssetsModel* hybridModel = dynamic_cast<ARM_MultiAssetsModel*>(&*GetPricingModel());
		ARM_QModel1F_Fx* fxModel = dynamic_cast<ARM_QModel1F_Fx*>(&*(*hybridModel->GetModelMap())[ARM_2IRFXModel::FxModel]->Model() );
		fxModel->SetIntegratedVersion(true);
	}

	return os.str();
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

