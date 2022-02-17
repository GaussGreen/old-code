
/*!
 /** Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file CRASpreadCalculator.cpp
 *  \brief Calculaltor to valuate Callable Range Accrual
 *  on spread 
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date August 2005
 *
 */

/// this header comes first as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalculators/craspreadcalculator.h"

/// gpbase
#include "gpbase/datestrip.h"
#include "gpbase/singleton.h"
#include "gpbase/stringconvert.h"
#include "gpbase/stringmanip.h"
#include "gpbase/curveconvert.h"
#include "gpbase/datestripcombiner.h"
#include "gpbase/gpmatrixlinalg.h"
#include "gpbase/gplinalgconvert.h"
#include "gpbase/autocleaner.h"
#include "gpbase/utilityport.h"

/// gpinfra
#include "gpinfra/gensecurity.h"
#include "gpinfra/pricingadviser.h"
#include "gpinfra/pricingadviser.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/modelnamemap.h"
#include "gpinfra/numeraire.h"
#include "gpinfra/numerairefactory.h"
#include "gpinfra/argconvdefault.h"
#include "gpinfra/genpricer.h"
#include "gpinfra/pricerinfo.h"

/// gpcalib
#include "gpcalib/modelparamsfactory.h"
#include "gpcalib/vanillaspreadoption.h"
#include "gpcalib/vanillaswaption.h"
#include "gpcalib/calibmethod.h"
#include "gpcalib/modelfitterdes.h"
#include "gpcalib/kerneltogp.h"
#include "gpcalib/stripper.h"


/// gpmodels
#include "gpmodels/modelparamshw1f.h"
#include "gpmodels/modelparamshw2f.h"
#include "gpmodels/hw1f.h"
#include "gpmodels/hw2f.h"
#include "gpmodels/local_normal_model.h"
#include "gpmodels/local_normal_modelparams.h"
#include "gpmodels/multiassets.h"
#include "gpmodels/ForwardMarginBasis.h"

/// gpclosedforms
#include "gpclosedforms/vanilla_normal.h"

/// gpnummethods
#include "gpnummethods/treebase.h"
#include "gpnummethods/treefactory.h"
#include "gpnummethods/scheduler.h"
#include "gpnummethods/sampler.h"
#include "gpnummethods/truncator.h"
#include "gpnummethods/reconnector.h"
#include "gpnummethods/smoother.h"


// Kernerl inst
#include <crv/zerocurv.h>
#include <inst/spreadoption.h>
#include <inst/corridordblcondition.h>
#include <inst/fixleg.h>
#include <util/fromto.h>
#include <mod/bsconvadjust.h>

CC_BEGIN_NAMESPACE( ARM )

const double SIGMA_LOWER_BOUND				= 0.0005;
const double SIGMA_UPPER_BOUND				= 0.05;
const double SIGMA_DEFAULT_VALUE			= 0.005;
const double SIGMA_DEFAULT_VALUE_FOR_VEGA	= 0.005;
const double SIGMA_SHIFT_FOR_VEGA			= 0.0001;
const double MRS_DEFAULT_VALUE_FOR_VEGA		= 0.;

const double VOLRATIO_DEFAULT_VALUE			= 1.0;
const double CORREL_DEFAULT_VALUE			= 0.0;


const double DEFAULT_PRICE			= 1.0e+100;
const double DEFAULT_PRECISION		= 1.0;
const double DEFAULT_WEIGHT			= 1.0;

const double SENSI_CMS				= 0.00001;

/// Tree using 15 time steps per year in 1D and 10 in 2D
const size_t TREE_NBSTEPS_BEFORE_1ST		= 4;
const size_t TREE_NBSTEPS_PER_YEAR_1D		= 15;
const size_t TREE_NBSTEPS_PER_YEAR_2D		= 10;

//const size_t TREE_NBSTEPS_SHORT_TERM		= 2.0*K_YEAR_LEN ;
//const size_t TREE_NBSTEPS_PER_YEAR_ST_2D	= 35 ;

const double TREE_NBSTEPS_OPTIMAL_YF		= 2.0 ;
const size_t TREE_NBSTEPS_PER_YEAR_BEF_OPT	= 16;
const size_t TREE_NBSTEPS_PER_YEAR_AFT_OPT	= 8;
const double TREE_NBSTEPS_LONGTERM_YF		= 5.0 ;
const size_t TREE_NBSTEPS_PER_YEAR_AFT_LT	= 4;


const string LOCAL_CORRIDOR_MODEL_NAME  = "LOCALCAPFLOOR";
const string LOCAL_CORRIDOR2_MODEL_NAME  = "LOCALCAPFLOOR2";
const string LOCAL_CORRIDOR3_MODEL_NAME  = "LOCALCAPFLOOR3";

const unsigned int CALL_CRASPREAD_SCHED	= 0;
const double NON_CALL_FEE			= 1.0e15;

/// strike spread used for ARM_SpreadOption
const double CORRIDOR_STRIKE_SPREAD = 0.01;

// Vega to select swaption in the calibration portfolio
const double VEGA_MIN_TO_SELECT = 1e-12;


const double MRS_DEFAULT_VALUE      = 0.0;

const bool SYNCHRONIZE = true;
const bool FLOWBYFLOW = false;

const int  CORREL_UNSQUEEZER			= 0;
const bool CORREL_UNSQUEEZER_DEFAULT	= false;
const int  VOL_UNSQUEEZER				= 1;
const bool VOL_UNSQUEEZER_DEFAULT		= false;
const int  PV_ADJUSTER					= 2;
const bool PV_ADJUSTER_DEFAULT			= false;
const int  BOOTSTRAP_OPTIMIZER			= 3;
const bool BOOTSTRAP_OPTIMIZER_DEFAULT	= false;
const int  OLD_CALIBRATION				= 4;
const bool OLD_CALIBRATION_DEFAULT		= false;


const size_t OPTIM_RESET						= 0;
const size_t OPTIM_RESET_DAILYLIMIT				= 1;
const size_t OPTIM_RESET_WEEKLYLIMIT			= 2;
const size_t OPTIM_RESET_MONTHLYLIMIT			= 3;
const size_t OPTIM_RESET_NBDATA					= 4;

const double OPTIM_RESET_AVE_YEARLEN			= 365.25;
const double OPTIM_RESET_DEFAULT				= 0.0; // not activated
const double OPTIM_RESET_DAILYLIMIT_DEFAULT		= 3.0/52.0*OPTIM_RESET_AVE_YEARLEN; // 3W
const double OPTIM_RESET_WEEKLYLIMIT_DEFAULT	= 3.0/12.0*OPTIM_RESET_AVE_YEARLEN; // 3M
const double OPTIM_RESET_MONTHLYLIMIT_DEFAULT	= 10.0*OPTIM_RESET_AVE_YEARLEN;		// 10Y


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator 
///	Routine: Constructor 
///	Returns: void 
///	Action : builds the object
/////////////////////////////////////////////////////////////////
ARM_CRASpreadCalculator::ARM_CRASpreadCalculator( 
		const ARM_Currency& ccy,
		const ARM_Date&	startDate,
		const ARM_Date&	endDate,
		int payRec,
		const ARM_ReferenceValue&  notional,
		int callFreq,
		int callNotice,
		const string& callCal,
		const ARM_ReferenceValue&  callFees,
		int fundFreq,
		const ARM_ReferenceValue& fundSpread,
		int fundDayCount,
		int cpnDayCount,
		int cpnPayFreq,
		const string& cpnResetCal,
		const string& cpnPayCal,
		int payIndex,
		int payIndexResetTiming,
		const ARM_ReferenceValue& boostedFixRate,
		const ARM_ReferenceValue& payIndexMult,
		const ARM_ReferenceValue& cpnBarrierDown,
		const ARM_ReferenceValue& cpnBarrierUp,
		int cpnResetFreq,
		int cpnResetTiming,
		int cpnResetGap,
		int refIndex1,
		double refCoeff1,
		int refIndex2,
		double refCoeff2,
		int refIndex3,
		const ARM_ReferenceValue& barrierDown3,
		const ARM_ReferenceValue& barrierUp3,
		bool triple,
		bool vms,
		const ARM_ReferenceValue& boostedFixRate2,
		const ARM_ReferenceValue& cpnBarrierDown2,
		const ARM_ReferenceValue& cpnBarrierUp2,
		const ARM_ReferenceValue& boostedFixRate3,
		const ARM_ReferenceValue& cpnBarrierDown3,
		const ARM_ReferenceValue& cpnBarrierUp3,
		const ARM_ReferenceValue& refTenor1,
		ARM_ModelType modelType,
		CalibrationType calibType,
		vector<CalibStrikeType> calibStrikeType,
		ARM_MarketIRModel::VnsPricingMethod vnsPricingMethod,
		const ARM_StringVector& mdmKeys,
		const ARM_MarketData_ManagerRep& mktDataManager,
		const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice,
		const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& localCalibFlags,
		const ARM_GP_Vector& optimResetData,
		bool isExerciseProbas,
		size_t exerciseProbaOffset,
		ARM_Currency* fundccy,
		const ARM_ReferenceValue& fundNominal)   

: ARM_CRALocalCalculator( ccy,
						  startDate,
						  endDate,
						  payRec,
						  notional,
						  callFreq,
						  callNotice,
						  callCal,
						  callFees,
						  fundFreq,
						  fundSpread,
						  fundDayCount,	
						  ARM_ReferenceValue(),
						  cpnPayFreq,
						  cpnResetCal,
						  cpnPayCal,
						  K_FIXED,
						  boostedFixRate,
						  "",
						  0,
						  payIndexResetTiming,
						  cpnDayCount,
						  K_MOD_FOLLOWING,
						  K_ADJUSTED,
						  cpnBarrierDown,
						  cpnBarrierUp,
						  cpnResetFreq,
						  cpnResetTiming,
						  cpnResetGap,
						  K_FIXED,
						  "",
						  KACTUAL_360,
						  refCoeff1,
						  false, // TOMP false,
						  LocalDownUp,
						  mdmKeys,
						  mktDataManager,
						  productsToPrice,
						  -1,
						  -1,
						  true,
						  isExerciseProbas,
						  exerciseProbaOffset,
						  fundccy,
						  fundNominal),

	itsModelType(modelType),
	itsCalibType(calibType),
	itsCalibStrikeType(calibStrikeType),
	itsVnsPricingMethod(vnsPricingMethod),
	itsCSOPF(NULL),
	itsCSOPF2(NULL),
	itsCSOPF3(NULL),
	itsCoeff1(refCoeff1),
	itsCoeff2(refCoeff2),
	itsPayIndex(payIndex),
	itsRefIndex1(refIndex1),
	itsRefIndex2(refIndex2),
	itsPayIndexMult(payIndexMult),
	itsRefIndexType2(K_FIXED),
	itsRefTerm2(""),
	itsRefDayCount2(KACTUAL_360),
	itsRefCoeff2(refCoeff2),
	itsSwoptCalib(true),
	itsvExerFees(0),
	itsvIsExerDate(0),
	itsExerSize(0),
	itsvFundIndex(0),
	itsvFundSpread(0),
	itsvFundNominal(0),
	itsvRealFundSpread(0),
	itsvRealFundNominal(0),
	itsFundSize(0),
	itsvCpnIndex(0),
	itsvCpnNominal(0),
	itsCpnSize(0),
	itsVanillaArgCSO(ARM_VanillaArgPtr(NULL)),
	itsFixedCpnLongSensi(0),
	itsFixedCpnShortSensi(0),
	itsVarCpnLongSensi(0),
	itsVarCpnShortSensi(0),
	itsPayIndexFwd(0),
	itsFixedCpnValue(0),
	itsVarCpnValue(0),
	itsStructDateStrip(ARM_DateStripPtr(NULL)),
	itsFundDateStrip(ARM_DateStripPtr(NULL)),
	itsExerDateStripUnadj(ARM_DateStripPtr(NULL)),
	itsDateStrip1(ARM_DateStripPtr(NULL)),
	itsDateStrip2(ARM_DateStripPtr(NULL)),
	itsDateStrip3(ARM_DateStripPtr(NULL)),
	itsIndexVector(ARM_GP_VectorPtr(NULL)),
	itsRefIndex3(refIndex3),
	itsRateBarrierDown(barrierDown3),
	itsRateBarrierUp(barrierUp3),
	itsRefIndexType3(K_FIXED),
	itsRefTerm3(""),
	itsIsTripleRange(triple),
	itsIsVms(vms),
	itsBoostedFixRate2(boostedFixRate2),
	itsBoostedFixRate3(boostedFixRate3),
	itsCpnBarrierDown2(cpnBarrierDown2),
	itsCpnBarrierDown3(cpnBarrierDown3),
	itsCpnBarrierUp2(cpnBarrierUp2),
	itsCpnBarrierUp3(cpnBarrierUp3),
	itsTenor(refTenor1),
	itsFixedCpnLongSensi2(0),
	itsFixedCpnShortSensi2(0),
	itsVarCpnLongSensi2(0),
	itsVarCpnShortSensi2(0),
	itsPayIndexFwd2(0),
	itsFixedCpnValue2(0),
	itsVarCpnValue2(0),
	itsFixedValues2(0),
	itsFixedCpnLongSensi3(0),
	itsFixedCpnShortSensi3(0),
	itsVarCpnLongSensi3(0),
	itsVarCpnShortSensi3(0),
	itsPayIndexFwd3(0),
	itsFixedCpnValue3(0),
	itsVarCpnValue3(0),
	itsFixedValues3(0),
	itsCorrelUnSqueezer(localCalibFlags[CORREL_UNSQUEEZER]),
	itsVolUnSqueezer(localCalibFlags[VOL_UNSQUEEZER]),
	itsPVAdjuster(localCalibFlags[PV_ADJUSTER]),
	itsBootstrapOptimizer(localCalibFlags[BOOTSTRAP_OPTIMIZER]),
	itsSwitchOldCalibration(localCalibFlags[OLD_CALIBRATION]),
	itsOptimResetData(optimResetData)
{
	// Ugly, fake enum, the enum will be change in order to supporte the basis.
	if(IsBasis()){
		string fundCcyName(GetFundingCcy().GetCcyName());
		string ccyName(GetCurrencyUnit()->GetCcyName());

		ARM_StringVector keys(NbKeys);
		keys[YcKey]						= "YC_"				+ ccyName;
		keys[OswModelKey]				= "OSWMOD_"			+ ccyName;
		keys[CfModelKey]				= "CFMOD_"			+ ccyName;
		keys[MrsKey]					= "MRS_"			+ ccyName;
		keys[VolRatioKey]				= "VOLRATIO_"		+ ccyName;
		keys[MrsSpreadKey]				= "MRSSPREAD_"		+ ccyName;
		keys[CorrelKey]					= "CORREL_"			+ ccyName;
		keys[CorrelRateSpreadKey]		= "CORRELSPREAD_"	+ ccyName;
		keys[YcFundKey]					= "YC_"				+ fundCcyName;
		keys[YcBasis]					= "YC_BASIS_"		+ ccyName;
		keys[YcFundBasis]				= "YC_BASIS_"		+ fundCcyName;
		keys[ForexKey]					= "FOREX_"			+ fundCcyName+ccyName;

		SetKeys(keys);
	}
		
	if( itsOptimResetData.size()<1 ||
		(itsOptimResetData[OPTIM_RESET] && itsOptimResetData.size() < OPTIM_RESET_NBDATA) ||
		!itsOptimResetData[OPTIM_RESET] )
	{
		bool isDefaultOptim = itsOptimResetData.size()<1;
		double inputOptim = OPTIM_RESET_DEFAULT;
		if(!isDefaultOptim)
			inputOptim = itsOptimResetData[OPTIM_RESET];

		itsOptimResetData.resize(OPTIM_RESET_NBDATA);
		itsOptimResetData[OPTIM_RESET] = inputOptim;
		itsOptimResetData[OPTIM_RESET_DAILYLIMIT]	= OPTIM_RESET_DAILYLIMIT_DEFAULT;
		itsOptimResetData[OPTIM_RESET_WEEKLYLIMIT]	= OPTIM_RESET_WEEKLYLIMIT_DEFAULT;
		itsOptimResetData[OPTIM_RESET_MONTHLYLIMIT]	= OPTIM_RESET_MONTHLYLIMIT_DEFAULT;
	}

	Init();

	if( mdmKeys.size() == 0 )
		SetName(ARM_CALLABLE_CORRIDOR_SPREADOPTION);

}


/////////////////////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator 
///	Routine: Constructor 
///	Returns: void 
///	Action : builds the object from a portfolio
///          Warning : this constructor is used only to build a
///                    standard CRA (refIndex2 should not be taken into account)
/////////////////////////////////////////////////////////////////////////////////
ARM_CRASpreadCalculator::ARM_CRASpreadCalculator(ARM_OptionPortfolio* optionPortfolio,
												 ARM_StringVector& mdmKeys,
												 ARM_ModelType modelType,
												 CalibrationType	calibType,
												 vector<CalibStrikeType>	calibStrikeType,
												 ARM_MarketIRModel::VnsPricingMethod vnsPricingMethod,
												 ARM_ReferenceValue& payIndexMult,
												 const ARM_MarketData_ManagerRep& mktDataManager,
												 const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice)
:ARM_CRALocalCalculator(optionPortfolio,
						mdmKeys,
						mktDataManager,
						productsToPrice),
	itsModelType(modelType),
	itsCalibType(calibType),
	itsCalibStrikeType(calibStrikeType),
	itsVnsPricingMethod(vnsPricingMethod),
	itsCSOPF(NULL),
	itsCSOPF2(NULL),
	itsCSOPF3(NULL),
	itsCoeff1(1),
	itsCoeff2(0),
	itsRefCoeff2(0),
	itsRefIndex1(-1),
	itsRefIndex2(-1),
	itsRefIndex3(0),
	itsPayIndex(-1),
	itsPayIndexMult(payIndexMult),
	itsRefIndexType2(K_FIXED),
	itsRefTerm2(""),
	itsRefIndexType3(K_FIXED),
	itsRefTerm3(""),
	itsRefDayCount2(KACTUAL_360),
	itsSwoptCalib(true),
	itsvExerFees(0),
	itsvIsExerDate(0),
	itsExerSize(0),
	itsvFundIndex(0),
	itsvFundSpread(0),
	itsvFundNominal(0),
	itsvRealFundSpread(0),
	itsvRealFundNominal(0),
	itsFundSize(0),
	itsvCpnIndex(0),
	itsvCpnNominal(0),
	itsCpnSize(0),
	itsVanillaArgCSO(ARM_VanillaArgPtr(NULL)),
	itsFixedCpnLongSensi(0),
	itsFixedCpnShortSensi(0),
	itsVarCpnLongSensi(0),
	itsVarCpnShortSensi(0),
	itsPayIndexFwd(0),
	itsFixedCpnValue(0),
	itsVarCpnValue(0),
	itsStructDateStrip(ARM_DateStripPtr(NULL)),
	itsFundDateStrip(ARM_DateStripPtr(NULL)),
	itsExerDateStripUnadj(ARM_DateStripPtr(NULL)),
	itsDateStrip1(ARM_DateStripPtr(NULL)),
	itsDateStrip2(ARM_DateStripPtr(NULL)),
	itsDateStrip3(ARM_DateStripPtr(NULL)),
	itsIndexVector(ARM_GP_VectorPtr(NULL)),
	itsIsTripleRange(false),
	itsIsVms(false),
	itsFixedCpnLongSensi2(0),
	itsFixedCpnShortSensi2(0),
	itsVarCpnLongSensi2(0),
	itsVarCpnShortSensi2(0),
	itsPayIndexFwd2(0),
	itsFixedCpnValue2(0),
	itsVarCpnValue2(0),
	itsFixedValues2(0),
	itsFixedCpnLongSensi3(0),
	itsFixedCpnShortSensi3(0),
	itsVarCpnLongSensi3(0),
	itsVarCpnShortSensi3(0),
	itsPayIndexFwd3(0),
	itsFixedCpnValue3(0),
	itsVarCpnValue3(0),
	itsFixedValues3(0),
	itsCorrelUnSqueezer(CORREL_UNSQUEEZER_DEFAULT),
	itsVolUnSqueezer(VOL_UNSQUEEZER_DEFAULT),
	itsPVAdjuster(PV_ADJUSTER_DEFAULT),
	itsBootstrapOptimizer(BOOTSTRAP_OPTIMIZER_DEFAULT),
	itsSwitchOldCalibration(OLD_CALIBRATION_DEFAULT),
	itsOptimResetData(ARM_GP_Vector(OPTIM_RESET_NBDATA))
{
	itsOptimResetData[OPTIM_RESET]				= OPTIM_RESET_DEFAULT;
	itsOptimResetData[OPTIM_RESET_DAILYLIMIT]	= OPTIM_RESET_DAILYLIMIT_DEFAULT;
	itsOptimResetData[OPTIM_RESET_WEEKLYLIMIT]	= OPTIM_RESET_WEEKLYLIMIT_DEFAULT;
	itsOptimResetData[OPTIM_RESET_MONTHLYLIMIT]	= OPTIM_RESET_MONTHLYLIMIT_DEFAULT;

	/// Default is NOT basis
	SetFundingCcy(*(GetCurrencyUnit()));

	// if no reference index were specified, we deal with a standard CRA !!
	// RefIndex1 = corridor's reference index
	// RefIndex2 not relevant, we'll set refCoeff2 = 0
	if ((itsRefIndex1 == -1) && (itsRefIndex2 == -1))
	{
		itsRefIndex1 = optionPortfolio->GetCorridorLeg()->GetRefIndex()->GetIndexType();
		itsCoeff1 = 1;
		itsRefCoeff1 = 1;
		itsRefIndex2 = K_CMS2; 
		itsCoeff2 = 0;
		itsRefCoeff2 = 0;
	}
	if (itsPayIndex == -1)
	{
		itsPayIndex = optionPortfolio->GetCorridorLeg()->GetPaymentIndex()->GetIndexType();
	}

    Init();
}

/////////////////////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator 
///	Routine: Constructor 
///	Returns: void 
///	Action : builds the object from a portfolio, but without market data
///          model and calibration methods are created by InitCRASpreadFromSummit
/////////////////////////////////////////////////////////////////////////////////
ARM_CRASpreadCalculator::ARM_CRASpreadCalculator(const ARM_Date& asOfDate,
												 ARM_OptionPortfolio* optionPortfolio,
												 ARM_ReferenceValue& payIndexMult)
:ARM_CRALocalCalculator(asOfDate,
						(ARM_Security*)optionPortfolio),
	itsCSOPF(NULL),
	itsCSOPF2(NULL),
	itsCSOPF3(NULL),
	itsCoeff1(1),
	itsCoeff2(0),
	itsRefCoeff2(0),
	itsRefIndex1(-1),
	itsRefIndex2(-1),
	itsRefIndex3(0),
	itsPayIndex(-1),
	itsPayIndexMult(payIndexMult),
	itsRefIndexType2(K_FIXED),
	itsRefTerm2(""),
	itsRefIndexType3(K_FIXED),
	itsRefTerm3(""),
	itsRefDayCount2(KACTUAL_360),
	itsSwoptCalib(true),
	itsvExerFees(0),
	itsvIsExerDate(0),
	itsExerSize(0),
	itsvFundIndex(0),
	itsvFundSpread(0),
	itsvFundNominal(0),
	itsvRealFundSpread(0),
	itsvRealFundNominal(0),
	itsFundSize(0),
	itsvCpnIndex(0),
	itsvCpnNominal(0),
	itsCpnSize(0),
	itsVanillaArgCSO(ARM_VanillaArgPtr(NULL)),
	itsFixedCpnLongSensi(0),
	itsFixedCpnShortSensi(0),
	itsVarCpnLongSensi(0),
	itsVarCpnShortSensi(0),
	itsPayIndexFwd(0),
	itsFixedCpnValue(0),
	itsVarCpnValue(0),
	itsStructDateStrip(ARM_DateStripPtr(NULL)),
	itsFundDateStrip(ARM_DateStripPtr(NULL)),
	itsExerDateStripUnadj(ARM_DateStripPtr(NULL)),
	itsDateStrip1(ARM_DateStripPtr(NULL)),
	itsDateStrip2(ARM_DateStripPtr(NULL)),
	itsDateStrip3(ARM_DateStripPtr(NULL)),
	itsIndexVector(ARM_GP_VectorPtr(NULL)),
	itsIsTripleRange(false),
	itsIsVms(false),
	itsFixedCpnLongSensi2(0),
	itsFixedCpnShortSensi2(0),
	itsVarCpnLongSensi2(0),
	itsVarCpnShortSensi2(0),
	itsPayIndexFwd2(0),
	itsFixedCpnValue2(0),
	itsVarCpnValue2(0),
	itsFixedValues2(0),
	itsFixedCpnLongSensi3(0),
	itsFixedCpnShortSensi3(0),
	itsVarCpnLongSensi3(0),
	itsVarCpnShortSensi3(0),
	itsPayIndexFwd3(0),
	itsFixedCpnValue3(0),
	itsVarCpnValue3(0),
	itsFixedValues3(0),
	itsCorrelUnSqueezer(CORREL_UNSQUEEZER_DEFAULT),
	itsVolUnSqueezer(VOL_UNSQUEEZER_DEFAULT),
	itsPVAdjuster(PV_ADJUSTER_DEFAULT),
	itsBootstrapOptimizer(BOOTSTRAP_OPTIMIZER_DEFAULT),
	itsSwitchOldCalibration(OLD_CALIBRATION_DEFAULT),
	itsOptimResetData(ARM_GP_Vector(OPTIM_RESET_NBDATA))
{
	itsOptimResetData[OPTIM_RESET]				= OPTIM_RESET_DEFAULT;
	itsOptimResetData[OPTIM_RESET_DAILYLIMIT]	= OPTIM_RESET_DAILYLIMIT_DEFAULT;
	itsOptimResetData[OPTIM_RESET_WEEKLYLIMIT]	= OPTIM_RESET_WEEKLYLIMIT_DEFAULT;
	itsOptimResetData[OPTIM_RESET_MONTHLYLIMIT]	= OPTIM_RESET_MONTHLYLIMIT_DEFAULT;

	/// Default is NOT basis
	SetFundingCcy(*(GetCurrencyUnit()));

	SetName(ARM_CALLABLE_CORRIDOR_SPREADOPTION);

	// if no reference index were specified, we deal with a standard CRA !!
	// RefIndex1 = corridor's reference index
	// RefIndex2 not relevant, we'll set refCoeff2 = 0
	if ((itsRefIndex1 == -1) && (itsRefIndex2 == -1))
	{
		SetDegeneratedCalculator(CRAFromCCSO);

		itsRefIndex1 = optionPortfolio->GetCorridorLeg()->GetRefIndex()->GetIndexType();
		itsCoeff1 = 1;
		itsRefCoeff1 = 1;
		itsRefIndex2 = K_CMS2; 
		itsCoeff2 = 0;
		itsRefCoeff2 = 0;
	}
	if (itsPayIndex == -1)
	{
		itsPayIndex = optionPortfolio->GetCorridorLeg()->GetPaymentIndex()->GetIndexType();
	}
}

/////////////////////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator 
///	Routine: InitCRASpreadFromSummit 
///	Returns: void 
///	Action : initializes an empty CRA with market data
/////////////////////////////////////////////////////////////////////////////////
// FIXMEFRED: mig.vc8 (25/05/2007 14:53:26):missing return type
void ARM_CRASpreadCalculator::InitCRASpreadFromSummit(ARM_ZeroCurve* zc,
												 ARM_VolCurve* capVol,
												 ARM_VolCurve* rhoCap,
												 ARM_VolCurve* nuCap,
												 ARM_VolCurve* betaCap,
												 ARM_VolCurve* swoptVol,
												 ARM_VolCurve* rhoSwopt,
												 ARM_VolCurve* nuSwopt,
												 ARM_VolCurve* betaSwopt,
												 int SABRSigmaOrAlpha,
												 ARM_VolCurve* convAdjustVolCap,
												 ARM_VolCurve* convAdjustVolSwopt,
												 long convAdjustType,
												 ARM_CorrelManager* correlCorr,
												 ARM_VolCurve* correlDiagCap,
												 ARM_VolCurve* correlDiagSwopt,
												 double mrs,
												 double volRatio,
												 double mrsSpread,
												 double correl,
												 ARM_ModelType modelType,
												 CalibrationType calibType,
												 vector<CalibStrikeType> calibStrikeType,
												 ARM_MarketIRModel::VnsPricingMethod vnsPricingMethod,
												 const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice,
												 const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& localCalibFlags)
{
	itsCorrelUnSqueezer = localCalibFlags[CORREL_UNSQUEEZER];
	itsVolUnSqueezer = localCalibFlags[VOL_UNSQUEEZER];
	itsPVAdjuster = localCalibFlags[PV_ADJUSTER];
	itsBootstrapOptimizer = localCalibFlags[BOOTSTRAP_OPTIMIZER];
	itsSwitchOldCalibration = localCalibFlags[OLD_CALIBRATION];

	itsOptimResetData.resize(OPTIM_RESET_NBDATA);
	itsOptimResetData[OPTIM_RESET]				= OPTIM_RESET_DEFAULT;
	itsOptimResetData[OPTIM_RESET_DAILYLIMIT]	= OPTIM_RESET_DAILYLIMIT_DEFAULT;
	itsOptimResetData[OPTIM_RESET_WEEKLYLIMIT]	= OPTIM_RESET_WEEKLYLIMIT_DEFAULT;
	itsOptimResetData[OPTIM_RESET_MONTHLYLIMIT]	= OPTIM_RESET_MONTHLYLIMIT_DEFAULT;

	
	itsModelType = modelType;
	itsCalibType = calibType;
	itsCalibStrikeType = calibStrikeType;
	itsVnsPricingMethod = vnsPricingMethod;

	ARM_StringVector keys(NbKeys);
	string ccyName(zc->GetCurrencyUnit()->GetCcyName());

	keys[YcKey]			= "YC_" + ccyName;
    keys[OswModelKey]	= "OSWMOD_"	+ ccyName;
    keys[CfModelKey]	= "CFMOD_"	+ ccyName;
    keys[MrsKey]		= "MRS_"	+ ccyName;
	keys[VolRatioKey]	= "VOLRATIO_"	+ ccyName;
	keys[MrsSpreadKey]	= "MRSSPREAD_"	+ ccyName;
	keys[CorrelKey]		= "CORREL_"	+ ccyName;
	SetKeys(keys);

	ARM_BSConvAdjust* convAdjustManager = new ARM_BSConvAdjust(convAdjustType);

	ARM_Date asof = zc->GetAsOfDate();

	//Market Data Manager ----------------------------------------------------------->
	ARM_MarketData_ManagerRep* marketDataManager = &*GetMktDataManager();
	marketDataManager->SetAsOfDate(asof);

	//Cap Model
	ARM_BSSmiledModel* capSabrModel = NULL;
	ARM_VolCurve* ATMCapVol = (capVol->GetName() == ARM_VOL_CUBE) ? ((ARM_VolCube*)capVol)->GetATMVol() : capVol;
	if (rhoCap && nuCap)
	{
		int methodType = K_SABR_ARITH;
		if (betaCap)
		{
			if (betaCap->IsEqualToOne())
				methodType = K_SABR_ARITH;
			else
				methodType = K_SABR_IMPLNVOL;
		}

		capSabrModel = new ARM_BSSmiledModel(asof, 
											 0.0, 
											 zc, 
											 zc, 
											 ATMCapVol, 
											 K_YIELD, 
											 rhoCap, 
											 nuCap, 
											 methodType, 
											 betaCap,
                                             0.5,  // SABR Weight: Irrelevant
                                             SABRSigmaOrAlpha);
	}

	ARM_BSModel* capModel = new ARM_BSModel(asof, 
											zc, 
											NULL, //spreadlock
											convAdjustVolCap,
											capVol,
											correlCorr,
											convAdjustManager,
											zc, //discount
											correlDiagCap,
											NULL, //cash vol
											NULL, //spread vol
											K_2LOG,
											K_COMPUTED,
											capSabrModel); 

	//Swopt Model
	ARM_BSSmiledModel* swoptSabrModel = NULL;
	ARM_VolCurve* ATMSwoptVol = (swoptVol->GetName() == ARM_VOL_CUBE) ? ((ARM_VolCube*)swoptVol)->GetATMVol() : swoptVol;
	if (rhoSwopt && nuSwopt)
	{
		int methodType = K_SABR_ARITH;
		if (betaSwopt)
		{
			if (betaSwopt->IsEqualToOne())
				methodType = K_SABR_ARITH;
			else
				methodType = K_SABR_IMPLNVOL;
		}

		swoptSabrModel = new ARM_BSSmiledModel( asof, 
												0.0, 
												zc, 
												zc, 
												ATMSwoptVol, 
												K_YIELD, 
												rhoSwopt, 
												nuSwopt, 
												methodType, 
												betaSwopt,
                                                0.5,  // SABR Weight: Irrelevant
                                                SABRSigmaOrAlpha);
	}

	ARM_BSModel* swoptModel = new ARM_BSModel(asof, 
											  zc, 
											  NULL, //spreadlock
											  convAdjustVolSwopt,
											  swoptVol,			///AL ??
											  correlCorr,
											  NULL, //conv adjust manager
											  NULL, //discount
											  correlDiagSwopt,
											  NULL, //cash vol
											  NULL, //spread vol
											  K_2LOG,
											  K_COMPUTED,
											  swoptSabrModel);

	ARM_GP_Vector breakPointTimesExt(2);breakPointTimesExt[0]=0;breakPointTimesExt[1]=1;
	ARM_GP_Vector valuesExt(2);valuesExt[0]=correl;valuesExt[1]=correl;


	ARM_GP_Vector breakPointTimes(1, 0.0 );

	ARM_GP_Vector values(1, mrs );
	ARM_CurveModelParam* mrsCurve = new ARM_CurveModelParam(ARM_ModelParamType::MeanReversion, &values, &breakPointTimes);

	values[0] = volRatio;
	ARM_CurveModelParam* volRatioCurve = new ARM_CurveModelParam(ARM_ModelParamType::VolatilityRatio, &values, &breakPointTimes);

	values[0] = mrsSpread;
	ARM_CurveModelParam* mrsSpreadCurve = new ARM_CurveModelParam(ARM_ModelParamType::MeanReversionSpread, &values, &breakPointTimes);

	values[0] = correl;

	ARM_CurveModelParam* correlCurve;
	if (modelType==ARM_PricingModelType::HWM1F)
		correlCurve = new ARM_CurveModelParam(ARM_ModelParamType::Correlation, &values, &breakPointTimes);
	else
		correlCurve = new ARM_CurveModelParam(ARM_ModelParamType::Correlation, &valuesExt, &breakPointTimesExt);

	//Fill mktData Manager with models and curves
	marketDataManager->RegisterData(keys[YcKey], zc);
	marketDataManager->RegisterData(keys[OswModelKey], swoptModel);
	marketDataManager->RegisterData(keys[CfModelKey], capModel);
	marketDataManager->RegisterData(keys[MrsKey], mrsCurve); 
	marketDataManager->RegisterData(keys[MrsSpreadKey], mrsSpreadCurve);
	marketDataManager->RegisterData(keys[VolRatioKey], volRatioCurve);
	marketDataManager->RegisterData(keys[CorrelKey], correlCurve);

	if (capSabrModel)
		delete capSabrModel;
	capSabrModel = NULL;

	if (swoptSabrModel)
		delete swoptSabrModel;
	swoptSabrModel = NULL;

	delete swoptModel;
	swoptModel = NULL;

	delete capModel;
	capModel = NULL;

	delete mrsCurve;
	mrsCurve = NULL;

	delete volRatioCurve;
	volRatioCurve = NULL;

	delete mrsSpreadCurve;
	mrsSpreadCurve = NULL;

	delete correlCurve;
	correlCurve = NULL;

	itsProductsToPrice = productsToPrice;

	CheckData();
	//initialize calibration methods:

	// Pour les CCSO, il ne faut pas construire d'OptionPortfolio
	if(GetName() != ARM_CALLABLE_CORRIDOR_SPREADOPTION)
		Initialize();

    Init();
}

/////////////////////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator 
///	Routine: Constructor 
///	Returns: void 
///	Action : builds the object from a Swaption, but without market data
///          model and calibration methods are created by InitCRASpreadFromSummit
/////////////////////////////////////////////////////////////////////////////////
ARM_CRASpreadCalculator::ARM_CRASpreadCalculator(const ARM_Date& asOfDate,
												 ARM_Swaption* swaption,
												 ARM_ReferenceValue& payIndexMult)
:ARM_CRALocalCalculator(asOfDate,
						(ARM_Security*)swaption),
	itsCSOPF(NULL),
	itsCSOPF2(NULL),
	itsCSOPF3(NULL),
	itsCoeff1(1),
	itsCoeff2(0),
	itsRefCoeff2(0),
	itsRefIndex1(-1),
	itsRefIndex2(-1),
	itsRefIndex3(0),
	itsPayIndexMult(payIndexMult),
	itsRefIndexType2(K_FIXED),
	itsRefTerm2(""),
	itsRefIndexType3(K_FIXED),
	itsRefTerm3(""),
	itsIsTripleRange(false),
	itsIsVms(false),
	itsRefDayCount2(KACTUAL_360),
	itsSwoptCalib(true),
	itsvExerFees(0),
	itsvIsExerDate(0),
	itsExerSize(0),
	itsvFundIndex(0),
	itsvFundSpread(0),
	itsvFundNominal(0),
	itsvRealFundSpread(0),
	itsvRealFundNominal(0),
	itsFundSize(0),
	itsvCpnIndex(0),
	itsvCpnNominal(0),
	itsCpnSize(0),
	itsVanillaArgCSO(ARM_VanillaArgPtr(NULL)),
	itsFixedCpnLongSensi(0),
	itsFixedCpnShortSensi(0),
	itsVarCpnLongSensi(0),
	itsVarCpnShortSensi(0),
	itsPayIndexFwd(0),
	itsFixedCpnValue(0),
	itsVarCpnValue(0),
	itsStructDateStrip(ARM_DateStripPtr(NULL)),
	itsFundDateStrip(ARM_DateStripPtr(NULL)),
	itsExerDateStripUnadj(ARM_DateStripPtr(NULL)),
	itsDateStrip1(ARM_DateStripPtr(NULL)),
	itsDateStrip2(ARM_DateStripPtr(NULL)),
	itsDateStrip3(ARM_DateStripPtr(NULL)),
	itsIndexVector(ARM_GP_VectorPtr(NULL)),
	itsCorrelUnSqueezer(CORREL_UNSQUEEZER_DEFAULT),
	itsVolUnSqueezer(VOL_UNSQUEEZER_DEFAULT),
	itsPVAdjuster(PV_ADJUSTER_DEFAULT),
	itsBootstrapOptimizer(BOOTSTRAP_OPTIMIZER_DEFAULT),
	itsSwitchOldCalibration(OLD_CALIBRATION_DEFAULT),
	itsOptimResetData(ARM_GP_Vector(OPTIM_RESET_NBDATA))
{
	itsOptimResetData[OPTIM_RESET]				= OPTIM_RESET_DEFAULT;
	itsOptimResetData[OPTIM_RESET_DAILYLIMIT]	= OPTIM_RESET_DAILYLIMIT_DEFAULT;
	itsOptimResetData[OPTIM_RESET_WEEKLYLIMIT]	= OPTIM_RESET_WEEKLYLIMIT_DEFAULT;
	itsOptimResetData[OPTIM_RESET_MONTHLYLIMIT]	= OPTIM_RESET_MONTHLYLIMIT_DEFAULT;

	SetName(ARM_CALLABLE_CORRIDOR_SPREADOPTION);

	// if no reference index were specified, we deal with a swaption !!
	// RefIndex1 = swaption's floating index
	// RefIndex2 not relevant, we'll set refCoeff2 = 0
	if ((itsRefIndex1 == -1) && (itsRefIndex2 == -1))
	{
		SetDegeneratedCalculator(SwaptionFromCCSO);

		itsRefIndex1 = swaption->GetFloatLeg()->GetIRIndex()->GetIndexType();
		itsCoeff1 = 1;
		itsRefCoeff1 = 1;
		itsRefIndex2 = K_CMS10; 
		itsCoeff2 = 0;
		itsRefCoeff2 = 0;
	}

	itsPayIndex = swaption->GetFixedLeg()->GetIRIndex()->GetIndexType();
}

void ARM_CRASpreadCalculator::Init(void)
{
	int payIndexType, indexType1, indexType2;
	string payIndexTerm = FromIndexTypeToTermAndType(itsPayIndex, payIndexType);
	string indexTerm1   = FromIndexTypeToTermAndType(itsRefIndex1, indexType1);
	string indexTerm2   = FromIndexTypeToTermAndType(itsRefIndex2, indexType2);

	SetBoostedVarTerm(payIndexTerm);
	SetBoostedIndexType(payIndexType);

	SetRefTerm1(indexTerm1);
	SetRefIndexType1(indexType1);

	SetRefTerm2(indexTerm2);
	SetRefIndexType2(indexType2);

	itsRefTerm3   = FromIndexTypeToTermAndType(itsRefIndex3, itsRefIndexType3);

	if(GetKeys().size() != 0)
	{
		ARM_StringVector pricedColumns = PricedColumnNames();

		DatesStructure();

		CreateCSOPortfolio_FakeForDateStripOnly();

		ARM_CstManagerPtr cstManagerPtr = CreateCstManager();


		ComputeProductVectorsFromCurves();

		char* ccyStr	= GetCcy().GetCcyName();
		string modelName(IsBasis() ?  GetKeys()[YcBasis] : string(GetCcy().GetCcyName()));

		ARM_StringVector::const_iterator found = find( pricedColumns.begin(), pricedColumns.end(),
			LocalCapFloorCRAColNamesTable[ LocalCapFloorCRAProductToPriceColumns[ARM_CRALocalCalculator::localCfCorridorLegPrices] ] );
		bool isOtherPayoff = (found!=pricedColumns.end());
		CreateAndSetDealDescriptionAndTimeIt(modelName, pricedColumns, cstManagerPtr,false,isOtherPayoff);

		CreateAndSetModel();

		CreateAndSetCalibration();

		/// Skip for CRA2 because it is quite long
		if(!(IsDbleCorridor()))
			Calibrate();
	}

	// else : les market data ne sont pas encore settées.
	// Init() sera rappelée par InitCRASpreadFromSummit() après le set des market data
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: Copy constructor
///	Returns: a built object
///	Action : builds the object
/////////////////////////////////////////////////////////////////
ARM_CRASpreadCalculator::ARM_CRASpreadCalculator(const ARM_CRASpreadCalculator& rhs)
:	ARM_CRALocalCalculator(rhs),
	itsCSOPF(rhs.itsCSOPF),
	itsCSOPF2(rhs.itsCSOPF2),
	itsCSOPF3(rhs.itsCSOPF3),
	itsModelType(rhs.itsModelType),
	itsCalibType(rhs.itsCalibType),
	itsCalibStrikeType(rhs.itsCalibStrikeType),
	itsVnsPricingMethod (rhs.itsVnsPricingMethod),
	itsPayIndex(rhs.itsPayIndex),
	itsCoeff1(rhs.itsCoeff1),
	itsCoeff2(rhs.itsCoeff2),
	itsRefIndex1(rhs.itsRefIndex1),
	itsRefIndex2(rhs.itsRefIndex2),
	itsPayIndexMult(rhs.itsPayIndexMult),
	itsRefIndexType2(rhs.itsRefIndexType2),
	itsRefTerm2(rhs.itsRefTerm2),
	itsRefDayCount2(rhs.itsRefDayCount2),
	itsRefCoeff2(rhs.itsRefCoeff2),
	itsSwoptCalib(rhs.itsSwoptCalib),
	itsvExerFees(rhs.itsvExerFees),
	itsvIsExerDate(rhs.itsvIsExerDate),
	itsExerSize(rhs.itsExerSize),
	itsvFundIndex(rhs.itsvFundIndex),
	itsvFundSpread(rhs.itsvFundSpread),
	itsvFundNominal(rhs.itsvFundNominal),
	itsvRealFundSpread(rhs.itsvRealFundSpread),
	itsvRealFundNominal(rhs.itsvRealFundNominal),
	itsFundSize(rhs.itsFundSize),
	itsvCpnIndex(rhs.itsvCpnIndex),
	itsvCpnNominal(rhs.itsvCpnNominal),
	itsCpnSize(rhs.itsCpnSize),
	itsVanillaArgCSO(ARM_VanillaArgPtr((!rhs.itsVanillaArgCSO.IsNull()?static_cast<ARM_VanillaArg*>(rhs.itsVanillaArgCSO->Clone()):NULL))),
	itsFixedCpnLongSensi(rhs.itsFixedCpnLongSensi),
	itsFixedCpnShortSensi(rhs.itsFixedCpnShortSensi),
	itsVarCpnLongSensi(rhs.itsVarCpnLongSensi),
	itsVarCpnShortSensi(rhs.itsVarCpnShortSensi),
	itsPayIndexFwd(rhs.itsPayIndexFwd),
	itsFixedCpnValue(rhs.itsFixedCpnValue),
	itsVarCpnValue(rhs.itsVarCpnValue),
	itsFundDateStrip(ARM_DateStripPtr((!rhs.itsFundDateStrip.IsNull()?static_cast<ARM_DateStrip*>(rhs.itsFundDateStrip->Clone()):NULL))),
	itsStructDateStrip(ARM_DateStripPtr((!rhs.itsStructDateStrip.IsNull()?static_cast<ARM_DateStrip*>(rhs.itsStructDateStrip->Clone()):NULL))),
	itsExerDateStripUnadj(ARM_DateStripPtr((!rhs.itsExerDateStripUnadj.IsNull()?static_cast<ARM_DateStrip*>(rhs.itsExerDateStripUnadj->Clone()):NULL))),
	itsDateStrip1(ARM_DateStripPtr((!rhs.itsDateStrip1.IsNull()?static_cast<ARM_DateStrip*>(rhs.itsDateStrip1->Clone()):NULL))),
	itsDateStrip2(ARM_DateStripPtr((!rhs.itsDateStrip2.IsNull()?static_cast<ARM_DateStrip*>(rhs.itsDateStrip2->Clone()):NULL))),
	itsDateStrip3(ARM_DateStripPtr((!rhs.itsDateStrip3.IsNull()?static_cast<ARM_DateStrip*>(rhs.itsDateStrip3->Clone()):NULL))),
	itsIndexVector(ARM_GP_VectorPtr((!rhs.itsIndexVector.IsNull()?static_cast<ARM_GP_Vector*>(rhs.itsIndexVector->Clone()):NULL))),
	itsRefIndex3(rhs.itsRefIndex3),
	itsRefIndexType3(rhs.itsRefIndexType3),
	itsRefTerm3(rhs.itsRefTerm3),
	itsRateBarrierDown(rhs.itsRateBarrierDown),
	itsRateBarrierUp(rhs.itsRateBarrierUp),
	itsIsTripleRange(rhs.itsIsTripleRange),
	itsIsVms(rhs.itsIsVms),
	itsBoostedFixRate2(rhs.itsBoostedFixRate2),
	itsBoostedFixRate3(rhs.itsBoostedFixRate3),
	itsCpnBarrierDown2(rhs.itsCpnBarrierDown2),
	itsCpnBarrierDown3(rhs.itsCpnBarrierDown3),
	itsCpnBarrierUp2(rhs.itsCpnBarrierUp2),
	itsCpnBarrierUp3(rhs.itsCpnBarrierUp3),
	itsTenor(rhs.itsTenor),
	itsCorrelUnSqueezer(rhs.itsCorrelUnSqueezer),
	itsVolUnSqueezer(rhs.itsVolUnSqueezer),
	itsBootstrapOptimizer(rhs.itsBootstrapOptimizer),
	itsSwitchOldCalibration(rhs.itsSwitchOldCalibration),
	itsPVAdjuster(rhs.itsPVAdjuster),
	itsOptimResetData(rhs.itsOptimResetData)

{
	SetExerDateStrip(ARM_DateStripPtr(!rhs.itsExerDateStrip.IsNull()?static_cast<ARM_DateStrip*>(rhs.itsExerDateStrip->Clone()):NULL));
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: Destructor
///	Returns: void
///	Action : destroys the object
/////////////////////////////////////////////////////////////////
ARM_CRASpreadCalculator::~ARM_CRASpreadCalculator()
{
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: assignment operator
///	Returns: 
///	Action : 
/////////////////////////////////////////////////////////////////
ARM_CRASpreadCalculator& ARM_CRASpreadCalculator::operator=(const ARM_CRASpreadCalculator& rhs)
{
	if( this != & rhs )
	{
		this->~ARM_CRASpreadCalculator();
        new (this) ARM_CRASpreadCalculator(rhs);
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class   : ARM_CRASpreadCalculator
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
///           Call copy constructor
////////////////////////////////////////////////////
ARM_Object* ARM_CRASpreadCalculator::Clone() const
{
	return(new ARM_CRASpreadCalculator(*this));
}


////////////////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: PricedColumnNames
///	Returns: ARM_StringVector
///	Action : create the priced column names of the deal description
//////////////////////////////////////////////////////////////////////////////
ARM_StringVector ARM_CRASpreadCalculator::PricedColumnNames() const
{
	ARM_StringVector pricedColumns(ARM_CRALocalCalculator::PricedColumnNames());

	//if(IsDbleCorridor() && itsCalibType == DIAG_CALIBRATION && itsCalibStrikeType[0] == FRONTIER)
	if(itsCalibType == DIAG_CALIBRATION && itsCalibStrikeType[0] == FRONTIER)
		pricedColumns.push_back(LocalCapFloorCRAColNamesTable[localCfFrontier]);

	return pricedColumns;
}


/////////////////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: DatesStructure
///	Returns: ARM_DateStripVector
///	Action : create the list of all event dates of the CRA. 
///			 The DateStripCombiner merges event dates of each leg
/////////////////////////////////////////////////////////////////////////////
ARM_DateStripCombiner ARM_CRASpreadCalculator::DatesStructure() const
{
	//General datas
	ARM_Currency ccy = GetCcy();
	
	ARM_DateStripVector SchedVect(1,NULL);

	double startDate = itsStartDate.GetJulian();
	double endDate = itsEndDate.GetJulian();

    ARM_DateStrip theCallSched;



    if ((itsExerDateStrip.IsNull())
        &&
        (!(itsIsPortfolioNoticeDays))
       )
    {
	    ARM_DateStrip CallSched(itsStartDate.GetJulian(), 
		                        itsEndDate.GetJulian(),
		                        itsCallFreq, 
		                        itsBoostedDayCount,
		                        itsCallCal.c_str(), 
		                        itsBoostedAdjRule, 
		                        itsBoostedIntRule, 
		                        K_SHORTSTART, 
		                        itsCallNotice,
		                        itsCallFreq, 
		                        GETDEFAULTVALUE, 
		                        itsCallCal.c_str(), 
		                        K_ADVANCE, 
		                        K_ARREARS, 
		                        true);

        theCallSched = CallSched;

	    SchedVect[CALL_CRASPREAD_SCHED]    = &theCallSched;

	    itsExerDateStrip = ARM_DateStripPtr(static_cast<ARM_DateStrip*>(CallSched.Clone()));

		ARM_DateStrip CallSchedUnadj(itsStartDate.GetJulian(), 
		                        itsEndDate.GetJulian(),
		                        itsCallFreq, 
		                        itsBoostedDayCount,
		                        itsCallCal.c_str(), 
		                        itsBoostedAdjRule, 
		                        K_UNADJUSTED, 
		                        K_SHORTSTART, 
		                        itsCallNotice,
		                        itsCallFreq, 
		                        GETDEFAULTVALUE, 
		                        itsCallCal.c_str(), 
		                        K_ADVANCE, 
		                        K_ARREARS, 
		                        true);

	    itsExerDateStripUnadj = ARM_DateStripPtr(static_cast<ARM_DateStrip*>(CallSchedUnadj.Clone()));
    }
    else
    {
        SchedVect[CALL_CRASPREAD_SCHED] = &(*itsExerDateStrip);
    }

	ARM_INDEX_TYPE indexType = ((ARM_Currency*)GetCurrencyUnit())->GetVanillaIndexType();
	int resetGap			 = - GetCurrencyUnit()->GetSpotDays();

	// funding schedule
    // N.B: take the real funding leg in the case of Bermuda
    itsFundDateStrip = ARM_DateStripPtr(
								new ARM_DateStrip(
		                                    itsStartDate.GetJulian(),
		                                    itsEndDate.GetJulian(),
		                                    itsFundFreq,
		                                    itsFundDayCount,
		                                    itsCpnResetCal.c_str(),
		                                    itsFundAdjRule,
		                                    itsFundIntRule,
		                                    itsFundStubRule,
		                                    resetGap,
		                                    itsFundFreq, 
		                                    GETDEFAULTVALUE, 
		                                    itsCpnPayCal.c_str(),
		                                    K_ADVANCE, 
		                                    K_ARREARS,
											true,
											GETDEFAULTVALUESTR,
											GETDEFAULTVALUE,
											GETDEFAULTVALUE,
											K_ACCRUED,
											K_FOLLOWING)); //adjust 1st date

	/// To avoid schedule inconsistency, if corridor keyword is based the additionnal dateStrip
	/// the "Structure DateStrip" must be based on the same schedule
	if(itsDateStrip3 == ARM_DateStripPtr(NULL))
	{
		itsStructDateStrip = ARM_DateStripPtr(
									new ARM_DateStrip(
												itsStartDate.GetJulian(),
												itsEndDate.GetJulian(),
												itsCpnResetFreq,
												itsBoostedDayCount,
												itsCpnResetCal.c_str(),
												itsCpnAdjRule,
												itsCpnIntRule,
												itsCpnStubRule,
												resetGap,
												itsCpnPayFreq,
												GETDEFAULTVALUE, 
												itsCpnPayCal.c_str(),
												K_ADVANCE, 
												K_ARREARS,
												true,
												GETDEFAULTVALUESTR,
												GETDEFAULTVALUE,
												GETDEFAULTVALUE,
												K_ACCRUED,
												K_FOLLOWING)); //adjust 1st date

		/// Optimise cpn reset dates
		if(itsOptimResetData[OPTIM_RESET])
			OptimiseResetDates(NULL,itsStructDateStrip);
	}
	else
		itsStructDateStrip = itsDateStrip3;

	//OB: correction of notice dates, to be between ResetFundDate et StartFund
	//

	ARM_GP_Vector* fundResetDates = GetFundDateStrip()->GetResetDates() ;
	ARM_GP_Vector* fundStartDates = GetFundDateStrip()->GetFlowStartDates();
	ARM_GP_Vector* exerciseDates = GetExerDateStrip()->GetResetDates() ;
	ARM_GP_Vector* startDates = GetExerDateStrip()->GetFlowStartDates();
	size_t ExerSize = exerciseDates->size();

	int i = 0;
	size_t k;
	for (k=0; k<ExerSize; k++)
	{
		while ((*startDates)[k]>(*fundStartDates)[i]) 
			i++;
		if((*exerciseDates)[k]>(*fundResetDates)[i])
			(*exerciseDates)[k] = (*fundResetDates)[i];
	}

	
    return(SchedVect);
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: CreateCstManager
///	Returns: void
///	Action : create the const manager (static data).
/////////////////////////////////////////////////////////////////
ARM_CstManagerPtr ARM_CRASpreadCalculator::CreateCstManager()
{
	try
	{
		double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

		vector<string>	cstNames; 
		cstNames.push_back("FixRateCst");
		cstNames.push_back("CpnSpreadCst");
		cstNames.push_back("PayIdxMultCst");
		cstNames.push_back("Coef1Cst");
		cstNames.push_back("Coef2Cst");

		if(IsDbleCorridor())
		{
			cstNames.push_back("SDownCst");
			cstNames.push_back("SUpCst");
			cstNames.push_back("RDownCst");
			cstNames.push_back("RUpCst");
		}
		else
		{
			cstNames.push_back("DownCst");
			cstNames.push_back("UpCst");
			if (IsTripleRange())
			{
				cstNames.push_back("FixRate2Cst");
				cstNames.push_back("FixRate3Cst");
				cstNames.push_back("Down2Cst");
				cstNames.push_back("Down3Cst");
				cstNames.push_back("Up2Cst");
				cstNames.push_back("Up3Cst");
			}
		}

		cstNames.push_back("NotioCst");
		cstNames.push_back("FundSpreadCst");

//		if (SYNCHRONIZE && !IsDbleCorridor())
		if (SYNCHRONIZE)
		{
			cstNames.push_back("dsStruct");
			cstNames.push_back("dsPay");
			cstNames.push_back("dsFix");
			cstNames.push_back("idx");
		}

		if (IsVms())
			cstNames.push_back("TenorCst");
		
		vector<ARM_GramFctorArg> cstVector;
		
		//BoostedFixRateCst
		ARM_Curve* tmpBoosted;
		tmpBoosted = RefValueToCurve(GetBoostedFixRate(), asOfDate, ARM_Constants::rateBase);
		cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpBoosted))));
	
		//CpnSpreadCst (always set to 0)
		ARM_GP_Vector dates(1, 0.0);
		ARM_GP_Vector values(1, 0.0);
		ARM_Curve* tmpCpnSpread = new ARM_Curve(dates, values,new ARM_LinInterpCstExtrapolDble);
		cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpCpnSpread))));

		//CpnPayIndexMultCst
		ARM_Curve* tmpPayIndexMult = RefValueToCurve(GetPayIndexMult(), asOfDate);
		cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpPayIndexMult))));

		//CpnCoeff1Cst
		ARM_Curve* tmpCoeff1 = RefValueToCurve(GetRefCoeff1(), ARM_Constants::rateBase);
		cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpCoeff1))));
		
		//CpnCoeff2Cst
		ARM_Curve* tmpCoeff2 = RefValueToCurve(GetRefCoeff2(), asOfDate);
		cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpCoeff2))));
		
		//BarrierDownCst (Stepup left)
		ARM_Curve* tmpBarrierDown = RefValueToCurve(itsCpnBarrierDown, asOfDate, ARM_Constants::rateBase);
		cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpBarrierDown))));

		//BarrierUpCst (Stepup left)
		ARM_Curve* tmpBarrierUp = RefValueToCurve(itsCpnBarrierUp, asOfDate, ARM_Constants::rateBase);
		cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpBarrierUp))));

		if(IsDbleCorridor())
		{
			//RateBarrierDownCst (Stepup left)
			ARM_Curve* tmpRateBarrierDown = RefValueToCurve(itsRateBarrierDown, asOfDate, ARM_Constants::rateBase);
			cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpRateBarrierDown))));

			//RateBarrierUpCst (Stepup left)
			ARM_Curve* tmpRateBarrierUp = RefValueToCurve(itsRateBarrierUp, asOfDate, ARM_Constants::rateBase);
			cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpRateBarrierUp))));
		}
		else
		{
			if(IsTripleRange())
			{
				//BoostedFixRateCst
				ARM_Curve* tmpBoosted2;
				tmpBoosted2 = RefValueToCurve(itsBoostedFixRate2, asOfDate, ARM_Constants::rateBase);
				cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpBoosted2))));

				ARM_Curve* tmpBoosted3;
				tmpBoosted3 = RefValueToCurve(itsBoostedFixRate3, asOfDate, ARM_Constants::rateBase);
				cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpBoosted3))));
			
				//BarrierDownCst (Stepup left)
				ARM_Curve* tmpBarrierDown2 = RefValueToCurve(itsCpnBarrierDown2, asOfDate, ARM_Constants::rateBase);
				cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpBarrierDown2))));

				ARM_Curve* tmpBarrierDown3 = RefValueToCurve(itsCpnBarrierDown3, asOfDate, ARM_Constants::rateBase);
				cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpBarrierDown3))));

				//BarrierUpCst (Stepup left)
				ARM_Curve* tmpBarrierUp2 = RefValueToCurve(itsCpnBarrierUp2, asOfDate, ARM_Constants::rateBase);
				cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpBarrierUp2))));

				ARM_Curve* tmpBarrierUp3 = RefValueToCurve(itsCpnBarrierUp3, asOfDate, ARM_Constants::rateBase);
				cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpBarrierUp3))));
			}

		}

		//NotionalCst (Stepup right)
		ARM_Curve* tmpNotional = RefValueToCurve(itsNotional, asOfDate);
		cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpNotional))));

		//FundSpreadCst (Stepup left)
		ARM_Curve* tmpFundSpread = RefValueToCurve(itsFundSpread, asOfDate, ARM_Constants::rateBase);
		cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpFundSpread))));

//		if (SYNCHRONIZE && !IsDbleCorridor())
		if (SYNCHRONIZE)
		{
			cstVector.push_back(ARM_GramFctorArg(GetDateStrip1()));
			cstVector.push_back(ARM_GramFctorArg(GetDateStrip2()));
			cstVector.push_back(ARM_GramFctorArg(GetDateStrip3()));
			cstVector.push_back(ARM_GramFctorArg(GetIndexVector()));
		}

		if (IsVms())
		{
			ARM_Curve* tmpTenor = RefValueToCurve(itsTenor, asOfDate,1.);
			cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpTenor))));
		}

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

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: IsVariableCorridor
///	Returns: boolean
///	Action : checks if the payments are floating
///			 this is determined by pay index mults
/////////////////////////////////////////////////////////////////
bool ARM_CRASpreadCalculator::IsVariableCorridor() const
{
	/// check if corridor is variable
	/// corridor is variable if pay index mult is set to 0
	ARM_Vector* mults = itsPayIndexMult.GetDiscreteValues();
	size_t size = mults->size();
	bool isVariable = false;
	for (size_t i=0; i<size; i++)
	{
		if (fabs((*mults)[i])>1e-12)
		{
			isVariable = true;
			break;
		}
	}
	return isVariable;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: IsVariableCorridor
///	Returns: boolean
///	Action : checks if the payments are floating
///			 this is determined by pay index mults
/////////////////////////////////////////////////////////////////
bool ARM_CRASpreadCalculator::IsVariableCpn(size_t cpnIndex) const
{
	if (cpnIndex<0 || cpnIndex>=itsCpnSize)
		ARM_THROW( ERR_INVALID_ARGUMENT, "CRA Spread calculator::IsVariableCpn : invalid cpn index");

	if (itsPayIndex == K_FIXED) 
		return false;

	/// to be checked....
	/// a priori c'est sur la start de la periode de payment correspondante qu'il faut interpoler
	///
	double t = GetStructDateStrip()->GetFlowStartDates()->Elt(cpnIndex);
	if (fabs(((ARM_ReferenceValue)itsPayIndexMult).Interpolate(t))>1e-12) 
		return true;
	else
		return false ;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: GetIndexType
///	Returns: the index type
///	Action : get the ARM index type of the coupon index
/////////////////////////////////////////////////////////////////
ARM_INDEX_TYPE ARM_CRASpreadCalculator::GetIndexType()
{
	return GetCurrencyUnit()->GetVanillaIndexType();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: CreateSwaptionPortfolio
///	Returns: a portfolio
///	Action : create the list of diagonal swaptions
/////////////////////////////////////////////////////////////////
ARM_StdPortfolioPtr ARM_CRASpreadCalculator::CreateSwaptionPortfolio()
{	
	ARM_GP_Vector* resetDates = itsExerDateStrip->GetResetDates();
	ARM_GP_Vector* startDates = itsExerDateStrip->GetFlowStartDates();

	ARM_BSModel* CFBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );

	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

	int i, nbFlows;

	list< ARM_Security* > swaptionList;

	nbFlows = startDates->size();

	if (	itsCalibType == BASKET_CALIBRATION || 
			itsCalibType == BASKET_CALIBRATION_SIMPLIFIED ||
			(itsCalibStrikeType[0] == EQUIVALENT &&
				(itsCalibType == DIAG_CALIBRATION ||
				 itsCalibType == DIAG_SPREAD_LONG ||
				 itsCalibType == DIAG_SPREAD_INDEX||
				 itsCalibType == DIAG_SPREAD_SHORT))	)
	{
		/// ComputeCSOPrices(); already computed in Init()
		if(IsDbleCorridor())
			ComputeCRA2NormalSensitivities();
		else
			ComputeCSONormalSensitivities();
	}

	// End Date of all the calibrated swaptions
	ARM_Date endDate = GetEndDate();

	double fees;
	
	for (i = 0; i < nbFlows; ++i)
	{
		ARM_Date startDate ((*startDates)[i]);

		fees = itsCallFees.Interpolate((*resetDates)[i]);

		if (((*resetDates)[i] > asOfDate) && (fees < NON_CALL_FEE))
		{
			if (itsCalibType == DIAG_CALIBRATION || itsCalibType == DIAG_SPREAD_LONG || itsCalibType == DIAG_SPREAD_SHORT || itsCalibType == DIAG_SPREAD_INDEX)
			{
				ARM_Swap stdSwap(startDate,
								endDate,
								GetIndexType(),
								0.0,
								1.0,
								GetPayRec(),
								K_DEF_FREQ,
								K_DEF_FREQ,
								GetCurrencyUnit());

				ARM_Date expiryDate((*(stdSwap.GetFloatLeg()->GetResetDates()))[0]);
				stdSwap.SetModel(CFBSModel);

				double residualPv=0.0;
				if(itsCalibStrikeType[0] == EQUIVALENT)
				{
					double floatLegPv,soLegPv;
					residualPv = ComputeResidualUnderlying(i,floatLegPv,soLegPv);

					ARM_SwapLeg* stdFloatLeg = stdSwap.Get1stLeg();
					if(stdFloatLeg->GetLegType() == K_FIXED_LEG)
						stdFloatLeg = stdSwap.Get2ndLeg();
					stdFloatLeg->SetModel(CFBSModel);
					double stdFloatNotio = floatLegPv / fabs(stdFloatLeg->ComputePrice());
					residualPv = GetPayRec() * (soLegPv-floatLegPv)/stdFloatNotio;
				}
				double equivStrike = stdSwap.PriceToRate((ARM_Date) asOfDate, residualPv);
					

				int RecOrPay = K_RCV;
				ARM_Swaption swaption((&stdSwap),RecOrPay,K_EUROPEAN,equivStrike,expiryDate);
					
				swaptionList.push_back(static_cast<ARM_Swaption*>(swaption.Clone()));
			}
			else if (itsCalibType == BASKET_CALIBRATION || itsCalibType == BASKET_CALIBRATION_SIMPLIFIED)
			{
				ARM_Swaption* swaption = NULL;

				swaption = CreateVarNotionalSwaptionAtExer(i);
				
				if (swaption)
					swaptionList.push_back(swaption);
			}
			else
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Bad calibration type." );
			}
		}
	}

	itsSwoptCalib = (swaptionList.size() != 0);

	ARM_StdPortfolioPtr swaptionPF(new ARM_StdPortfolio(swaptionList));
    	
	for(i=0;i<swaptionPF->size();++i)
	{
		swaptionPF->SetPrecision(DEFAULT_PRECISION,i);
        swaptionPF->SetWeight(DEFAULT_WEIGHT,i);
		swaptionPF->SetPrice(DEFAULT_PRICE,i);
	}
	return swaptionPF;

}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: CreateSOPortfolio
///	Returns: a portfolio
///	Action : create the list of CMS spread options. Only one
///			 option is created for each call date
/////////////////////////////////////////////////////////////////
ARM_StdPortfolioPtr ARM_CRASpreadCalculator::CreateSOPortfolio(bool isSOCond,bool isDouble,bool isLong)
{
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

	ARM_GP_Vector* resetDates = itsExerDateStrip->GetResetDates();
	ARM_GP_Vector* startDates = itsExerDateStripUnadj->GetFlowStartDates();
	ARM_GP_Vector* endDates = itsExerDateStripUnadj->GetFlowEndDates();
	size_t i,nbFlows = startDates->size();

	list< ARM_Security* > soList;
	ARM_SpreadOption* soSec;

	ARM_ReferenceValue strikes(0.0); // strike will be computed later
	ARM_ReferenceValue weight1(itsCoeff1);
	ARM_ReferenceValue weight2(itsCoeff2);
	ARM_INDEX_TYPE index1 = (ARM_INDEX_TYPE)itsRefIndex1;
	ARM_INDEX_TYPE index2 = (ARM_INDEX_TYPE)itsRefIndex2;
	ARM_ReferenceValue* unUsedFixing=NULL;
	int defaultResetGap = 10000;
	if(!isSOCond)
	{
		/// Second index is fictive... a CMS caplet list will be created !
		weight2=ARM_ReferenceValue(0.0);
		if (isDouble)
		{
			if (itsCalibType ==  DIAG_SPREAD_INDEX)
				index1 = (ARM_INDEX_TYPE)itsPayIndex;
			else
			{
				index1 = (ARM_INDEX_TYPE)itsRefIndex3;
				if(itsRefIndex1 != itsRefIndex3)
					index2 = (ARM_INDEX_TYPE)itsRefIndex1;
			}
		}
		else
		{
			if (isLong)
			{
				index2 = (ARM_INDEX_TYPE)itsRefIndex2;
				index1 = (ARM_INDEX_TYPE)itsRefIndex1;
			}
			else
			{
				index2 = (ARM_INDEX_TYPE)itsRefIndex1;
				index1 = (ARM_INDEX_TYPE)itsRefIndex2;
			}
		}
	}
	//int cpnPayFreqNbMonths = (int)(12*itsCpnPayFreq);
	for(i=0;i<nbFlows;++i)
	{
		double fees = itsCallFees.Interpolate((*resetDates)[i]);

		if (((*resetDates)[i] > asOfDate) && (fees < NON_CALL_FEE))
		{
			/// Compute unadjusted flow end date
			ARM_Date startDate((*startDates)[i]);
			ARM_Date unadjEndDate((*endDates)[i]);
			/*ARM_Date unadjEndDate((*startDates)[i]);
			unadjEndDate.AddMonths(cpnPayFreqNbMonths);*/

			/// Create the CMS spread option with a single flow at call frequency
			soSec = new ARM_SpreadOption(startDate, unadjEndDate, K_CAP, &strikes,
				index2,index1,&weight2,&weight1,
				itsBoostedDayCount, itsCallFreq, itsCallFreq, K_ADVANCE, K_ARREARS, 
				&itsCcy,defaultResetGap,unUsedFixing, unUsedFixing, itsCpnIntRule, K_LONGSTART, 1,
				const_cast<char*>(itsCpnResetCal.c_str()),const_cast<char*>(itsCpnPayCal.c_str()),1,itsCpnAdjRule);

			soList.push_back(soSec);
		}
	}

	ARM_StdPortfolioPtr soPF(new ARM_StdPortfolio(soList));
    	
	for(i=0;i<soPF->size();++i)
	{
		soPF->SetPrecision(DEFAULT_PRECISION,i);
        soPF->SetWeight(DEFAULT_WEIGHT,i);
		soPF->SetPrice(DEFAULT_PRICE,i);
	}
	return soPF;
}

////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: GetSOPF
///	Returns: ARM_StfPortfolioPtr
///	Action : Return the SO portfolio associated to the CMS spread
///			 or the rate condition of a double condition CRA
/////////////////////////////////////////////////////////////////
ARM_StdPortfolioPtr ARM_CRASpreadCalculator::GetSOPF(bool isSOCond) const
{
	if (IsDbleCorridor())
	{
		if(GetCalibMethod() != ARM_CalibMethodPtr(NULL))
		{
			if(isSOCond && GetCalibMethod()->GetlinkedMethod())
				return GetCalibMethod()->GetlinkedMethod()->GetPortfolio();
			else if(!isSOCond && GetCalibMethod()->GetlinkedMethod()->GetlinkedMethod())
				return GetCalibMethod()->GetlinkedMethod()->GetlinkedMethod()->GetPortfolio();
			else
				return ARM_StdPortfolioPtr(NULL);
		}
		else
			return ARM_StdPortfolioPtr(NULL);
	}
	else if (itsCalibType==DIAG_SPREAD_LONG || itsCalibType==DIAG_SPREAD_SHORT || itsCalibType==DIAG_SPREAD_INDEX)
	{
		return GetCalibMethod()->GetlinkedMethod()->GetPortfolio();
	}
	else if (itsCalibType==SHORT_LONG_SPREAD)
	{
		return GetCalibMethod()->GetlinkedMethod()->GetlinkedMethod()->GetPortfolio();
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, " ARM_CRASpreadCalculator:: no SO portfolio available");

}

////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: GetCMSLONGPortolio()
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_StdPortfolioPtr ARM_CRASpreadCalculator::GetCMSLONGPortfolio() const
{ 

	if (itsCalibType==DIAG_SPREAD_LONG)
	{
		return GetCalibMethod()->GetlinkedMethod()->GetlinkedMethod()->GetPortfolio();
	}
	else if (itsCalibType==SHORT_LONG_SPREAD)
	{
		return GetCalibMethod()->GetlinkedMethod()->GetPortfolio();
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, " ARM_LocalCSOCalculator:: no CMSLONG portfolio available");

    
	return GetCalibMethod()->GetPortfolio();
}

////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: GetCMSLONGPortolio()
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_StdPortfolioPtr ARM_CRASpreadCalculator::GetCMSSHORTPortfolio() const
{ 

	if (itsCalibType==DIAG_SPREAD_SHORT || itsCalibType==DIAG_SPREAD_INDEX)
	{
		return GetCalibMethod()->GetlinkedMethod()->GetlinkedMethod()->GetPortfolio();
	}
	else if (itsCalibType==SHORT_LONG_SPREAD)
	{
		return GetCalibMethod()->GetPortfolio();
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, " ARM_LocalCSOCalculator:: no CMSLONG portfolio available");

    
	return GetCalibMethod()->GetPortfolio();
}

////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: ComputeSOPrices
///	Returns: nothing
///	Action : compute market target prices of the SO portfolio
///			 (CSM Spread or CMS caplet)
/////////////////////////////////////////////////////////////////
void ARM_CRASpreadCalculator::ComputeSOPrices(ARM_StdPortfolioPtr& pf,CalibStrikeType cst,bool isSOBarrier)
{
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

	ARM_BSModel* SOBSModel;
	if (itsSwitchOldCalibration)
		SOBSModel = dynamic_cast< ARM_BSModel* >(GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );
	else
		SOBSModel = dynamic_cast< ARM_BSModel* >(GetMktDataManager()->GetData(GetKeys()[OswModelKey]) );

	ARM_SpreadOption* spreadOption;
	double price,strike;
	
	double soFwdRate1,soFwdRate2,weight1,weight2,soFwd;
	ARM_ReferenceValue strikes;

	for (size_t i = 0; i < pf->GetSize(); ++i)
	{
		spreadOption = static_cast< ARM_SpreadOption* >(pf->GetAsset(i));
		if(spreadOption->GetNumFlows()!=1)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : SO in calibration portfolio must have a single flow" );
		spreadOption->SetModel(SOBSModel);

		/// Compute ATM strike
		price=spreadOption->ComputePrice();
		soFwdRate1 = (*(spreadOption->GetSpreadLeg()->GetFirstLeg()->GetFwdRates()))[0];
		soFwdRate2 = (*(spreadOption->GetSpreadLeg()->GetSecondLeg()->GetFwdRates()))[0];
		weight1 = spreadOption->GetWeight1();
		weight2 = spreadOption->GetWeight2();
		soFwd = weight2*soFwdRate2 - weight1*soFwdRate1;

		if (cst==ZERO)
			strike=0.;
		else if (cst==ATM)
			strike=soFwd;
		else if (cst==CAP)
		{
			if(isSOBarrier)
				strike=GetCpnBarrierUp().Interpolate((*(spreadOption->GetSpreadLeg()->GetFirstLeg()->GetResetDates()))[0]-asOfDate);
			else
				strike=itsRateBarrierDown.Interpolate((*(spreadOption->GetSpreadLeg()->GetFirstLeg()->GetResetDates()))[0]-asOfDate);
		}
		else if (cst==FLOOR)
		{
			if(isSOBarrier)
				strike=GetCpnBarrierDown().Interpolate((*(spreadOption->GetSpreadLeg()->GetFirstLeg()->GetResetDates()))[0]-asOfDate);
			else
				strike=itsRateBarrierUp.Interpolate((*(spreadOption->GetSpreadLeg()->GetFirstLeg()->GetResetDates()))[0]-asOfDate);
		}
		else
			strike=0.;

		spreadOption->SetStrike(strike);
		strikes = ARM_ReferenceValue(strike);
		spreadOption->SetStrikes(&strikes);

		/// Compute & save price
		price=spreadOption->ComputePrice();
		pf->SetPrice(price, i);
	}
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: ComputeProductVectorsFromCurves
///	Returns: void
///	Action : compute vector attributes from ARM_Curve 
/////////////////////////////////////////////////////////////////
void ARM_CRASpreadCalculator::ComputeProductVectorsFromCurves()
{
	size_t i;
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

	/// Convert ARM_Curves into relevant vectors
	/// this will be much easier for computations
	ARM_GP_Vector* fundResetDates = GetFundDateStrip()->GetResetDates() ;
	itsFundSize = fundResetDates->size();

	ARM_GP_Vector* cpnResetDates = GetStructDateStrip()->GetResetDates() ;
	itsCpnSize = cpnResetDates->size();

	
	ARM_GP_Vector* exerciseDates = GetExerDateStrip()->GetResetDates() ;
	ARM_GP_Vector* startDates = GetExerDateStrip()->GetFlowStartDates();
	itsExerSize = exerciseDates->size();

	/// past  no call
	size_t nbPastNoCall=0;
	while(nbPastNoCall < itsExerSize && (*exerciseDates)[nbPastNoCall] < asOfDate)
		++nbPastNoCall;

	for (i=0; i<itsExerSize; i++)
	{
		double lag = (*exerciseDates)[i] - asOfDate;
		double fee =  i < nbPastNoCall ? NON_CALL_FEE : itsCallFees.Interpolate((*exerciseDates)[i]);
		itsvExerFees.push_back(fee);

		( fabs((fee-NON_CALL_FEE))<K_NEW_DOUBLE_TOL )? itsvIsExerDate.push_back(false) : itsvIsExerDate.push_back(true);
	}

	/// last: compute indexes to relate exer - funding - cpn
	ARM_GP_Vector* cpnStartDates  = GetStructDateStrip()->GetFlowStartDates();
	ARM_GP_Vector* cpnPayDates  = GetStructDateStrip()->GetPaymentDates();
	double cpnEndDate = GetStructDateStrip()->GetFlowEndDates()->Elt(itsCpnSize-1);
	ARM_GP_Vector* fundStartDates = GetFundDateStrip()->GetFlowStartDates();
	ARM_GP_Vector* fundPayDates = GetFundDateStrip()->GetPaymentDates();


	for (i=0; i<itsCpnSize; i++)
	{
		itsvCpnNominal.push_back(itsNotional.Interpolate(cpnPayDates->Elt(i)));
	}

	for (i=0; i<itsFundSize; i++)
	{
		itsvFundNominal.push_back(itsNotional.Interpolate(fundPayDates->Elt(i)));
		itsvFundSpread.push_back(itsFundSpread.Interpolate(fundStartDates->Elt(i))/ARM_Constants::rateBase);
		if(IsBasis()){
			itsvRealFundSpread.push_back(itsFundSpread.Interpolate(fundStartDates->Elt(i))/ARM_Constants::rateBase);
			itsvRealFundNominal.push_back(itsRealFundNotional.Interpolate(fundPayDates->Elt(i)));
		}
	}

	if(IsBasis())
		itsvFundSpread = ComputeDomesticBasis();
	
	i = 0;
	size_t k;
	for (k=0; k<itsExerSize; k++)
	{
		while ((*startDates)[k]>(*fundStartDates)[i]) 
			i++;
		if((*exerciseDates)[k]>(*fundResetDates)[i])
		{
			char msg[500];
			sprintf(msg, "CRA Spread calculator : call date[%d] (%.0f) between fundingReset[%d] (%.0f) and fundingStart[%d] (%.0f). Not allowed !", 
					k, (*exerciseDates)[k], i, (*fundResetDates)[i], i, (*fundStartDates)[i]);
			ARM_THROW( ERR_INVALID_ARGUMENT, msg );
		}
		itsvFundIndex.push_back(i);
	}
	i = 0;
	for (k=0; k<itsExerSize; k++)
	{
		while ((*startDates)[k]>(*cpnStartDates)[i]) 
			i++;
		if((*exerciseDates)[k]>(*cpnResetDates)[i])
		{
			char msg[500];
			sprintf(msg, "CRA Spread calculator : call date[%d] (%.0f) between couponReset[%d] (%.0f) and couponStart[%d] (%.0f). Not allowed !", 
					k, (*exerciseDates)[k], i, (*cpnResetDates)[i], i, (*cpnStartDates)[i]);
			ARM_THROW( ERR_INVALID_ARGUMENT, msg );
		}
		itsvCpnIndex.push_back(i);
	}
	itsvFundIndex.push_back(fundStartDates->size());
	itsvCpnIndex.push_back(cpnStartDates->size());
		
	/* 2 call dates cannot belong to the same period */
	/* this may happen when freq (call) > freq (exotic) or freq (funding) */
	for (k=0; k<itsExerSize; k++)
	{
		if (itsvFundIndex[k]>=itsvFundIndex[k+1])
		{
			char msg[500];
			sprintf(msg, "CRA Spread calculator : call date[%d] and [%d] found within the same funding leg period[%d]. Freq (call) > Freq (Funding) is forbidden", k, k+1, k);
			ARM_THROW( ERR_INVALID_ARGUMENT, msg);
		}
		if (itsvCpnIndex[k]>=itsvCpnIndex[k+1])
		{
			char msg[500];
			sprintf(msg, "CRA Spread calculator : call date[%d] and [%d] found within the same coupon leg period[%d]. Freq (call) > Freq (Coupon) is forbidden", k, k+1, k);
			ARM_THROW( ERR_INVALID_ARGUMENT, msg);
		}
	}

	/* after each call call date, we must have same start dates for exo and funding legs */
	for (k=0; k<itsExerSize; k++) 
	{		
		//if (fabs((*fundStartDates)[itsvFundIndex[k]]-(*cpnStartDates)[itsvCpnIndex[k]])>0)
		if (fabs((*fundStartDates)[itsvFundIndex[k]]-(*startDates)[k])>0)
			ARM_THROW( ERR_INVALID_ARGUMENT, "CSO calculator : call date  defines a swap with mismatching exo and floating start dates (no broken period allowed)");
	}
	
}

////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: GetSwoptPF
///	Returns: ARM_StfPortfolioPtr
///	Action : Return the swaption portfolio
/////////////////////////////////////////////////////////////////
ARM_StdPortfolioPtr ARM_CRASpreadCalculator::GetSwoptPF() const
{
	if(GetCalibMethod() != ARM_CalibMethodPtr(NULL))
		return GetCalibMethod()->GetPortfolio();
	else
		return ARM_StdPortfolioPtr(NULL);
}

////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: ComputeSwaptionPrices
///	Returns: nothing
///	Action : compute market target prices of the swaption portfolio
/////////////////////////////////////////////////////////////////
void ARM_CRASpreadCalculator::ComputeSwaptionPrices()
{
    double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
//	ARM_BSModel* oswBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswModelKey]) );
	ARM_ZeroCurve* pCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));
	ARM_CurveModelParam* mrsParam = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(GetKeys()[MrsKey]));
	double mrs = mrsParam->GetCurve()->Interpolate(0);
	double price,vega,volratio,correl,mrs_spread;

	if (itsModelType == ARM_PricingModelType::HWM2F)
	{
		ARM_CurveModelParam* volRatioParam = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(GetKeys()[VolRatioKey]));
		ARM_CurveModelParam* mrsSpreadParam = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(GetKeys()[MrsSpreadKey]));
		ARM_CurveModelParam* correlParam = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(GetKeys()[CorrelKey]));
		volratio = volRatioParam->GetCurve()->Interpolate(0);
		correl = correlParam->GetCurve()->Interpolate(0);
		mrs_spread = mrsSpreadParam->GetCurve()->Interpolate(0);
	}

	ARM_StdPortfolioPtr pf = GetCalibMethod()->GetPortfolio();
	for(int i=0;i<pf->GetSize();++i)
	{
		if (itsCalibType == DIAG_CALIBRATION || itsCalibType == DIAG_SPREAD_LONG || itsCalibType == DIAG_SPREAD_SHORT || itsCalibType == DIAG_SPREAD_INDEX)
		{
			ARM_BSModel* modelGen = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswModelKey]) );
			ARM_BSModel* oswBSModel = itsSwitchOldCalibration?modelGen:modelGen->GetSabrModel();

			ARM_Swaption* swaption = static_cast< ARM_Swaption* >(pf->GetAsset(i));
			swaption->SetModel(oswBSModel);
			price = swaption->ComputePrice();
		}
		else if (itsCalibType == BASKET_CALIBRATION || itsCalibType == BASKET_CALIBRATION_SIMPLIFIED)
		{
			ARM_BSModel* modelGen = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswModelKey]) );
			ARM_BSModel* oswBSModel = modelGen;

			ARM_MarketIRModel analyticModel(*GetMktDataManager(), GetKeys(), itsVnsPricingMethod);
			analyticModel.SetZeroCurveKey(GetKeys()[YcKey]);
			analyticModel.SetBsModelKey(GetKeys()[OswModelKey]);
			ARM_VanillaSwaptionArg* swaptionGP = static_cast<ARM_VanillaSwaptionArg*>(ARM_ConverterFromKernel::ConvertSecuritytoArgObject(pf->GetAsset(i), asOfDate, GetKeys()[YcKey]));
			price = swaptionGP->Price(&analyticModel);

			ARM_PricingModel* pModel;
			ARM_PricingModel* pModelShifted;

			if (itsModelType == ARM_PricingModelType::HWM1F)
			{
				//mrs = MRS_DEFAULT_VALUE_FOR_VEGA;		//pour ne pas impacter la sensi meanrev
				ARM_GP_Vector defaultTimes(1,0.);
				ARM_GP_Vector defaultSigmas(1,SIGMA_DEFAULT_VALUE_FOR_VEGA);
				ARM_GP_Vector defaultSigmasShifted(1,SIGMA_DEFAULT_VALUE_FOR_VEGA+SIGMA_SHIFT_FOR_VEGA);
				ARM_GP_Vector defaultMeanRevs(1,mrs);
				ARM_CurveModelParam* volParam = new ARM_CurveModelParam(ARM_ModelParamType::Volatility, &defaultSigmas, &defaultTimes );
				ARM_CurveModelParam* volParamShifted = new ARM_CurveModelParam(ARM_ModelParamType::Volatility, &defaultSigmasShifted, &defaultTimes );
				ARM_CurveModelParam* mrsParam = new ARM_CurveModelParam(ARM_ModelParamType::MeanReversion, &defaultMeanRevs, &defaultTimes );
				ARM_CurveModelParam* mrsParamShifted = new ARM_CurveModelParam(ARM_ModelParamType::MeanReversion, &defaultMeanRevs, &defaultTimes );

				ARM_ModelParamVector modelParams(2);
				modelParams[0] = volParam;
				modelParams[1] = mrsParam;

				ARM_ModelParamVector modelParamsShifted(2);
				modelParamsShifted[0] = volParamShifted;
				modelParamsShifted[1] = mrsParamShifted;
	   
				pModel = new ARM_HullWhite1F( CreateClonedPtr( pCurve ) );
				pModelShifted = new ARM_HullWhite1F( CreateClonedPtr( pCurve ) );

				pModel->SetModelParams(ARM_ModelParamsHW1FStd(modelParams));
				pModelShifted->SetModelParams(ARM_ModelParamsHW1FStd(modelParamsShifted));

				vega = (swaptionGP->Price(pModelShifted)-swaptionGP->Price(pModel))/SIGMA_SHIFT_FOR_VEGA;

				delete volParam;
				delete volParamShifted;
				delete mrsParam;
				delete mrsParamShifted;

				if (fabs(vega) < VEGA_MIN_TO_SELECT)
					pf->SetWeight (0.0, i);

			}
			else
			{
				ARM_GP_Vector defaultTimes(1,0.);
				ARM_GP_Vector defaultSigmas(1,SIGMA_DEFAULT_VALUE_FOR_VEGA);
				ARM_GP_Vector defaultSigmasShifted(1,SIGMA_DEFAULT_VALUE_FOR_VEGA+SIGMA_SHIFT_FOR_VEGA);
				ARM_GP_Vector defaultMeanRevs(1,mrs);
				ARM_GP_Vector defaultMRspread(1,mrs_spread);
				ARM_GP_Vector defaultCorrel(1,correl);
				ARM_GP_Vector defaultVolRatio(1,volratio);

				ARM_CurveModelParam* volParam = new ARM_CurveModelParam(ARM_ModelParamType::Volatility, &defaultSigmas, &defaultTimes );
				ARM_CurveModelParam* volParamShifted = new ARM_CurveModelParam(ARM_ModelParamType::Volatility, &defaultSigmasShifted, &defaultTimes );
				ARM_CurveModelParam* mrsParam = new ARM_CurveModelParam(ARM_ModelParamType::MeanReversion, &defaultMeanRevs, &defaultTimes );
				ARM_CurveModelParam* mrsParamShifted = new ARM_CurveModelParam(ARM_ModelParamType::MeanReversion, &defaultMeanRevs, &defaultTimes );
				ARM_CurveModelParam* mrssParam = new ARM_CurveModelParam(ARM_ModelParamType::MeanReversionSpread, &defaultMRspread, &defaultTimes );
				ARM_CurveModelParam* mrssParamShifted = new ARM_CurveModelParam(ARM_ModelParamType::MeanReversionSpread, &defaultMRspread, &defaultTimes );
				ARM_CurveModelParam* correlParam = new ARM_CurveModelParam(ARM_ModelParamType::Correlation, &defaultCorrel, &defaultTimes );
				ARM_CurveModelParam* correlParamShifted = new ARM_CurveModelParam(ARM_ModelParamType::Correlation, &defaultCorrel, &defaultTimes );
				ARM_CurveModelParam* volratioParam = new ARM_CurveModelParam(ARM_ModelParamType::VolatilityRatio, &defaultVolRatio, &defaultTimes );
				ARM_CurveModelParam* volratioParamShifted = new ARM_CurveModelParam(ARM_ModelParamType::VolatilityRatio, &defaultVolRatio, &defaultTimes );

				ARM_ModelParamVector modelParams(5);
				modelParams[0] = volParam;
				modelParams[1] = mrsParam;
				modelParams[2] = mrssParam;
				modelParams[3] = correlParam;
				modelParams[4] = volratioParam;

				ARM_ModelParamVector modelParamsShifted(5);
				modelParamsShifted[0] = volParamShifted;
				modelParamsShifted[1] = mrsParamShifted;
				modelParamsShifted[2] = mrssParamShifted;
				modelParamsShifted[3] = correlParamShifted;
				modelParamsShifted[4] = volratioParamShifted;

				pModel = new ARM_HullWhite2F( CreateClonedPtr( pCurve ) );
				pModelShifted = new ARM_HullWhite2F( CreateClonedPtr( pCurve ) );

				pModel->SetModelParams(ARM_ModelParamsHW2FStd(modelParams));
				pModelShifted->SetModelParams(ARM_ModelParamsHW2FStd(modelParamsShifted));

				vega = (swaptionGP->Price(pModelShifted)-swaptionGP->Price(pModel))/SIGMA_SHIFT_FOR_VEGA;

				delete volParam;
				delete volParamShifted;
				delete mrsParam;
				delete mrsParamShifted;
				delete mrssParam;
				delete mrssParamShifted;
				delete correlParam;
				delete correlParamShifted;
				delete volratioParam;
				delete volratioParamShifted;

				if (fabs(vega) < VEGA_MIN_TO_SELECT)
					pf->SetWeight (0.0, i);
			}

			delete pModel;
			delete pModelShifted;
			delete swaptionGP;
			
		}

		pf->SetPrice (price, i);
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: SkipResetDates
///	Returns: void
///	Action : Skip reset dates w.r.t. to required size and limit
/////////////////////////////////////////////////////////////////
void ARM_CRASpreadCalculator::SkipResetDates(const ARM_GP_Vector& resets,const ARM_GP_Vector& payments,
											double refDate, double skipSize, double limit,
											size_t &resetIdx, ARM_IntVector& resetsToErase) const
{
	size_t nbResets = resets.size();
	if(payments.size() != nbResets)
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : schedule mismatch in reset dates optimisation" );
	}
	double lastReset,lastPay;
	while(resetIdx < nbResets && resets[resetIdx]-refDate < limit)
	{
		if(resetIdx>0 && payments[resetIdx-1] == payments[resetIdx])
		{
			lastReset = resets[resetIdx];
			lastPay = payments[resetIdx];
			while(resetIdx < nbResets && resets[resetIdx]-lastReset < skipSize)
			{
				if(lastPay == payments[resetIdx])
					resetsToErase.push_back(resetIdx);
				else
				{
					/// Keep the first reset of a new payment flow
					/// This allows at least one reset per payment
					lastReset	= resets[resetIdx];
					lastPay		= payments[resetIdx];
				}

				++resetIdx;
			}
		}

		/// Keep the current reset
		++resetIdx;
	}
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: OptimiseResetDates
///	Returns: void
///	Action : Given an input spread corridor leg, it reduces reset dates
///			 to optimise further computational evaluation time
///			 in local models. Input date vectors of "sec" are modified
///			 accordingly
/////////////////////////////////////////////////////////////////
void ARM_CRASpreadCalculator::OptimiseResetDates(ARM_SpreadOption* sec,ARM_DateStripPtr& sched) const
{
	double weekLen		= 7.0;
	double monthLen		= OPTIM_RESET_AVE_YEARLEN/12.0;
	double quaterLen	= OPTIM_RESET_AVE_YEARLEN/4.0;

	ARM_GP_Vector resets(sec ? To_ARM_GP_Vector(*(sec->GetSwapLeg()->GetResetDates())) : *(sched->GetResetDates()));
	size_t nbResets = resets.size();
	if(nbResets <= 1)
		return;

	/// Reference date for reset optimisation is the very first reset date
	double refDate = resets[0];

	ARM_GP_Vector payments(sec ?  To_ARM_GP_Vector(*(sec->GetSwapLeg()->GetPaymentDates())) : *(sched->GetPaymentDates()));
	size_t resetIdx = 0;
	double lastPay = payments[resetIdx];

	ARM_IntVector resetsToErase(0);
	if(itsCpnResetFreq==K_DAILY)
	{
		/// Keep daily reset up to "dailyLimit"
		while(resetIdx < nbResets && resets[resetIdx]-refDate <= itsOptimResetData[OPTIM_RESET_DAILYLIMIT])
			++resetIdx;

		/// Switch to weekly then monthly then quaterly reset up to proper limit
		SkipResetDates(resets,payments,refDate,weekLen,itsOptimResetData[OPTIM_RESET_WEEKLYLIMIT],resetIdx,resetsToErase);
		SkipResetDates(resets,payments,refDate,monthLen,itsOptimResetData[OPTIM_RESET_MONTHLYLIMIT],resetIdx,resetsToErase);
		SkipResetDates(resets,payments,refDate,quaterLen,ARM_NumericConstants::ARM_INFINITY,resetIdx,resetsToErase);

	} // Daily reset

	else if(itsCpnResetFreq==K_WEEKLY)
	{
		/// Keep weekly reset up to "weeklyLimit"
		while(resetIdx < nbResets && resets[resetIdx]-refDate <= itsOptimResetData[OPTIM_RESET_WEEKLYLIMIT])
			++resetIdx;

		/// Switch to monthly then quaterly reset up to proper limit
		SkipResetDates(resets,payments,refDate,monthLen,itsOptimResetData[OPTIM_RESET_MONTHLYLIMIT],resetIdx,resetsToErase);
		SkipResetDates(resets,payments,refDate,quaterLen,ARM_NumericConstants::ARM_INFINITY,resetIdx,resetsToErase);

	} // Weekly reset

	else if(itsCpnResetFreq==K_MONTHLY || itsCpnResetFreq==K_BIMONTHLY)
	{
		/// Keep (bi)monthly reset up to "monthlyLimit"
		while(resetIdx < nbResets && resets[resetIdx]-refDate <= itsOptimResetData[OPTIM_RESET_MONTHLYLIMIT])
			++resetIdx;

		/// Switch to quaterly reset beyond "monthlyLimit"
		SkipResetDates(resets,payments,refDate,quaterLen,ARM_NumericConstants::ARM_INFINITY,resetIdx,resetsToErase);

	} // (Bi)montly reset


	size_t nbResetsToErase = resetsToErase.size();
	size_t nbNewResets = nbResets - nbResetsToErase;

	resetsToErase.push_back(nbResets+1); // +1 to copy non erased resets after last erased one

	size_t i,newResetIdx=0,lastNewResetIdx=0;
	ARM_IntVector newToOldReset(nbNewResets);
	ARM_IntVector newToOldFlowEndDate(nbNewResets);
	resetIdx=0;
	for(i=0;i<nbResetsToErase+1;++i)
	{
		while(resetIdx < nbResets && resetIdx < resetsToErase[i])
		{
			newToOldReset[newResetIdx]		= resetIdx;
			newToOldFlowEndDate[newResetIdx]= resetIdx;

			lastNewResetIdx = newResetIdx;
			++newResetIdx;
			++resetIdx;
		}
		if(resetIdx == resetsToErase[i])
		{
			/// Here is a reset to erase => save flow end date of previous non erased reset
			/// for furher update
			newToOldFlowEndDate[lastNewResetIdx] = resetIdx;
			++resetIdx;
		}
	}

	if(sec)
	{
		UpdateLegSchedules(newToOldReset,newToOldFlowEndDate,sec->GetSwapLeg());
		UpdateLegSchedules(newToOldReset,newToOldFlowEndDate,sec->GetSpreadLeg()->GetFirstLeg());
		UpdateLegSchedules(newToOldReset,newToOldFlowEndDate,sec->GetSpreadLeg()->GetSecondLeg());
		sec->SetNumFlows(nbNewResets);
	}
	else
		UpdateLegSchedules(newToOldReset,newToOldFlowEndDate,NULL,sched);
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: UpdateLegSchedules
///	Returns: void
///	Action : Resize and fill schedules of a swapleg or a DateStrip
///			 according to resets to keep
/////////////////////////////////////////////////////////////////
void ARM_CRASpreadCalculator::UpdateLegSchedules(const ARM_IntVector& newToOldReset,
												 const ARM_IntVector& newToOldFlowEndDate,
												 ARM_SwapLeg* leg,ARM_DateStripPtr& sched) const
{
	size_t i,nbNewResets = newToOldReset.size();

	bool isResetDates,isPaymentDates,isTheoPayDates;
	bool isFwdRateStartDates,isFwdRateEndDates,isFlowStartDates,isFlowEndDates;
	bool isInterestDays,isInterestTerms,isRateMaturities,isFwdRates;

	ARM_GP_Vector newResets,newPayments,newFwdRateStarts,newFwdRateEnds;
	ARM_GP_Vector newFlowStarts,newFlowEnds,newInterestDays,newInterestTerms;
	ARM_Vector newTheoPays,newRateMaturities,newFwdRates;

	bool isSched = sched != ARM_DateStripPtr(NULL);

	if(isResetDates = ((leg && leg->GetResetDates()) || (isSched && sched->GetResetDates())))
		newResets.resize(nbNewResets);
	if(isPaymentDates = ((leg && leg->GetPaymentDates()) || (isSched && sched->GetPaymentDates())))
		newPayments.resize(nbNewResets);
	if(isTheoPayDates = (leg && leg->GetTheoPayDates()))
		newTheoPays.Resize(nbNewResets);
	if(isFwdRateStartDates = ((leg && leg->GetFwdRateStartDates()) || (isSched && sched->GetFwdRateStartDates())))
		newFwdRateStarts.resize(nbNewResets);
	if(isFwdRateEndDates = ((leg && leg->GetFwdRateEndDates()) || (isSched && sched->GetFwdRateEndDates())))
		newFwdRateEnds.resize(nbNewResets);
	if(isFlowStartDates = ((leg && leg->GetFlowStartDates()) || (isSched && sched->GetFlowStartDates())))
		newFlowStarts.resize(nbNewResets);
	if(isFlowEndDates = ((leg && leg->GetFlowEndDates()) || (isSched && sched->GetFlowEndDates())))
		newFlowEnds.resize(nbNewResets);
	if(isInterestDays = ((leg && leg->GetInterestDays()) || (isSched && sched->GetInterestDays())))
		newInterestDays.resize(nbNewResets);
	if(isInterestTerms = ((leg && leg->GetInterestTerms()) || (isSched && sched->GetInterestTerms())))
		newInterestTerms.resize(nbNewResets);
	if(isRateMaturities = (leg && leg->GetRateMaturities()))
		newRateMaturities.Resize(nbNewResets);
	if(isFwdRates = (leg && leg->GetFwdRates()))
		newFwdRates.Resize(nbNewResets);

	for(i=0;i<nbNewResets;++i)
	{
		if(isResetDates)		newResets[i]		= leg ? (*(leg->GetResetDates()))[newToOldReset[i]] : (*(sched->GetResetDates()))[newToOldReset[i]];
		if(isPaymentDates)		newPayments[i]		= leg ? (*(leg->GetPaymentDates()))[newToOldReset[i]] : (*(sched->GetPaymentDates()))[newToOldReset[i]];
		if(isTheoPayDates)		newTheoPays[i]		= (*(leg->GetTheoPayDates()))[newToOldReset[i]];
		if(isFwdRateStartDates) newFwdRateStarts[i]	= leg ? (*(leg->GetFwdRateStartDates()))[newToOldReset[i]] : (*(sched->GetFwdRateStartDates()))[newToOldReset[i]];
		if(isFwdRateEndDates)	newFwdRateEnds[i]	= leg ? (*(leg->GetFwdRateEndDates()))[newToOldReset[i]] : (*(sched->GetFwdRateEndDates()))[newToOldReset[i]];
		if(isFlowStartDates)	newFlowStarts[i]	= leg ? (*(leg->GetFlowStartDates()))[newToOldReset[i]] : (*(sched->GetFlowStartDates()))[newToOldReset[i]];

		if(isFlowEndDates)		newFlowEnds[i]		= leg ? (*(leg->GetFlowEndDates()))[newToOldFlowEndDate[i]] : (*(sched->GetFlowEndDates()))[newToOldFlowEndDate[i]];

		if(isInterestDays)		newInterestDays[i]	= newFlowEnds[i] - newFlowStarts[i];
		if(isInterestTerms)		newInterestTerms[i] = newInterestDays[i]
			* (leg ? (*(leg->GetInterestTerms()))[newToOldReset[i]] : (*(sched->GetInterestTerms()))[newToOldReset[i]])
			/ (leg ? (*(leg->GetInterestDays()))[newToOldReset[i]] : (*(sched->GetInterestDays()))[newToOldReset[i]]);

		if(isRateMaturities)	newRateMaturities[i]= (*(leg->GetRateMaturities()))[newToOldReset[i]];
		if(isFwdRates)			newFwdRates[i]= (*(leg->GetFwdRates()))[newToOldReset[i]];
	}

	if(isResetDates)		leg ? leg->SetResetDates(&To_ARM_Vector(newResets))					: sched->SetResetDates(&newResets);

	/// Careful : payments, attribute of ARM_Security, is not cloned by SetPaymentDates()...
	/// ...but each set of SwapLeg clones date vectors !!!
	if(isPaymentDates)		leg ? leg->SetPaymentDates(To_pARM_Vector(&newPayments))			: sched->SetPaymentDates(&newPayments);

	if(isTheoPayDates)		leg->SetTheoPayDates(&newTheoPays);
	if(isFwdRateStartDates) leg ? leg->SetFwdRateStartDates(&To_ARM_Vector(newFwdRateStarts))	: sched->SetFwdRateStartDates(&newFwdRateStarts);
	if(isFwdRateEndDates)	leg ? leg->SetFwdRateEndDates(&To_ARM_Vector(newFwdRateEnds))		: sched->SetFwdRateEndDates(&newFwdRateEnds);
	if(isFlowStartDates)	leg ? leg->SetFlowStartDates(&To_ARM_Vector(newFlowStarts))			: sched->SetFlowStartDates(&newFlowStarts);
	if(isFlowEndDates)		leg ? leg->SetFlowEndDates(&To_ARM_Vector(newFlowEnds))				: sched->SetFlowEndDates(&newFlowEnds);
	if(isInterestDays)		leg ? leg->SetInterestDays(&To_ARM_Vector(newInterestDays))			: sched->SetInterestDays(&newInterestDays);
	if(isInterestTerms)		leg ? leg->SetInterestTerms(&To_ARM_Vector(newInterestTerms))		: sched->SetInterestTerms(&newInterestTerms);
	if(isRateMaturities)	leg->SetRateMaturities(&newRateMaturities);
	if(isFwdRates)			leg->SetFwdRates(&newFwdRates);

	if(leg)
		leg->SetNumFlows(nbNewResets);
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: AgregateSameResetDates
///	Returns: void
///	Action : Find same reset dates and merge associated flows
/////////////////////////////////////////////////////////////////
void ARM_CRASpreadCalculator::AgregateSameResetDates(ARM_SpreadOption* sec)
{
	ARM_GP_Vector resets(To_ARM_GP_Vector(*(sec->GetSwapLeg()->GetResetDates())));
	size_t nbResets = resets.size();
	if(nbResets <= 1)
		return;

	/// Find doublons
	size_t resetIdx = 0;
	double lastReset = resets[0];
	ARM_IntVector resetsToErase(0);
	for(resetIdx=1;resetIdx<nbResets;++resetIdx)
	{
		if(resets[resetIdx] == lastReset)
			resetsToErase.push_back(resetIdx);
		lastReset = resets[resetIdx];
	}

	/// Agregate flows with same reset
	size_t nbResetsToErase = resetsToErase.size();
	size_t nbNewResets = nbResets - nbResetsToErase;

	resetsToErase.push_back(nbResets+1); // +1 to copy non erased resets after last erased one

	size_t i,newResetIdx=0,lastNewResetIdx=0;
	ARM_IntVector newToOldReset(nbNewResets);
	ARM_IntVector newToOldFlowEndDate(nbNewResets);
	resetIdx=0;
	for(i=0;i<nbResetsToErase+1;++i)
	{
		while(resetIdx < nbResets && resetIdx < resetsToErase[i])
		{
			newToOldReset[newResetIdx]		= resetIdx;
			newToOldFlowEndDate[newResetIdx]= resetIdx;

			lastNewResetIdx = newResetIdx;
			++newResetIdx;
			++resetIdx;
		}
		if(resetIdx == resetsToErase[i])
		{
			/// Here is a reset to erase => save flow end date of previous non erased reset
			/// for furher update
			newToOldFlowEndDate[lastNewResetIdx] = resetIdx;
			++resetIdx;
		}
	}

	/// Update legs
	UpdateLegSchedules(newToOldReset,newToOldFlowEndDate,sec->GetSwapLeg());
	UpdateLegSchedules(newToOldReset,newToOldFlowEndDate,sec->GetSpreadLeg()->GetFirstLeg());
	UpdateLegSchedules(newToOldReset,newToOldFlowEndDate,sec->GetSpreadLeg()->GetSecondLeg());
	sec->SetNumFlows(nbNewResets);
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: CreateCSOPortfolio_FakeForDateStripOnly
///	Returns: void
///	Action : 
/////////////////////////////////////////////////////////////////
void ARM_CRASpreadCalculator::CreateCSOPortfolio_FakeForDateStripOnly()
{
	bool flowByFlow = IsVms()?itsBootstrapOptimizer:FLOWBYFLOW;
	
	if (IsVms() && flowByFlow)
	{
		CreateCSOPortfolio_FakeForDateStripOnly_ForVMS();
	}
	else
	{
		ARM_IRIndex payIndex;

		if (IsLiborIndex((ARM_INDEX_TYPE)itsPayIndex))
		{					
			payIndex = ARM_IRIndex((ARM_INDEX_TYPE)itsPayIndex, K_DEF_FREQ, K_DEF_FREQ, &itsCcy);
			payIndex.SetIntRule(K_ADJUSTED);
			
			payIndex.SetResetTiming(GetBoostedResetTiming());
		}
		else
		{
			int dayCount	= itsCcy.GetFixedDayCount();
			int resetFreq	= itsCcy.GetFixedPayFreq();
			int payFreq		= itsCcy.GetFixedPayFreq();
			double term		= -1;
			int compMeth	= 0;
			int fwdRule		= itsCcy.GetFwdRule();
			int intRule		= K_ADJUSTED;
			int resetTiming = GetBoostedResetTiming();
			int resetGap	= -itsCcy.GetSpotDays();
			int payTiming	= K_ARREARS;
			int payGap		= 0;
					
			payIndex = ARM_IRIndex( dayCount, resetFreq, payFreq, term, 
								   compMeth, fwdRule, itsBoostedIntRule/*intRule*/, resetTiming, 
								   resetGap, payTiming, payGap, 
								   &itsCcy,
								   (ARM_INDEX_TYPE)itsPayIndex);
		}

		ARM_ReferenceValue payIndexRate(0.0);

		if (itsPayIndexMult.GetSize() != 0)
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_CRASpreadCalculator::CreateCSOPortfolio :  step-up pay index mult not supported");

		/// spread option does not handle a vector of pay index mults
		double payIndexMult = itsPayIndexMult.Interpolate(GetMktDataManager()->GetAsOfDate().GetJulian());

		ARM_ReferenceValue aWeight1(itsCoeff1);
		ARM_ReferenceValue aWeight2(itsCoeff2);
		ARM_ReferenceValue aNotUsed(1.);

		/// CMS spread corridor
		ARM_SpreadOption spreadOption(
			itsStartDate,
			itsEndDate,
			K_CAP, 
			&aNotUsed,
			&payIndex,
			&payIndexRate,
			(ARM_INDEX_TYPE)itsRefIndex2,  
			(ARM_INDEX_TYPE)itsRefIndex1,
			&aWeight2, 
			&aWeight1,
			itsBoostedDayCount, 
			itsCpnResetFreq, 
			itsCpnPayFreq, 
			itsCpnResetTiming, 
			K_ARREARS,
			&itsCcy,
			itsCpnResetGap,
			-CORRIDOR_STRIKE_SPREAD,
			CORRIDOR_STRIKE_SPREAD,
			(ARM_Vector*) NULL,
			(ARM_Vector*) NULL,
			itsBoostedIntRule,
			K_SHORTSTART,
			const_cast<char*>(itsCpnResetCal.c_str()),
			const_cast<char*>(itsCpnPayCal.c_str()),
			1,
			1,
			&aNotUsed,
			payIndexMult, 
			NULL,
			0);

		/// Agregate flow with same resets (artifac of weekly reset case)
		AgregateSameResetDates(&spreadOption);

		if(itsOptimResetData[OPTIM_RESET])
			OptimiseResetDates(&spreadOption);

		size_t nbPeriods = spreadOption.GetPayIndexLeg()->GetFlowStartDates()->size();
		size_t nbFlows = spreadOption.GetNumFlows();
		
		ARM_GP_Vector resetTimesForStruct(nbPeriods);
		ARM_GP_Vector startTimesForStruct(nbPeriods);
		ARM_GP_Vector endTimesForStruct(nbPeriods);
		ARM_GP_Vector payTimesForStruct(nbPeriods);
		ARM_GP_Vector payPeriodsForStruct(nbPeriods);

		ARM_GP_Vector resetTimesForPay(nbPeriods);
		ARM_GP_Vector startTimesForPay(nbPeriods);
		ARM_GP_Vector endTimesForPay(nbPeriods);
		ARM_GP_Vector payTimesForPay(nbPeriods);
		ARM_GP_Vector payPeriodsForPay(nbPeriods);

		ARM_GP_Vector resetTimesForFix(nbFlows);
		ARM_GP_Vector startTimesForFix(nbFlows);
		ARM_GP_Vector endTimesForFix(nbFlows);
		ARM_GP_Vector fwdStartTimesForFix(nbFlows);
		ARM_GP_Vector fwdEndTimesForFix(nbFlows);
		ARM_GP_Vector payTimesForFix(nbFlows);
		ARM_GP_Vector itForFix(nbFlows);
		ARM_GP_Vector weightForFix(nbFlows);

		size_t i;
		double asof = GetMktDataManager()->GetAsOfDate().GetJulian();
		for (i=0;i<nbPeriods;i++)
		{
		// DATESTRIP structure
			resetTimesForStruct[i]		= asof;
			// récupéré par le gramfunctor mais pas utilisé dans le cas corridor CMS spread (pour le modèle local)
			startTimesForStruct[i]		= spreadOption.GetPayIndexLeg()->GetFlowStartDates()->Elt(i);
			endTimesForStruct[i]		= spreadOption.GetPayIndexLeg()->GetFlowEndDates()->Elt(i);

			// récupéré par le gramfunctor et vraiment utilisé par le modèle local
			payTimesForStruct[i]		= spreadOption.GetPayIndexLeg()->GetPaymentDates()->Elt(i);
			payPeriodsForStruct[i]		= 0.;

		// DATESTRIP index paiement
			resetTimesForPay[i]		= spreadOption.GetPayIndexLeg()->GetResetDates()->Elt(i);
			startTimesForPay[i]		= spreadOption.GetPayIndexLeg()->GetFwdRateStartDates()->Elt(i);
			endTimesForPay[i]		= asof;
			payTimesForPay[i]		= asof;
			payPeriodsForPay[i]		= (*spreadOption.GetPayIndexLeg()->GetInterestTerms())[i];
		}

		for (i=0;i<nbFlows;i++)
		{
		// DATESTRIP fixings
			//récupéré par le gram functor, et utilisé par le modèle local
			resetTimesForFix[i]		= spreadOption.GetSwapLeg()->GetResetDates()->Elt(i);
			startTimesForFix[i]		= spreadOption.GetSwapLeg()->GetFlowStartDates()->Elt(i);
			endTimesForFix[i]		= spreadOption.GetSwapLeg()->GetFlowEndDates()->Elt(i);

			fwdStartTimesForFix[i]	= spreadOption.GetSpreadLeg()->GetFirstLeg()->GetFwdRateStartDates()->Elt(i);

			//récupéré par le gram functor pour first index = LIBOR
			fwdEndTimesForFix[i]	= spreadOption.GetSpreadLeg()->GetFirstLeg()->GetFwdRateEndDates()->Elt(i);
			itForFix[i]				= (*spreadOption.GetSpreadLeg()->GetFirstLeg()->GetInterestTerms())[i];
			
			//seulement utilisé pour le calcul des idx ci-dessous
			payTimesForFix[i]		= spreadOption.GetSwapLeg()->GetPaymentDates()->Elt(i);
		
			//récupéré par le gram functor, et utilisé par le modèle local
			//double start					= spreadOption.GetSpreadLeg()->GetFirstLeg()->GetFlowStartDates()->Elt(i);
			//double end					= spreadOption.GetSpreadLeg()->GetFirstLeg()->GetFlowEndDates()->Elt(i);
			//weightForFix[i]		= CountYears(itsBoostedDayCount, start, end);
			weightForFix[i]			= CountYears(itsBoostedDayCount, startTimesForFix[i], endTimesForFix[i]);
		}

		int nb =0;
		double payDate = payTimesForStruct[0];
		ARM_GP_Vector idx(0);
		for (i=0;i<nbFlows;i++)
		{
			if (fabs(payTimesForFix[i] - payDate) < ARM_VanillaSpreadOptionArg::PayDateToleranceInDays)
				nb++;
			else
			{
				idx.push_back(nb);
				payDate=payTimesForFix[i];
				nb=1;
			}
		}
		idx.push_back(nb);

		itsDateStrip1 = ARM_DateStripPtr(new ARM_DateStrip(
			&startTimesForStruct,
			&endTimesForStruct,
			&startTimesForStruct,
			&endTimesForStruct,
			&resetTimesForStruct,
			&payTimesForStruct,
			&payPeriodsForStruct,
			&payPeriodsForStruct));

		itsDateStrip2 = ARM_DateStripPtr(new ARM_DateStrip(
			&startTimesForPay,
			&endTimesForPay,
			&startTimesForPay,
			&endTimesForPay,
			&resetTimesForPay,
			&payTimesForPay,
			&payPeriodsForPay,
			&payPeriodsForPay));

		itsDateStrip3 = ARM_DateStripPtr(new ARM_DateStrip(
			&startTimesForFix,
			&endTimesForFix,
			&fwdStartTimesForFix,
			&fwdEndTimesForFix,
			&resetTimesForFix,
			&payTimesForFix,
			&weightForFix,
			&itForFix));

		itsIndexVector = ARM_GP_VectorPtr(new ARM_GP_Vector(idx));
	}

}

void ARM_CRASpreadCalculator::CreateCSOPortfolio_FakeForDateStripOnly_ForVMS()
{
	int forcedIntRule = K_UNADJUSTED;
	ARM_IRIndex payIndex;

	if (IsLiborIndex((ARM_INDEX_TYPE)itsPayIndex))
	{					
		payIndex = ARM_IRIndex((ARM_INDEX_TYPE)itsPayIndex, K_DEF_FREQ, K_DEF_FREQ, &itsCcy);
		payIndex.SetIntRule(K_ADJUSTED);
		
		payIndex.SetResetTiming(GetBoostedResetTiming());
	}
	else
	{
		int dayCount	= itsCcy.GetFixedDayCount();
		int resetFreq	= itsCcy.GetFixedPayFreq();
		int payFreq		= itsCcy.GetFixedPayFreq();
		double term		= -1;
		int compMeth	= 0;
		int fwdRule		= itsCcy.GetFwdRule();
		int intRule		= K_ADJUSTED;
		int resetTiming = GetBoostedResetTiming();
		int resetGap	= -itsCcy.GetSpotDays();
		int payTiming	= K_ARREARS;
		int payGap		= 0;
				
		payIndex = ARM_IRIndex( dayCount, resetFreq, payFreq, term, 
							   compMeth, fwdRule, itsBoostedIntRule/*intRule*/, resetTiming, 
							   resetGap, payTiming, payGap, 
							   &itsCcy,
							   (ARM_INDEX_TYPE)itsPayIndex);
	}

	ARM_ReferenceValue payIndexRate(0.0);

	double payIndexMult = 1.;

	ARM_ReferenceValue aWeight1(itsCoeff1);
	ARM_ReferenceValue aWeight2(itsCoeff2);
	ARM_ReferenceValue aNotUsed(1.);

	ARM_GP_Vector* startDatesAdj = GetExerDateStrip()->GetFlowStartDates();
	ARM_GP_Vector* resetDatesAdj = GetExerDateStrip()->GetResetDates();

	ARM_GP_Vector* startDates = GetExerDateStripUnadj()->GetFlowStartDates();
	ARM_GP_Vector* resetDates = GetExerDateStripUnadj()->GetResetDates();

	ARM_GP_Vector resetTimesForStruct;
	ARM_GP_Vector startTimesForStruct;
	ARM_GP_Vector endTimesForStruct;
	ARM_GP_Vector payTimesForStruct;
	ARM_GP_Vector payPeriodsForStruct;

	ARM_GP_Vector resetTimesForPay;
	ARM_GP_Vector startTimesForPay;
	ARM_GP_Vector endTimesForPay;
	ARM_GP_Vector payTimesForPay;
	ARM_GP_Vector payPeriodsForPay;

	ARM_GP_Vector resetTimesForFix;
	ARM_GP_Vector startTimesForFix;
	ARM_GP_Vector endTimesForFix;
	ARM_GP_Vector fwdStartTimesForFix;
	ARM_GP_Vector fwdEndTimesForFix;
	ARM_GP_Vector payTimesForFix;
	ARM_GP_Vector itForFix;
	ARM_GP_Vector weightForFix;

	ARM_GP_Vector idx(0);
	for ( int k = 0; k < startDates->size(); k++ )
	{
		ARM_Date startDate((*startDates)[k]);
		ARM_Date endDate((*startDates)[k]);
		endDate.AddMonths(12/itsCpnPayFreq);

		ARM_ReferenceValue barrierDown(itsCpnBarrierDown.Interpolate((*startDatesAdj)[k]));
		ARM_ReferenceValue barrierUp(itsCpnBarrierUp.Interpolate((*startDatesAdj)[k]));

		int refIndex = (int) itsTenor.Interpolate((*startDatesAdj)[k]);
		ARM_ReferenceValue refNotional(itsvCpnNominal[k]);

		ARM_SpreadOption spreadOption(
			startDate,
			endDate,
			K_CAP, 
			&barrierDown,
			&payIndex,
			&payIndexRate,
			(ARM_INDEX_TYPE)K_CMS2,  
			(ARM_INDEX_TYPE)refIndex,
			&aWeight2, 
			&aWeight1,
			itsBoostedDayCount, 
			itsCpnResetFreq, 
			itsCpnPayFreq, 
			itsCpnResetTiming, 
			K_ARREARS,
			&itsCcy,
			itsCpnResetGap,
			-CORRIDOR_STRIKE_SPREAD,
			CORRIDOR_STRIKE_SPREAD,
			(ARM_Vector*) NULL,
			(ARM_Vector*) NULL,
			forcedIntRule,
			K_SHORTSTART,
			const_cast<char*>(itsCpnResetCal.c_str()),
			const_cast<char*>(itsCpnPayCal.c_str()),
			1,
			1,
			&itsBoostedFixRate,
			payIndexMult, 
			NULL,
			0);

		size_t nbPeriods = spreadOption.GetPayIndexLeg()->GetFlowStartDates()->size();
		size_t nbFlows = spreadOption.GetNumFlows();
	
		size_t i;
		double asof = GetMktDataManager()->GetAsOfDate().GetJulian();
		for (i=0;i<nbPeriods;i++)
		{
			resetTimesForStruct.push_back(asof);
			startTimesForStruct.push_back(spreadOption.GetPayIndexLeg()->GetFlowStartDates()->Elt(i));
			endTimesForStruct.push_back(spreadOption.GetPayIndexLeg()->GetFlowEndDates()->Elt(i));
			payTimesForStruct.push_back(spreadOption.GetPayIndexLeg()->GetPaymentDates()->Elt(i));
			payPeriodsForStruct.push_back(0.);

			resetTimesForPay.push_back(spreadOption.GetPayIndexLeg()->GetResetDates()->Elt(i));
			startTimesForPay.push_back(spreadOption.GetPayIndexLeg()->GetFwdRateStartDates()->Elt(i));
			endTimesForPay.push_back(asof);
			payTimesForPay.push_back(asof);
			payPeriodsForPay.push_back((*spreadOption.GetPayIndexLeg()->GetInterestTerms())[i]);
		}

		double payDate = payTimesForStruct[0];
		int nb =0;
		for (i=0;i<nbFlows;i++)
		{
			resetTimesForFix.push_back(spreadOption.GetSwapLeg()->GetResetDates()->Elt(i));
			startTimesForFix.push_back(spreadOption.GetSwapLeg()->GetFlowStartDates()->Elt(i));
			endTimesForFix.push_back(spreadOption.GetSwapLeg()->GetFlowEndDates()->Elt(i));
			fwdStartTimesForFix.push_back(spreadOption.GetSpreadLeg()->GetFirstLeg()->GetFwdRateStartDates()->Elt(i));
			fwdEndTimesForFix.push_back(spreadOption.GetSpreadLeg()->GetFirstLeg()->GetFwdRateEndDates()->Elt(i));
			itForFix.push_back((*spreadOption.GetSpreadLeg()->GetFirstLeg()->GetInterestTerms())[i]);
			payTimesForFix.push_back(spreadOption.GetSwapLeg()->GetPaymentDates()->Elt(i));
			weightForFix.push_back((*spreadOption.GetSwapLeg()->GetInterestTerms())[i]);
			//weightForFix.push_back(CountYears(itsBoostedDayCount, startTimesForFix[i], endTimesForFix[i]));

			if (fabs(payTimesForFix[i] - payDate) < ARM_VanillaSpreadOptionArg::PayDateToleranceInDays)
				nb++;
			else
			{
				idx.push_back(nb);
				payDate=payTimesForFix[i];
				nb=1;
			}
		}
		idx.push_back(nb);
	}

	itsDateStrip1 = ARM_DateStripPtr(new ARM_DateStrip(
		&startTimesForStruct,
		&endTimesForStruct,
		&startTimesForStruct,
		&endTimesForStruct,
		&resetTimesForStruct,
		&payTimesForStruct,
		&payPeriodsForStruct,
		&payPeriodsForStruct));

	itsDateStrip2 = ARM_DateStripPtr(new ARM_DateStrip(
		&startTimesForPay,
		&endTimesForPay,
		&startTimesForPay,
		&endTimesForPay,
		&resetTimesForPay,
		&payTimesForPay,
		&payPeriodsForPay,
		&payPeriodsForPay));

	itsDateStrip3 = ARM_DateStripPtr(new ARM_DateStrip(
		&startTimesForFix,
		&endTimesForFix,
		&fwdStartTimesForFix,
		&fwdEndTimesForFix,
		&resetTimesForFix,
		&payTimesForFix,
		&weightForFix,
		&itForFix));

	itsIndexVector = ARM_GP_VectorPtr(new ARM_GP_Vector(idx));
}
	

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: CreateCorridorSpreadOption
///	Returns: ARM_StdPortfolioPtr
///	Action : Create the corridor the spread option
/// portfolio. 
/////////////////////////////////////////////////////////////////
void ARM_CRASpreadCalculator::CreateCSOPortfolio()
{
	if (IsVms())
	{
		CreateCSOPortfolio_ForVMS();
	}
	else
	{
		// The portfolio contains just a single spread option used for 
		// the local model calibration

		ARM_IRIndex payIndex;

		if (IsLiborIndex((ARM_INDEX_TYPE)itsPayIndex))
		{					
			payIndex = ARM_IRIndex((ARM_INDEX_TYPE)itsPayIndex, K_DEF_FREQ, K_DEF_FREQ, &itsCcy);
			payIndex.SetIntRule(K_ADJUSTED);
			
			payIndex.SetResetTiming(GetBoostedResetTiming());
		}
		else
		{
			int dayCount	= itsCcy.GetFixedDayCount();
			int resetFreq	= itsCcy.GetFixedPayFreq();
			int payFreq		= itsCcy.GetFixedPayFreq();
			double term		= -1;
			int compMeth	= 0;
			int fwdRule		= itsCcy.GetFwdRule();
			int intRule		= K_ADJUSTED;
			int resetTiming = GetBoostedResetTiming();
			int resetGap	= -itsCcy.GetSpotDays();
			int payTiming	= K_ARREARS;
			int payGap		= 0;
					
			payIndex = ARM_IRIndex( dayCount, resetFreq, payFreq, term, 
								   compMeth, fwdRule, itsBoostedIntRule/*intRule*/, resetTiming, 
								   resetGap, payTiming, payGap, 
								   &itsCcy,
								   (ARM_INDEX_TYPE)itsPayIndex);
		}


		
					
		ARM_ReferenceValue payIndexRate(0.0);

		if (itsPayIndexMult.GetSize() != 0)
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_CRASpreadCalculator::CreateCSOPortfolio :  step-up pay index mult not supported");


		/// spread option does not handle a vector of pay index mults
		double payIndexMult = itsPayIndexMult.Interpolate(GetMktDataManager()->GetAsOfDate().GetJulian());

		ARM_ReferenceValue aWeight1(itsCoeff1);
		ARM_ReferenceValue aWeight2(itsCoeff2);

		size_t nbSec = itsRefIndexType3 == K_FIXED ? 2 : ARM_Local_Normal_Model::DBLECOR_NB_BARRIERS;

		vector<double> mktPrices(nbSec,DEFAULT_PRICE);
		vector<double> weights(nbSec,DEFAULT_WEIGHT);
		vector<double> precisions(nbSec,DEFAULT_WEIGHT);
		vector<ARM_Security*> securities(nbSec);
		if(!IsDbleCorridor())
		{
			/// CMS spread corridor L <= S
			ARM_SpreadOption spreadOptionDOWN(
				itsStartDate,
				itsEndDate,
				K_CAP, 
				&itsCpnBarrierDown,
				&payIndex,
				&payIndexRate,
				(ARM_INDEX_TYPE)itsRefIndex2,  
				(ARM_INDEX_TYPE)itsRefIndex1,
				&aWeight2, 
				&aWeight1,
				itsBoostedDayCount, 
				itsCpnResetFreq, 
				itsCpnPayFreq, 
				itsCpnResetTiming, 
				K_ARREARS,
				&itsCcy,
				itsCpnResetGap,
				-CORRIDOR_STRIKE_SPREAD,
				CORRIDOR_STRIKE_SPREAD,
				(ARM_Vector*) NULL,
				(ARM_Vector*) NULL,
				itsBoostedIntRule,
				K_SHORTSTART,
				const_cast<char*>(itsCpnResetCal.c_str()),
				const_cast<char*>(itsCpnPayCal.c_str()),
				1,
				1,
				&itsBoostedFixRate,
				payIndexMult, 
				NULL,
				0);

			spreadOptionDOWN.SetAmount(&itsNotional);

			/// CMS spread corridor S <= U
			ARM_SpreadOption spreadOptionUP(
				itsStartDate,
				itsEndDate,
				K_CAP, 
				&itsCpnBarrierUp,
				&payIndex,
				&payIndexRate,
				(ARM_INDEX_TYPE)itsRefIndex2,  
				(ARM_INDEX_TYPE)itsRefIndex1,
				&aWeight2, 
				&aWeight1,
				itsBoostedDayCount, 
				itsCpnResetFreq, 
				itsCpnPayFreq, 
				itsCpnResetTiming, 
				K_ARREARS,
				&itsCcy,
				itsCpnResetGap,
				-CORRIDOR_STRIKE_SPREAD,
				CORRIDOR_STRIKE_SPREAD,
				(ARM_Vector*) NULL,
				(ARM_Vector*) NULL,
				itsBoostedIntRule,
				K_SHORTSTART,
				const_cast<char*>(itsCpnResetCal.c_str()),
				const_cast<char*>(itsCpnPayCal.c_str()),
				1,
				1,
				&itsBoostedFixRate,
				payIndexMult,
				NULL,
				0);
				 
			spreadOptionUP.SetAmount(&itsNotional);

			if(itsOptimResetData[OPTIM_RESET])
			{
				OptimiseResetDates(&spreadOptionDOWN);
				OptimiseResetDates(&spreadOptionUP);
			}


			securities[0] = &spreadOptionDOWN;
			securities[1] = &spreadOptionUP;

			itsVanillaArgCSO = ARM_VanillaArgPtr(CreateVanillaArgCSO(&spreadOptionDOWN));

			itsCSOPF = ARM_StdPortfolioPtr(new ARM_StdPortfolio(securities, weights, mktPrices, precisions));

			if (IsTripleRange())
			{
				vector<ARM_Security*> securities2(2);
				vector<ARM_Security*> securities3(2);
		
				/// CMS spread corridor L2 <= S
				ARM_SpreadOption spreadOption2DOWN(
					itsStartDate,itsEndDate,K_CAP, 
					&itsCpnBarrierDown2,
					&payIndex,&payIndexRate,(ARM_INDEX_TYPE)itsRefIndex2,(ARM_INDEX_TYPE)itsRefIndex1,&aWeight2,&aWeight1,itsBoostedDayCount,itsCpnResetFreq,itsCpnPayFreq,itsCpnResetTiming,K_ARREARS,&itsCcy,itsCpnResetGap,-CORRIDOR_STRIKE_SPREAD,CORRIDOR_STRIKE_SPREAD,(ARM_Vector*) NULL,(ARM_Vector*) NULL,itsBoostedIntRule,K_SHORTSTART,const_cast<char*>(itsCpnResetCal.c_str()),const_cast<char*>(itsCpnPayCal.c_str()),1,1,
					&itsBoostedFixRate2,
					payIndexMult,NULL,0);

				spreadOption2DOWN.SetAmount(&itsNotional);

				/// CMS spread corridor S <= U2
				ARM_SpreadOption spreadOption2UP(
					itsStartDate,itsEndDate,K_CAP, 
					&itsCpnBarrierUp2,
					&payIndex,&payIndexRate,(ARM_INDEX_TYPE)itsRefIndex2,(ARM_INDEX_TYPE)itsRefIndex1,&aWeight2,&aWeight1,itsBoostedDayCount,itsCpnResetFreq,itsCpnPayFreq,itsCpnResetTiming,K_ARREARS,&itsCcy,itsCpnResetGap,-CORRIDOR_STRIKE_SPREAD,CORRIDOR_STRIKE_SPREAD,(ARM_Vector*) NULL,(ARM_Vector*) NULL,itsBoostedIntRule,K_SHORTSTART,const_cast<char*>(itsCpnResetCal.c_str()),const_cast<char*>(itsCpnPayCal.c_str()),1,1,
					&itsBoostedFixRate2,
					payIndexMult,NULL,0);
					 
				spreadOption2UP.SetAmount(&itsNotional);

				if(itsOptimResetData[OPTIM_RESET])
				{
					OptimiseResetDates(&spreadOption2DOWN);
					OptimiseResetDates(&spreadOption2UP);
				}

				securities2[0] = &spreadOption2DOWN;
				securities2[1] = &spreadOption2UP;

				itsCSOPF2 = ARM_StdPortfolioPtr(new ARM_StdPortfolio(securities2, weights, mktPrices, precisions));
				
				/// CMS spread corridor L2 <= S
				ARM_SpreadOption spreadOption3DOWN(
					itsStartDate,itsEndDate,K_CAP, 
					&itsCpnBarrierDown3,
					&payIndex,&payIndexRate,(ARM_INDEX_TYPE)itsRefIndex2,(ARM_INDEX_TYPE)itsRefIndex1,&aWeight2,&aWeight1,itsBoostedDayCount,itsCpnResetFreq,itsCpnPayFreq,itsCpnResetTiming,K_ARREARS,&itsCcy,itsCpnResetGap,-CORRIDOR_STRIKE_SPREAD,CORRIDOR_STRIKE_SPREAD,(ARM_Vector*) NULL,(ARM_Vector*) NULL,itsBoostedIntRule,K_SHORTSTART,const_cast<char*>(itsCpnResetCal.c_str()),const_cast<char*>(itsCpnPayCal.c_str()),1,1,
					&itsBoostedFixRate3,
					payIndexMult,NULL,0);

				spreadOption3DOWN.SetAmount(&itsNotional);

				/// CMS spread corridor S <= U3
				ARM_SpreadOption spreadOption3UP(
					itsStartDate,itsEndDate,K_CAP, 
					&itsCpnBarrierUp3,
					&payIndex,&payIndexRate,(ARM_INDEX_TYPE)itsRefIndex2,(ARM_INDEX_TYPE)itsRefIndex1,&aWeight2,&aWeight1,itsBoostedDayCount,itsCpnResetFreq,itsCpnPayFreq,itsCpnResetTiming,K_ARREARS,&itsCcy,itsCpnResetGap,-CORRIDOR_STRIKE_SPREAD,CORRIDOR_STRIKE_SPREAD,(ARM_Vector*) NULL,(ARM_Vector*) NULL,itsBoostedIntRule,K_SHORTSTART,const_cast<char*>(itsCpnResetCal.c_str()),const_cast<char*>(itsCpnPayCal.c_str()),1,1,
					&itsBoostedFixRate3,
					payIndexMult,NULL,0);
					 
				spreadOption3UP.SetAmount(&itsNotional);

				if(itsOptimResetData[OPTIM_RESET])
				{
					OptimiseResetDates(&spreadOption3DOWN);
					OptimiseResetDates(&spreadOption3UP);
				}

				securities3[0] = &spreadOption3DOWN;
				securities3[1] = &spreadOption3UP;

				itsCSOPF3 = ARM_StdPortfolioPtr(new ARM_StdPortfolio(securities3, weights, mktPrices, precisions));
			}
		}
		else
		{
			ARM_ReferenceValue aWeight4(0.0);  // temp
			int itsRefIndex4 = (itsRefIndex3 != itsRefIndex2 ? itsRefIndex2 : itsRefIndex1); // temp
			ARM_ReferenceValue aWeight3(1.0);  // temp
			ARM_VolLInterpol* spreadRateCorrelCurve=NULL;
			if(GetKeys().size() > CorrelRateSpreadKey && !GetMktDataManager()->TestIfKeyMissing(GetKeys()[CorrelRateSpreadKey]))
			{
				spreadRateCorrelCurve = dynamic_cast< ARM_VolLInterpol* >( GetMktDataManager()->GetData(GetKeys()[CorrelRateSpreadKey]) );
				if(!spreadRateCorrelCurve)
					ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : CMS Spread/Rate correl curve for key=" + GetKeys()[CorrelRateSpreadKey] + " is expected in the Market Data Manager");
			}

			/// Test activation of barriers
			bool isRateBarDownOff	=	itsRateBarrierDown.IsConstant() && (*(itsRateBarrierDown.GetDiscreteValues()))[0] < -ARM_CorridorDblCondition::BarrierLimit;
			bool isRateBarUpOff		=	itsRateBarrierUp.IsConstant() && (*(itsRateBarrierUp.GetDiscreteValues()))[0] > ARM_CorridorDblCondition::BarrierLimit;
			bool isSpreadBarDownOff	=	itsCpnBarrierDown.IsConstant() && (*(itsCpnBarrierDown.GetDiscreteValues()))[0] < -ARM_CorridorDblCondition::BarrierLimit;
			bool isSpreadBarUpOff	=	itsCpnBarrierUp.IsConstant() && (*(itsCpnBarrierUp.GetDiscreteValues()))[0] > ARM_CorridorDblCondition::BarrierLimit;

			/// Option type on spread or rate condition is always K_CAP
			/// because in further use we need to get the value of the whole
			/// double condition corridor. But the rank in the portfolio
			/// give the type of barriers
			ARM_CorridorDblCondition* SDownRDownCorridor=NULL;
			ARM_CorridorDblCondition* SDownRUpCorridor=NULL;
			ARM_CorridorDblCondition* SUpRDownCorridor=NULL;
			ARM_CorridorDblCondition* SUpRUpCorridor=NULL;

			bool isSDownRDown	= false;
			bool isSDownRUp		= false;
			bool isSUpRDown		= false;
			bool isSUpRUp		= false;
			int rateBarType		= K_CAP;
			int spreadBarType	= K_CAP;

			if((isRateBarDownOff || isRateBarUpOff) && (isSpreadBarDownOff || isSpreadBarUpOff))
			{
				/// The RA2 can be described by a single double condition
				/// by using barrier types K_CAP/K_FLOOR
				if(isRateBarDownOff)
				{
					if(isSpreadBarDownOff)
					{
						/// Double corridor -Inf<=R<=Ur & -Inf<=S<=Us => Floor on both conditions
						isSUpRUp		= true;
						rateBarType		= K_FLOOR;
						spreadBarType	= K_FLOOR;
					}
					else
					{
						/// Double corridor -Inf<=R<=Ur & Ls<=S<=+Inf => Floor on Rate, Cap on Spread
						isSDownRUp		= true;
						rateBarType		= K_FLOOR;
					}
				}
				else
				{
					if(isSpreadBarDownOff)
					{
						/// Double corridor Lr<=R<=+Inf & -Inf<=S<=Us => Cap on Rate, Floor on Spread
						isSUpRDown		= true;
						spreadBarType	= K_FLOOR;
					}
					else
					{
						/// Double corridor Lr<=R<=+Inf & Ls<=S<=+Inf => Cap on both conditions
						isSDownRDown	= true;
					}
				}
			}
			else
			{
				/// Four RA2s are used to described actual ranges on both conditions
				/// and all barrier types are K_CAP
				isSDownRDown	= true;
				isSDownRUp		= true;
				isSUpRDown		= true;
				isSUpRUp		= true;
			}

			ARM_CorridorDblCondition* argCRA2=NULL;

			/// Double corridor Lr<=R & Ls<=S
			if(isSDownRDown)
			{
				SDownRDownCorridor = new ARM_CorridorDblCondition(
					itsStartDate,itsEndDate,
					rateBarType,
					spreadBarType,
					&itsRateBarrierDown,
					&itsCpnBarrierDown,
					&payIndex,
					&payIndexRate,
					&itsBoostedFixRate, // payIndexMargin
					(ARM_INDEX_TYPE)itsRefIndex4,
					(ARM_INDEX_TYPE)itsRefIndex3,
					(ARM_INDEX_TYPE)itsRefIndex2,  
					(ARM_INDEX_TYPE)itsRefIndex1,
					&aWeight4,
					&aWeight3,
					&aWeight2, 
					&aWeight1,
					NULL,NULL,NULL,NULL,NULL,
					itsBoostedDayCount, 
					itsCpnResetFreq, 
					itsCpnPayFreq, 
					itsCpnResetTiming, 
					K_ARREARS,
					&itsCcy,
					itsCpnResetGap,
					CORRIDOR_STRIKE_SPREAD,CORRIDOR_STRIKE_SPREAD,
					itsBoostedIntRule,
					K_SHORTSTART,
					const_cast<char*>(itsCpnResetCal.c_str()),
					const_cast<char*>(itsCpnPayCal.c_str()),
					payIndexMult,
					0,
					10000,
					spreadRateCorrelCurve);

				SDownRDownCorridor->SetAmount(&itsNotional);

				if(itsOptimResetData[OPTIM_RESET])
				{
					OptimiseResetDates(SDownRDownCorridor);
					OptimiseResetDates(SDownRDownCorridor->GetSpreadDigital());
				}

				argCRA2 = SDownRDownCorridor;
			}

			/// Double corridor R<=Ur & Ls<=S
			if(isSDownRUp)
			{
				SDownRUpCorridor = new ARM_CorridorDblCondition(
					itsStartDate,itsEndDate,
					rateBarType,
					spreadBarType,
					&itsRateBarrierUp,
					&itsCpnBarrierDown,
					&payIndex,
					&payIndexRate,
					&itsBoostedFixRate, // payIndexMargin
					(ARM_INDEX_TYPE)itsRefIndex4,
					(ARM_INDEX_TYPE)itsRefIndex3,
					(ARM_INDEX_TYPE)itsRefIndex2,  
					(ARM_INDEX_TYPE)itsRefIndex1,
					&aWeight4,
					&aWeight3,
					&aWeight2, 
					&aWeight1,
					NULL,NULL,NULL,NULL,NULL,
					itsBoostedDayCount, 
					itsCpnResetFreq, 
					itsCpnPayFreq, 
					itsCpnResetTiming, 
					K_ARREARS,
					&itsCcy,
					itsCpnResetGap,
					CORRIDOR_STRIKE_SPREAD,CORRIDOR_STRIKE_SPREAD,
					itsBoostedIntRule,
					K_SHORTSTART,
					const_cast<char*>(itsCpnResetCal.c_str()),
					const_cast<char*>(itsCpnPayCal.c_str()),
					payIndexMult,
					0,
					10000,
					spreadRateCorrelCurve);

				SDownRUpCorridor->SetAmount(&itsNotional);

				if(itsOptimResetData[OPTIM_RESET])
				{
					OptimiseResetDates(SDownRUpCorridor);
					OptimiseResetDates(SDownRUpCorridor->GetSpreadDigital());
				}

				if(!argCRA2)
					argCRA2 = SDownRUpCorridor;
			}

			/// Double corridor Lr<=R & S<=Ls
			if(isSUpRDown)
			{
				SUpRDownCorridor = new ARM_CorridorDblCondition(
					itsStartDate,itsEndDate,
					rateBarType,
					spreadBarType,
					&itsRateBarrierDown,
					&itsCpnBarrierUp,
					&payIndex,
					&payIndexRate,
					&itsBoostedFixRate, // payIndexMargin
					(ARM_INDEX_TYPE)itsRefIndex4,
					(ARM_INDEX_TYPE)itsRefIndex3,
					(ARM_INDEX_TYPE)itsRefIndex2,  
					(ARM_INDEX_TYPE)itsRefIndex1,
					&aWeight4,
					&aWeight3,
					&aWeight2, 
					&aWeight1,
					NULL,NULL,NULL,NULL,NULL,
					itsBoostedDayCount, 
					itsCpnResetFreq, 
					itsCpnPayFreq, 
					itsCpnResetTiming, 
					K_ARREARS,
					&itsCcy,
					itsCpnResetGap,
					CORRIDOR_STRIKE_SPREAD,CORRIDOR_STRIKE_SPREAD,
					itsBoostedIntRule,
					K_SHORTSTART,
					const_cast<char*>(itsCpnResetCal.c_str()),
					const_cast<char*>(itsCpnPayCal.c_str()),
					payIndexMult,
					0,
					10000,
					spreadRateCorrelCurve);

				SUpRDownCorridor->SetAmount(&itsNotional);

				if(itsOptimResetData[OPTIM_RESET])
				{
					OptimiseResetDates(SUpRDownCorridor);
					OptimiseResetDates(SUpRDownCorridor->GetSpreadDigital());
				}

				if(!argCRA2)
					argCRA2 = SUpRDownCorridor;
			}


			/// Double corridor R<=Ur & S<=Us
			if(isSUpRUp)
			{
				SUpRUpCorridor = new ARM_CorridorDblCondition(
					itsStartDate,itsEndDate,
					rateBarType,
					spreadBarType,
					&itsRateBarrierUp,
					&itsCpnBarrierUp,
					&payIndex,
					&payIndexRate,
					&itsBoostedFixRate, // payIndexMargin
					(ARM_INDEX_TYPE)itsRefIndex4,
					(ARM_INDEX_TYPE)itsRefIndex3,
					(ARM_INDEX_TYPE)itsRefIndex2,  
					(ARM_INDEX_TYPE)itsRefIndex1,
					&aWeight4,
					&aWeight3,
					&aWeight2, 
					&aWeight1,
					NULL,NULL,NULL,NULL,NULL,
					itsBoostedDayCount, 
					itsCpnResetFreq, 
					itsCpnPayFreq, 
					itsCpnResetTiming, 
					K_ARREARS,
					&itsCcy,
					itsCpnResetGap,
					CORRIDOR_STRIKE_SPREAD,CORRIDOR_STRIKE_SPREAD,
					itsBoostedIntRule,
					K_SHORTSTART,
					const_cast<char*>(itsCpnResetCal.c_str()),
					const_cast<char*>(itsCpnPayCal.c_str()),
					payIndexMult,
					0,
					10000,
					spreadRateCorrelCurve);

				SUpRUpCorridor->SetAmount(&itsNotional);

				if(itsOptimResetData[OPTIM_RESET])
				{
					OptimiseResetDates(SUpRUpCorridor);
					OptimiseResetDates(SUpRUpCorridor->GetSpreadDigital());
				}

				if(!argCRA2)
					argCRA2 = SUpRUpCorridor;
			}

			securities[ARM_Local_Normal_Model::DBLECOR_SDOWN_RDOWN]	= SDownRDownCorridor;
			securities[ARM_Local_Normal_Model::DBLECOR_SDOWN_RUP]	= SDownRUpCorridor;
			securities[ARM_Local_Normal_Model::DBLECOR_SUP_RDOWN]	= SUpRDownCorridor;
			securities[ARM_Local_Normal_Model::DBLECOR_SUP_RUP]		= SUpRUpCorridor;

			itsVanillaArgCSO = ARM_VanillaArgPtr(CreateVanillaArgCSO(argCRA2));

			itsCSOPF = ARM_StdPortfolioPtr(new ARM_StdPortfolio(securities, weights, mktPrices, precisions));

			delete SDownRDownCorridor;
			delete SDownRUpCorridor;
			delete SUpRDownCorridor;
			delete SUpRUpCorridor;
		}
	}
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: CreateCSOPortfolio_ForVMS
///	Returns: ARM_StdPortfolioPtr
///	Action : Create the corridor the spread option
/// portfolio. 
/////////////////////////////////////////////////////////////////
void ARM_CRASpreadCalculator::CreateCSOPortfolio_ForVMS()
{
	int forcedIntRule = K_UNADJUSTED;
	ARM_IRIndex payIndex;

	if (IsLiborIndex((ARM_INDEX_TYPE)itsPayIndex))
	{					
		payIndex = ARM_IRIndex((ARM_INDEX_TYPE)itsPayIndex, K_DEF_FREQ, K_DEF_FREQ, &itsCcy);
		payIndex.SetIntRule(K_ADJUSTED);
			
		payIndex.SetResetTiming(GetBoostedResetTiming());
	}
	else
	{
		int dayCount	= itsCcy.GetFixedDayCount();
		int resetFreq	= itsCcy.GetFixedPayFreq();
		int payFreq		= itsCcy.GetFixedPayFreq();
		double term		= -1;
		int compMeth	= 0;
		int fwdRule		= itsCcy.GetFwdRule();
		int intRule		= K_ADJUSTED;
		int resetTiming = GetBoostedResetTiming();
		int resetGap	= -itsCcy.GetSpotDays();
		int payTiming	= K_ARREARS;
		int payGap		= 0;
					
		payIndex = ARM_IRIndex( dayCount, resetFreq, payFreq, term, 
							   compMeth, fwdRule, intRule, resetTiming, 
							   resetGap, payTiming, payGap, 
							   &itsCcy,
							   (ARM_INDEX_TYPE)itsPayIndex);
	}
					
	ARM_ReferenceValue payIndexRate(0.0);

	if (itsPayIndexMult.GetSize() != 0)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_CRASpreadCalculator::CreateCSOPortfolio :  step-up pay index mult not supported");


	/// spread option does not handle a vector of pay index mults
	double payIndexMult = itsPayIndexMult.Interpolate(GetMktDataManager()->GetAsOfDate().GetJulian());

	ARM_ReferenceValue aWeight1(itsCoeff1);
	ARM_ReferenceValue aWeight2(0.);

	ARM_GP_Vector* startDatesAdj = GetExerDateStrip()->GetFlowStartDates();
	ARM_GP_Vector* resetDatesAdj = GetExerDateStrip()->GetResetDates();

	ARM_GP_Vector* startDates = GetExerDateStripUnadj()->GetFlowStartDates();
	ARM_GP_Vector* resetDates = GetExerDateStripUnadj()->GetResetDates();

	list< ARM_Security* > spreadOptionList;
	
	for ( int i = 0; i < startDates->size(); i++ )
	{
		ARM_Date startDate((*startDates)[i]);
		ARM_Date endDate((*startDates)[i]);
		endDate.AddMonths(12/itsCpnPayFreq);
	
		ARM_ReferenceValue barrierDown(itsCpnBarrierDown.Interpolate((*startDatesAdj)[i]));
		ARM_ReferenceValue barrierUp(itsCpnBarrierUp.Interpolate((*startDatesAdj)[i]));

		int refIndex = (int) itsTenor.Interpolate((*startDatesAdj)[i]);
		ARM_ReferenceValue refNotional(itsvCpnNominal[i]);

		ARM_SpreadOption spreadOptionDOWN(
			startDate,
			endDate,
			K_CAP, 
			&barrierDown,
			&payIndex,
			&payIndexRate,
			(ARM_INDEX_TYPE)K_CMS2,  
			(ARM_INDEX_TYPE)refIndex,
			&aWeight2, 
			&aWeight1,
			itsBoostedDayCount, 
			itsCpnResetFreq, 
			itsCpnPayFreq, 
			itsCpnResetTiming, 
			K_ARREARS,
			&itsCcy,
			itsCpnResetGap,
			-CORRIDOR_STRIKE_SPREAD,
			CORRIDOR_STRIKE_SPREAD,
			(ARM_Vector*) NULL,
			(ARM_Vector*) NULL,
			forcedIntRule,
			K_SHORTSTART,
			const_cast<char*>(itsCpnResetCal.c_str()),
			const_cast<char*>(itsCpnPayCal.c_str()),
			1,
			1,
			&itsBoostedFixRate,
			payIndexMult, 
			NULL,
			0);

		spreadOptionDOWN.SetAmount(&refNotional);

		spreadOptionList.push_back(static_cast<ARM_SpreadOption*>(spreadOptionDOWN.Clone()));

		ARM_SpreadOption spreadOptionUP(
			startDate,
			endDate,
			K_CAP, 
			&barrierUp,
			&payIndex,
			&payIndexRate,
			(ARM_INDEX_TYPE)24,  
			(ARM_INDEX_TYPE)refIndex,
			&aWeight2, 
			&aWeight1,
			itsBoostedDayCount, 
			itsCpnResetFreq, 
			itsCpnPayFreq, 
			itsCpnResetTiming, 
			K_ARREARS,
			&itsCcy,
			itsCpnResetGap,
			-CORRIDOR_STRIKE_SPREAD,
			CORRIDOR_STRIKE_SPREAD,
			(ARM_Vector*) NULL,
			(ARM_Vector*) NULL,
			forcedIntRule,
			K_SHORTSTART,
			const_cast<char*>(itsCpnResetCal.c_str()),
			const_cast<char*>(itsCpnPayCal.c_str()),
			1,
			1,
			&itsBoostedFixRate,
			payIndexMult,
			NULL,
			0);

		spreadOptionUP.SetAmount(&refNotional);

		spreadOptionList.push_back(static_cast<ARM_SpreadOption*>(spreadOptionUP.Clone()));
		
	}

	itsCSOPF = ARM_StdPortfolioPtr(new ARM_StdPortfolio(spreadOptionList));
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: CreateCorridorSpreadOption
///	Returns: ARM_StdPortfolioPtr
///	Action : Compute prices in the spread option
/// portfolio 
/////////////////////////////////////////////////////////////////
void ARM_CRASpreadCalculator::ComputeCSOPrices()
{
	double asOf = GetMktDataManager()->GetAsOfDate().GetJulian();
	ARM_Date asOfDate(asOf);

	/// Get the B&S model and vol curves from the market data manager
	ARM_BSModel* BSModel = dynamic_cast< ARM_BSModel* >(GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );

	for (size_t i = 0; i < itsCSOPF->GetSize(); ++i)
	{
		if(itsCSOPF->GetAsset(i))
		{
			ARM_SpreadOption* spreadOption=static_cast< ARM_SpreadOption* >(itsCSOPF->GetAsset(i));
			spreadOption->SetModel(BSModel);
			double price=spreadOption->ComputePrice();
			itsCSOPF->SetPrice(price, i);
		}
	}

	if (IsTripleRange())
	{
		for (size_t i = 0; i < itsCSOPF2->GetSize(); ++i)
		{
			if(itsCSOPF2->GetAsset(i))
			{
				ARM_SpreadOption* spreadOption=static_cast< ARM_SpreadOption* >(itsCSOPF2->GetAsset(i));
				spreadOption->SetModel(BSModel);
				double price=spreadOption->ComputePrice();
				itsCSOPF2->SetPrice(price, i);
			}
		}
		for (i = 0; i < itsCSOPF3->GetSize(); ++i)
		{
			if(itsCSOPF3->GetAsset(i))
			{
				ARM_SpreadOption* spreadOption=static_cast< ARM_SpreadOption* >(itsCSOPF3->GetAsset(i));
				spreadOption->SetModel(BSModel);
				double price=spreadOption->ComputePrice();
				itsCSOPF3->SetPrice(price, i);
			}
		}
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: UpdateCalibration
///	Returns: void
///	Action : update the calibration datas w.r.t. MarketDataManager
///          The context is an hedge ratios computation
/////////////////////////////////////////////////////////////////
void ARM_CRASpreadCalculator::UpdateCalibration(bool isUpdateStrike)
{
	ComputeCSOPrices();

	CreateEmptyCalibration();

	ComputeCalibPortfolioPrices();

	Calibrate();

	/// Reset any previous prices
	ResetHasBeenPriced();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: CreateAndSetCalibration
///	Returns: void
///	Action : create the calibration methods
/////////////////////////////////////////////////////////////////
void ARM_CRASpreadCalculator::CreateAndSetCalibration()
{
	// The corridor spread option need to be created and priced for
	// the basket calibration.

	CreateCSOPortfolio();

	ComputeCSOPrices();

	CreateEmptyCalibration();

	ComputeCalibPortfolioPrices();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: ComputeCalibPortfolioPrices
///	Returns: void
///	Action : create the calibration methods
/////////////////////////////////////////////////////////////////
void ARM_CRASpreadCalculator::ComputeCalibPortfolioPrices()
{

	if (itsCalibType == DIAG_SPREAD_LONG || itsCalibType == DIAG_SPREAD_SHORT || itsCalibType == DIAG_SPREAD_INDEX ||
		itsCalibType == SHORT_LONG_SPREAD ||
		(itsModelType == ARM_PricingModelType::HWM2F && IsDbleCorridor()))
	{
		if(!GetCalibMethod()->GetlinkedMethod() || !GetCalibMethod()->GetlinkedMethod()->GetlinkedMethod())
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Unable to find model calibration portfolios");
		}
		if (itsCalibType == SHORT_LONG_SPREAD)
		{
			ComputeSOPrices(GetCalibMethod()->GetPortfolio(),itsCalibStrikeType[0]);
			ComputeSOPrices(GetCalibMethod()->GetlinkedMethod()->GetPortfolio(),itsCalibStrikeType[2]);
			ComputeSOPrices(GetCalibMethod()->GetlinkedMethod()->GetlinkedMethod()->GetPortfolio(),itsCalibStrikeType[1]);
		}
		else
		{
			ComputeSwaptionPrices();
			ComputeSOPrices(GetCalibMethod()->GetlinkedMethod()->GetPortfolio(),itsCalibStrikeType[1],true);
			bool isSOBarrier = !(IsDbleCorridor());
			ComputeSOPrices(GetCalibMethod()->GetlinkedMethod()->GetlinkedMethod()->GetPortfolio(),itsCalibStrikeType[2],isSOBarrier);
		}
	}
	else
	{
		if (itsCalibStrikeType[0] != FRONTIER) 
			ComputeSwaptionPrices();
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: CreateEmptyCalibration
///	Returns: void
///	Action : create the calibration methods
/////////////////////////////////////////////////////////////////
void ARM_CRASpreadCalculator::CreateEmptyCalibration()
{
	if (	itsModelType == ARM_PricingModelType::HWM2F && IsDbleCorridor()
		||	(itsCalibType==DIAG_SPREAD_LONG)
		||	(itsCalibType==DIAG_SPREAD_SHORT)
		||	(itsCalibType==DIAG_SPREAD_INDEX)
		||	(itsCalibType==SHORT_LONG_SPREAD)
		)
	{
		CreateEmptyCalibrationFor2F();
	}
	else
	{
		ARM_StdPortfolioPtr swaptionPF = CreateSwaptionPortfolio();

		//if(IsDbleCorridor() && itsCalibType == DIAG_CALIBRATION && itsCalibStrikeType[0] == FRONTIER)
		if(itsCalibType == DIAG_CALIBRATION && itsCalibStrikeType[0] == FRONTIER)
		{
			SetDiagStrikeToFrontier(swaptionPF);
		}

		ARM_GP_Vector calibTimes(1, 0.00);
		ARM_GP_Vector volatility(1, SIGMA_DEFAULT_VALUE);
		ARM_GP_Vector lowerVolatility(1, SIGMA_LOWER_BOUND);
		ARM_GP_Vector upperVolatility(1, SIGMA_UPPER_BOUND);
		ARM_CurveModelParam* volParam = new ARM_CurveModelParam(ARM_ModelParamType::Volatility, &volatility, &calibTimes ,"","STEPUPRIGHT",&lowerVolatility,&upperVolatility);
		ARM_ModelParamVector volParamVec(1,volParam);

		ARM_CalibMethodPtr calibMethod(NULL);

		/// To be upgraded...
		ARM_ModelFitterDes modelfitter(ARM_ModelFitterSolverType::Brent);
		calibMethod = ARM_CalibMethodPtr( new ARM_CalibMethod(swaptionPF, volParamVec, ARM_CalibMethodType::Bootstrap1D,&modelfitter) );

		SetCalibMethod(calibMethod);

		delete volParam;
	}
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: SetDiagStrikeToFrontier
///	Returns: void
///	Action : compute exercise strikes and set them in the diagonal
///			 swaption portfolio
/////////////////////////////////////////////////////////////////
void ARM_CRASpreadCalculator::SetDiagStrikeToFrontier(ARM_StdPortfolioPtr& swaptionPF)
{
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

	ARM_BSModel* modelGen = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswModelKey]) );
	ARM_BSSmiledModel* oswBSModel = modelGen->GetSabrModel();
	
	if (oswBSModel)
	{
		/// Set ATM strike
		size_t i,nbSwaptions = swaptionPF->size();
		double atmStrike,price,vega;
		ARM_Swaption *swaption;
		for(i=0;i<nbSwaptions;++i)
		{
			swaption = static_cast< ARM_Swaption* >(swaptionPF->GetAsset(i));
			swaption->SetModel(oswBSModel);
			atmStrike = swaption->PriceToRate((ARM_Date) asOfDate, 0.0);
			swaption->UpdateStrike(atmStrike);
			price = swaption->ComputePrice();
			swaptionPF->SetPrice(price,i);
		}

		/// Calibrate ATM diagonal swaption
		ARM_CalibMethodPtr saveCalibMethod = GetCalibMethod();

		ARM_GP_Vector calibTimes(1, 0.00);
		ARM_GP_Vector volatility(1, SIGMA_DEFAULT_VALUE);
		ARM_GP_Vector lowerVolatility(1, SIGMA_LOWER_BOUND);
		ARM_GP_Vector upperVolatility(1, SIGMA_UPPER_BOUND);
		ARM_CurveModelParam volParam(ARM_ModelParamType::Volatility, &volatility, &calibTimes ,"","STEPUPRIGHT",&lowerVolatility,&upperVolatility);
		ARM_ModelParamVector volParamVect(1,&volParam);
		ARM_ModelFitterDes modelfitter(ARM_ModelFitterSolverType::Brent);
		ARM_CalibMethodPtr AtmVolCalib(new ARM_CalibMethod(swaptionPF,volParamVect,ARM_CalibMethodType::Bootstrap1D,&modelfitter));
		SetCalibMethod(AtmVolCalib);
		Calibrate(); /// calibrate model vol and local models !

		/// Compute exercise frontier projection w.r.t. swap rate
		ARM_GenSecurityPtr genSec = GetGenSecurity();
		bool savePayoffsFlag = genSec->GetOtherPayoffsFlag();
		genSec->SetOtherPayoffsFlag(true);
		ARM_GenPricer* genPricer = new ARM_GenPricer( &*genSec,&*GetPricingModel() );
		ARM_AutoCleaner<ARM_GenPricer> HoldGP(genPricer);
		genPricer->Price();

		ARM_VectorPtr frontier = genPricer->GetPricerInfo()->GetContents(LocalCapFloorCRAColNamesTable[localCfFrontier]).GetData("Intermediateprices").GetVector();
		ARM_GP_Vector* exerciseDates = itsExerDateStrip->GetResetDates() ;
		size_t sizeFrontier = frontier->size();
		size_t offset = exerciseDates->size()-sizeFrontier;

		//deleted in refValue destructor
		ARM_Vector* refValueX = new ARM_Vector(sizeFrontier);
		ARM_Vector* refValueY = new ARM_Vector(sizeFrontier);
		for(i=0;i<sizeFrontier;i++)
		{
			(*refValueX)[i]=(*exerciseDates)[i+offset];
			(*refValueY)[i]=(*frontier)[i];
		}
		ARM_ReferenceValue strikeFrontier(refValueX,refValueY);
		
		/// Set diagonal swaption strike to frontier swap rates
		double frontierStrike,strike,strikeStep;
		size_t j,nbStrikeSteps=20;
		for(i=0;i<nbSwaptions;++i)
		{
			swaption = static_cast< ARM_Swaption* >(swaptionPF->GetAsset(i));
			atmStrike = swaption->GetStrike();
			frontierStrike = strikeFrontier.Interpolate(swaption->GetExpiryDate().GetJulian())*100.;
			swaption->UpdateStrike(frontierStrike);
			price = swaption->ComputePrice();
			vega = swaption->ComputeSensitivity(K_VEGA)/100.;
			if(vega<VEGA_MIN_TO_SELECT)
			{
				/// Find a more convenient strike !
				strikeStep = (frontierStrike-atmStrike)/nbStrikeSteps;
				for(j=0,strike=frontierStrike;j<nbStrikeSteps;++j,strike-=strikeStep)
				{
					swaption->UpdateStrike(strike);
					price = swaption->ComputePrice();
					vega = swaption->ComputeSensitivity(K_VEGA)/100.;
					if(vega >= VEGA_MIN_TO_SELECT)
						break;
				}
				if(j>=nbStrikeSteps)
				{
					/// No chance back to ATM strike !
					swaption->UpdateStrike(atmStrike);
					price = swaption->ComputePrice();
				}
			}
			swaptionPF->SetPrice(price,i);
		}

		/// Restore previous values
		genSec->SetOtherPayoffsFlag(savePayoffsFlag);
		SetCalibMethod(saveCalibMethod);
	}
	else
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, " ARM_CRASpreadCalculator::SetDiagStrikeToFrontier: should contain BS smiled");
	}
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: CreateEmptyCalibrationFor2F
///	Returns: void
///	Action : create the calibration methods
/////////////////////////////////////////////////////////////////
void ARM_CRASpreadCalculator::CreateEmptyCalibrationFor2F()
{
	ARM_GP_Vector calibTimes(1, 0.00);
	ARM_GP_Vector volatility(1, SIGMA_DEFAULT_VALUE);
	ARM_GP_Vector lowerVolatility(1, SIGMA_LOWER_BOUND);
	ARM_GP_Vector upperVolatility(1, SIGMA_UPPER_BOUND);
	ARM_CurveModelParam* volParam = new ARM_CurveModelParam(ARM_ModelParamType::Volatility, &volatility, &calibTimes ,"","STEPUPRIGHT",&lowerVolatility,&upperVolatility);
	ARM_ModelParamVector volParamVec(1,volParam);

	calibTimes.push_back(1.0); // to switch H&W2F params to time dependent version
	ARM_GP_Vector volatilityRatio(2, VOLRATIO_DEFAULT_VALUE);
	ARM_CurveModelParam* volRatioParam = new ARM_CurveModelParam(ARM_ModelParamType::VolatilityRatio, &volatilityRatio, &calibTimes ,"","STEPUPRIGHT");
	ARM_ModelParamVector volRatioParamVec(1,volRatioParam);
	
	ARM_GP_Vector correl(2, CORREL_DEFAULT_VALUE);
	ARM_CurveModelParam* correlParam = new ARM_CurveModelParam(ARM_ModelParamType::Correlation, &correl, &calibTimes ,"","STEPUPRIGHT");
	ARM_ModelParamVector correlParamVec(1,correlParam);
	

	ARM_CalibMethodPtr calibMethod(NULL);
	size_t N=3; /// Newton-Raphson max iter is used here to input this number

	ARM_StdPortfolioPtr pf1; 
	ARM_StdPortfolioPtr pf2; 
	ARM_StdPortfolioPtr pf3; 
	
	/// Create the build-in H&W2F calibMethod : portfolio must
	/// be made of swaption (VNS or bullet) a spread option (actual
	/// or degenerated with coef2=0 to get CMS caplet)

	/// The calibMethod link order is used if calibration failed
	/// : the more internal PF is then skipped
	/// Link order must be pf3 <-L- pf2 <-L- pf1

	if (IsDbleCorridor() && (itsCalibType == DIAG_CALIBRATION ||
							 itsCalibType == BASKET_CALIBRATION_SIMPLIFIED) )
	{
		pf1 = CreateSwaptionPortfolio();	//swaptionPF
		pf2 = CreateSOPortfolio(true);		//so1stCondPF
		pf3 = CreateSOPortfolio(false);		//caplet2ndCondPF
	}
	else if (itsCalibType==DIAG_SPREAD_LONG)
	{
		pf1 = CreateSwaptionPortfolio();	//swaptionPF
		pf2 = CreateSOPortfolio(true);
		pf3 = CreateSOPortfolio(false,false,true);
	}
	else if (itsCalibType==DIAG_SPREAD_SHORT)
	{
		pf1 = CreateSwaptionPortfolio();	//swaptionPF
		pf2 = CreateSOPortfolio(true);
		pf3 = CreateSOPortfolio(false,false,false);
	}
	else if (itsCalibType==DIAG_SPREAD_INDEX)
	{
		pf1 = CreateSwaptionPortfolio();	//swaptionPF
		pf2 = CreateSOPortfolio(true);
		pf3 = CreateSOPortfolio(false);
	}
	else if (itsCalibType==SHORT_LONG_SPREAD)
	{
		pf1 = CreateSOPortfolio(false,false,false);
		pf2 = CreateSOPortfolio(false,false,true);
		pf3 = CreateSOPortfolio(true);
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, " ARM_LocalCSOCalculator::CreateEmptyCalibration: invalid calibration type");

	int aux = (itsBootstrapOptimizer?10:0);


	ARM_CalibMethod calibVolMethod(pf3,
						volParamVec,ARM_CalibMethodType::HW2FOnly,aux+N,ARM_CalibrationTarget::UnknownTarget,
						NULL,NULL,false,0,1,false);

	/// calibVolMethod is cloned by ARM_CalibMethod (not shared by default)
	ARM_CalibMethod calibVolRatioMethod(pf2,
						volRatioParamVec,ARM_CalibMethodType::HW2FOnly,aux+N,ARM_CalibrationTarget::UnknownTarget,
						&calibVolMethod,NULL,false,0,1,false);

	/// calibVolRatioMethod is cloned by ARM_CalibMethod (not shared by default)
	calibMethod = ARM_CalibMethodPtr ( new ARM_CalibMethod (pf1,
						correlParamVec,ARM_CalibMethodType::HW2FOnly,aux+N,ARM_CalibrationTarget::UnknownTarget,
						&calibVolRatioMethod,NULL,false,0,1,false) );
	
	SetCalibMethod(calibMethod);

	delete volRatioParam;
	delete correlParam;
	delete volParam;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: CreateAndSetModel
///	Returns: void
///	Action : create the encapsulated model
/////////////////////////////////////////////////////////////////
void ARM_CRASpreadCalculator::CreateAndSetModel()
{
	/// Get yield curve
	ARM_ZeroCurve* curve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

	///---------------------------------------------------------
	/// Build HW model
	///---------------------------------------------------------
	ARM_PricingModel* HWmodel = NULL;
	ARM_GP_Vector defaultTimes(1,0.0);
	ARM_GP_Vector defaultSigmas(1,SIGMA_DEFAULT_VALUE);

	ARM_CurveModelParam volParam( ARM_ModelParamType::Volatility, &defaultSigmas, &defaultTimes );
	ARM_CurveModelParam* mrsParam = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(GetKeys()[MrsKey]));

	if(!mrsParam || mrsParam->GetType() != ARM_ModelParamType::MeanReversion)
		ARM_THROW( ERR_INVALID_ARGUMENT, "CSO calculator : invalid MeanReversion parameter");

	if (itsModelType == ARM_PricingModelType::HWM1F)
	{
		ARM_ModelParamVector paramVector(0);
		paramVector.push_back( &volParam );
		paramVector.push_back( mrsParam );
		
		ARM_ModelParamsHW1FStd params(paramVector);
		/// params are cloned in ARM_PricingModel constructor...
		HWmodel = new ARM_HullWhite1F( CreateClonedPtr( curve ), &params ) ;

	}
	else if (itsModelType == ARM_PricingModelType::HWM2F)
	{
		ARM_CurveModelParam* volRatioParam = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(GetKeys()[VolRatioKey]));
		
		if(!volRatioParam || volRatioParam->GetType() != ARM_ModelParamType::VolatilityRatio)
			ARM_THROW( ERR_INVALID_ARGUMENT, "CSO calculator : invalid VolatilityRatio parameter");

		ARM_CurveModelParam* mrsSpreadParam = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(GetKeys()[MrsSpreadKey]));
		
		if(!mrsSpreadParam || mrsSpreadParam->GetType() != ARM_ModelParamType::MeanReversionSpread)
			ARM_THROW( ERR_INVALID_ARGUMENT, "CSO calculator : invalid MeanReversionSpread parameter");

		ARM_CurveModelParam* correlParam = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(GetKeys()[CorrelKey]));
		
		if(!correlParam || correlParam->GetType() != ARM_ModelParamType::Correlation)
			ARM_THROW( ERR_INVALID_ARGUMENT, "CSO calculator : invalid Correlation parameter");
		
		ARM_ModelParamVector paramVector(0);
		paramVector.push_back( &volParam );
		paramVector.push_back( mrsParam );
		paramVector.push_back( volRatioParam );
		paramVector.push_back( mrsSpreadParam );
		paramVector.push_back( correlParam );

        if(volRatioParam->GetCurve()->GetAbscisses().size() <= 1 && correlParam->GetCurve()->GetAbscisses().size() <= 1 )
		{
			/// Constant correlation and 2nd factor volatility
			ARM_ModelParamsHW2FStd params(paramVector);
			HWmodel = new ARM_HullWhite2F( CreateClonedPtr( curve ), &params );
		}
        else
		{
			/// Time dependent correlation and/or 2nd factor volatility
			ARM_ModelParamsHW2FExt params(paramVector);
			HWmodel = new ARM_HullWhite2F( CreateClonedPtr( curve ), &params );
		}
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "CorridorSpreadOption calculator : invalid Model Type");
	

	
	///---------------------------------------------------------
	/// Build the local models
	///---------------------------------------------------------
	ARM_IntVector paramTypes(2);
	ARM_IntVector paramNbSurfaces(2);
	paramTypes[0] = ARM_ModelParamType::ForwardAdjustment;
	paramTypes[1] = ARM_ModelParamType::Volatility;
	if(IsDbleCorridor())
	{
		paramTypes.push_back(ARM_ModelParamType::Shift);
		paramTypes.push_back(ARM_ModelParamType::Correlation);

		paramNbSurfaces[0] = ARM_Local_Normal_Model::DBLECOR_ADJ_SIZE;
		paramNbSurfaces[1] = ARM_Local_Normal_Model::DBLECOR_VOL_SIZE;
		paramNbSurfaces.push_back(ARM_Local_Normal_Model::DBLECOR_KADJ_SIZE);
		paramNbSurfaces.push_back(ARM_Local_Normal_Model::DBLECOR_CORREL_SIZE);
	}
	else
	{
		paramTypes.push_back(ARM_ModelParamType::Shift);

		paramNbSurfaces[0] = ARM_Local_Normal_Model::SOCOR_ADJ_SIZE;
		paramNbSurfaces[1] = ARM_Local_Normal_Model::SOCOR_VOL_SIZE;

		paramNbSurfaces.push_back(ARM_Local_Normal_Model::SOCOR_KADJ_SIZE);

	}
	ARM_Local_Normal_ModelParams* defaultParams = ARM_Local_Normal_Model::CreateDefaultModelParams(paramTypes,paramNbSurfaces);
	ARM_Local_Model* LocalCorridorModel = new ARM_Local_Normal_Model (ARM_ZeroCurvePtr((ARM_ZeroCurve*)curve->Clone()), *defaultParams,paramTypes,paramNbSurfaces);
	ARM_Local_Model* LocalCorridor2Model = NULL;
	ARM_Local_Model* LocalCorridor3Model= NULL;
	if (IsTripleRange())
	{
		LocalCorridor2Model = new ARM_Local_Normal_Model (ARM_ZeroCurvePtr((ARM_ZeroCurve*)curve->Clone()), *defaultParams,paramTypes,paramNbSurfaces);
		LocalCorridor3Model = new ARM_Local_Normal_Model (ARM_ZeroCurvePtr((ARM_ZeroCurve*)curve->Clone()), *defaultParams,paramTypes,paramNbSurfaces);
	}
	delete defaultParams;

	///---------------------------------------------------------
	/// Build MultiAssets
	///---------------------------------------------------------
	size_t mamoSize = IsTripleRange()?4:2;
	
	ARM_StringVector names (mamoSize);
	vector<ARM_PricingModelPtr> models (mamoSize);
	ARM_StringVectorVector depends(mamoSize);

	char* ccy	= GetCcy().GetCcyName();

	names[0] =  ccy;
	names[1] =	LOCAL_CORRIDOR_MODEL_NAME;

	models[0] = ARM_PricingModelPtr( HWmodel );
	models[1] = ARM_PricingModelPtr( LocalCorridorModel );
		
	depends[0] = ARM_StringVector (1,names[0]);
	depends[1] = ARM_StringVector (1,names[0]);

	if (IsTripleRange())
	{
		names[2] =	LOCAL_CORRIDOR2_MODEL_NAME;
		names[3] =	LOCAL_CORRIDOR3_MODEL_NAME;

		models[2] = ARM_PricingModelPtr( LocalCorridor2Model );
		models[3] = ARM_PricingModelPtr( LocalCorridor3Model );

		depends[2] = ARM_StringVector (1,names[0]);
		depends[3] = ARM_StringVector (1,names[0]);
	}
	
	if(IsBasis())
	{
		names.resize(mamoSize+1);
		models.resize(mamoSize+1);
        depends.resize(mamoSize+1);
        names[mamoSize] = GetKeys()[YcBasis];

		// Build the Basis forward margin model
	    ARM_ZeroCurve* basisCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcBasis]));
		models[mamoSize] = ARM_PricingModelPtr(static_cast< ARM_PricingModel* >(new ARM_ForwardMarginBasis( CreateClonedPtr(basisCurve))));
		depends[mamoSize] = ARM_StringVector (1,names[0]);
	}
	
	ARM_ModelNameMap modelMap(names, models, depends);

	ARM_PricingModelPtr refModel (new ARM_MultiAssetsModel ( &modelMap) );
	
				
	///---------------------------------------------------------
	/// Set numerical method
	///---------------------------------------------------------
	/// compute nb steps for tree
	ARM_GP_Vector* exerciseDates = itsExerDateStrip->GetResetDates();
	double firstEventTime = (*exerciseDates)[0] - asOfDate;
	double lastEventTime = (*exerciseDates)[exerciseDates->size()-1] - asOfDate;

	size_t nbStepPerYear = TREE_NBSTEPS_PER_YEAR_1D;
	size_t nbStepBefore1st = 1;
	size_t nbDim = 1;
	int smootherType=ARM_SmootherBase::DoNothing;
	int nbSteps = static_cast<int>(floor(nbStepPerYear*lastEventTime/K_YEAR_LEN));

	int schedulerType=ARM_SchedulerBase::ConstantVarianceMeanReverting;
	ARM_GP_Vector schedulerDatas;

	if(itsModelType == ARM_PricingModelType::HWM2F)
	{
		if(IsDbleCorridor())
		{
			nbStepBefore1st	= TREE_NBSTEPS_BEFORE_1ST;
/***
			int nbStepsST = static_cast<int>(floor(TREE_NBSTEPS_PER_YEAR_ST_2D*(lastEventTime<TREE_NBSTEPS_SHORT_TERM ? lastEventTime : TREE_NBSTEPS_SHORT_TERM)/K_YEAR_LEN));
			int nbStepsMT = static_cast<int>(lastEventTime>TREE_NBSTEPS_SHORT_TERM ? floor(TREE_NBSTEPS_PER_YEAR_2D*(lastEventTime-TREE_NBSTEPS_SHORT_TERM)/K_YEAR_LEN) : 0);

			nbSteps = nbStepBefore1st + nbStepsST + nbStepsMT ;
***/

			schedulerType=ARM_SchedulerBase::MultiRegime;
			schedulerDatas.resize(8);
			schedulerDatas[0] = nbStepBefore1st;
			schedulerDatas[1] = 1000;
			schedulerDatas[2] = TREE_NBSTEPS_PER_YEAR_BEF_OPT;
			schedulerDatas[3] = TREE_NBSTEPS_OPTIMAL_YF;
			schedulerDatas[4] = TREE_NBSTEPS_PER_YEAR_BEF_OPT;
			schedulerDatas[5] = TREE_NBSTEPS_PER_YEAR_AFT_OPT;
			schedulerDatas[6] = TREE_NBSTEPS_LONGTERM_YF;
			schedulerDatas[7] = TREE_NBSTEPS_PER_YEAR_AFT_LT;
		}
		else
		{
			nbStepBefore1st	= TREE_NBSTEPS_BEFORE_1ST;
			nbSteps = static_cast<int>(floor(TREE_NBSTEPS_PER_YEAR_2D*lastEventTime/K_YEAR_LEN));
		}
		smootherType	= ARM_SmootherBase::Linear;
		nbDim = 2;
	}

	if(schedulerType == ARM_SchedulerBase::ConstantVarianceMeanReverting)
	{
		schedulerDatas.resize(3);
		schedulerDatas[0] = nbSteps;
		schedulerDatas[1] = nbStepBefore1st;
		schedulerDatas[2] = 1.0e-3;
	}
	int samplerType=ARM_SamplerBase::MeanReverting;
	ARM_GP_Vector samplerDatas(1,1.0e-3);
	int truncatorType=ARM_TruncatorBase::StandardDeviation;
	ARM_GP_Vector truncatorDatas(1,5.0);
	int reconnectorType=ARM_ReconnectorBase::Mean;
	bool probaFlag = itsPVAdjuster; // to allow spot probabilities price computation (... and PV calibration under terminal ZC numeraire)
	ARM_TreeBase* tree = ARM_TreeFactory.Instance()->CreateTreeND(nbDim,schedulerType,schedulerDatas,
		samplerType,samplerDatas,truncatorType,truncatorDatas,probaFlag,reconnectorType,smootherType);
	refModel->SetNumMethod( ARM_NumMethodPtr( tree ) );

	///---------------------------------------------------------
	/// Create a Numeraire and set it
	///---------------------------------------------------------
    ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::TerminalZc ) );
    refModel->SetNumeraire(numeraire);

	///---------------------------------------------------------
	/// Set the model
	///---------------------------------------------------------
	SetPricingModel(refModel);
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: UpdateModel
///	Returns: void
///	Action : update the model datas w.r.t. MarketDataManager
/////////////////////////////////////////////////////////////////
void ARM_CRASpreadCalculator::UpdateModel()
{
	CreateAndSetModel();

	/// Reset any previous prices
	ResetHasBeenPriced();
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRALocalCalculator
///	Routine: MiddleRows
///	Returns: ARM_RowInfo
///	Action : create a row of a deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_CRASpreadCalculator::MiddleRows(size_t eventIdx, 
                                          const ARM_DateStripCombiner& datesStructure) const
{
	bool flowByFlow = IsVms()?itsBootstrapOptimizer:FLOWBYFLOW;
	bool synchronize = IsVms()?!itsSwitchOldCalibration:SYNCHRONIZE;

	try
	{
		for(size_t firstEventIdx(0); firstEventIdx<itsExerSize; firstEventIdx++)
		{
			if(itsvIsExerDate[firstEventIdx]) 
				break;
		}
		bool isFirstCall = (eventIdx==firstEventIdx);

		ARM_GP_Vector* exerciseDates = GetExerDateStrip()->GetResetDates();
		double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
		size_t nbPastNoCall=0;
		while(nbPastNoCall < itsExerSize && (*exerciseDates)[nbPastNoCall] < asOfDate)
			++nbPastNoCall;

		bool isFirstEvent = (eventIdx==nbPastNoCall);

		//Useful datas
		string modelName(IsBasis() ?  GetKeys()[YcBasis] : string(GetCcy().GetCcyName()));
		double maturityDate = GetEndDate().GetJulian();
		int intPayRec		= GetPayRec();
		string payRec; 
		string invPayRec;
		if ( intPayRec == -1 )
		{
			payRec		= "P";
			invPayRec	= "R";
		}
		else
		{	payRec		= "R";
			invPayRec	= "P";
		}
		int intFundFreq		= GetFundFreq();
		string fundFreq		= ARM_ArgConvReverse_MatFrequency.GetString(intFundFreq);
		int intFundDayCount = GetFundDayCount();
		string fundDayCount = ARM_ArgConvReverse_DayCount.GetString(intFundDayCount);
		int intBoosIndexType = GetBoostedIndexType();
		string boostedIndexType;
		if (intBoosIndexType == Fixed)
			boostedIndexType = "FIXED";
		else if (intBoosIndexType == Libor)
			boostedIndexType = "LIBOR";
		else if (intBoosIndexType == Cms)
			boostedIndexType = "CMS";

		/// IsVariableCorridor is a test on pay index mults
		/// this will save computation time
		if (!IsVariableCorridor())
			boostedIndexType = "FIXED";
			

		int intCpnPayFreq			= GetCpnPayFreq();
		string cpnPayFreq			= ARM_ArgConvReverse_StdFrequency.GetString(intCpnPayFreq);
		string boostedVarTerm		= GetBoostedVarIndexTerm();
		int intBoosResetTiming		= GetBoostedResetTiming();
		string boostedResetTiming	= ARM_ArgConvReverse_Timing.GetString(intBoosResetTiming);
		int intBoosPayDayCount		= GetBoostedDayCount(); 
		string boostedPayDayCount	= ARM_ArgConvReverse_DayCount.GetString(intBoosPayDayCount);
		int cpnResetGap				= GetCpnResetGap();
		int intCpnIntRule			= GetBoostedIntRule();
		string cpnIntRule			= ARM_ArgConvReverse_IntRule.GetString(intCpnIntRule);
		int intCpnResetFreq			= GetCpnResetFreq();
		string cpnResetFreq			= ARM_ArgConvReverse_StdFrequency.GetString(intCpnResetFreq);
		int intCpnResetTiming		= GetCpnResetTiming();
		string cpnResetTiming		= ARM_ArgConvReverse_Timing.GetString(intCpnResetTiming);
		int intCpnIndexType1		= GetRefIndexType1();
		string cpnIndexType1		= ARM_ArgConvReverse_IndexClass.GetString(intCpnIndexType1);
		string cpnTerm1				= GetRefTerm1();
		int intCpnIndexDayCount1	= GetRefDayCount1();  
		string cpnIndexDayCount1	= ARM_ArgConvReverse_DayCount.GetString(intCpnIndexDayCount1);
		int intCpnIndexType2		= GetRefIndexType2();
		string cpnIndexType2		= ARM_ArgConvReverse_IndexClass.GetString(intCpnIndexType2);
		string cpnTerm2				= GetRefTerm2();
		string cpnIndexDayCount2	= cpnIndexDayCount1;

		/// Double corridor rate index
		string cpnIndexType3		= ARM_ArgConvReverse_IndexClass.GetString(itsRefIndexType3);
		string cpnTerm3				= itsRefTerm3;

		//Number of columns
		size_t descSize = GetDescSize();

		size_t callEventSize = datesStructure.GetDateStrip(CALL_CRASPREAD_SCHED)->GetResetDates()->size();
		vector< string > rowDescVec(descSize);
		vector< ARM_GP_VALUE_TYPE > rowTypeVec(descSize, ARM_MISSING_TYPE); 
		
		InitPriceableColumns(rowDescVec,rowTypeVec);
		
		//Nb Non Call Fees
		int fwdStart;
		int j = 0;
		int feeSize = itsCallFees.size();
		double currFee = itsCallFees.GetDiscreteValues()->Elt(j);
		while ((currFee >= NON_CALL_FEE) && (j < feeSize-1))
		while ((j < feeSize) && (currFee >= NON_CALL_FEE))
		{
			j++;
			if (j < feeSize)
				currFee = itsCallFees.GetDiscreteValues()->Elt(j);
		}
		fwdStart = j;

		//Weights for exotic swap
		double wCorridor = 1.0;
		double wFunding  = 1.0;

		//NOTICE DATE
		double	noticeDate = (*(datesStructure.GetDateStrip(CALL_CRASPREAD_SCHED)->GetResetDates()))[eventIdx];
		CC_Ostringstream noticeDateDesc;
		noticeDateDesc << CC_NS(std,fixed) << noticeDate;
		rowDescVec[localCfResetDate] = noticeDateDesc.str();
		rowTypeVec[localCfResetDate] = ARM_DATE_TYPE;

		//START DATE / CALL DATE
		CC_Ostringstream startDateDesc;
		if ( flowByFlow )
		{
			double startDate = (*itsExerDateStripUnadj->GetFlowStartDates())[eventIdx] ;
			startDateDesc << CC_NS(std,fixed) << startDate;
			rowDescVec[localCfStartDate] = startDateDesc.str();
			rowTypeVec[localCfStartDate] = ARM_DATE_TYPE;
		}
		else
		{
			double startDate = (*(datesStructure.GetDateStrip(CALL_CRASPREAD_SCHED)->GetFlowStartDates()))[eventIdx] ;
			startDateDesc << CC_NS(std,fixed) << startDate;
			rowDescVec[localCfStartDate] = startDateDesc.str();
			rowTypeVec[localCfStartDate] = ARM_DATE_TYPE;
		}

		//MATURITY DATE
		CC_Ostringstream maturityDateDesc;
		if ( flowByFlow )
		{
			//double endDate = (*(datesStructure.GetDateStrip(CALL_CRASPREAD_SCHED)->GetFlowEndDates()))[eventIdx] ;
			double endDate = (*itsExerDateStripUnadj->GetFlowEndDates())[eventIdx] ;
			maturityDateDesc << CC_NS(std,fixed) << endDate;
		}
		else
		{
			maturityDateDesc << CC_NS(std,fixed) << maturityDate;
		}
		rowDescVec[localCfMaturityDate] = maturityDateDesc.str();
		rowTypeVec[localCfMaturityDate] = ARM_DATE_TYPE;

		//FEES
		double fees = const_cast< ARM_ReferenceValue& >(itsCallFees).Interpolate(noticeDate);
		CC_Ostringstream feesDesc;
		feesDesc << CC_NS(std,fixed) << fees;
		rowDescVec[localCfFees] = feesDesc.str();
		rowTypeVec[localCfFees] = ARM_DOUBLE;

		//FUNDING LEG
		CC_Ostringstream fundingLegDesc;
		if( flowByFlow )
		{
			fundingLegDesc << "SWAP(" << modelName << ","<< LocalCapFloorCRAColNamesTable[localCfStartDate] << "[i],"<< LocalCapFloorCRAColNamesTable[localCfMaturityDate] << "[i], 0," << payRec << ", ";
			fundingLegDesc << fundFreq << ", " << fundDayCount << ", " << fundFreq << ", " << fundDayCount << ", ";
			if(eventIdx == callEventSize-1)
				fundingLegDesc << "FundSpreadCst, NotioCst" << ")";
			else
				fundingLegDesc << "FundSpreadCst, NotioCst" << ")+PV(FundingLeg[i+1])";
			rowDescVec[localCfFundingLeg] = fundingLegDesc.str();
			rowTypeVec[localCfFundingLeg] = ARM_STRING;
		}
		else
		{
			fundingLegDesc << "SWAP(" << modelName << ","<< LocalCapFloorCRAColNamesTable[localCfStartDate] << "[i],"<< LocalCapFloorCRAColNamesTable[localCfMaturityDate] << "[i], 0," << payRec << ", ";
			fundingLegDesc << fundFreq << ", " << fundDayCount << ", " << fundFreq << ", " << fundDayCount << ", ";
			fundingLegDesc << "FundSpreadCst, NotioCst" << ")";
			rowDescVec[localCfFundingLeg] = fundingLegDesc.str();
			rowTypeVec[localCfFundingLeg] = ARM_STRING;
		}	

		//CORRIDOR LEG1
		CC_Ostringstream corridorLeg1Desc;
		if( flowByFlow )
		{
			double sd = (*(datesStructure.GetDateStrip(CALL_CRASPREAD_SCHED)->GetFlowStartDates()))[eventIdx];

			int vmsIndex;
			if (IsVms())
				vmsIndex = (int) const_cast< ARM_ReferenceValue& >(itsTenor).Interpolate(sd);
			else
				vmsIndex = itsRefIndex1;

			int intVmsType;
			string vmsTerm = FromIndexTypeToTermAndType(vmsIndex, intVmsType);
			string vmsType = ARM_ArgConvReverse_IndexClass.GetString(intVmsType);
			 
			corridorLeg1Desc << "CORRIDOR(" << "LOCALCAPFLOOR" << ",";			// Model
			corridorLeg1Desc << LocalCapFloorCRAColNamesTable[localCfStartDate] << "[i],";				// StartDate
			corridorLeg1Desc << LocalCapFloorCRAColNamesTable[localCfMaturityDate] << "[i],";			// EndDate
			corridorLeg1Desc << payRec << ",";						// PayRec
			corridorLeg1Desc << boostedIndexType << ",";			// PayIndexType
			corridorLeg1Desc << "FixRateCst" << ",";				// Fix Value
			corridorLeg1Desc << "PayIdxMultCst" << ",";				// PayIndexMult Value
			corridorLeg1Desc << cpnPayFreq << ",";					// Payment Freq
			corridorLeg1Desc << boostedVarTerm << ",";				// Pay Index Term
			corridorLeg1Desc << boostedResetTiming << ",";			// Pay Index Timing
			corridorLeg1Desc << boostedPayDayCount << ",";			// Pay Day Count
			corridorLeg1Desc << cpnResetGap << ",";					// Pay Gap
			corridorLeg1Desc << cpnIntRule << ",";					// Int Rule
			corridorLeg1Desc << "CpnSpreadCst" << ",";				// Spread
			corridorLeg1Desc <<  cpnResetFreq << ",";				// Fixing Freq
			corridorLeg1Desc << cpnResetTiming << ",";				// Fixing Timing
			corridorLeg1Desc << vmsType << ",";						// VMS
			corridorLeg1Desc << vmsTerm << ",";						// VMS
			corridorLeg1Desc << "Coef1Cst" << ",";					// Coeff 1
			corridorLeg1Desc << "DownCst" << ",";			// Barrier Down
			corridorLeg1Desc << "UpCst" << ",";				// Barrier Up
			corridorLeg1Desc << "NotioCst" << ",";					// Notional
			corridorLeg1Desc << cpnIndexType2 << ",";				// Fixing Index Type 2
			corridorLeg1Desc << cpnTerm2 << ",";					// Fixing Index Term 2
			corridorLeg1Desc << "Coef2Cst";							// Coeff 2
			corridorLeg1Desc << ",,,,";

			if (synchronize)
			{
				int pastPeriods=0;
				int pastFixings=0;
				while (itsvCpnIndex[eventIdx]>pastFixings)
				{
					pastFixings+=GetIndexVector()->Elt(pastPeriods);
					pastPeriods++;
				}

				corridorLeg1Desc << ",dsStruct,dsPay,dsFix,idx,"<<pastPeriods<<",1";
			}
			if(eventIdx == callEventSize-1)
				corridorLeg1Desc << ")";
			else
				corridorLeg1Desc << ")+PV(CorridorLeg1[i+1])";
		}
		else
		{
			string term1;
			if (IsVms())
				term1 = "TenorCst";
			else
				term1 = cpnTerm1;

			corridorLeg1Desc << "CORRIDOR(" << "LOCALCAPFLOOR" << ",";			// Model
			corridorLeg1Desc << LocalCapFloorCRAColNamesTable[localCfStartDate] << "[i],";				// StartDate
			corridorLeg1Desc << LocalCapFloorCRAColNamesTable[localCfMaturityDate] << "[i],";			// EndDate
			corridorLeg1Desc << payRec << ",";						// PayRec
			corridorLeg1Desc << boostedIndexType << ",";			// PayIndexType
			corridorLeg1Desc << "FixRateCst" << ",";				// Fix Value
			corridorLeg1Desc << "PayIdxMultCst" << ",";				// PayIndexMult Value
			corridorLeg1Desc << cpnPayFreq << ",";					// Payment Freq
			corridorLeg1Desc << boostedVarTerm << ",";				// Pay Index Term
			corridorLeg1Desc << boostedResetTiming << ",";			// Pay Index Timing
			corridorLeg1Desc << boostedPayDayCount << ",";			// Pay Day Count
			corridorLeg1Desc << cpnResetGap << ",";					// Pay Gap
			corridorLeg1Desc << cpnIntRule << ",";					// Int Rule
			corridorLeg1Desc << "CpnSpreadCst" << ",";				// Spread
			corridorLeg1Desc <<  cpnResetFreq << ",";				// Fixing Freq
			corridorLeg1Desc << cpnResetTiming << ",";				// Fixing Timing
			corridorLeg1Desc << cpnIndexType1 << ",";				// Fixing Index Type 1
			corridorLeg1Desc << term1 << ",";					// Fixing Index Term 1
			corridorLeg1Desc << "Coef1Cst" << ",";					// Coeff 1
			if(IsDbleCorridor())
			{
				corridorLeg1Desc << "SDownCst" << ",";				// Lower barrier for the spread condition
				corridorLeg1Desc << "SUpCst" << ",";				// Upper barrier for the spread condition
			}
			else
			{
				corridorLeg1Desc << "DownCst" << ",";		// Barrier Down
				corridorLeg1Desc << "UpCst" << ",";			// Barrier Up
			}
			corridorLeg1Desc << "NotioCst" << ",";					// Notional
			corridorLeg1Desc << cpnIndexType2 << ",";				// Fixing Index Type 2
			corridorLeg1Desc << cpnTerm2 << ",";					// Fixing Index Term 2
			corridorLeg1Desc << "Coef2Cst";							// Coeff 2

			if(IsDbleCorridor())
			{
				corridorLeg1Desc << "," << cpnIndexType3 << ",";	// Fixing Index Type 3
				corridorLeg1Desc << cpnTerm3 << ",";				// Fixing Index Term 3
				corridorLeg1Desc << "RDownCst,";					// Lower barrier for the rate condition
				corridorLeg1Desc << "RUpCst";						// Upper barrier for the rate condition
			}
			else
			{
				corridorLeg1Desc << ",,,,";
			}

			if (synchronize)
			{
				int pastPeriods=0;
				int pastFixings=0;
				while (itsvCpnIndex[eventIdx]>pastFixings)
				{
					pastFixings+=GetIndexVector()->Elt(pastPeriods);
					pastPeriods++;
				}

				corridorLeg1Desc << ",dsStruct,dsPay,dsFix,idx,"<<pastPeriods<<")";
			}
			else
				corridorLeg1Desc << ")";
		}

		rowDescVec[localCfCorridorLeg1] = corridorLeg1Desc.str();
		rowTypeVec[localCfCorridorLeg1] = ARM_STRING;

		if (IsTripleRange())
		{

			//CORRIDOR LEG2
			CC_Ostringstream corridorLeg2Desc;
			corridorLeg2Desc << "CORRIDOR(" << "LOCALCAPFLOOR2" << ",";			// Model
			corridorLeg2Desc << LocalCapFloorCRAColNamesTable[localCfStartDate] << "[i],";				// StartDate
			corridorLeg2Desc << LocalCapFloorCRAColNamesTable[localCfMaturityDate] << "[i],";			// EndDate
			corridorLeg2Desc << payRec << ",";						// PayRec
			corridorLeg2Desc << boostedIndexType << ",";			// PayIndexType
			corridorLeg2Desc << "FixRate2Cst" << ",";				// Fix Value
			corridorLeg2Desc << "PayIdxMultCst" << ",";				// PayIndexMult Value
			corridorLeg2Desc << cpnPayFreq << ",";					// Payment Freq
			corridorLeg2Desc << boostedVarTerm << ",";				// Pay Index Term
			corridorLeg2Desc << boostedResetTiming << ",";			// Pay Index Timing
			corridorLeg2Desc << boostedPayDayCount << ",";			// Pay Day Count
			corridorLeg2Desc << cpnResetGap << ",";					// Pay Gap
			corridorLeg2Desc << cpnIntRule << ",";					// Int Rule
			corridorLeg2Desc << "CpnSpreadCst" << ",";				// Spread
			corridorLeg2Desc <<  cpnResetFreq << ",";				// Fixing Freq
			corridorLeg2Desc << cpnResetTiming << ",";				// Fixing Timing
			corridorLeg2Desc << cpnIndexType1 << ",";				// Fixing Index Type 1
			corridorLeg2Desc << cpnTerm1 << ",";					// Fixing Index Term 1
			corridorLeg2Desc << "Coef1Cst" << ",";					// Coeff 1
			if(IsDbleCorridor())
			{
				corridorLeg2Desc << "SDownCst" << ",";				// Lower barrier for the spread condition
				corridorLeg2Desc << "SUpCst" << ",";				// Upper barrier for the spread condition
			}
			else
			{
				corridorLeg2Desc << "Down2Cst" << ",";		// Barrier Down
				corridorLeg2Desc << "Up2Cst" << ",";			// Barrier Up
			}
			corridorLeg2Desc << "NotioCst" << ",";					// Notional
			corridorLeg2Desc << cpnIndexType2 << ",";				// Fixing Index Type 2
			corridorLeg2Desc << cpnTerm2 << ",";					// Fixing Index Term 2
			corridorLeg2Desc << "Coef2Cst";							// Coeff 2

			if(IsDbleCorridor())
			{
				corridorLeg2Desc << "," << cpnIndexType3 << ",";	// Fixing Index Type 3
				corridorLeg2Desc << cpnTerm3 << ",";				// Fixing Index Term 3
				corridorLeg2Desc << "RDownCst,";					// Lower barrier for the rate condition
				corridorLeg2Desc << "RUpCst" << ")";				// Upper barrier for the rate condition
			}
			else
			{
				corridorLeg2Desc << ")";
			}

			
			rowDescVec[localCfCorridorLeg2] = corridorLeg2Desc.str();
			rowTypeVec[localCfCorridorLeg2] = ARM_STRING;

			//CORRIDOR LEG3
			CC_Ostringstream corridorLeg3Desc;
			corridorLeg3Desc << "CORRIDOR(" << "LOCALCAPFLOOR3" << ",";			// Model
			corridorLeg3Desc << LocalCapFloorCRAColNamesTable[localCfStartDate] << "[i],";				// StartDate
			corridorLeg3Desc << LocalCapFloorCRAColNamesTable[localCfMaturityDate] << "[i],";			// EndDate
			corridorLeg3Desc << payRec << ",";						// PayRec
			corridorLeg3Desc << boostedIndexType << ",";			// PayIndexType
			corridorLeg3Desc << "FixRate3Cst" << ",";				// Fix Value
			corridorLeg3Desc << "PayIdxMultCst" << ",";				// PayIndexMult Value
			corridorLeg3Desc << cpnPayFreq << ",";					// Payment Freq
			corridorLeg3Desc << boostedVarTerm << ",";				// Pay Index Term
			corridorLeg3Desc << boostedResetTiming << ",";			// Pay Index Timing
			corridorLeg3Desc << boostedPayDayCount << ",";			// Pay Day Count
			corridorLeg3Desc << cpnResetGap << ",";					// Pay Gap
			corridorLeg3Desc << cpnIntRule << ",";					// Int Rule
			corridorLeg3Desc << "CpnSpreadCst" << ",";				// Spread
			corridorLeg3Desc <<  cpnResetFreq << ",";				// Fixing Freq
			corridorLeg3Desc << cpnResetTiming << ",";				// Fixing Timing
			corridorLeg3Desc << cpnIndexType1 << ",";				// Fixing Index Type 1
			corridorLeg3Desc << cpnTerm1 << ",";					// Fixing Index Term 1
			corridorLeg3Desc << "Coef1Cst" << ",";					// Coeff 1
			if(IsDbleCorridor())
			{
				corridorLeg3Desc << "SDownCst" << ",";				// Lower barrier for the spread condition
				corridorLeg3Desc << "SUpCst" << ",";				// Upper barrier for the spread condition
			}
			else
			{
				corridorLeg3Desc << "Down3Cst" << ",";		// Barrier Down
				corridorLeg3Desc << "Up3Cst" << ",";			// Barrier Up
			}
			corridorLeg3Desc << "NotioCst" << ",";					// Notional
			corridorLeg3Desc << cpnIndexType2 << ",";				// Fixing Index Type 2
			corridorLeg3Desc << cpnTerm2 << ",";					// Fixing Index Term 2
			corridorLeg3Desc << "Coef2Cst";							// Coeff 2

			if(IsDbleCorridor())
			{
				corridorLeg3Desc << "," << cpnIndexType3 << ",";	// Fixing Index Type 3
				corridorLeg3Desc << cpnTerm3 << ",";				// Fixing Index Term 3
				corridorLeg3Desc << "RDownCst,";					// Lower barrier for the rate condition
				corridorLeg3Desc << "RUpCst" << ")";				// Upper barrier for the rate condition
			}
			else
			{
				corridorLeg3Desc << ")";
			}

			
			rowDescVec[localCfCorridorLeg3] = corridorLeg3Desc.str();
			rowTypeVec[localCfCorridorLeg3] = ARM_STRING;

			CC_Ostringstream corridorLegDesc;
			corridorLegDesc << "CorridorLeg1[i]+CorridorLeg2[i]+CorridorLeg3[i]";
			rowDescVec[localCfCorridorLeg] = corridorLegDesc.str();
			rowTypeVec[localCfCorridorLeg] = ARM_STRING;
		}
		else
		{
			CC_Ostringstream corridorLeg2Desc;
			corridorLeg2Desc << "0";
			rowDescVec[localCfCorridorLeg2] = corridorLeg2Desc.str();
			rowTypeVec[localCfCorridorLeg2] = ARM_STRING;

			CC_Ostringstream corridorLeg3Desc;
			corridorLeg3Desc << "0";
			rowDescVec[localCfCorridorLeg3] = corridorLeg3Desc.str();
			rowTypeVec[localCfCorridorLeg3] = ARM_STRING;

			CC_Ostringstream corridorLegDesc;
			corridorLegDesc << "CorridorLeg1[i]";
			rowDescVec[localCfCorridorLeg] = corridorLegDesc.str();
			rowTypeVec[localCfCorridorLeg] = ARM_STRING;
		}

		//EXOTIC SWAP
		CC_Ostringstream ExoticSwapDesc;
		ExoticSwapDesc << wCorridor << "*CorridorLeg[i]+" << wFunding << "*FundingLeg[i]";
		rowDescVec[localCfExoticSwap] = ExoticSwapDesc.str();
		rowTypeVec[localCfExoticSwap] = ARM_STRING;   

		//OPTION & EXERCISE INDICATOR
		CC_Ostringstream optionDesc;
		CC_Ostringstream indicExerDesc;
		if(eventIdx == callEventSize-1)
		{
			optionDesc << "MAX(ExoticSwap[i]-Fees[i],0)";
			indicExerDesc << "IF(ExoticSwap[i]-Fees[i]>0,1,0)";
		}
		else
		{
			optionDesc << "Exercise(0,ExoticSwap[i]-Fees[i],Option[i+1])";
			indicExerDesc << "IF(ExoticSwap[i]-Fees[i]>PV(Option[i+1]),1,0)";
		}
		rowDescVec[localCfOption] = optionDesc.str();
		rowTypeVec[localCfOption] = ARM_STRING;

		rowDescVec[localCfIndicExer] = indicExerDesc.str();
		rowTypeVec[localCfIndicExer] = ARM_STRING;

		//BERMUDA
		CC_Ostringstream bermudaDesc;
		if(isFirstCall)
		{
			if (itsOptionPortfolio)
			{
				// this information is available only when the calculator is created from summit portfolio !
				if (itsOptionPortfolio->GetPorS() == K_PAY)
					bermudaDesc << "-Option[i]";
				else
					bermudaDesc << "Option[i]";
			}
			else if (itsSwaption)
			{
				if (itsSwaption->GetPorS() == K_PAY)
					bermudaDesc << "-Option[i]";
				else
					bermudaDesc << "Option[i]";
			}
			else
			{
				bermudaDesc << "Option[i]";
			}
			rowDescVec[localCfBermuda] = bermudaDesc.str();
			rowTypeVec[localCfBermuda] = ARM_STRING;
		}
		else
		{
			bermudaDesc << "0";
			rowDescVec[localCfBermuda] = bermudaDesc.str();
			rowTypeVec[localCfBermuda] = ARM_DOUBLE;
		}

		size_t pIdx,probaOffset,nbProbas = CC_Min(callEventSize-itsExerciseProbaOffset,LocalCapFloorCRAExerciseProbaNb);
		if(itsIsExerciseProbas && eventIdx+1 <= itsExerciseProbaOffset + nbProbas)
		{
			//NEXT NOTICE DATE
			double	nextNoticeDate = (*(datesStructure.GetDateStrip(CALL_CRASPREAD_SCHED)->GetResetDates()))[CC_Min(eventIdx+1,callEventSize-1)];
			CC_Ostringstream nextNoticeDateDesc;
			nextNoticeDateDesc << CC_NS(std,fixed) << nextNoticeDate;
			rowDescVec[localCfNextResetDate] = nextNoticeDateDesc.str();
			rowTypeVec[localCfNextResetDate] = ARM_DATE_TYPE;

			/// EXERCISE PROBABILITY
			if(eventIdx==0)
			{
				probaOffset = localCfProba1;
				for(pIdx=0;pIdx<nbProbas;++pIdx,probaOffset+=2)
				{
					CC_Ostringstream probaDesc;
					probaDesc << LocalCapFloorCRALocalExerciseProbaNamesTable[pIdx] << "[i]*DF("
							  << modelName << "," << LocalCapFloorCRAColNamesTable[localCfMaturityDate] << "[i])";
					rowDescVec[probaOffset] = probaDesc.str();
					rowTypeVec[probaOffset] = ARM_STRING;
				}
			}

			/// LOCAL PROBABILITY
			probaOffset = localCfLocalProba1;
			for(pIdx=0;pIdx < nbProbas;++pIdx,probaOffset+=2)
			{
				CC_Ostringstream localProbaDesc;
				if(itsExerciseProbaOffset+pIdx==eventIdx)
				{
					/// Proba computation initialisation : exercise or not at the current line
					localProbaDesc << LocalCapFloorCRAColNamesTable[localCfIndicExer] << "[i]";
				}
				else if(eventIdx < itsExerciseProbaOffset+pIdx)
				{
					localProbaDesc	<< "MAX(" << LocalCapFloorCRAColNamesTable[localCfIndicExer] << "[i]"
									<< ",PV(" << LocalCapFloorCRALocalExerciseProbaNamesTable[pIdx] << "[i+1])/DF("
									<< modelName << "," << LocalCapFloorCRAColNamesTable[localCfNextResetDate] << "[i]))";
				}
				rowDescVec[probaOffset]	= localProbaDesc.str();
				rowTypeVec[probaOffset]	= ARM_STRING;
			}
		}

		CC_Ostringstream corridorDesc;
		CC_Ostringstream fundingDesc;
		if(isFirstEvent)
		{
			//CORRIDOR 
			corridorDesc << "CorridorLeg[i]";
			rowDescVec[localCfCorridor] = corridorDesc.str();
			rowTypeVec[localCfCorridor] = ARM_STRING;

			//FUNDING 
			fundingDesc << "FundingLeg[i]";
			rowDescVec[localCfFunding] = fundingDesc.str();
			rowTypeVec[localCfFunding] = ARM_STRING;
		}
		else
		{
			corridorDesc << "0";
			rowDescVec[localCfCorridor] = corridorDesc.str();
			rowTypeVec[localCfCorridor] = ARM_DOUBLE;
		
			fundingDesc << "0";
			rowDescVec[localCfFunding] = fundingDesc.str();
			rowTypeVec[localCfFunding] = ARM_DOUBLE;
		}

		//FWD CORRIDOR / FUNDING
		CC_Ostringstream fwdCorridorDesc;
		CC_Ostringstream fwdFundingDesc;
		if (isFirstCall)
		{
			fwdCorridorDesc << "CorridorLeg[i]";
			rowDescVec[localCfFwdCorridor] = fwdCorridorDesc.str();
			rowTypeVec[localCfFwdCorridor] = ARM_STRING;

			fwdFundingDesc << "FundingLeg[i]";
			rowDescVec[localCfFwdFunding] = fwdFundingDesc.str();
			rowTypeVec[localCfFwdFunding] = ARM_STRING;
		}
		else
		{
			fwdCorridorDesc << "0";
			rowDescVec[localCfFwdCorridor] = fwdCorridorDesc.str();
			rowTypeVec[localCfFwdCorridor] = ARM_DOUBLE;

			CC_Ostringstream fwdFundingDesc;
			fwdFundingDesc << "0";
			rowDescVec[localCfFwdFunding] = fwdFundingDesc.str();
			rowTypeVec[localCfFwdFunding] = ARM_DOUBLE;
		}


		// SWAPRATE for exercise frontier projection
		CC_Ostringstream swaprateDesc;
		swaprateDesc << "SWAPRATE(" << modelName << ",";
		swaprateDesc << LocalCapFloorCRAColNamesTable[localCfStartDate] << "[i],";
		swaprateDesc << LocalCapFloorCRAColNamesTable[localCfMaturityDate] << "[i])";
		
		rowDescVec[localCfExerSwapRate] = swaprateDesc.str();
		rowTypeVec[localCfExerSwapRate] = ARM_STRING;

		// FRONTIER
		CC_Ostringstream frontierDesc;
		if(eventIdx == callEventSize-1)
		{
			frontierDesc << "FRONTIER(" << LocalCapFloorCRAColNamesTable[localCfExoticSwap] << "[i]";
			frontierDesc << "-" << LocalCapFloorCRAColNamesTable[localCfFees] << "[i],";
			frontierDesc << "0," << LocalCapFloorCRAColNamesTable[localCfExerSwapRate] << "[i])";
		}
		else
		{
			frontierDesc << "FRONTIER(" << LocalCapFloorCRAColNamesTable[localCfExoticSwap] << "[i]";
			frontierDesc << "-" << LocalCapFloorCRAColNamesTable[localCfFees] << "[i],";
			frontierDesc << LocalCapFloorCRAColNamesTable[localCfOption] << "[i+1],";
			frontierDesc << LocalCapFloorCRAColNamesTable[localCfExerSwapRate] << "[i])";
		}
		rowDescVec[localCfFrontier] = frontierDesc.str();
		rowTypeVec[localCfFrontier] = ARM_STRING;

	
		return ARM_RowInfo(rowDescVec,rowTypeVec);
	}
	catch (Exception& x)
	{
	    x.DebugPrint();

		throw x;
	}
	catch (...)
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in CraLocalCalculator::MiddleRows" );
	}
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: Calibrate
///	Returns: void
///	Action : cfalibrate the model 
/////////////////////////////////////////////////////////////////
void ARM_CRASpreadCalculator::Calibrate()
{
//ARM_Timer timer;
//timer.ClockStartTime();
	///-------------------------------------
	/// HW calibration
	///-------------------------------------
	ARM_MultiAssetsModel*  refModel = static_cast<ARM_MultiAssetsModel*>(&*GetPricingModel());
	char* ccy	= GetCcy().GetCcyName();
	
	ARM_PricingModelPtr HWmodel = (*refModel->GetModelMap())[ccy]->Model();
	if(itsSwoptCalib)
	{
		GetCalibMethod()->Calibrate(&*HWmodel);
	}

	///-------------------------------------
	/// Local models calibration
	///-------------------------------------
	ARM_PricingModelPtr model;
	model = (*refModel->GetModelMap())[LOCAL_CORRIDOR_MODEL_NAME]->Model();
	ARM_Local_Model* LocalCorridorModel = static_cast<ARM_Local_Model*> (&*model);

	/// eval times
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	ARM_GP_Vector* exerciseDates = itsExerDateStrip->GetResetDates() ;
	ARM_GP_Vector evalTimes(0);
	
	for (size_t i(0); i<exerciseDates->size(); i++)
			evalTimes.push_back((*exerciseDates)[i] - asOfDate);
	
	LocalCorridorModel->SetVolUnSqueezer(itsVolUnSqueezer);
	LocalCorridorModel->SetCorrelUnSqueezer(itsCorrelUnSqueezer); // only used if double condition
	LocalCorridorModel->SetPVAdjuster(itsPVAdjuster);

	LocalCorridorModel->SetResetCalib(true);
	LocalCorridorModel->CalibrateLocalModel(*itsCSOPF, evalTimes);

	if (IsTripleRange())
	{
		ARM_PricingModelPtr model2;
		model2 = (*refModel->GetModelMap())[LOCAL_CORRIDOR2_MODEL_NAME]->Model();
		ARM_Local_Model* LocalCorridor2Model = static_cast<ARM_Local_Model*> (&*model2);
		LocalCorridor2Model->SetResetCalib(true);
		LocalCorridor2Model->CalibrateLocalModel(*itsCSOPF2, evalTimes);

		ARM_PricingModelPtr model3;
		model3 = (*refModel->GetModelMap())[LOCAL_CORRIDOR3_MODEL_NAME]->Model();
		ARM_Local_Model* LocalCorridor3Model = static_cast<ARM_Local_Model*> (&*model3);
		LocalCorridor3Model->SetResetCalib(true);
		LocalCorridor3Model->CalibrateLocalModel(*itsCSOPF3, evalTimes);
	}
//timer.ClockEndTime();
//FILE* f=fopen("c:\\temp\\dumpCRA2Timer.txt","a");
//fprintf(f,"Calib Duration = %10.5lf ms\n",timer.GetDuration()*1000.0);
//fclose(f);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: CreateVanillaArgCSO
///	Returns: ARM_VanillaSpreadOptionArg
///	Action : create  a VanillaArg from an ARM_SpreadOption and generate the fix leg and float leg ofthe CMSs 
/////////////////////////////////////////////////////////////////
ARM_VanillaSpreadOptionArg* ARM_CRASpreadCalculator::CreateVanillaArgCSO(ARM_SpreadOption* spreadOption) 
{
	double asOfDate	  = GetMktDataManager()->GetAsOfDate().GetJulian();
	ARM_Currency* ccy = GetCurrencyUnit();
	ARM_VanillaSpreadOptionArg* soVanillaArg = (ARM_VanillaSpreadOptionArg*)ARM_ConverterFromKernel::ConvertSecuritytoArgObject(spreadOption,asOfDate);
	soVanillaArg->ComputeIndexSchedulesAndAdjustDates(ccy, asOfDate);
	return soVanillaArg;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: ComputeCRA2NormalSensitivities
///	Returns: void
///	Action : simplified version of sensitivities computation
///			 for double corridor underlying : only fixed
///			 payment flows are taking into account
///			 If rate paid is variable, forward rates are
///			 used (but not sensitivities of the index to the forward
///			 yield curve : next step...)
/////////////////////////////////////////////////////////////////
void ARM_CRASpreadCalculator::ComputeCRA2NormalSensitivities()
{
	if( !(itsCalibType == BASKET_CALIBRATION_SIMPLIFIED ||
			(itsCalibStrikeType[0] == EQUIVALENT &&
				(itsCalibType == DIAG_CALIBRATION ||
				 itsCalibType == DIAG_SPREAD_LONG ||
				 itsCalibType == DIAG_SPREAD_INDEX||
				 itsCalibType == DIAG_SPREAD_SHORT)))	)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : only simple basket or diag based + equiv strike calibration type allowed" );

	itsFixedCpnLongSensi.resize(itsCpnSize);
	itsFixedCpnShortSensi.resize(itsCpnSize);
	itsFixedCpnValue.resize(itsCpnSize);
	itsVarCpnLongSensi.resize(itsCpnSize);
	itsVarCpnShortSensi.resize(itsCpnSize);
	itsVarCpnValue.resize(itsCpnSize);
	itsPayIndexFwd.resize(itsCpnSize);

	ARM_CorridorDblCondition* raSDownRDown	= static_cast< ARM_CorridorDblCondition* >(itsCSOPF->GetAsset(ARM_Local_Normal_Model::DBLECOR_SDOWN_RDOWN));
	ARM_CorridorDblCondition* raSDownRUp	= static_cast< ARM_CorridorDblCondition* >(itsCSOPF->GetAsset(ARM_Local_Normal_Model::DBLECOR_SDOWN_RUP));
	ARM_CorridorDblCondition* raSUpRDown	= static_cast< ARM_CorridorDblCondition* >(itsCSOPF->GetAsset(ARM_Local_Normal_Model::DBLECOR_SUP_RDOWN));
	ARM_CorridorDblCondition* raSUpRUp		= static_cast< ARM_CorridorDblCondition* >(itsCSOPF->GetAsset(ARM_Local_Normal_Model::DBLECOR_SUP_RUP));
	bool isSDownRDown	= raSDownRDown && !(raSDownRUp || raSUpRDown || raSUpRUp);
	bool isSDownRUp		= raSDownRUp && !(raSDownRDown || raSUpRDown || raSUpRUp);
	bool isSUpRDown		= raSUpRDown && !(raSDownRDown || raSDownRUp || raSUpRUp);
	bool isSUpRUp		= raSUpRUp && !(raSDownRDown || raSDownRUp || raSUpRDown);
	for(size_t i=0;i<itsCpnSize;++i)
	{
		itsFixedCpnLongSensi[i]		= 0.0;
		itsFixedCpnShortSensi[i]	= 0.0;
		itsVarCpnLongSensi[i]		= 0.0;
		itsVarCpnShortSensi[i]		= 0.0;

		if(isSUpRUp)
		{
			itsFixedCpnValue[i]		= raSUpRUp->GetSpreadRateProba()[i];
			itsVarCpnValue[i]		= raSUpRUp->GetSpreadRateProbaFLT()[i];
			itsPayIndexFwd[i]		= 0.01 * (*(raSUpRUp->GetPayFwdRateFLT()))[i];
		}
		else if(isSDownRUp)
		{
			itsFixedCpnValue[i]		= raSDownRUp->GetRateProba()[i] - raSDownRUp->GetSpreadRateProba()[i];
			itsVarCpnValue[i]		= raSDownRUp->GetRateProbaFLT()[i] - raSDownRUp->GetSpreadRateProbaFLT()[i];
			itsPayIndexFwd[i]		= 0.01 * (*(raSDownRUp->GetPayFwdRateFLT()))[i];
		}
		else if(isSUpRDown)
		{
			itsFixedCpnValue[i]		= raSUpRDown->GetSpreadProba()[i] - raSUpRDown->GetSpreadRateProba()[i];
			itsVarCpnValue[i]		= raSUpRDown->GetSpreadProbaFLT()[i] - raSUpRDown->GetSpreadRateProbaFLT()[i];
			itsPayIndexFwd[i]		= 0.01 * (*(raSUpRDown->GetPayFwdRateFLT()))[i];
		}
		else if(isSDownRDown)
		{
			itsFixedCpnValue[i]		= 1 - raSDownRDown->GetSpreadProba()[i]
									- raSDownRDown->GetRateProba()[i] + raSDownRDown->GetSpreadRateProba()[i];
			itsVarCpnValue[i]		= 1 - raSDownRDown->GetSpreadProbaFLT()[i]
									- raSDownRDown->GetRateProbaFLT()[i] + raSDownRDown->GetSpreadRateProbaFLT()[i];
			itsPayIndexFwd[i]		= 0.01 * (*(raSDownRDown->GetPayFwdRateFLT()))[i];
		}
		else
		{
			/// General case with at least one non degenerated range
			itsFixedCpnValue[i]			= raSUpRUp->GetSpreadRateProba()[i] - raSUpRDown->GetSpreadRateProba()[i]
										- raSDownRUp->GetSpreadRateProba()[i] + raSDownRDown->GetSpreadRateProba()[i];

			itsVarCpnValue[i]			= raSUpRUp->GetSpreadRateProbaFLT()[i] - raSUpRDown->GetSpreadRateProbaFLT()[i]
										- raSDownRUp->GetSpreadRateProbaFLT()[i] + raSDownRDown->GetSpreadRateProbaFLT()[i];

			itsPayIndexFwd[i]			= 0.01 * (*(raSDownRDown->GetPayFwdRateFLT()))[i];
		}

	}
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: ComputeSONormalSensitivities
///	Returns: void
///	Action : compute so cap/floor cpn sensis w.r.t. long & short CMS
///			 under the NORMAL ASSUMPTION on spread dynamics
///			 NB: So portfolio prices are supposed to be already computed
/////////////////////////////////////////////////////////////////
void ARM_CRASpreadCalculator::ComputeCSONormalSensitivities()
{
	/// ---- NOTE ------
	/// itsFixedCpnLongSensi[i] and itsFixedCpnShortSensi[i] are the sensitivities
	/// of coupon[i] with respect to long & short CMS rates. 
	/// >>> IT DOES NOT TAKE INTO ACCOUNT notional, interest period and payment discount factor !!!
	/// Same thing for itsFixedCpnValue
	/// 
	ARM_BSModel* SOBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );
	ARM_ZeroCurve* curve   = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));

	ARM_VanillaSpreadOptionArg* soVanillaArg = static_cast<ARM_VanillaSpreadOptionArg*>(&*itsVanillaArgCSO);
	ARM_GP_Vector* payPeriods = soVanillaArg->GetPayPeriods();
	
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

	itsFixedCpnLongSensi.resize(itsCpnSize);
	itsFixedCpnShortSensi.resize(itsCpnSize);
	itsFixedCpnValue.resize(itsCpnSize);
	itsVarCpnLongSensi.resize(itsCpnSize);
	itsVarCpnShortSensi.resize(itsCpnSize);
	itsVarCpnValue.resize(itsCpnSize);
	itsPayIndexFwd.resize(itsCpnSize);

	// We take the first spread option in the porfolio which should contain 2 products
	// A CAP on the down barrier
	ARM_SpreadOption* csoDown = static_cast< ARM_SpreadOption* >(itsCSOPF->GetAsset(0));
	// A FLOOR on the up barrier
	ARM_SpreadOption* csoUp = static_cast< ARM_SpreadOption* >(itsCSOPF->GetAsset(1));

/// TRIPLE --------------------------------------------------------------------------
	ARM_SpreadOption* csoDown2;
	ARM_SpreadOption* csoUp2;
	ARM_SpreadOption* csoDown3;
	ARM_SpreadOption* csoUp3;

	if (IsTripleRange())
	{
		csoDown2 = static_cast< ARM_SpreadOption* >(itsCSOPF2->GetAsset(0));
		csoUp2 = static_cast< ARM_SpreadOption* >(itsCSOPF2->GetAsset(1));

		csoDown3 = static_cast< ARM_SpreadOption* >(itsCSOPF3->GetAsset(0));
		csoUp3 = static_cast< ARM_SpreadOption* >(itsCSOPF3->GetAsset(1));

		itsFixedCpnLongSensi2.resize(itsCpnSize);
		itsFixedCpnShortSensi2.resize(itsCpnSize);
		itsFixedCpnValue2.resize(itsCpnSize);
		itsVarCpnLongSensi2.resize(itsCpnSize);
		itsVarCpnShortSensi2.resize(itsCpnSize);
		itsVarCpnValue2.resize(itsCpnSize);
		itsPayIndexFwd2.resize(itsCpnSize);
		itsFixedValues2.resize(itsCpnSize);

		itsFixedCpnLongSensi3.resize(itsCpnSize);
		itsFixedCpnShortSensi3.resize(itsCpnSize);
		itsFixedCpnValue3.resize(itsCpnSize);
		itsVarCpnLongSensi3.resize(itsCpnSize);
		itsVarCpnShortSensi3.resize(itsCpnSize);
		itsVarCpnValue3.resize(itsCpnSize);
		itsPayIndexFwd3.resize(itsCpnSize);
		itsFixedValues3.resize(itsCpnSize);
	}
/// TRIPLE --------------------------------------------------------------------------
	
	ARM_Vector* periodResetDates =  csoDown->GetPayIndexLeg()->GetResetDates();
	ARM_Vector* periodPayDates   = csoDown->GetPayIndexLeg()->GetPaymentDates();
	ARM_Vector* periodStartDates =  csoDown->GetPayIndexLeg()->GetFlowStartDates();
	ARM_Vector* flowPayDates = csoDown->GetSpreadLeg()->GetFirstLeg()->GetPaymentDates();

	const double NSTDEV_NO_IMPLIED_VOL = 6.0;
			
	for (size_t i (0); i<itsCpnSize; i++)
	{
		double payDate=(*flowPayDates)[i];
		int periodIndex			= periodPayDates->find(payDate);

		double shortFwd	 = 0.01 * csoDown->GetSpreadLeg()->GetFirstLeg()->GetFwdRates()->Elt(i);
		double longFwd   = 0.01 * csoDown->GetSpreadLeg()->GetSecondLeg()->GetFwdRates()->Elt(i);
		double correl	 = 0.01 * csoDown->GetCorrelVector()->Elt(i);
		double coeffShort= csoDown->GetWeight1();
		double coeffLong = csoDown->GetWeight2();
		double optMat    = (GetStructDateStrip()->GetResetDates()->Elt(i) - asOfDate)/K_YEAR_LEN;

		if (optMat < 0)
			continue;
		
		// Warning : SO oplet price is computed with Strike and Margin interpolated on RESET dates
		double strikeDown = 0.01 * csoDown->GetStrikes()->CptReferenceValue((*periodResetDates)[periodIndex]);
		double strikeUp   = 0.01 * csoUp->GetStrikes()->CptReferenceValue((*periodResetDates)[periodIndex]);
		double fixedRate  = csoDown->GetPayIndexMargins()->CptReferenceValue((*periodResetDates)[periodIndex]);

/// TRIPLE --------------------------------------------------------------------------
		double strikeDown2,strikeUp2,fixedRate2,strikeDown3,strikeUp3,fixedRate3;
		if (IsTripleRange())
		{
			strikeDown2 = 0.01 * csoDown2->GetStrikes()->CptReferenceValue((*periodResetDates)[periodIndex]);
			strikeUp2   = 0.01 * csoUp2->GetStrikes()->CptReferenceValue((*periodResetDates)[periodIndex]);
			fixedRate2  = csoDown2->GetPayIndexMargins()->CptReferenceValue((*periodResetDates)[periodIndex]);
			itsFixedValues2[i] = fixedRate2;

			strikeDown3 = 0.01 * csoDown3->GetStrikes()->CptReferenceValue((*periodResetDates)[periodIndex]);
			strikeUp3   = 0.01 * csoUp3->GetStrikes()->CptReferenceValue((*periodResetDates)[periodIndex]);
			fixedRate3  = csoDown3->GetPayIndexMargins()->CptReferenceValue((*periodResetDates)[periodIndex]);
			itsFixedValues3[i] = fixedRate3;
		}
/// TRIPLE --------------------------------------------------------------------------

		double spread1 = 0.01*csoDown->GetSpread1();
		double spread2 = 0.01*csoDown->GetSpread2();

		double longTenor  = (soVanillaArg->GetSwapLongFloatEndTime()->Elt(i) - soVanillaArg->GetSwapLongFloatStartTime()->Elt(i)) / K_YEAR_LEN;
		double shortTenor = (soVanillaArg->GetSwapShortFloatEndTime()->Elt(i) - soVanillaArg->GetSwapShortFloatStartTime()->Elt(i)) / K_YEAR_LEN;
					
		double atmShortVol  =  0.01 * csoDown->GetVol1ATM(optMat, shortTenor);
		double atmLongVol   =  0.01 * csoDown->GetVol2ATM(optMat, longTenor);
		atmShortVol	*= coeffShort * shortFwd ;
		atmLongVol	*= coeffLong  * longFwd;
		double atmSpreadVol =  sqrt( atmShortVol * atmShortVol + atmLongVol * atmLongVol - 2.* correl * atmShortVol * atmLongVol );
		double atmSpreadStdev = atmSpreadVol * sqrt(optMat);
		
		double forwardSpread = coeffLong * longFwd - coeffShort * shortFwd;

		double payDf		= curve->DiscountPrice((csoDown->GetSwapLeg()->GetPaymentDates()->Elt(i)-asOfDate)/K_YEAR_LEN);
		double payNotional	= csoDown->GetPayIndexLeg()->GetAmount()->CptReferenceValue(csoDown->GetSwapLeg()->GetPaymentDates()->Elt(i));
		double payPeriod	= (*payPeriods)[i];

		/// take vols ATM
		double downVolLeft  = atmSpreadVol;
		double downVolRight = atmSpreadVol;
		double upVolLeft	= atmSpreadVol;
		double upVolRight	= atmSpreadVol;

		double downVolLeftVar  = atmSpreadVol;
		double downVolRightVar = atmSpreadVol;
		double upVolLeftVar	   = atmSpreadVol;
		double upVolRightVar   = atmSpreadVol;
				
		bool isVariable = IsVariableCpn(i);

		double forwardSpreadVarUp   (0.0);
		double forwardSpreadVarDown (0.0);
		
		if (isVariable)
		{
			double shortFwdVarDown	= 0.01 * csoDown->GetFirstFwdRateFLT()->Elt(i);
			double longFwdVarDown	= 0.01 * csoDown->GetSecondFwdRateFLT()->Elt(i);
			double shortFwdVarUp	= 0.01 * csoUp->GetFirstFwdRateFLT()->Elt(i);
			double longFwdVarUp		= 0.01 * csoUp->GetSecondFwdRateFLT()->Elt(i);
			forwardSpreadVarUp		= coeffLong * longFwdVarUp - coeffShort * shortFwdVarUp;
			forwardSpreadVarDown	= coeffLong * longFwdVarDown - coeffShort * shortFwdVarDown;
		}

/// TRIPLE ---------------------------------------------------
		double forwardSpreadVarUp2   (0.0);
		double forwardSpreadVarDown2 (0.0);
		double forwardSpreadVarUp3   (0.0);
		double forwardSpreadVarDown3 (0.0);
		
		if (IsTripleRange())
		{
			if (isVariable)
			{
				double shortFwdVarDown	= 0.01 * csoDown2->GetFirstFwdRateFLT()->Elt(i);
				double longFwdVarDown	= 0.01 * csoDown2->GetSecondFwdRateFLT()->Elt(i);
				double shortFwdVarUp	= 0.01 * csoUp2->GetFirstFwdRateFLT()->Elt(i);
				double longFwdVarUp		= 0.01 * csoUp2->GetSecondFwdRateFLT()->Elt(i);
				forwardSpreadVarUp2		= coeffLong * longFwdVarUp - coeffShort * shortFwdVarUp;
				forwardSpreadVarDown2	= coeffLong * longFwdVarDown - coeffShort * shortFwdVarDown;

				shortFwdVarDown			= 0.01 * csoDown3->GetFirstFwdRateFLT()->Elt(i);
				longFwdVarDown			= 0.01 * csoDown3->GetSecondFwdRateFLT()->Elt(i);
				shortFwdVarUp			= 0.01 * csoUp3->GetFirstFwdRateFLT()->Elt(i);
				longFwdVarUp			= 0.01 * csoUp3->GetSecondFwdRateFLT()->Elt(i);
				forwardSpreadVarUp3		= coeffLong * longFwdVarUp - coeffShort * shortFwdVarUp;
				forwardSpreadVarDown3	= coeffLong * longFwdVarDown - coeffShort * shortFwdVarDown;
			}
		}
/// FIN TRIPLE ---------------------------------------------------
		
							
		/// overwrite only if not too far from ATM
		if (fabs(strikeDown - forwardSpread)<NSTDEV_NO_IMPLIED_VOL * atmSpreadStdev)
		{
			int callPut = csoDown->IsCap()-csoDown->IsFloor();
			double targetPriceLeft_N = csoDown->GetCap1Prices()->Elt(i);
			downVolLeft = VanillaImpliedVol_N (forwardSpread, targetPriceLeft_N, strikeDown+spread1, optMat, callPut, &atmSpreadVol);
			double targetPriceRight_N = csoDown->GetCap2Prices()->Elt(i);
			downVolRight = VanillaImpliedVol_N (forwardSpread, targetPriceRight_N, strikeDown+spread2, optMat, callPut, &atmSpreadVol);

			if (isVariable)
			{
				targetPriceLeft_N = csoDown->GetCap1PricesFLT()->Elt(i);
				downVolLeftVar = VanillaImpliedVol_N (forwardSpreadVarDown, targetPriceLeft_N,  strikeDown+spread1, optMat, callPut, &atmSpreadVol);
				targetPriceRight_N = csoDown->GetCap2PricesFLT()->Elt(i);
				downVolRightVar = VanillaImpliedVol_N (forwardSpreadVarDown, targetPriceRight_N, strikeDown+spread2, optMat, callPut, &atmSpreadVol);
			}
		}

		/// overwrite only if not too far from ATM
		if (fabs(strikeUp - forwardSpread)<NSTDEV_NO_IMPLIED_VOL * atmSpreadStdev)
		{
			int callPut = csoUp->IsCap()-csoUp->IsFloor();
			double targetPriceLeft_N = csoUp->GetCap1Prices()->Elt(i);
			upVolLeft = VanillaImpliedVol_N (forwardSpread, targetPriceLeft_N, strikeUp+spread1, optMat, callPut, &atmSpreadVol);
			double targetPriceRight_N = csoUp->GetCap2Prices()->Elt(i);
			upVolRight = VanillaImpliedVol_N (forwardSpread, targetPriceRight_N, strikeUp+spread2, optMat, callPut, &atmSpreadVol);

			if (isVariable)
			{
				targetPriceLeft_N = csoUp->GetCap1PricesFLT()->Elt(i);
				upVolLeftVar = VanillaImpliedVol_N (forwardSpreadVarUp, targetPriceLeft_N, strikeUp+spread1, optMat, callPut, &atmSpreadVol);
				targetPriceRight_N = csoUp->GetCap2PricesFLT()->Elt(i);
				upVolRightVar = VanillaImpliedVol_N (forwardSpreadVarUp, targetPriceRight_N, strikeUp+spread2, optMat, callPut, &atmSpreadVol);
			}
		}

/// TRIPLE -----------------------------------------------------------------
		double downVolLeft2  = atmSpreadVol;
		double downVolRight2 = atmSpreadVol;
		double upVolLeft2	= atmSpreadVol;
		double upVolRight2	= atmSpreadVol;
		double downVolLeftVar2  = atmSpreadVol;
		double downVolRightVar2 = atmSpreadVol;
		double upVolLeftVar2	   = atmSpreadVol;
		double upVolRightVar2  = atmSpreadVol;

		double downVolLeft3  = atmSpreadVol;
		double downVolRight3 = atmSpreadVol;
		double upVolLeft3	= atmSpreadVol;
		double upVolRight3	= atmSpreadVol;
		double downVolLeftVar3  = atmSpreadVol;
		double downVolRightVar3 = atmSpreadVol;
		double upVolLeftVar3	   = atmSpreadVol;
		double upVolRightVar3  = atmSpreadVol;

		if (IsTripleRange())
		{
			if (fabs(strikeDown2 - forwardSpread)<NSTDEV_NO_IMPLIED_VOL * atmSpreadStdev)
			{
				int callPut = csoDown2->IsCap()-csoDown2->IsFloor();
				double targetPriceLeft_N = csoDown2->GetCap1Prices()->Elt(i);
				downVolLeft2 = VanillaImpliedVol_N (forwardSpread, targetPriceLeft_N, strikeDown2+spread1, optMat, callPut, &atmSpreadVol);
				double targetPriceRight_N = csoDown2->GetCap2Prices()->Elt(i);
				downVolRight2 = VanillaImpliedVol_N (forwardSpread, targetPriceRight_N, strikeDown2+spread2, optMat, callPut, &atmSpreadVol);

				if (isVariable)
				{
					targetPriceLeft_N = csoDown2->GetCap1PricesFLT()->Elt(i);
					downVolLeftVar2 = VanillaImpliedVol_N (forwardSpreadVarDown2, targetPriceLeft_N,  strikeDown2+spread1, optMat, callPut, &atmSpreadVol);
					targetPriceRight_N = csoDown2->GetCap2PricesFLT()->Elt(i);
					downVolRightVar2 = VanillaImpliedVol_N (forwardSpreadVarDown2, targetPriceRight_N, strikeDown2+spread2, optMat, callPut, &atmSpreadVol);
				}
			}

			if (fabs(strikeUp2 - forwardSpread)<NSTDEV_NO_IMPLIED_VOL * atmSpreadStdev)
			{
				int callPut = csoUp2->IsCap()-csoUp2->IsFloor();
				double targetPriceLeft_N = csoUp2->GetCap1Prices()->Elt(i);
				upVolLeft2 = VanillaImpliedVol_N (forwardSpread, targetPriceLeft_N, strikeUp2+spread1, optMat, callPut, &atmSpreadVol);
				double targetPriceRight_N = csoUp2->GetCap2Prices()->Elt(i);
				upVolRight2 = VanillaImpliedVol_N (forwardSpread, targetPriceRight_N, strikeUp2+spread2, optMat, callPut, &atmSpreadVol);

				if (isVariable)
				{
					targetPriceLeft_N = csoUp2->GetCap1PricesFLT()->Elt(i);
					upVolLeftVar2 = VanillaImpliedVol_N (forwardSpreadVarUp2, targetPriceLeft_N, strikeUp2+spread1, optMat, callPut, &atmSpreadVol);
					targetPriceRight_N = csoUp2->GetCap2PricesFLT()->Elt(i);
					upVolRightVar2 = VanillaImpliedVol_N (forwardSpreadVarUp2, targetPriceRight_N, strikeUp2+spread2, optMat, callPut, &atmSpreadVol);
				}
			}

			if (fabs(strikeDown3 - forwardSpread)<NSTDEV_NO_IMPLIED_VOL * atmSpreadStdev)
			{
				int callPut = csoDown3->IsCap()-csoDown3->IsFloor();
				double targetPriceLeft_N = csoDown3->GetCap1Prices()->Elt(i);
				downVolLeft3 = VanillaImpliedVol_N (forwardSpread, targetPriceLeft_N, strikeDown3+spread1, optMat, callPut, &atmSpreadVol);
				double targetPriceRight_N = csoDown3->GetCap2Prices()->Elt(i);
				downVolRight3 = VanillaImpliedVol_N (forwardSpread, targetPriceRight_N, strikeDown3+spread2, optMat, callPut, &atmSpreadVol);

				if (isVariable)
				{
					targetPriceLeft_N = csoDown3->GetCap1PricesFLT()->Elt(i);
					downVolLeftVar3 = VanillaImpliedVol_N (forwardSpreadVarDown3, targetPriceLeft_N,  strikeDown3+spread1, optMat, callPut, &atmSpreadVol);
					targetPriceRight_N = csoDown3->GetCap2PricesFLT()->Elt(i);
					downVolRightVar3 = VanillaImpliedVol_N (forwardSpreadVarDown3, targetPriceRight_N, strikeDown3+spread2, optMat, callPut, &atmSpreadVol);
				}
			}

			if (fabs(strikeUp3 - forwardSpread)<NSTDEV_NO_IMPLIED_VOL * atmSpreadStdev)
			{
				int callPut = csoUp3->IsCap()-csoUp3->IsFloor();
				double targetPriceLeft_N = csoUp3->GetCap1Prices()->Elt(i);
				upVolLeft3 = VanillaImpliedVol_N (forwardSpread, targetPriceLeft_N, strikeUp3+spread1, optMat, callPut, &atmSpreadVol);
				double targetPriceRight_N = csoUp3->GetCap2Prices()->Elt(i);
				upVolRight3 = VanillaImpliedVol_N (forwardSpread, targetPriceRight_N, strikeUp3+spread2, optMat, callPut, &atmSpreadVol);

				if (isVariable)
				{
					targetPriceLeft_N = csoUp3->GetCap1PricesFLT()->Elt(i);
					upVolLeftVar3 = VanillaImpliedVol_N (forwardSpreadVarUp3, targetPriceLeft_N, strikeUp3+spread1, optMat, callPut, &atmSpreadVol);
					targetPriceRight_N = csoUp3->GetCap2PricesFLT()->Elt(i);
					upVolRightVar3 = VanillaImpliedVol_N (forwardSpreadVarUp3, targetPriceRight_N, strikeUp3+spread2, optMat, callPut, &atmSpreadVol);
				}
			}
	
		}

/// FIN TRIPLE -----------------------------------------------------------------


		///---------------------------------
		/// initial prices
		///---------------------------------
		/// it is safer to recompute them instead of taking the price
		/// already computed in the portfolio
		double sensi;
		double price;
		
		price  = (VanillaOption_N (forwardSpread, downVolLeft, strikeDown+spread1, optMat, K_CALL)-VanillaOption_N (forwardSpread, downVolRight, strikeDown+spread2, optMat, K_CALL))/(spread2-spread1);
		price += (VanillaOption_N (forwardSpread, upVolRight, strikeUp+spread2, optMat, K_CALL)-VanillaOption_N (forwardSpread, upVolLeft, strikeUp+spread1, optMat, K_CALL))/(spread2-spread1);
		
		
		///---------------------------------
		/// sensi w.r.t short CMS rate
		///---------------------------------
		sensi  = (VanillaOption_N (forwardSpread - coeffShort * SENSI_CMS, downVolLeft, strikeDown+spread1, optMat, K_CALL)-VanillaOption_N (forwardSpread - coeffShort * SENSI_CMS, downVolRight, strikeDown+spread2, optMat, K_CALL))/(spread2-spread1);
		sensi += (VanillaOption_N (forwardSpread - coeffShort * SENSI_CMS, upVolRight,  strikeUp+spread2,   optMat, K_CALL)-VanillaOption_N (forwardSpread - coeffShort * SENSI_CMS, upVolLeft,  strikeUp+spread1,   optMat, K_CALL))/(spread2-spread1);

		sensi-=price;

		sensi /= SENSI_CMS;

		itsFixedCpnShortSensi[i] = sensi;


		///---------------------------------
		/// sensi w.r.t long CMS rate
		///---------------------------------
		sensi  = (VanillaOption_N (forwardSpread + coeffLong * SENSI_CMS, downVolLeft, strikeDown+spread1, optMat, K_CALL)-VanillaOption_N (forwardSpread + coeffLong * SENSI_CMS, downVolRight, strikeDown+spread2, optMat, K_CALL))/(spread2-spread1);
		sensi += (VanillaOption_N (forwardSpread + coeffLong * SENSI_CMS, upVolRight,  strikeUp+spread2,   optMat, K_CALL)-VanillaOption_N (forwardSpread + coeffLong * SENSI_CMS, upVolLeft,  strikeUp+spread1,   optMat, K_CALL))/(spread2-spread1);

		sensi-=price;

		sensi /= SENSI_CMS;

		itsFixedCpnLongSensi[i] = sensi;


		///---------------------------------
		/// Cpn value
		///---------------------------------
		itsFixedCpnValue[i] = price;
				
		if (isVariable)
		{
			price  = (VanillaOption_N (forwardSpreadVarDown, downVolLeftVar, strikeDown+spread1, optMat, K_CALL)-VanillaOption_N (forwardSpreadVarDown, downVolRightVar, strikeDown+spread2, optMat, K_CALL))/(spread2-spread1);
			price += (VanillaOption_N (forwardSpreadVarUp, upVolRightVar, strikeUp+spread2, optMat, K_CALL)-VanillaOption_N (forwardSpreadVarUp, upVolLeftVar, strikeUp+spread1, optMat, K_CALL))/(spread2-spread1);
			
			
			///---------------------------------
			/// sensi w.r.t short CMS rate
			///---------------------------------
			sensi  = (VanillaOption_N (forwardSpreadVarDown - coeffShort * SENSI_CMS, downVolLeftVar, strikeDown+spread1, optMat, K_CALL)-VanillaOption_N (forwardSpreadVarDown - coeffShort * SENSI_CMS, downVolRightVar, strikeDown+spread2, optMat, K_CALL))/(spread2-spread1);
			sensi += (VanillaOption_N (forwardSpreadVarUp - coeffShort * SENSI_CMS, upVolRightVar,  strikeUp+spread2,   optMat, K_CALL)-VanillaOption_N (forwardSpreadVarUp - coeffShort * SENSI_CMS, upVolLeftVar,  strikeUp+spread1,   optMat, K_CALL))/(spread2-spread1);

			sensi-=price;

			sensi /= SENSI_CMS;

			itsVarCpnShortSensi[i] = sensi;


			///---------------------------------
			/// sensi w.r.t long CMS rate
			///---------------------------------
			sensi  = (VanillaOption_N (forwardSpreadVarDown + coeffLong * SENSI_CMS, downVolLeftVar, strikeDown+spread1, optMat, K_CALL)-VanillaOption_N (forwardSpreadVarDown + coeffLong * SENSI_CMS, downVolRightVar, strikeDown+spread2, optMat, K_CALL))/(spread2-spread1);
			sensi += (VanillaOption_N (forwardSpreadVarUp + coeffLong * SENSI_CMS, upVolRightVar,  strikeUp+spread2,   optMat, K_CALL)-VanillaOption_N (forwardSpreadVarUp + coeffLong * SENSI_CMS, upVolLeftVar,  strikeUp+spread1,   optMat, K_CALL))/(spread2-spread1);

			sensi-=price;

			sensi /= SENSI_CMS;

			itsVarCpnLongSensi[i] = sensi;


			///---------------------------------
			/// Cpn value
			///---------------------------------
			itsVarCpnValue[i] = price;


			///---------------------------------
			/// Fwd of payment index
			///---------------------------------
			itsPayIndexFwd[i] = 0.01 * csoDown->GetPayFwdRateFLT()->Elt(i);

		}
		else
		{
			itsVarCpnValue[i]      = 0.0;
			itsVarCpnLongSensi[i]  = 0.0;
			itsVarCpnShortSensi[i] = 0.0;
			itsPayIndexFwd[i]      = 0.0;
		}


/// TRIPLE ----------------------------------------------------------
		if (IsTripleRange())
		{
			price  = (VanillaOption_N (forwardSpread, downVolLeft2, strikeDown2+spread1, optMat, K_CALL)-VanillaOption_N (forwardSpread, downVolRight, strikeDown2+spread2, optMat, K_CALL))/(spread2-spread1);
			price += (VanillaOption_N (forwardSpread, upVolRight2, strikeUp2+spread2, optMat, K_CALL)-VanillaOption_N (forwardSpread, upVolLeft, strikeUp2+spread1, optMat, K_CALL))/(spread2-spread1);

			// short CMS
			sensi  = (VanillaOption_N (forwardSpread - coeffShort * SENSI_CMS, downVolLeft2, strikeDown2+spread1, optMat, K_CALL)-VanillaOption_N (forwardSpread - coeffShort * SENSI_CMS, downVolRight2, strikeDown2+spread2, optMat, K_CALL))/(spread2-spread1);
			sensi += (VanillaOption_N (forwardSpread - coeffShort * SENSI_CMS, upVolRight2,  strikeUp2+spread2,   optMat, K_CALL)-VanillaOption_N (forwardSpread - coeffShort * SENSI_CMS, upVolLeft2,  strikeUp2+spread1,   optMat, K_CALL))/(spread2-spread1);
			sensi-=price;
			sensi /= SENSI_CMS;
			itsFixedCpnShortSensi2[i] = sensi;

			// long CMS
			sensi  = (VanillaOption_N (forwardSpread + coeffLong * SENSI_CMS, downVolLeft2, strikeDown2+spread1, optMat, K_CALL)-VanillaOption_N (forwardSpread + coeffLong * SENSI_CMS, downVolRight2, strikeDown2+spread2, optMat, K_CALL))/(spread2-spread1);
			sensi += (VanillaOption_N (forwardSpread + coeffLong * SENSI_CMS, upVolRight2,  strikeUp2+spread2,   optMat, K_CALL)-VanillaOption_N (forwardSpread + coeffLong * SENSI_CMS, upVolLeft2,  strikeUp2+spread1,   optMat, K_CALL))/(spread2-spread1);
			sensi-=price;
			sensi /= SENSI_CMS;
			itsFixedCpnLongSensi2[i] = sensi;


			/// Cpn value
			itsFixedCpnValue2[i] = price;
					
			if (isVariable)
			{
				price  = (VanillaOption_N (forwardSpreadVarDown2, downVolLeftVar2, strikeDown2+spread1, optMat, K_CALL)-VanillaOption_N (forwardSpreadVarDown2, downVolRightVar2, strikeDown2+spread2, optMat, K_CALL))/(spread2-spread1);
				price += (VanillaOption_N (forwardSpreadVarUp2, upVolRightVar2, strikeUp2+spread2, optMat, K_CALL)-VanillaOption_N (forwardSpreadVarUp2, upVolLeftVar2, strikeUp2+spread1, optMat, K_CALL))/(spread2-spread1);
				
				/// short CMS
				sensi  = (VanillaOption_N (forwardSpreadVarDown2 - coeffShort * SENSI_CMS, downVolLeftVar2, strikeDown2+spread1, optMat, K_CALL)-VanillaOption_N (forwardSpreadVarDown2 - coeffShort * SENSI_CMS, downVolRightVar2, strikeDown2+spread2, optMat, K_CALL))/(spread2-spread1);
				sensi += (VanillaOption_N (forwardSpreadVarUp2 - coeffShort * SENSI_CMS, upVolRightVar2,  strikeUp2+spread2,   optMat, K_CALL)-VanillaOption_N (forwardSpreadVarUp2 - coeffShort * SENSI_CMS, upVolLeftVar2,  strikeUp2+spread1,   optMat, K_CALL))/(spread2-spread1);
				sensi-=price;
				sensi /= SENSI_CMS;
				itsVarCpnShortSensi2[i] = sensi;

				/// long CMS
				sensi  = (VanillaOption_N (forwardSpreadVarDown2 + coeffLong * SENSI_CMS, downVolLeftVar2, strikeDown2+spread1, optMat, K_CALL)-VanillaOption_N (forwardSpreadVarDown2 + coeffLong * SENSI_CMS, downVolRightVar2, strikeDown2+spread2, optMat, K_CALL))/(spread2-spread1);
				sensi += (VanillaOption_N (forwardSpreadVarUp2 + coeffLong * SENSI_CMS, upVolRightVar2,  strikeUp2+spread2,   optMat, K_CALL)-VanillaOption_N (forwardSpreadVarUp2 + coeffLong * SENSI_CMS, upVolLeftVar2,  strikeUp2+spread1,   optMat, K_CALL))/(spread2-spread1);
				sensi-=price;
				sensi /= SENSI_CMS;
				itsVarCpnLongSensi2[i] = sensi;

				itsVarCpnValue2[i] = price;
			}
			else
			{
				itsVarCpnValue2[i]      = 0.0;
				itsVarCpnLongSensi2[i]  = 0.0;
				itsVarCpnShortSensi2[i] = 0.0;
				itsPayIndexFwd2[i]      = 0.0;
			}

			price  = (VanillaOption_N (forwardSpread, downVolLeft3, strikeDown3+spread1, optMat, K_CALL)-VanillaOption_N (forwardSpread, downVolRight, strikeDown3+spread2, optMat, K_CALL))/(spread2-spread1);
			price += (VanillaOption_N (forwardSpread, upVolRight3, strikeUp3+spread2, optMat, K_CALL)-VanillaOption_N (forwardSpread, upVolLeft, strikeUp3+spread1, optMat, K_CALL))/(spread2-spread1);

			// short CMS
			sensi  = (VanillaOption_N (forwardSpread - coeffShort * SENSI_CMS, downVolLeft3, strikeDown3+spread1, optMat, K_CALL)-VanillaOption_N (forwardSpread - coeffShort * SENSI_CMS, downVolRight3, strikeDown3+spread2, optMat, K_CALL))/(spread2-spread1);
			sensi += (VanillaOption_N (forwardSpread - coeffShort * SENSI_CMS, upVolRight3,  strikeUp3+spread2,   optMat, K_CALL)-VanillaOption_N (forwardSpread - coeffShort * SENSI_CMS, upVolLeft3,  strikeUp3+spread1,   optMat, K_CALL))/(spread2-spread1);
			sensi-=price;
			sensi /= SENSI_CMS;
			itsFixedCpnShortSensi3[i] = sensi;

			// long CMS
			sensi  = (VanillaOption_N (forwardSpread + coeffLong * SENSI_CMS, downVolLeft3, strikeDown3+spread1, optMat, K_CALL)-VanillaOption_N (forwardSpread + coeffLong * SENSI_CMS, downVolRight3, strikeDown3+spread2, optMat, K_CALL))/(spread2-spread1);
			sensi += (VanillaOption_N (forwardSpread + coeffLong * SENSI_CMS, upVolRight3,  strikeUp3+spread2,   optMat, K_CALL)-VanillaOption_N (forwardSpread + coeffLong * SENSI_CMS, upVolLeft3,  strikeUp3+spread1,   optMat, K_CALL))/(spread2-spread1);
			sensi-=price;
			sensi /= SENSI_CMS;
			itsFixedCpnLongSensi3[i] = sensi;


			/// Cpn value
			itsFixedCpnValue3[i] = price;
					
			if (isVariable)
			{
				price  = (VanillaOption_N (forwardSpreadVarDown3, downVolLeftVar3, strikeDown3+spread1, optMat, K_CALL)-VanillaOption_N (forwardSpreadVarDown3, downVolRightVar3, strikeDown3+spread2, optMat, K_CALL))/(spread2-spread1);
				price += (VanillaOption_N (forwardSpreadVarUp3, upVolRightVar3, strikeUp3+spread2, optMat, K_CALL)-VanillaOption_N (forwardSpreadVarUp3, upVolLeftVar3, strikeUp3+spread1, optMat, K_CALL))/(spread2-spread1);
				
				/// short CMS
				sensi  = (VanillaOption_N (forwardSpreadVarDown3 - coeffShort * SENSI_CMS, downVolLeftVar3, strikeDown3+spread1, optMat, K_CALL)-VanillaOption_N (forwardSpreadVarDown3 - coeffShort * SENSI_CMS, downVolRightVar3, strikeDown3+spread2, optMat, K_CALL))/(spread2-spread1);
				sensi += (VanillaOption_N (forwardSpreadVarUp3 - coeffShort * SENSI_CMS, upVolRightVar3,  strikeUp3+spread2,   optMat, K_CALL)-VanillaOption_N (forwardSpreadVarUp3 - coeffShort * SENSI_CMS, upVolLeftVar3,  strikeUp3+spread1,   optMat, K_CALL))/(spread2-spread1);
				sensi-=price;
				sensi /= SENSI_CMS;
				itsVarCpnShortSensi3[i] = sensi;

				/// long CMS
				sensi  = (VanillaOption_N (forwardSpreadVarDown3 + coeffLong * SENSI_CMS, downVolLeftVar3, strikeDown3+spread1, optMat, K_CALL)-VanillaOption_N (forwardSpreadVarDown3 + coeffLong * SENSI_CMS, downVolRightVar3, strikeDown3+spread2, optMat, K_CALL))/(spread2-spread1);
				sensi += (VanillaOption_N (forwardSpreadVarUp3 + coeffLong * SENSI_CMS, upVolRightVar3,  strikeUp3+spread2,   optMat, K_CALL)-VanillaOption_N (forwardSpreadVarUp3 + coeffLong * SENSI_CMS, upVolLeftVar3,  strikeUp3+spread1,   optMat, K_CALL))/(spread2-spread1);
				sensi-=price;
				sensi /= SENSI_CMS;
				itsVarCpnLongSensi3[i] = sensi;

				itsVarCpnValue3[i] = price;
			}
			else
			{
				itsVarCpnValue3[i]      = 0.0;
				itsVarCpnLongSensi3[i]  = 0.0;
				itsVarCpnShortSensi3[i] = 0.0;
				itsPayIndexFwd3[i]      = 0.0;
			}
		}

/// FIN TRIPLE ------------------------------------------------------

	}
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: CreateVarNotionalSwaptionAtExer
///	Returns: void
///	Action : create the variable notional swaption
///			 
/////////////////////////////////////////////////////////////////
ARM_Swaption* ARM_CRASpreadCalculator::CreateVarNotionalSwaptionAtExer (int exerIdx)
{
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	ARM_ZeroCurve* curve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));

	/// get defaults swaption params for currency
	ARM_Currency* ccy		 = GetCurrencyUnit();
	int spotDays			 = ccy->GetSpotDays(); 
	ARM_INDEX_TYPE indexType = ccy->GetVanillaIndexType();
	char* resetCalendar		 = ccy->GetResetCalName(indexType);
	char* payCalendar		 = ccy->GetPayCalName(indexType);
	int stdFixFreq			 = ccy->GetFixedPayFreq();
	int stdFixDayCount		 = ccy->GetFixedDayCount();
	int resetGap			 = - spotDays;
	int fwdRule				 = K_MOD_FOLLOWING;	// for forward dates
    int intRule				 = K_ADJUSTED;			// for interest dates
    int stubRule			 = K_SHORTSTART;
	int resetTiming			 = K_ADVANCE;
	int payTiming			 = K_ARREARS;
		
	ARM_GP_Vector* exerDates = GetExerDateStrip()->GetResetDates();	
	double ExerciseDate = (*exerDates)[exerIdx];
	
	ARM_Date vnsExpiryDate(ExerciseDate);

	// compute start date
	ARM_Date startDate(ExerciseDate);
	startDate.GapBusinessDay(spotDays,resetCalendar);

	if(IsDbleCorridor() || itsCalibType==BASKET_CALIBRATION_SIMPLIFIED)
	{
		/// Forward yield curve & VNS will start at currend residual
		/// exotic leg start and reset spot days before
		startDate = ARM_Date((*(itsExerDateStrip->GetFlowStartDates()))[exerIdx]);
		vnsExpiryDate = startDate;
		vnsExpiryDate.GapBusinessDay(-spotDays,resetCalendar);
	}


	/// find (unadjusted) end date for variable notional swaption
	ARM_VanillaSpreadOptionArg* soVanillaArg = static_cast<ARM_VanillaSpreadOptionArg*>(&*itsVanillaArgCSO);
	ARM_GP_Vector* fixValues = soVanillaArg->GetFixValues();
	ARM_GP_Vector* payMults   = soVanillaArg->GetPayIndexLeverages();
	ARM_GP_Vector* cpnCoverages = soVanillaArg->GetPayPeriods();
	ARM_GP_Vector* cpnPaydates  = GetStructDateStrip()->GetPaymentDates();
	ARM_IntVector* periodIndexes = soVanillaArg->GetPeriodIndex();

	ARM_Date realEndDate((*cpnPaydates)[itsCpnSize-1]); /// last payment of the exotic leg
	bool IsStdCalibDbleCorridor = IsDbleCorridor() && itsCalibType==BASKET_CALIBRATION_SIMPLIFIED;
	if(!IsStdCalibDbleCorridor)
	{
		if (itsRefIndex1>itsRefIndex2)
			realEndDate = soVanillaArg->GetSwapLongFloatEndTime()->Elt(itsCpnSize-1) + asOfDate;
		else
			realEndDate = soVanillaArg->GetSwapShortFloatEndTime()->Elt(itsCpnSize-1) + asOfDate;
	}
	
	ARM_Date endDate(startDate), prevEndDate;

	while (endDate<realEndDate)
	{
		prevEndDate = endDate;
		endDate.AddMonths( 12/stdFixFreq ) ;
	}

	if( !IsStdCalibDbleCorridor && fabs(prevEndDate.GetJulian()-realEndDate.GetJulian()) < fabs(endDate.GetJulian()-realEndDate.GetJulian()) )
		endDate = prevEndDate;


	/// create date strip
	ARM_DateStrip dateStrip  (	startDate,endDate,stdFixFreq,stdFixDayCount,
								resetCalendar,
								fwdRule,intRule,stubRule,
								resetGap,stdFixFreq, GETDEFAULTVALUE, payCalendar, resetTiming, payTiming);	

	/// instanciate stripper
	ARM_STRIPPER stripper(asOfDate, ARM_DateStripPtr(new ARM_DateStrip(dateStrip)));

	ARM_GP_Vector* resetDates	 = dateStrip.GetResetDates();
	ARM_GP_Vector* startDates	 = dateStrip.GetFlowStartDates();
	ARM_GP_Vector* endDates		 = dateStrip.GetFlowEndDates();
	ARM_GP_Vector* payDates		 = dateStrip.GetPaymentDates();
	ARM_GP_Vector* interestTerms = dateStrip.GetInterestTerms();


	/// compute start df + fwd swap rates to feed stripper
	double startDf = curve->DiscountPrice(((*startDates)[0]-asOfDate)/K_YEAR_LEN);
	
	int i, j, sizeSched = startDates->size(); 

	double level (0.0);
	double df;
	ARM_GP_Vector forwards(sizeSched);

	stripper.setStartDf(startDf);

	for (j=0; j<sizeSched; j++)
	{
		df = curve->DiscountPrice(((*endDates)[j]-asOfDate)/K_YEAR_LEN);
		level += df * (*interestTerms)[j];
		forwards[j] = (startDf - df) / level;
		stripper.setSwapRate (j, forwards[j]);
	}
	
	// strip, given start df and forward swap rates
	stripper.strip();

	
	///
	/// Fixed Leg notios
	///
	/// create notional
	ARM_Vector fixAbs(sizeSched);
	ARM_Vector fixNotios(sizeSched);

	/// check that notionals are all > 0 (as snv requires so)
	ARM_Vector* notio = itsNotional.GetDiscreteValues();
	for (j=0; j<notio->size(); j++)
	{
		if (notio->Elt(j)<0)
			ARM_THROW( ERR_INVALID_ARGUMENT, "CreateVarNotionalSwaptionAtExer : negative notional" );
	}


	// compute fixIndex so that we are as close as possible from product end date
	size_t fixIndex (0);

	if ((*endDates)[sizeSched-1] < GetEndDate().GetJulian())
		ARM_THROW( ERR_INVALID_ARGUMENT, "CreateVarNotionalSwaptionAtExer : unresolved date problem" );

	while ((*endDates)[fixIndex]<GetEndDate().GetJulian())
		fixIndex ++;


	if (	(fixIndex>0)
		&&  ( (*endDates)[fixIndex] - GetEndDate().GetJulian() >= GetEndDate().GetJulian()  - (*endDates)[fixIndex-1] )  )
		fixIndex = fixIndex - 1;
	
	for (j=0; j<sizeSched; j++)
	{
		fixAbs[j] = (*endDates)[j];
		
		if (j>fixIndex)
			fixNotios[j] = 0.0;
		else
		{	/// rather approximative ...
			fixNotios[j] = itsNotional.Interpolate((*payDates)[j]); 
			
			/// take 0 if <0 (with lin interp, result can be sometimes <0 for small notio values)
			if (fixNotios[j]<0)
				fixNotios[j] = 0.0;
		}
	}

		
	/// if fixNotios are all 0, exclude this swaption of calib set
	bool AllZero = true;
	for (j=0; j<sizeSched; j++)
	{
		if (fixNotios[j])
		{
			AllZero = false;
			break;
		}
	}
	
	if (AllZero) 
		return NULL;
	
	// precomputation
	ARM_GP_Vector A(itsCpnSize), Bshort(itsCpnSize), Blong (itsCpnSize), Bpay (itsCpnSize);
		
	double soLegPv (0.0);

	double weight = (itsCalibType==BASKET_CALIBRATION_SIMPLIFIED?0.:1.);

	if (IsTripleRange())
	{
		for (i = itsvCpnIndex[exerIdx]; i<itsCpnSize; i++)
		{
			double payMult = (*payMults)[(*periodIndexes)[i]];
			df			= curve->DiscountPrice(((*cpnPaydates)[i]-asOfDate)/K_YEAR_LEN);
			A[i]		= (*cpnCoverages)[i] * itsvCpnNominal[i] * (
					itsFixedCpnValue[i] * (*fixValues)[i]
				+	itsFixedCpnValue2[i] * itsFixedValues2[i]
				+	itsFixedCpnValue3[i] * itsFixedValues3[i]	)/ARM_Constants::rateBase;
			A[i]	   += (*cpnCoverages)[i] * itsvCpnNominal[i] * (
					itsVarCpnValue[i]   
				+	itsVarCpnValue2[i]   
				+	itsVarCpnValue3[i]	) * payMult * itsPayIndexFwd[i];
			Blong[i]	= (*cpnCoverages)[i] * itsvCpnNominal[i] * df * (
					itsFixedCpnLongSensi[i]  * (*fixValues)[i] 
				+	itsFixedCpnLongSensi2[i] * itsFixedValues2[i]
				+	itsFixedCpnLongSensi3[i] * itsFixedValues3[i]	)/ARM_Constants::rateBase * weight;
			Bshort[i]	= (*cpnCoverages)[i] * itsvCpnNominal[i] * df * (
					itsFixedCpnShortSensi[i] * (*fixValues)[i]
				+	itsFixedCpnShortSensi2[i] * itsFixedValues2[i]
				+	itsFixedCpnShortSensi3[i] * itsFixedValues3[i]	)/ARM_Constants::rateBase * weight;
			Blong[i]   += (*cpnCoverages)[i] * itsvCpnNominal[i] * df * (
					itsVarCpnLongSensi[i]  
				+	itsVarCpnLongSensi2[i]
				+	itsVarCpnLongSensi3[i]	) * payMult * itsPayIndexFwd[i] * weight;
			Bshort[i]  += (*cpnCoverages)[i] * itsvCpnNominal[i] * df * (
					itsVarCpnShortSensi[i]
				+	itsVarCpnShortSensi2[i]
				+	itsVarCpnShortSensi3[i]	) * payMult * itsPayIndexFwd[i]	* weight;
			Bpay[i]		= (*cpnCoverages)[i] * itsvCpnNominal[i] * df * (
					itsVarCpnValue[i]
				+	itsVarCpnValue2[i]
				+	itsVarCpnValue3[i]	) * payMult * weight;
			soLegPv	   += A[i] * df;
		}
	}
	else
	{
		for (i = itsvCpnIndex[exerIdx]; i<itsCpnSize; i++)
		{
			double payMult = (*payMults)[(*periodIndexes)[i]];
			df			= curve->DiscountPrice(((*cpnPaydates)[i]-asOfDate)/K_YEAR_LEN);
			A[i]		= (*cpnCoverages)[i] * itsvCpnNominal[i] * itsFixedCpnValue[i] * (*fixValues)[i]/ARM_Constants::rateBase;
			A[i]	   += (*cpnCoverages)[i] * itsvCpnNominal[i] * itsVarCpnValue[i]   * payMult * itsPayIndexFwd[i];
			Blong[i]	= (*cpnCoverages)[i] * itsvCpnNominal[i] * df * itsFixedCpnLongSensi[i]  * (*fixValues)[i]/ARM_Constants::rateBase		* weight;
			Bshort[i]	= (*cpnCoverages)[i] * itsvCpnNominal[i] * df * itsFixedCpnShortSensi[i] * (*fixValues)[i]/ARM_Constants::rateBase		* weight;
			Blong[i]   += (*cpnCoverages)[i] * itsvCpnNominal[i] * df * itsVarCpnLongSensi[i]  * payMult * itsPayIndexFwd[i]					* weight;
			Bshort[i]  += (*cpnCoverages)[i] * itsvCpnNominal[i] * df * itsVarCpnShortSensi[i] * payMult * itsPayIndexFwd[i]					* weight;
			Bpay[i]		= (*cpnCoverages)[i] * itsvCpnNominal[i] * df * itsVarCpnValue[i] * payMult												* weight;
			soLegPv	   += A[i] * df;
		}
	}

	if ( GetPayRec()==K_RCV )
		soLegPv -= itsvExerFees[exerIdx] * curve->DiscountPrice((ExerciseDate-asOfDate)/K_YEAR_LEN);
	else
		soLegPv += itsvExerFees[exerIdx] * curve->DiscountPrice((ExerciseDate-asOfDate)/K_YEAR_LEN);

	ARM_GP_Vector df0(itsCpnSize), Slong0(itsCpnSize), Sshort0(itsCpnSize), Spay0(itsCpnSize);
	double dfFee0;

	dfFee0 =  stripper.df(ExerciseDate-asOfDate);

	/// Notional used for normalization
	double fixNotional (1.0);
	for (j=0; j<sizeSched; j++)
	{
		if (fixNotios[j])
		{
			fixNotional = fixNotios[j];
			break;
		}
	}

	double Numeraire0 (0.0);
	ARM_GP_Vector NumerizedLiborFlow0 (sizeSched, 0.0);

	for (j=0; j<sizeSched; j++)
		Numeraire0 += (*interestTerms)[j] * (fixNotios[j]/fixNotional) * stripper.df((*endDates)[j]-asOfDate) ;
	
	for (j=0; j<sizeSched; j++)
		NumerizedLiborFlow0[j] =  (stripper.df((*startDates)[j]-asOfDate) -  stripper.df((*endDates)[j]-asOfDate)) / Numeraire0;


	/// initial values for dfs and swaprates
	for (i = itsvCpnIndex[exerIdx]; i<itsCpnSize; i++)
	{
		df0[i] = stripper.df((*cpnPaydates)[i]-asOfDate);
				
		double longStartTime = soVanillaArg->GetSwapLongFloatStartTime()->Elt(i);
		double longEndTime   = soVanillaArg->GetSwapLongFloatEndTime()->Elt(i);
		ARM_GP_Vector* longPayTimes    = soVanillaArg->GetSwapLongFixPayTimes()[i];
		ARM_GP_Vector* longPayPeriods  = soVanillaArg->GetSwapLongFixPayPeriods()[i];
		Slong0[i] = stripper.swapRate(longStartTime, longEndTime, *longPayTimes, *longPayPeriods);

		double shortStartTime = soVanillaArg->GetSwapShortFloatStartTime()->Elt(i);
		double shortEndTime   = soVanillaArg->GetSwapShortFloatEndTime()->Elt(i);
		ARM_GP_Vector* shortPayTimes    = soVanillaArg->GetSwapShortFixPayTimes()[i];
		ARM_GP_Vector* shortPayPeriods  = soVanillaArg->GetSwapShortFixPayPeriods()[i];
		Sshort0[i] = stripper.swapRate(shortStartTime, shortEndTime, *shortPayTimes, *shortPayPeriods);

		Spay0[i] = 0.0;

		if (IsVariableCpn(i))
		{
			/// only if new period
			if (i==0 || ( (i>0) && ((*periodIndexes)[i]!=(*periodIndexes)[i-1]) ) )
			{
				int periodIndex = (*periodIndexes)[i];
				double payStartTime = soVanillaArg->GetSwapPayFloatStartTime()->Elt(periodIndex);
				double payEndTime   = soVanillaArg->GetSwapPayFloatEndTime()->Elt(periodIndex);
				ARM_GP_Vector* payTimes    = soVanillaArg->GetSwapPayFixPayTimes()[periodIndex];
				ARM_GP_Vector* payPeriods  = soVanillaArg->GetSwapPayFixPayPeriods()[periodIndex];
				Spay0[i] = stripper.swapRate(payStartTime, payEndTime, *payTimes, *payPeriods);
			}
			else
				Spay0[i] = 0.0;
		}
	}

	/// initial value for funding leg
	ARM_GP_Vector fundingStarts (0) ;
	ARM_GP_Vector fundingEnds	(0) ;
	ARM_GP_Vector fundingPeriods(0) ;
	ARM_GP_Vector fundingNotios (0) ;
	ARM_GP_Vector fundingMargin (0) ;

	for (i=itsvFundIndex[exerIdx]; i<itsFundSize; i++)
	{
		fundingStarts.push_back ( GetFundDateStrip()->GetFlowStartDates()->Elt(i) - asOfDate );
		fundingEnds.push_back   ( GetFundDateStrip()->GetFlowEndDates()->Elt(i) - asOfDate );
		fundingPeriods.push_back( GetFundDateStrip()->GetInterestTerms()->Elt(i));
		fundingNotios.push_back ( itsvFundNominal[i] );
		fundingMargin.push_back ( itsvFundSpread[i] );
	}
	
	double floatLeg0 = stripper.stdFundingLeg(fundingStarts, fundingEnds, fundingPeriods, fundingNotios, fundingMargin);

	/// donc forget to store underlying PV for later strike determination
	double underlyingPv = floatLeg0 - soLegPv;
	
	/// compute deal sensitivities w.r.t fwd swap rates
	ARM_GP_Vector sensitivities (sizeSched, 0.0);
	double Sshort, Slong, floatLeg, Spay, SpaySensi(0.0), Numeraire;
	ARM_GP_Vector invNumeraireSensis (sizeSched, 0.0);
	ARM_GP_Matrix toInvert(sizeSched, sizeSched);
		
	/// main loop on bumped forwards
	for (j=0; j<sizeSched; j++)
	{
		/// shift fwd swap rate #j
		stripper.setSwapRate(j, forwards[j] + SENSI_CMS);
		stripper.strip();

		/// (1) for spread option leg
		for (i = itsvCpnIndex[exerIdx]; i<itsCpnSize; i++)
		{
			df = stripper.df((*cpnPaydates)[i]-asOfDate);
			
			sensitivities[j] -= A[i] * (df - df0[i]) / SENSI_CMS;

			double longStartTime = soVanillaArg->GetSwapLongFloatStartTime()->Elt(i);
			double longEndTime   = soVanillaArg->GetSwapLongFloatEndTime()->Elt(i);
			ARM_GP_Vector* longPayTimes    = soVanillaArg->GetSwapLongFixPayTimes()[i];
			ARM_GP_Vector* longPayPeriods  = soVanillaArg->GetSwapLongFixPayPeriods()[i];
			Slong = stripper.swapRate(longStartTime, longEndTime, *longPayTimes, *longPayPeriods);

			double shortStartTime = soVanillaArg->GetSwapShortFloatStartTime()->Elt(i);
			double shortEndTime   = soVanillaArg->GetSwapShortFloatEndTime()->Elt(i);
			ARM_GP_Vector* shortPayTimes    = soVanillaArg->GetSwapShortFixPayTimes()[i];
			ARM_GP_Vector* shortPayPeriods  = soVanillaArg->GetSwapShortFixPayPeriods()[i];
			Sshort = stripper.swapRate(shortStartTime, shortEndTime, *shortPayTimes, *shortPayPeriods);

			sensitivities[j] -= Blong[i]  * (Slong  - Slong0[i])  / SENSI_CMS;
			sensitivities[j] -= Bshort[i] * (Sshort - Sshort0[i]) / SENSI_CMS;

			if (IsVariableCpn(i))
			{
				/// compute SpaySensi only if new period
				if (i==0 || ( (i>0) && ((*periodIndexes)[i]!=(*periodIndexes)[i-1]) ) )
				{
					int periodIndex = (*periodIndexes)[i];
					double payStartTime = soVanillaArg->GetSwapPayFloatStartTime()->Elt(periodIndex);
					double payEndTime   = soVanillaArg->GetSwapPayFloatEndTime()->Elt(periodIndex);
					ARM_GP_Vector* payTimes    = soVanillaArg->GetSwapPayFixPayTimes()[periodIndex];
					ARM_GP_Vector* payPeriods  = soVanillaArg->GetSwapPayFixPayPeriods()[periodIndex];
					Spay = stripper.swapRate(payStartTime, payEndTime, *payTimes, *payPeriods);
					SpaySensi = (Spay - Spay0[i]) / SENSI_CMS;

				}

				/// SpaySensi remains the same for along the payment period
				sensitivities[j] -= Bpay[i] * SpaySensi;
			}
		}

		/// take fee into account
		df = stripper.df(ExerciseDate-asOfDate);
		if ( GetPayRec()==K_RCV )
			sensitivities[j] += itsvExerFees[exerIdx] * (df - dfFee0);
		else
			sensitivities[j] -= itsvExerFees[exerIdx] * (df - dfFee0);
		


		/// (2) for funding leg + margin
		floatLeg = stripper.stdFundingLeg(fundingStarts, fundingEnds, fundingPeriods, fundingNotios, fundingMargin);

		sensitivities[j] += (floatLeg - floatLeg0) / SENSI_CMS;

		/// (3) for 1/Numeraire
		Numeraire = 0.0;
		for (size_t k=0; k<sizeSched; k++)
		{	
			df = stripper.df((*endDates)[k]-asOfDate);
			Numeraire += (*interestTerms)[k] * (fixNotios[k]/fixNotional) * df;
		}

		for (k=0; k<sizeSched; k++)
		{
			toInvert(j,k)  =  (stripper.df((*startDates)[k]-asOfDate) -  stripper.df((*endDates)[k]-asOfDate)) / Numeraire;
			toInvert(j,k) -= NumerizedLiborFlow0[k];
			toInvert(j,k) /=  SENSI_CMS;
		}
		
		invNumeraireSensis[j] = (1./Numeraire - 1./Numeraire0) / SENSI_CMS;


		/// RAZ
		stripper.setSwapRate(j, forwards[j]);
	}

	/// Reset stripper
	stripper.strip();

	///
	/// build matrix to invert
	/// for the moment its contains the sensis of the Underlying w.r.t swap rates
	/// now, we want it to contain the sensi of Underlying/Numeraire w.r.t. swap rates
	for (j=0; j<sizeSched; j++)
	{
		sensitivities[j] *= 1./Numeraire0;
		sensitivities[j] += underlyingPv * invNumeraireSensis[j];
	}
	
	/// invert matrix
	LinSolve(&toInvert, &sensitivities);
	ARM_Vector floatNotios (sizeSched);
	for (j=0; j<sizeSched; j++)
		floatNotios[j] = sensitivities[j];
	
		
	/// compute strike
	double strike;
	double varFloatLeg = 0.0;

	for (j=0; j<sizeSched; j++)
	{			
		varFloatLeg += floatNotios[j] * (stripper.df((*startDates)[j]-asOfDate) - stripper.df((*endDates)[j]-asOfDate));
	}

	if (itsCalibStrikeType[0] == ATM)
		underlyingPv = 0.0;

	strike = (varFloatLeg - underlyingPv) / (Numeraire0 * fixNotional) ;
	

	///  vector need to be cloned (deleted in ARM_ReferenceValue destructor)
	ARM_ReferenceValue floatNotionalRefVal ((ARM_Vector*)fixAbs.Clone(), (ARM_Vector*)floatNotios.Clone());
	floatNotionalRefVal.SetCalcMethod(K_STEPUP_RIGHT);

		///  vectors need to be cloned (deleted in ARM_ReferenceValue destructor)
	ARM_ReferenceValue fixNotionalRefVal((ARM_Vector*)fixAbs.Clone(), (ARM_Vector*)fixNotios.Clone());
	fixNotionalRefVal.SetCalcMethod(K_STEPUP_RIGHT);


		
	// index type: same frequency as fixed leg
	
	ARM_IRIndex irIndex (indexType, stdFixFreq, stdFixFreq, ccy);
	irIndex.SetTerm(stdFixFreq);
	irIndex.SetYearTerm(1.0/stdFixFreq);

	ARM_SwapLeg armFloatLeg( startDate, 
							 endDate, 
							 &irIndex, 
							 K_PAY,
							 0.0, 
							 K_SHORTSTART, 
							 K_COMP_PROP,
							 ccy,
							 ccy->GetLiborIndexDayCount());

	/// floatNotionalRefVal is cloned in SetAmount
	armFloatLeg.SetAmount(&floatNotionalRefVal);


	// create fixed leg with good strike
	ARM_FixLeg armFixLeg ( startDate,
						   endDate, 
						   strike * 100.0,
						   K_RCV, 
						   stdFixFreq,
						   stdFixDayCount,
						   K_COMP_PROP,
						   K_ARREARS, 
						   K_ADJUSTED,
						   K_SHORTSTART,
						   ccy);
	
	/// fixNotionalRefVal is cloned in SetAmount
	armFixLeg.SetAmount(&fixNotionalRefVal);
	
	
	/// create swap
	ARM_Swap swap(&armFixLeg, &armFloatLeg);
			
	/// create swaptionj
	ARM_Swaption* swaption = new ARM_Swaption((&swap), K_RCV,K_EUROPEAN, strike * 100.0, vnsExpiryDate);

	if (resetCalendar)
		delete resetCalendar;

	if (payCalendar)
		delete payCalendar;

	return swaption;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: ComputeResidualUnderlying
///	Returns: void
///	Action : compute PV of the resiudal underlying
///			 at current the exercise date
/////////////////////////////////////////////////////////////////
double ARM_CRASpreadCalculator::ComputeResidualUnderlying(int exerIdx, double& floatLegPv, double& exoLegPv)
{
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	ARM_ZeroCurve* curve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));
		
	double ExerciseDate = (*(GetExerDateStrip()->GetResetDates()))[exerIdx];
	
	/// Restore coupon datas
	ARM_VanillaSpreadOptionArg* soVanillaArg = static_cast<ARM_VanillaSpreadOptionArg*>(&*itsVanillaArgCSO);
	ARM_GP_Vector* fixValues = soVanillaArg->GetFixValues();
	ARM_GP_Vector* payMults   = soVanillaArg->GetPayIndexLeverages();
	ARM_GP_Vector* cpnCoverages = soVanillaArg->GetPayPeriods();
	ARM_GP_Vector* cpnPaydates  = GetStructDateStrip()->GetPaymentDates();
	ARM_IntVector* periodIndexes = soVanillaArg->GetPeriodIndex();

	/// Compute exotic coupon leg
	size_t i;
	double df,flow;
	exoLegPv = 0.0;
	if (IsTripleRange())
	{
		for (i = itsvCpnIndex[exerIdx]; i<itsCpnSize; i++)
		{
			double payMult = (*payMults)[(*periodIndexes)[i]];
			df			= curve->DiscountPrice(((*cpnPaydates)[i]-asOfDate)/K_YEAR_LEN);
			flow		= (*cpnCoverages)[i] * itsvCpnNominal[i] * (
					itsFixedCpnValue[i] * (*fixValues)[i]
				+	itsFixedCpnValue2[i] * itsFixedValues2[i]
				+	itsFixedCpnValue3[i] * itsFixedValues3[i]	)/ARM_Constants::rateBase;
			flow	   += (*cpnCoverages)[i] * itsvCpnNominal[i] * (
					itsVarCpnValue[i]   
				+	itsVarCpnValue2[i]   
				+	itsVarCpnValue3[i]	) * payMult * itsPayIndexFwd[i];
			exoLegPv	   += flow * df;
		}
	}
	else
	{
		for (i = itsvCpnIndex[exerIdx]; i<itsCpnSize; i++)
		{
			double payMult = (*payMults)[(*periodIndexes)[i]];
			df			= curve->DiscountPrice(((*cpnPaydates)[i]-asOfDate)/K_YEAR_LEN);
			flow		= (*cpnCoverages)[i] * itsvCpnNominal[i] * itsFixedCpnValue[i] * (*fixValues)[i]/ARM_Constants::rateBase;
			flow	   += (*cpnCoverages)[i] * itsvCpnNominal[i] * itsVarCpnValue[i]   * payMult * itsPayIndexFwd[i];
			exoLegPv	   += flow * df;
		}
	}

	if ( GetPayRec()==K_RCV )
		exoLegPv -= itsvExerFees[exerIdx] * curve->DiscountPrice((ExerciseDate-asOfDate)/K_YEAR_LEN);
	else
		exoLegPv += itsvExerFees[exerIdx] * curve->DiscountPrice((ExerciseDate-asOfDate)/K_YEAR_LEN);


	/// Compute floating and fixed part of the funding leg
	double t,it,dfStart,dfEnd;
	floatLegPv = 0.0;
	for (i=itsvFundIndex[exerIdx]; i<itsFundSize; i++)
	{
		t = (*(GetFundDateStrip()->GetFlowStartDates()))[i];
		dfStart = curve->DiscountPrice((t-asOfDate)/K_YEAR_LEN);
		t = (*(GetFundDateStrip()->GetFlowEndDates()))[i];
		dfEnd = curve->DiscountPrice((t-asOfDate)/K_YEAR_LEN);

		it = (*(GetFundDateStrip()->GetInterestTerms()))[i];

		exoLegPv -= it * itsvFundSpread[i]*dfEnd * itsvFundNominal[i];
		floatLegPv += (dfStart-dfEnd) * itsvFundNominal[i];
	}
	
	return GetPayRec() * (exoLegPv - floatLegPv);
}


////////////////////////////////////////////////////
///	Class   : ARM_CRASpreadCalculator
///	Routines: toString
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
string ARM_CRASpreadCalculator::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream craData;

    return craData.str();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: View
///	Returns: 
///	Action : .
/////////////////////////////////////////////////////////////////
void ARM_CRASpreadCalculator::View(char* id, FILE* ficOut) const
{
	FILE* fOut;
    char fOutName[200];
	
    if ( ficOut == NULL )
    {
		ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
		(void) unlink(fOutName);
		fOut = fopen(fOutName, "w");
    }
    else
		fOut = ficOut;

    /// Callable Range Accrual Calculator specific datas viewing
    fprintf(fOut,"\n\n ==================================================================== \n");
	fprintf(fOut,"\n\n =======>      CALLABLE RANGE ACCRUAL CALCULATOR SPREAD	 <========= \n");
    fprintf(fOut,"\n\n ==================================================================== \n");

	fprintf(fOut,"%s", DealDesDataToString().c_str());

	fprintf(fOut,"\n\n ============>    Calib Flags    <=========== \n");
	fprintf(fOut,(itsSwoptCalib			?		"Volatility bootstrapping     : ON\n" : "Volatility bootstrapping     : OFF\n"));
	fprintf(fOut,(itsVolUnSqueezer		?	"Local variance unsqueezer : ON\n" : "Local variance unsqueezer : OFF\n"));
	if(IsDbleCorridor())
		fprintf(fOut,(itsCorrelUnSqueezer	?	"Local correlation unsqueezer : ON\n" : "Local correlation unsqueezer : OFF\n"));
	fprintf(fOut,(itsPVAdjuster			?	"PV Adjuster : ON\n" : "PV Adjuster : OFF\n"));

	fprintf(fOut,"\n\n ============>    Reset Optimiser Params    <=========== \n");
	if(itsOptimResetData[OPTIM_RESET])
	{
		fprintf(fOut,"Reset Optimizer : ON\n");
		fprintf(fOut,"	Daily Limit	  : %8.1lf days\n",itsOptimResetData[OPTIM_RESET_DAILYLIMIT]);
		fprintf(fOut,"	Weekly Limit  : %8.1lf days\n",itsOptimResetData[OPTIM_RESET_WEEKLYLIMIT]);
		fprintf(fOut,"	Monthly Limit : %8.1lf days\n",itsOptimResetData[OPTIM_RESET_MONTHLYLIMIT]);
	}
	else
		fprintf(fOut,"Reset Optimizer : OFF\n");

	CC_Ostringstream genData;
	fprintf(fOut,"\n\n ============>    GEN CALCULATOR    <=========== \n");
	genData << "\n\nCommon Calculator Part\n" << ARM_GenCalculator::toString() << "\n\n";
	fprintf(fOut,"%s", genData.str().c_str());
	fprintf(fOut,"\n\n ============>  EOF GEN CALCULATOR  <=========== \n");

	fprintf(fOut,"\n\n ============> SWAPTION PORTFOLIO    <=========== \n");
	if (!GetSwoptPF().IsNull())
		fprintf(fOut,"%s", (&*GetSwoptPF())->toString().c_str());
	fprintf(fOut,"\n\n ============>  SWAPTION PORTFOLIO  <=========== \n");

	fprintf(fOut,"\n\n ============> CORRIDOR SPREAD OPTION PORTFOLIO    <=========== \n");
	if (!itsCSOPF.IsNull())
		fprintf(fOut,"%s", (&*itsCSOPF)->toString().c_str());
	fprintf(fOut,"\n\n ============>  CORRIDOR SPREAD OPTION PORTFOLIO  <=========== \n");

	fprintf(fOut,"\n\n ============> CORRIDOR SPREAD OPTION PORTFOLIO 2   <=========== \n");
	if (!itsCSOPF2.IsNull())
		fprintf(fOut,"%s", (&*itsCSOPF2)->toString().c_str());
	fprintf(fOut,"\n\n ============>  CORRIDOR SPREAD OPTION PORTFOLIO 2 <=========== \n");

	fprintf(fOut,"\n\n ============> CORRIDOR SPREAD OPTION PORTFOLIO 3   <=========== \n");
	if (!itsCSOPF3.IsNull())
		fprintf(fOut,"%s", (&*itsCSOPF3)->toString().c_str());
	fprintf(fOut,"\n\n ============>  CORRIDOR SPREAD OPTION PORTFOLIO 3 <=========== \n");

	fprintf(fOut,"\n\n ============>  STRUCT DATE STRIP <=========== \n");
	if(itsStructDateStrip != ARM_DateStripPtr(NULL))
		itsStructDateStrip->View(id,fOut);

	if ( ficOut == NULL )
       fclose(fOut);
}



// In the case of CRA from CCSO 
double ARM_CRASpreadCalculator::InterpolMeanRevParam(ARM_OptionPortfolio* optionPortfolio, 
                                                     ARM_VolLInterpol* meanRevParam, 
                                                     ARM_Date asOfDate)
{
    double	mrs = MRS_DEFAULT_VALUE;

    if ( meanRevParam == NULL )
    {
       return(mrs);
    }

    ARM_Date startDate = optionPortfolio->GetCorridorLeg()->GetStartDate();

    ARM_Date endDate   = optionPortfolio->GetCorridorLeg()->GetEndDate();

    double   maturity = 0.;


    if ( startDate > asOfDate )
       maturity = (endDate-startDate)/365.0;
    else
       maturity = (endDate-asOfDate)/365.0;

 

    ARM_Date nextNoticeDate = optionPortfolio->GetNextNoticeDate(asOfDate);

    double callability      = (nextNoticeDate-asOfDate)/365.;


    mrs = meanRevParam->ComputeVol(maturity, callability);

    return(mrs);
}



// In the case of Bermuda from CCSO 
double ARM_CRASpreadCalculator::InterpolMeanRevParam(ARM_Swaption* bermuda, 
                                                     ARM_VolLInterpol* meanRevParam, 
                                                     ARM_Date asOfDate)
{
    double	mrs = MRS_DEFAULT_VALUE;

    if ( meanRevParam == NULL )
    {
       return(mrs);
    }

    ARM_Date startDate = bermuda->GetFloatLeg()->GetStartDate();

    ARM_Date endDate   = bermuda->GetFloatLeg()->GetEndDate();

    double   maturity = 0.;


    if ( startDate > asOfDate )
       maturity = (endDate-startDate)/365.0;
    else
       maturity = (endDate-asOfDate)/365.0;

 

    ARM_Date nextNoticeDate = bermuda->GetNextNoticeDate(asOfDate);

    double callability      = (nextNoticeDate-asOfDate)/365.;


    mrs = meanRevParam->ComputeVol(maturity, callability);

    return(mrs);
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRASpreadCalculator
///	Routine: ConvertToVarNotinalSwaption
///	Returns: ARM_VanillaSwaptionArg
///	Action : convert a general swaption (variable funding spreads,
///			 strikes and notionals) to a fixed strike & unspreaded
///			 funding leg variable notional standard swaption
/////////////////////////////////////////////////////////////////
ARM_Swaption* ARM_CRASpreadCalculator::ConvertToVarNotionalSwaption(ARM_ZeroCurve* curve, const ARM_VanillaSwaptionArg& oswVanillaArg)
{
	double asOfDate = curve->GetAsOfDate().GetJulian();
	double expiryDate = oswVanillaArg.GetExpiry() + asOfDate;	

	/// Get defaults swaption params for currency
	ARM_Currency* ccy		 = curve->GetCurrencyUnit();
	int spotDays			 = ccy->GetSpotDays(); 
	ARM_INDEX_TYPE indexType = ccy->GetVanillaIndexType();
	char* resetCalendar		 = ccy->GetResetCalName(indexType);
	char* payCalendar		 = ccy->GetPayCalName(indexType);
	int stdFixFreq			 = ccy->GetFixedPayFreq();
	int stdFixDayCount		 = ccy->GetFixedDayCount();
	int resetGap			 = - spotDays;
	int fwdRule				 = K_MOD_FOLLOWING;	// for forward dates
    int intRule				 = K_ADJUSTED;		// for interest dates
    int stubRule			 = K_SHORTSTART;
	int resetTiming			 = K_ADVANCE;
	int payTiming			 = K_ARREARS;
		
	/// Standard start date
	ARM_Date startDate(expiryDate);
	startDate.GapBusinessDay(spotDays,resetCalendar);
		
	/// Find the standard (unadjusted) end date for variable notional
	/// swaption which ends just after the very last swaption flow
	size_t nbFloatFlows = oswVanillaArg.GetFloatEndTimes()->size();
	size_t nbFixFlows = oswVanillaArg.GetFixPayTimes()->size();
	double floatEndTime = (*(oswVanillaArg.GetFloatEndTimes()))[nbFloatFlows-1];
	double fixEndTime = (*(oswVanillaArg.GetFixPayTimes()))[nbFixFlows-1];
	double oswEndDate = floatEndTime;
	if(oswEndDate < fixEndTime)
		oswEndDate = fixEndTime;
	oswEndDate += asOfDate;
	
	ARM_Date endDate(startDate);
	while(endDate.GetJulian() < oswEndDate)
		endDate.AddMonths( 12/stdFixFreq ) ;

	/// Create a standard date strip (w.r.t. yield curve currency)
	ARM_DateStrip dateStrip  (	startDate,endDate,stdFixFreq,stdFixDayCount,
								resetCalendar,
								fwdRule,intRule,stubRule,
								resetGap,stdFixFreq, GETDEFAULTVALUE, payCalendar, resetTiming, payTiming);	

	/// Create the stripper object
	ARM_STRIPPER stripper(asOfDate, ARM_DateStripPtr(new ARM_DateStrip(dateStrip)));

	ARM_GP_Vector* resetDates	 = dateStrip.GetResetDates();
	ARM_GP_Vector* startDates	 = dateStrip.GetFlowStartDates();
	ARM_GP_Vector* endDates		 = dateStrip.GetFlowEndDates();
	ARM_GP_Vector* payDates		 = dateStrip.GetPaymentDates();
	ARM_GP_Vector* interestTerms = dateStrip.GetInterestTerms();


	/// Compute a yield curve starting at start date to feed the stripper
	double startDf = curve->DiscountPrice(((*startDates)[0]-asOfDate)/K_YEAR_LEN);
	
	int i, j, sizeSched = startDates->size(); 

	double level = 0.0;
	double df;
	ARM_GP_Vector forwards(sizeSched);

	stripper.setStartDf(startDf);

	for (j=0; j<sizeSched; j++)
	{
		df = curve->DiscountPrice(((*endDates)[j]-asOfDate)/K_YEAR_LEN);
		level += df * (*interestTerms)[j];
		forwards[j] = (startDf - df) / level;
		stripper.setSwapRate (j, forwards[j]);
	}
	
	/// Strip in zero coupon bonds
	stripper.strip();

	
	///
	/// Fixed Leg notios
	///
	/// Create the standard notional profile interpolator is not cloned)
	ARM_StepUpRightOpenCstExtrapolDble* interpolator = new ARM_StepUpRightOpenCstExtrapolDble(ARM_StepUpRightOpenCstExtrapolDble());
	ARM_Curve oswFixNotionals(*(oswVanillaArg.GetFixPayTimes()),*(oswVanillaArg.GetFixNotional()),interpolator);

	ARM_Vector stdSwapDates(sizeSched);
	ARM_Vector stdFixNotionals(sizeSched);


	/// Compute the index that refer to the closest standard end date to the very last swaption flow
	/// Standard fixed leg notional profile is non null up to this index
	size_t stdFixIndex=0;

	if((*endDates)[sizeSched-1] < oswEndDate)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ConvertToVarNotionalSwaption : sdt dates must cover swaption end date" );

	while((*endDates)[stdFixIndex] < oswEndDate)
		++stdFixIndex;


	if(stdFixIndex > 0 && (*endDates)[stdFixIndex] - oswEndDate >= oswEndDate  - (*endDates)[stdFixIndex-1] )
		--stdFixIndex;
	
	bool isZero = true;
	bool isNegative = false;
	double refStdFixNotional;
	for (j=0; j<sizeSched; j++)
	{
		stdSwapDates[j] = (*endDates)[j];
		
		if (j>stdFixIndex)
			stdFixNotionals[j] = 0.0;
		else
			stdFixNotionals[j] = oswFixNotionals.Interpolate((*payDates)[j]); 
		if(stdFixNotionals[j] > 0.0)
		{
			isZero = false;
			refStdFixNotional=stdFixNotionals[j];
		}
		else if(stdFixNotionals[j] < 0.0)
			isNegative = true;
	}
	
	if(isZero || isNegative) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "ConvertToVarNotionalSwaption : all fixed leg notionals are all null or one is negative" );
	
	/// Fixed leg & standard numeraire reference PVs
	const ARM_GP_Vector* fixPayRates	= oswVanillaArg.GetStrikes();
	const ARM_GP_Vector* fixPayPeriods	= oswVanillaArg.GetFixPayPeriods();
	const ARM_GP_Vector* fixPayTimes	= oswVanillaArg.GetFixPayTimes();
	const ARM_GP_Vector* fixNotionals	= oswVanillaArg.GetFixNotional();
	
	double fixLeg0=0.0;

	ARM_GP_Vector A(nbFixFlows);
	for(i = 0; i<nbFixFlows; ++i)
	{
		df			= curve->DiscountPrice((*fixPayTimes)[i]/K_YEAR_LEN);
		A[i]		= (*fixPayPeriods)[i] * (*fixNotionals)[i] * (*fixPayRates)[i];
		fixLeg0		+= A[i] * df;
	}


	ARM_GP_Vector df0(nbFixFlows);

	double numeraire0 = 0.0;
	for (j=0; j<sizeSched; j++)
		numeraire0 += (*interestTerms)[j] * (stdFixNotionals[j]/refStdFixNotional) * stripper.df((*endDates)[j]-asOfDate) ;
	
	ARM_GP_Vector numerizedLiborFlow0(sizeSched, 0.0);
	for (j=0; j<sizeSched; j++)
		numerizedLiborFlow0[j] =  (stripper.df((*startDates)[j]-asOfDate) -  stripper.df((*endDates)[j]-asOfDate)) / numeraire0;


	/// Dfs, float leg & swap reference PVs
	for(i = 0; i<nbFixFlows; i++)
		df0[i] = stripper.df((*fixPayTimes)[i]);

	double floatLeg0 = stripper.stdFundingLeg(*(oswVanillaArg.GetFloatStartTimes()),
		*(oswVanillaArg.GetFloatEndTimes()), *(oswVanillaArg.GetFloatIntTerms()),
		*(oswVanillaArg.GetFloatNotional()), *(oswVanillaArg.GetFundingSpread()));

	double swapValue = floatLeg0 - fixLeg0;
	
	ARM_GP_Vector sensitivities(sizeSched, 0.0);
	double floatLeg,numeraire,numerizedLiborFlow;
	ARM_GP_Vector invNumeraireSensis (sizeSched, 0.0);
	ARM_GP_Matrix toInvert(sizeSched, sizeSched);
		
	/// Bump yield curve swap rates and compute sensitivities
	/// of underlying swap of input swaption 
	for (j=0; j<sizeSched; j++)
	{
		/// shift fwd swap rate #j
		stripper.setSwapRate(j, forwards[j] + SENSI_CMS);
		stripper.strip();

		/// (1) for spread option leg
		for (i = 0; i<nbFixFlows; i++)
		{
			df = stripper.df((*fixPayTimes)[i]);
			
			sensitivities[j] -= A[i] * (df - df0[i]) / SENSI_CMS;
		}

		/// (2) for funding leg + margin
		floatLeg = stripper.stdFundingLeg(*(oswVanillaArg.GetFloatStartTimes()),
		*(oswVanillaArg.GetFloatEndTimes()), *(oswVanillaArg.GetFloatIntTerms()),
		*(oswVanillaArg.GetFloatNotional()), *(oswVanillaArg.GetFundingSpread()));

		sensitivities[j] += (floatLeg - floatLeg0) / SENSI_CMS;

		/// (3) for 1/Numeraire
		numeraire = 0.0;
		for (size_t k=0; k<sizeSched; k++)
		{	
			df = stripper.df((*endDates)[k]-asOfDate);
			numeraire += (*interestTerms)[k] * (stdFixNotionals[k]/refStdFixNotional) * df;
		}
		invNumeraireSensis[j] = (1./numeraire - 1./numeraire0) / SENSI_CMS;

		for (k=0; k<sizeSched; k++)
		{
			numerizedLiborFlow = (stripper.df((*startDates)[k]-asOfDate) -  stripper.df((*endDates)[k]-asOfDate)) / numeraire; 
			toInvert(j,k) = (numerizedLiborFlow-numerizedLiborFlow0[k]) / SENSI_CMS;
		}
		

		/// RAZ
		stripper.setSwapRate(j, forwards[j]);
	}

	/// Reset stripper
	stripper.strip();

	///
	/// build matrix to invert
	/// for the moment its contains the sensis of the Underlying w.r.t swap rates
	/// now, we want it to contain the sensi of Underlying/Numeraire w.r.t. swap rates
	for (j=0; j<sizeSched; j++)
	{
		sensitivities[j] *= 1./numeraire0;
		sensitivities[j] += swapValue * invNumeraireSensis[j];
	}
	
	/// invert matrix
	LinSolve(&toInvert, &sensitivities);
	ARM_Vector stdFloatNotionals(sizeSched);
	for (j=0; j<sizeSched; j++)
		stdFloatNotionals[j] = sensitivities[j];
	
		
	/// Compute strike assuming the same PV for VNS and underlying swap
	double strike;
	double stdFloatLeg = 0.0;

	for (j=0; j<sizeSched; j++)
	{			
		stdFloatLeg += stdFloatNotionals[j] * (stripper.df((*startDates)[j]-asOfDate) - stripper.df((*endDates)[j]-asOfDate));
	}

	strike = (stdFloatLeg - swapValue) / (numeraire0 * refStdFixNotional) ;
	

	///  vector need to be cloned (deleted in ARM_ReferenceValue destructor)
	ARM_ReferenceValue armSdtFloatNotionals((ARM_Vector*)stdSwapDates.Clone(), (ARM_Vector*)stdFloatNotionals.Clone());
	armSdtFloatNotionals.SetCalcMethod(K_STEPUP_RIGHT);

	///  vectors need to be cloned (deleted in ARM_ReferenceValue destructor)
	ARM_ReferenceValue armStdFixNotionals((ARM_Vector*)stdSwapDates.Clone(), (ARM_Vector*)stdFixNotionals.Clone());
	armStdFixNotionals.SetCalcMethod(K_STEPUP_RIGHT);


	bool isFixPayer = oswVanillaArg.GetCallPut() == K_CALL;
		
	/// index type: same frequency as fixed leg
	ARM_IRIndex irIndex (indexType, stdFixFreq, stdFixFreq, ccy);
	irIndex.SetTerm(stdFixFreq);
	irIndex.SetYearTerm(1.0/stdFixFreq);

	ARM_SwapLeg armFloatLeg( startDate, 
							 endDate, 
							 &irIndex, 
							 isFixPayer ? K_RCV : K_PAY,
							 0.0, 
							 K_SHORTSTART, 
							 K_COMP_PROP,
							 ccy,
							 ccy->GetLiborIndexDayCount());

	/// RefVal is cloned in SetAmount
	armFloatLeg.SetAmount(&armSdtFloatNotionals);


	// create fixed leg with good strike
	ARM_FixLeg armFixLeg ( startDate,
						   endDate, 
						   strike * 100.0,
						   isFixPayer ? K_PAY : K_RCV, 
						   stdFixFreq,
						   stdFixDayCount,
						   K_COMP_PROP,
						   K_ARREARS, 
						   K_ADJUSTED,
						   K_SHORTSTART,
						   ccy);
	
	/// RefVal is cloned in SetAmount
	armFixLeg.SetAmount(&armStdFixNotionals);
	
	
	/// create swap
	ARM_Swap armSwap(&armFixLeg, &armFloatLeg);
			
	/// create swaption
	ARM_Date exerDate(expiryDate);
	ARM_Swaption* armSwaption = new ARM_Swaption((&armSwap), isFixPayer ? K_PAY : K_RCV,K_EUROPEAN, strike * 100.0,exerDate);

	if (resetCalendar)
		delete resetCalendar;

	if (payCalendar)
		delete payCalendar;

	return armSwaption;
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
