/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file maturitycapcalculator.cpp
 *
 *  \brief calculator for the maturity cap
 *
 *	\author  R. GUILLEMOT
 *	\version 1.0
 *	\date March 2004
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalculators/maturitycapcalculator.h"

/// gpbase
#include "gpbase/env.h"
#include "gpbase/datestripcombiner.h"
#include "gpbase/autocleaner.h"
#include "gpbase/singleton.h"
#include "gpbase/datestrip.h"
#include "gpbase/functor.h"
#include "gpbase\gptrigomatrix.h"

/// gpinfra
#include "gpinfra/argconvdefault.h"
#include "gpinfra/gensecurity.h"
#include "gpinfra/cstmanager.h"
#include "gpinfra/dealdescription.h"
#include "gpinfra/mktdatamanagerrep.h"
#include "gpinfra/pricingadviser.h"
#include "gpinfra/genpricer.h"
#include "gpinfra/numeraire.h"
#include "gpinfra/numerairefactory.h"
#include "gpinfra/pricerinfo.h"
#include "gpinfra/gramfunctorargdict.h"
#include "gpinfra/correlmatparam.h"
#include "gpinfra/curvemodelparam.h"

/// gpcalib
#include "gpcalib/calibmethod.h"
#include "gpcalib/vanillaarg.h"
#include "gpcalib/kerneltogp.h"
#include "gpcalib\modelparamsfactory.h"

/// gpmodels
#include "gpmodels/modelparamssfrm.h"
#include "gpmodels/modelparamssfrmdiag.h"
#include "gpmodels/modelparamssfrmfactory.h"
#include "gpmodels/sfrm.h"
#include "gpmodels/irfwdmod.h"

/// gpnumlib
#include "gpnumlib/randomgenfactory.h"
#include "gpnumlib/compositegen.h"
#include "gpnumlib/antitheticgen.h"
#include "gpnumlib/solver.h"
#include "gpnumlib/newtonraphson.h"
#include "gpnumlib/numfunction.h"

/// gpnummethods
#include "gpnummethods/mcmethod.h"
#include "gpnummethods/finummethod.h"
#include "gpnummethods/normalcentredsampler.h"
#include "gpnummethods/scheduler.h"

/// kernel
#include <inst/portfolio.h>
#include <inst/forex.h>
#include <mod/bssmiled.h>
#include <glob/paramview.h>

/// STL
#include <iomanip> /// for setprecision()
#include <memory>


CC_BEGIN_NAMESPACE( ARM )

const string ARM_MaturityCapCalculator::MaturityCapColNamesTable [] =
{
    "EventDate",
    "StartDate",
    "EndDate",
	"FirstCouponDate",
	"UnderlyingEndDate",
    "PaymentDate",
	"IndexStartDate",
    "IT",
	"Coeff",
	"DF",
	"Annuity",
	"Yield0",
	"TRINominal",
	"INFINENominal",
	"INFINENominalForCalib",
	"RefStdCapNominal",
	"Amortizing",
    "Index",
    "Spread",
    "Yield",
	"INFINEYield",
	"TRIInterest",
	"INFINEInterest",
	"Time",
	"Volatility",
	"Shock",
	"EstimTRIYield",
	"EstimTRINominal",
	"EstimINFINEYield",
	"EstimINFINENominal",
	"TRIMaturityCap",
	"INFINEMaturityCap",
	"INFINENominalFinal",
	"RefStdCap",
};

// SFRM default factors number
const int SFRM_NB_FACTORS			= 2;
const int SFRM_VOL_TYPE				= K_DIAG;

// MC default pricing values
const double MC_NB_INTER_STEPS		= 1;

/// SFRM sigma range [1bp,10000bp] with a 500bp default value
const double SIGMA_LOWER_BOUND			= 0.001;
const double SIGMA_UPPER_BOUND			= 1.0;
const double SIGMA_DEFAULT_VALUE		= 0.2;
const double MRS_DEFAULT_VALUE			= 0.0;
const double BETA_DEFAULT_VALUE			= 1.0;
const double CORREL_MAT_DEFAULT_VALUE   = 20;
const double CORREL_THETA_DEFAULT_VALUE   = 0.5;

const double CF_DEFAULT_PRICE=1.0e+100;
const double CF_DEFAULT_WEIGHT=1.0;

/// Default MDM key names
const string YC_KEY_NAME            = "YC_";
const string CFMODEL_KEY_NAME       = "CFMOD_";
const string MRS_KEY_NAME           = "MRS_";
const string BETA_KEY_NAME          = "BETA_";
const string CORREL_KEY_NAME        = "CORREL_";

/// Reference schedules for Maturity Cap date structure
const unsigned int LOANRESET_SCHED =0;
const unsigned int LOANPAY_SCHED =1;
const unsigned int NB_MATCAP_SCHED =2;

// Shock calculation constant
const double MINSHOCK			= 0.0;
const double MAXSHOCK			= 3.0;
const double SHOCKPRECISION		= 1e-3;
const double SHOCKINC			= 0.5;
const string SHOCKCST			= "SHOCK";

const double MaxYield			= 5;

const double DAY_TOLERANCE		= 7;
const double DEFAULT_PREC		= 1e-8;


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MatruityCap
///	Routine: Constructor
///	Returns: void
///	Action : builds the object (stand alone version)
/////////////////////////////////////////////////////////////////

ARM_MaturityCapCalculator::ARM_MaturityCapCalculator(
	const ARM_Date& startDate,
    const ARM_Date& endDate,
	const ARM_Date& underlyingEndDate,
    int longShort,
    int capFloor,
    int resetFreq,
    int payFreq,
    const string& indexTerm,
    int dayCount,
    int intRule,
    double spread,
    double initNominal,
    double initTRI,
	double annuity,
    ARM_MaturityCapCalculator::ProductMode productMode,
    double coeff,
    const ARM_Curve& amortizing,
	int resetGap,
    const string& payCal,
	const string& resetCal,
	CalibrationMode calibrationMode,
	int nbIterations,
	std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
    const ARM_MarketData_ManagerRep& mktDataManager,
    const ARM_StringVector& mdmKeys
	) : ARM_GenCalculator(mktDataManager),
	itsStartDate(startDate),
	itsEndDate(endDate),
	itsUnderlyingEndDate(underlyingEndDate),
	itsLongShort(longShort),
	itsCapFloor(capFloor),
	itsResetFreq(resetFreq),
	itsPayFreq(payFreq),
	itsIndexTerm(indexTerm),
	itsDayCount(dayCount),
	itsIntRule(intRule),
	itsSpread(spread),
	itsInitNominal(initNominal),
	itsInitTRI(initTRI),
	itsAnnuity(annuity),
	itsProductMode(productMode),
	itsCoeff(coeff),
	itsAmortizing(amortizing),
	itsResetGap(resetGap),
	itsResetCal(resetCal),
	itsPayCal(payCal),
	itsCalibrationMode(calibrationMode),
	itsNbIterations(nbIterations),
	itsProductsToPrice(productsToPrice),
	itsHasBeenPriced(false),
	itsMaturityCapPrice(0.0),
	itsMaturityCapStdDev(0.0),
	itsRefStdCapPrice(0.0),
	itsRefStdCapStdDev(0.0)	
{
	SetName(ARM_MATURITYCAP);

	SetKeys(mdmKeys);

	// Check consitency
	CheckDataAndTimeIt();


    /// Set currency
    SetCurrencyUnit( (static_cast< ARM_ZeroCurve* >(GetMktDataManager()->GetData(GetKeys()[YcKey]))->GetCurrencyUnit()) );
  
    /// Set the coupon & funding model name alias i.e. keys to access to
    /// the right model for each GP keyword (1st argument)
    SetModelKeys();

    /// Create the Generic Security
    CreateAndSetDealDescriptionAndTimeIt();

    /// Create the SFRM 2F pricing model with its default parameters
    CreateAndSetModelAndTimeIt();

    /// Create the calibration set for volatility bootstapping
    CreateAndSetCalibrationAndTimeIt();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_MatruityCap
///	Routine: SUMMIT Constructor
///	Returns: void
///	Action : builds the object (stand alone version)
/////////////////////////////////////////////////////////////////

ARM_MaturityCapCalculator::ARM_MaturityCapCalculator(
		const ARM_Date& asOfDate,
		const ARM_Date& startDate,
		const ARM_Date& endDate,
		const ARM_Date& underlyingEndDate,
		int longShort,
		int capFloor,
		int resetFreq,
		int payFreq,
		const string& indexTerm,
		int dayCount,
		int intRule,
		double spread,
		double initNominal,
		double initTRI,
		double annuity,
		ARM_MaturityCapCalculator::ProductMode productMode,
		double coeff,
		const ARM_Curve& amortizing,
		int resetGap,
		const string& payCal,
		const string& resetCal,
		CalibrationMode calibrationMode,
		int nbIterations,
		std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
		const ARM_Currency& ccy
		)  : ARM_GenCalculator(asOfDate),
	itsStartDate(startDate),
	itsEndDate(endDate),
	itsUnderlyingEndDate(underlyingEndDate),
	itsLongShort(longShort),
	itsCapFloor(capFloor),
	itsResetFreq(resetFreq),
	itsPayFreq(payFreq),
	itsIndexTerm(indexTerm),
	itsDayCount(dayCount),
	itsIntRule(intRule),
	itsSpread(spread),
	itsInitNominal(initNominal),
	itsInitTRI(initTRI),
	itsAnnuity(annuity),
	itsProductMode(productMode),
	itsCoeff(coeff),
	itsAmortizing(amortizing),
	itsResetGap(resetGap),
	itsResetCal(resetCal),
	itsPayCal(payCal),
	itsCalibrationMode(calibrationMode),
	itsNbIterations(nbIterations),
	itsProductsToPrice(productsToPrice),
	itsHasBeenPriced(false),
	itsMaturityCapPrice(0.0),
	itsMaturityCapStdDev(0.0),
	itsRefStdCapPrice(0.0),
	itsRefStdCapStdDev(0.0)	
{
	SetName(ARM_MATURITYCAP);

    /// Set the coupon/payment currency (inherited from ARM_Security)
    SetCurrencyUnit(const_cast< ARM_Currency* >(&ccy));
    string ccyName(ccy.GetCcyName());


    /// Set default keys for MDM data access
    ARM_StringVector keys(NbKeys);
    keys[YcKey]             = YC_KEY_NAME       + ccyName;
    keys[CfModelKey]        = CFMODEL_KEY_NAME  + ccyName;
    keys[MrsKey]            = MRS_KEY_NAME      + ccyName;
	keys[BetaKey]           = BETA_KEY_NAME     + ccyName;
	keys[CorrelKey]         = CORREL_KEY_NAME   + ccyName;

    SetKeys(keys);
 
	// Check consitency
	CheckDataAndTimeIt();

    /// Set the coupon & funding model name alias i.e. keys to access to
    /// the right model for each GP keyword (1st argument)
    SetModelKeys();

	/// Create the Generic Security
    CreateAndSetDealDescriptionAndTimeIt();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_MaturityCapCalculator
///	Routine: Copy constructor
///	Returns: a built object
///	Action : builds the object
/////////////////////////////////////////////////////////////////
ARM_MaturityCapCalculator::ARM_MaturityCapCalculator(const ARM_MaturityCapCalculator& rhs)
: ARM_GenCalculator(rhs),
itsStartDate(rhs.itsStartDate),
itsEndDate(rhs.itsEndDate),
itsUnderlyingEndDate(rhs.itsUnderlyingEndDate),
itsLongShort(rhs.itsLongShort),
itsCapFloor(rhs.itsCapFloor),
itsResetFreq(rhs.itsResetFreq),
itsPayFreq(rhs.itsPayFreq),
itsIndexTerm(rhs.itsIndexTerm),
itsDayCount(rhs.itsDayCount),
itsIntRule(rhs.itsIntRule),
itsSpread(rhs.itsSpread),
itsInitNominal(rhs.itsInitNominal),
itsInitTRI(rhs.itsInitTRI),
itsProductMode(rhs.itsProductMode),
itsCoeff(rhs.itsCoeff),
itsAmortizing(rhs.itsAmortizing),
itsResetGap(rhs.itsResetGap),
itsResetCal(rhs.itsResetCal),
itsPayCal(rhs.itsPayCal),
itsCalibrationMode(rhs.itsCalibrationMode),
itsNbIterations(rhs.itsNbIterations),
itsProductsToPrice(rhs.itsProductsToPrice),
itsHasBeenPriced(false),
itsMaturityCapPrice(0.0),
itsMaturityCapStdDev(0.0),
itsRefStdCapPrice(0.0),
itsRefStdCapStdDev(0.0)
{
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_MaturityCapCalculator
///	Routine: assignment operator
///	Returns: 
///	Action : 
/////////////////////////////////////////////////////////////////
ARM_MaturityCapCalculator& ARM_MaturityCapCalculator::operator=(const ARM_MaturityCapCalculator& rhs)
{
	if (this != &rhs)
	{
		ARM_GenCalculator::operator=(rhs);

		itsStartDate = rhs.itsStartDate;
		itsEndDate = rhs.itsEndDate;
		itsUnderlyingEndDate = rhs.itsUnderlyingEndDate;
		itsLongShort = rhs.itsLongShort;
		itsCapFloor = rhs.itsCapFloor;
		itsResetFreq = rhs.itsResetFreq;
		itsPayFreq = rhs.itsPayFreq;
		itsIndexTerm = rhs.itsIndexTerm;
		itsDayCount = rhs.itsDayCount;
		itsIntRule = rhs.itsIntRule;
		itsSpread = rhs.itsSpread;
		itsInitNominal = rhs.itsInitNominal;
		itsInitTRI = rhs.itsInitTRI;
		itsProductMode = rhs.itsProductMode;
		itsCoeff = rhs.itsCoeff;
		itsAmortizing = rhs.itsAmortizing;
		itsResetGap = rhs.itsResetGap;
		itsResetCal = rhs.itsResetCal;
		itsPayCal = rhs.itsPayCal;
		itsCalibrationMode = rhs.itsCalibrationMode;
		itsNbIterations = rhs.itsNbIterations;
		itsProductsToPrice = rhs.itsProductsToPrice;
		itsHasBeenPriced = false;
		itsMaturityCapPrice = 0.0;
		itsMaturityCapStdDev = 0.0;
		itsRefStdCapPrice = 0.0;
		itsRefStdCapStdDev = 0.0;
	}

	return *this;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MaturityCapCalculator
///	Routine: Destructor
///	Returns: void
///	Action : destroys the object
/////////////////////////////////////////////////////////////////
ARM_MaturityCapCalculator::~ARM_MaturityCapCalculator()
{
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_MaturityCapCalculator
///	Routine: Set???ToPrice() and Is???ToPrice()
///	Returns: void/boolean
///	Action : Set the flag to know which product to price.
///          Test the product currently able to be priced
/////////////////////////////////////////////////////////////////
void ARM_MaturityCapCalculator::SetMaturityCapToPrice(bool toPrice)			{itsProductsToPrice[MaturityCapPrice]				=   toPrice;}
bool ARM_MaturityCapCalculator::IsMaturityCapToPrice() const					{return itsProductsToPrice[MaturityCapPrice];}

void ARM_MaturityCapCalculator::SetRefStdCapToPrice(bool toPrice)			{itsProductsToPrice[RefStdCapPrice]				=   toPrice;}
bool ARM_MaturityCapCalculator::IsRefStdCapToPrice() const					{return itsProductsToPrice[RefStdCapPrice];}

void ARM_MaturityCapCalculator::SetEstimatedTRIToPrice(bool toPrice)	{itsProductsToPrice[EstimatedTRI]						=   toPrice;}
bool ARM_MaturityCapCalculator::IsEstimatedTRIToPrice() const			{return itsProductsToPrice[EstimatedTRI];}

void ARM_MaturityCapCalculator::SetEstimatedNominalToPrice(bool toPrice)	{itsProductsToPrice[EstimatedNominal]				=   toPrice;}
bool ARM_MaturityCapCalculator::IsEstimatedNominalToPrice() const			{return itsProductsToPrice[EstimatedNominal];}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MaturityCapCalculator
///	Routine: SetModelKeys
///	Returns: void
///	Action : set the model name
/////////////////////////////////////////////////////////////////
void ARM_MaturityCapCalculator::SetModelKeys()
{
	itsModelKey = YcKey;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_MaturityCapCalculator
///	Routine: GetMRS & SetMRS
///	Returns: 
///	Action : get & set the MRS param for SUMMIT interface
/////////////////////////////////////////////////////////////////
const ARM_CurveModelParam& ARM_MaturityCapCalculator::GetMRS() const
{
    return *(static_cast< ARM_CurveModelParam* >(GetMktDataManager()->GetData(GetKeys()[MrsKey])));
}

void ARM_MaturityCapCalculator::SetMRS(ARM_CurveModelParam* mrsParam)
{
    if(!mrsParam || mrsParam->GetType() != ARM_ModelParamType::MeanReversion)
		    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : an MRS Param is expected !" );

    GetMktDataManager()->RegisterData(GetKeys()[MrsKey],static_cast< ARM_Object* >(mrsParam));
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_MaturityCapCalculator
///	Routine: GetMRS & SetMRS
///	Returns: 
///	Action : get & set the Beta param for SUMMIT interface
/////////////////////////////////////////////////////////////////
const ARM_CurveModelParam& ARM_MaturityCapCalculator::GetBeta() const
{
    return *(static_cast< ARM_CurveModelParam* >(GetMktDataManager()->GetData(GetKeys()[BetaKey])));
}

void ARM_MaturityCapCalculator::SetBeta(ARM_CurveModelParam* betaParam)
{
    if(!betaParam || betaParam->GetType() != ARM_ModelParamType::Beta)
		    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : an Beta Param is expected !" );

    GetMktDataManager()->RegisterData(GetKeys()[BetaKey],static_cast< ARM_Object* >(betaParam));
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_MaturityCapCalculator
///	Routine: GetCorrel & SetCorrel
///	Returns: 
///	Action : get & set the Correlation param for SUMMIT interface
/////////////////////////////////////////////////////////////////
const ARM_CorrelMatParam& ARM_MaturityCapCalculator::GetCorrel() const
{
    return *(static_cast< ARM_CorrelMatParam* >(GetMktDataManager()->GetData(GetKeys()[CorrelKey])));
}

void ARM_MaturityCapCalculator::SetCorrel(ARM_CorrelMatParam* correlParam)
{
    if(!correlParam || correlParam->GetType() != ARM_ModelParamType::Correlation)
		    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : a Correl Param is expected !" );

    GetMktDataManager()->RegisterData(GetKeys()[CorrelKey],static_cast< ARM_Object* >(correlParam));
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_MaturityCapCalculator
///	Routine: GetCFCalibMethod
///	Returns: ARM_CalibMethod
///	Action : get calibration method for cap floor
/////////////////////////////////////////////////////////////////
ARM_CalibMethod* ARM_MaturityCapCalculator::GetCFCalibMethod()
{
#ifdef __GP_STRICT_VALIDATION
    if( GetCalibMethod() == ARM_CalibMethodPtr(NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Bootstrap calib method not found");
#endif

        return &(*GetCalibMethod());
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MaturityCapCalculator
///	Routine: GetCFPortfolio
///	Returns: ARM_Portfolio
///	Action : get the cap floor calibration portfolio
/////////////////////////////////////////////////////////////////
const ARM_StdPortfolioPtr ARM_MaturityCapCalculator::GetCFPortfolio(void)
{
#ifdef __GP_STRICT_VALIDATION
    if( (GetCalibMethod() == ARM_CalibMethodPtr(NULL)) ||
		(GetCalibMethod()->GetPortfolio() == ARM_StdPortfolioPtr(NULL)))
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Bootstrap calib method or portfolio not found");
#endif
	
	return GetCalibMethod()->GetPortfolio();
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MaturityCapCalculator
///	Routine: SetCFPortfolio
///	Returns: void
///	Action : affectation of a calibration portfolio made of cap floors
/////////////////////////////////////////////////////////////////
void ARM_MaturityCapCalculator::SetCFPortfolio(const ARM_StdPortfolio& port)
{
    int pfSize=port.GetSize();
    if(pfSize<1)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Not any asset in the given portfolio");


    /// Test for portfolio consistency
    for(int i=0;i<pfSize;++i)
    {
        if(port.GetAsset(i)->GetName() != ARM_CAPFLOOR)
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Diagonal swaption portfolio is not only made of swaptions");
    }

#ifdef __GP_STRICT_VALIDATION
    if( GetCalibMethod() == ARM_CalibMethodPtr(NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Bootstrap calib method not found");
#endif

	GetCalibMethod()->SetPortfolio( ARM_StdPortfolioPtr((ARM_StdPortfolio*)(const_cast< ARM_StdPortfolio& >(port).Clone())) );
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_MaturityCapCalculator
///	Routine: InitMaturityCapForSummit
///	Returns: nothing
///	Action : Feed the MDM with input market objects
/////////////////////////////////////////////////////////////////
void ARM_MaturityCapCalculator::InitMaturityCapForSummit(ARM_ZeroCurve* zc, 
                                                         ARM_VolCurve* capVC, 
                                                         ARM_VolLInterpol* capRo, 
                                                         ARM_VolLInterpol* capNu, 
                                                         ARM_VolLInterpol* capBeta,
                                                         int SABRSigmaOrAlpha) // SABR Sigma or Alpha
                                                                               // By default=1(Sigma)
{
    /// Give a default value to MRS if not set.
    ARM_CurveModelParam mr;
    if(GetMktDataManager()->TestIfKeyMissing(GetKeys()[MrsKey]))
    {
        ARM_GP_Vector values(1,MRS_DEFAULT_VALUE);
        ARM_GP_Vector times(1,0.0);
        mr = ARM_CurveModelParam(ARM_ModelParamType::MeanReversion,&values,&times,"MRS");
    }
    else
    {
	    mr = GetMRS();
    }
	/// Give a default value to Beta if not set.
	ARM_CurveModelParam beta;
	if(GetMktDataManager()->TestIfKeyMissing(GetKeys()[BetaKey]))
    {
        ARM_GP_Vector values(1,BETA_DEFAULT_VALUE);
        ARM_GP_Vector times(1,0.0);
        beta = ARM_CurveModelParam(ARM_ModelParamType::Beta,&values,&times,"Beta");
    }
    else
    {
	    beta = GetBeta();
    }
	/// Give a default value to Correl if not set.
	bool localCorrel = false;
	ARM_CorrelMatParam* correl;
	if(GetMktDataManager()->TestIfKeyMissing(GetKeys()[CorrelKey]))
    {
		localCorrel = true;
		int n = CORREL_MAT_DEFAULT_VALUE;
		
		// Buid a DateStrip to have business days in times
		ARM_Date maxDate(zc->GetAsOfDate());
		maxDate.AddYears(n);
	    ARM_DateStrip resetLegSched(zc->GetAsOfDate(), 
									maxDate, 
									K_ANNUAL,
									KACTUAL_ACTUAL,
									"DEFAULT",
									K_MOD_FOLLOWING,
									K_ADJUSTED,
									K_SHORTSTART,
									GETDEFAULTVALUE,
									K_ANNUAL); // obliger d'aller jusqu'au PayFreq car valeur par defaut non valide...

		// Conversion from ARM_Vector To ARM_RealVector
		ARM_GP_Vector* times = NULL;
		times = resetLegSched.GetResetDates();
		*times -= zc->GetAsOfDate().GetJulian(); // times will be deleted in DateStrip destructor

		// Conversion from ARM_GP_Matrix To ARM_RealVector
		ARM_GP_Matrix* trigoMatrix = TrigoMatrix(n,CORREL_THETA_DEFAULT_VALUE);
        ARM_GP_Vector* values = new ARM_GP_Vector(n*n);

		size_t i = 0,j = 0;
		for (i = 0; i < n; ++i)
			for (j = 0; j < n; ++j)
				(*values)[i*n+j] = (*trigoMatrix)(i,j);
        

        correl = static_cast<ARM_CorrelMatParam*>(ARM_ModelParamFactory.Instance()->CreateModelParam(
			ARM_ModelParamType::Correlation,
            values,
            times));

		if (trigoMatrix)
			delete trigoMatrix;
			trigoMatrix = NULL;
		if (values)
			delete values;
			values = NULL;
    }
    else
    {
	    correl = const_cast<ARM_CorrelMatParam*>(&GetCorrel());
    }

	vector <ARM_Object*> marketDatas;

	marketDatas.push_back(zc);

	ARM_BSModel* capBSmod = NULL;
	
	//SABR with Beta = 1: 
	if (capRo && capNu && !capBeta)
	{
		capBSmod = new ARM_BSSmiledModel(zc->GetAsOfDate(), 0,
										 zc,
										 zc,
										 capVC,
										 K_YIELD,
										 capRo,
										 capNu,
										 K_SABR_ARITH,
                                         NULL,
                                         0.5, // SABR Weight: Irrelevant
                                         SABRSigmaOrAlpha);
	}
	//Complete SABR :
	else if (capRo && capNu && capBeta)
	{
		int methodType;

		if (capBeta->IsEqualToOne())
			methodType = K_SABR_ARITH;
		else
			methodType = K_SABR_IMPLNVOL;

		capBSmod = new ARM_BSSmiledModel(zc->GetAsOfDate(),
                                         0,
										 zc,
										 zc,
										 capVC,
										 K_YIELD,
										 capRo,
										 capNu,
										 methodType,
										 capBeta,
                                         0.5, // SABR Weight: Irrelevant
                                         SABRSigmaOrAlpha);


	}
	//No SABR:  --> Vol Cube
	else if (!capRo && !capNu && !capBeta)
	{
		capBSmod = new ARM_BSModel(zc->GetAsOfDate(),
								   0.0,
								   zc,
								   zc,
								   capVC,
								   K_YIELD);

	}

	marketDatas.push_back(capBSmod);
	marketDatas.push_back(&mr);
	marketDatas.push_back(&beta);
	marketDatas.push_back(correl);

	Init(marketDatas);

	if (capBSmod)
	{
		delete capBSmod;
		capBSmod = NULL;
	}

	if (localCorrel)
	{
		delete correl;
		correl = NULL;
	}
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MaturityCapCalculator
///	Routine: CheckData
///	Returns: void
///	Action : check if Maturity Cap datas are consistent
/////////////////////////////////////////////////////////////////
void ARM_MaturityCapCalculator::CheckData()
{
	// MaturityCap parameters checking

	if ( itsInitNominal < K_NEW_DOUBLE_TOL )
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Init nominal has to be strictly positive.");
	}

	if ( itsInitTRI < K_NEW_DOUBLE_TOL )
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Init TRI has to be strictly positive.");
	}

	if ( itsEndDate < itsStartDate )
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The Start Date of the deal has to be before the end date.");
	}

	if ( itsUnderlyingEndDate < itsEndDate )
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The End Date of the deal has to be before the underlying end date.");
	}

	if ( itsPayFreq < itsResetFreq )
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The payment frequency has to be greater  .");
	}

	if (!(itsAmortizing  > -K_NEW_DOUBLE_TOL))
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Amortizing has to be strictly positive.");
	}

	if ( itsNbIterations <= 0 )
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Nb iterations has to be strictly positive.");
	}

	ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();

	ARM_Date firstResetDate(itsStartDate);

	firstResetDate.GapBusinessDay(itsResetGap, const_cast<char *>(itsResetCal.c_str()));

	if (asOfDate > firstResetDate)
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The first reset date should be after the asofdate (the calculator cannot price in the past).");
	}

	if (itsIndexTerm == "12M")
	{
		itsIndexTerm = "1Y";
	}
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MaturityCapCalculator
///	Routine: CheckMktData
///	Returns: void
///	Action : check if Maturity Cap mkt datas are consistent
/////////////////////////////////////////////////////////////////
void ARM_MaturityCapCalculator::CheckMktData()
{
	/// MdM datas checking
	ARM_ZeroCurve* cpnCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));
    if(!cpnCurve)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Yield Curve for key=" + GetKeys()[YcKey] + " is expected in the Market Data Manager");
    string cpnCcy(cpnCurve->GetCurrencyUnit()->GetCcyName());


	ARM_BSModel* cfBSModel = dynamic_cast< ARM_BSModel* >(GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );
    if(!cfBSModel)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : cap/floor B&S model for key=" + GetKeys()[CfModelKey] + " is expected in the Market Data Manager");



	ARM_ModelParam* mrsParam = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[MrsKey]));
    if(!mrsParam || mrsParam->GetType() != ARM_ModelParamType::MeanReversion)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : MRS Param for key=" + GetKeys()[MrsKey] + " is expected in the Market Data Manager");

	ARM_ModelParam* betaParam = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[BetaKey]));
    if(!betaParam || betaParam->GetType() != ARM_ModelParamType::Beta)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Beta Param for key=" + GetKeys()[BetaKey] + " is expected in the Market Data Manager");

	ARM_ModelParam* correlParam = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[CorrelKey]));
    if(!correlParam || correlParam->GetType() != ARM_ModelParamType::Correlation)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Correl Param for key=" + GetKeys()[CorrelKey] + " is expected in the Market Data Manager");
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_MaturityCapCalculator
///	Routine: ColumnNames
///	Returns: ARM_RowInfo
///	Action : create the column names of the deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_MaturityCapCalculator::ColumnNames() const
{
	size_t colNamesSize = sizeof(MaturityCapColNamesTable)/sizeof(MaturityCapColNamesTable[0]);
    vector< string > colNamesVec(colNamesSize);
    vector< ARM_GP_VALUE_TYPE > colTypeVec(colNamesSize, ARM_STRING); 

    for(size_t i=0;i<colNamesSize; ++i)
        colNamesVec[i] = MaturityCapColNamesTable[i];

    ARM_RowInfo rowInfo(colNamesVec,colTypeVec);

    return rowInfo;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_MaturityCapCalculator
///	Routine: MiddleRows
///	Returns: ARM_RowInfo
///	Action : create a row of a deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_MaturityCapCalculator::MiddleRows( size_t eventIdx, const ARM_DateStripCombiner& datesStructure ) const
{
    /// Here only exercise events are generated. At each date the coupon
    /// value is known with its analytical formula because
    /// keywords LIBOR & CAPLET were extended to handle the in arrears case


    size_t eventSize = datesStructure.GetDateStrip(0)->GetResetDates()->size();
    size_t descSize = sizeof(MaturityCapColNamesTable)/sizeof(MaturityCapColNamesTable[0]);

    vector< string > rowDescVec(descSize);
    vector< ARM_GP_VALUE_TYPE > rowTypeVec(descSize, ARM_MISSING_TYPE); 

	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

	string modelName = GetKeys()[itsModelKey];


	string payFreq=ARM_ArgConvReverse_MatFrequency.GetString(itsPayFreq);
    string dayCount=ARM_ArgConvReverse_DayCount.GetString(itsDayCount);

	bool isFirstEvent = (eventIdx == 0);
	bool isLastEvent = (eventIdx == eventSize-1);

	double endDate=(*(datesStructure.GetDateStrip(LOANPAY_SCHED)->GetFlowEndDates()))[eventIdx];

	bool beforeEndDate = (endDate <= itsEndDate.GetJulian()+DAY_TOLERANCE);

	/// Set default 0 value for each column to be able to sum it
    InitPriceableColumns(rowDescVec,rowTypeVec);

    /// EventDate description
    double eventDate=(*(datesStructure.GetMergeData()))[eventIdx];
    CC_Ostringstream eventDateDesc;
    eventDateDesc << CC_NS(std,fixed) << eventDate;
    rowDescVec[EventDate] = eventDateDesc.str();
    rowTypeVec[EventDate] = ARM_DATE_TYPE;

	// StartDate description
	double startDate=(*(datesStructure.GetDateStrip(LOANPAY_SCHED)->GetFlowStartDates()))[eventIdx];
	CC_Ostringstream startDateDesc;
    startDateDesc << CC_NS(std,fixed) << startDate;
    rowDescVec[StartDate] = startDateDesc.str();
    rowTypeVec[StartDate] = ARM_DATE_TYPE;

	// EndDate description
	CC_Ostringstream endDateDesc;
    endDateDesc << CC_NS(std,fixed) << endDate;
    rowDescVec[EndDate] = endDateDesc.str();
    rowTypeVec[EndDate] = ARM_DATE_TYPE;

	// FirstCouponDate description
	CC_Ostringstream firstCouponDateDesc;
	ARM_Date firstCouponDate(startDate);
	firstCouponDate.AddPeriod(payFreq);
    firstCouponDateDesc << CC_NS(std,fixed) << firstCouponDate.GetJulian();
    rowDescVec[FirstCouponDate] = firstCouponDateDesc.str();
    rowTypeVec[FirstCouponDate] = ARM_DATE_TYPE;

	// UnderlyingEndDate description
	double underlyingEndDate=itsUnderlyingEndDate.GetJulian();
	CC_Ostringstream underlyingEndDateDesc;
    underlyingEndDateDesc << CC_NS(std,fixed) << underlyingEndDate;
    rowDescVec[UnderlyingEndDate] = underlyingEndDateDesc.str();
    rowTypeVec[UnderlyingEndDate] = ARM_DATE_TYPE;

	// PaymentDate description
	double paymentDate=(*(datesStructure.GetDateStrip(LOANPAY_SCHED)->GetPaymentDates()))[eventIdx];
	CC_Ostringstream paymentDateDesc;
    paymentDateDesc << CC_NS(std,fixed) << paymentDate;
    rowDescVec[PaymentDate] = paymentDateDesc.str();
    rowTypeVec[PaymentDate] = ARM_DATE_TYPE;

	/// Cpn Index Date and interest term descriptions
	CC_Ostringstream indexStartDateDesc;
	
	ARM_Date indexStartDate(eventDate);
	
	int defaultResetGap = GetCurrencyUnit()->GetSpotDays();
	indexStartDate.GapBusinessDay(defaultResetGap, const_cast<char *>(itsResetCal.c_str()));
	
	indexStartDateDesc << CC_NS(std,fixed) << indexStartDate.GetJulian();
	rowDescVec[IndexStartDate] = indexStartDateDesc.str();
	rowTypeVec[IndexStartDate] = ARM_DATE_TYPE;

	// IT description
	CC_Ostringstream itDesc;
	itDesc << CC_NS(std,fixed) << CC_NS(std,setprecision)(8) << (*(datesStructure.GetDateStrip(LOANPAY_SCHED)->GetInterestTerms()))[eventIdx];
    rowDescVec[IT] = itDesc.str();
    rowTypeVec[IT] = ARM_DOUBLE_TYPE;

	// Coeff description
	CC_Ostringstream itCoeff;
	itCoeff << CC_NS(std,fixed) << CC_NS(std,setprecision)(8) << itsCoeff;
    rowDescVec[Coeff] = itCoeff.str();
    rowTypeVec[Coeff] = ARM_DOUBLE_TYPE;

	// DF description
	CC_Ostringstream dfDesc;
	dfDesc << "DF(" << modelName << ",";
	dfDesc << MaturityCapColNamesTable[PaymentDate] << "[i])";
    rowDescVec[DF] = dfDesc.str();
    rowTypeVec[DF] = ARM_STRING_TYPE;

	/// Annuity description
	CC_Ostringstream annuityDesc;
	annuityDesc << CC_NS(std,fixed) << itsAnnuity;
	rowDescVec[Annuity] = annuityDesc.str();
	rowTypeVec[Annuity] = ARM_DOUBLE;

	/// INFINE Nominal for pricing description
	CC_Ostringstream INFINEnominalDesc;

	if (isFirstEvent)
	{
		INFINEnominalDesc << "(" << CC_NS(std,fixed) << itsInitNominal << ")";
	}
	else
	{
		INFINEnominalDesc << "MAX(";
		INFINEnominalDesc << MaturityCapColNamesTable[INFINENominal] << "[i-1]";
		INFINEnominalDesc << "*(1+" << MaturityCapColNamesTable[INFINEInterest] << "[i-1])-";
		INFINEnominalDesc << MaturityCapColNamesTable[Annuity] << "[i-1]";
		INFINEnominalDesc << ",0)";
	}

	rowDescVec[INFINENominal] = INFINEnominalDesc.str();
	rowTypeVec[INFINENominal] = ARM_STRING;

	/// INFINE Nominal for calibration description
	CC_Ostringstream INFINEnominalForCalibDesc;

	if (isFirstEvent)
	{
		INFINEnominalForCalibDesc << "(" << CC_NS(std,fixed) << itsInitNominal << ")";
	}
	else
	{
		INFINEnominalForCalibDesc << "MAX(";
		INFINEnominalForCalibDesc << MaturityCapColNamesTable[INFINENominalForCalib] << "[i-1]";
		INFINEnominalForCalibDesc << "*(1+" << MaturityCapColNamesTable[INFINEInterest] << "[i-1])-";
		INFINEnominalForCalibDesc << MaturityCapColNamesTable[Annuity] << "[i-1]";
		INFINEnominalForCalibDesc << ",";
		INFINEnominalForCalibDesc << "YieldToPrice(" << MaturityCapColNamesTable[StartDate] << "[i],";
		INFINEnominalForCalibDesc << MaturityCapColNamesTable[FirstCouponDate] << "[i],";
		INFINEnominalForCalibDesc << MaturityCapColNamesTable[UnderlyingEndDate] << "[i],";
		INFINEnominalForCalibDesc << MaturityCapColNamesTable[Annuity] << "[i]*";
		INFINEnominalForCalibDesc << itsPayFreq << ",0,";
		INFINEnominalForCalibDesc << CC_NS(std,fixed) << CC_NS(std,setprecision)(8) << MaxYield << ",";
		INFINEnominalForCalibDesc << payFreq << ",";
		INFINEnominalForCalibDesc << "30/360))";

	}

	rowDescVec[INFINENominalForCalib] = INFINEnominalForCalibDesc.str();
	rowTypeVec[INFINENominalForCalib] = ARM_STRING;

	CC_Ostringstream RefStdCapnominalDesc;

	if (isFirstEvent)
	{
		RefStdCapnominalDesc << "(" << CC_NS(std,fixed) << itsInitNominal << ")";
	}
	else
	{
		RefStdCapnominalDesc << "MAX(";
		RefStdCapnominalDesc << MaturityCapColNamesTable[RefStdCapNominal] << "[i-1]";
		RefStdCapnominalDesc << "*(1+" << MaturityCapColNamesTable[Yield0] << "[i-1]*";
		RefStdCapnominalDesc << MaturityCapColNamesTable[Coeff] << "[i])-";
		RefStdCapnominalDesc << MaturityCapColNamesTable[Annuity] << "[i-1]";
		RefStdCapnominalDesc << ",0)";
	}

	rowDescVec[RefStdCapNominal] = RefStdCapnominalDesc.str();
	rowTypeVec[RefStdCapNominal] = ARM_STRING;

	/// INFINE Nominal Final description
	CC_Ostringstream INFINEnominalFinalDesc;
	if (isLastEvent)
	{
		INFINEnominalFinalDesc << MaturityCapColNamesTable[INFINENominalForCalib] << "[i]";
		INFINEnominalFinalDesc << "*(1+" << MaturityCapColNamesTable[INFINEInterest] << "[i])-";
		INFINEnominalFinalDesc << MaturityCapColNamesTable[Annuity] << "[i]";
		rowDescVec[INFINENominalFinal] = INFINEnominalFinalDesc.str();
		rowTypeVec[INFINENominalFinal] = ARM_STRING;
	}

	/// INFINE Yield description
	CC_Ostringstream INFINEYieldDesc;
	if (isFirstEvent)
	{
		INFINEYieldDesc << CC_NS(std,fixed) << itsInitTRI;

		rowDescVec[INFINEYield] = INFINEYieldDesc.str();
		rowTypeVec[INFINEYield] = ARM_DOUBLE;
	}
	else
	{
		INFINEYieldDesc << "PriceToYield(" << MaturityCapColNamesTable[StartDate] << "[i],";
		INFINEYieldDesc << MaturityCapColNamesTable[FirstCouponDate] << "[i],";
		INFINEYieldDesc << MaturityCapColNamesTable[UnderlyingEndDate] << "[i],";
		INFINEYieldDesc << MaturityCapColNamesTable[Annuity] << "[i]*";
		INFINEYieldDesc << itsPayFreq << ",0,";
		INFINEYieldDesc << MaturityCapColNamesTable[INFINENominalForCalib] << "[i],";
		INFINEYieldDesc << payFreq << ",";
		INFINEYieldDesc << "30/360)";

		rowDescVec[INFINEYield] = INFINEYieldDesc.str();
		rowTypeVec[INFINEYield] = ARM_STRING;
	}

	/// Amortizing
	CC_Ostringstream amortizingDesc;
	
	double amortizing = const_cast< ARM_Curve& >(itsAmortizing).Interpolate(paymentDate-asOfDate);

	amortizingDesc << CC_NS(std,fixed) << CC_NS(std,setprecision)(8) << amortizing;

	rowDescVec[Amortizing] = amortizingDesc.str();
	rowTypeVec[Amortizing] = ARM_STRING;

	/// Index description

	CC_Ostringstream indexDesc;
	if((*(datesStructure.GetDateStrip(LOANRESET_SCHED)->GetResetDates()))[eventIdx] != ARM_DateStripCombiner::DateStripCombiner_BlankData)
	{	
		indexDesc << "LIBOR(" << modelName << ",";
		indexDesc << MaturityCapColNamesTable[IndexStartDate] << "[i],";
		indexDesc << itsIndexTerm << ")";
	}
	else
	{
		indexDesc << MaturityCapColNamesTable[Index] << "[i-1]";
	}

	rowDescVec[Index] = indexDesc.str();
	rowTypeVec[Index] = ARM_STRING;

	/// Spread description
	CC_Ostringstream spreadDesc;
	spreadDesc << CC_NS(std,fixed) << itsSpread;

	rowDescVec[Spread] = spreadDesc.str();
	rowTypeVec[Spread] = ARM_STRING;

	/// Volatility description
	CC_Ostringstream volatilityStr;

	volatilityStr << 0.0;

	rowDescVec[Volatility] = volatilityStr.str();
	rowTypeVec[Volatility] = ARM_DOUBLE;

	/// Time description
	CC_Ostringstream timeStr;
	double time = (eventDate - asOfDate) / K_YEAR_LEN;

	timeStr << time;

	rowDescVec[Time] = timeStr.str();
	rowTypeVec[Time] = ARM_DOUBLE;

	/// Shock description
	CC_Ostringstream shockStr;
	shockStr << SHOCKCST;
	rowDescVec[Shock] = shockStr.str();
	rowTypeVec[Shock] = ARM_STRING;

	/// INFINE Interest description
	CC_Ostringstream INFINEInterestDesc;

	INFINEInterestDesc << "(" << MaturityCapColNamesTable[Index] << "[i]+";
	INFINEInterestDesc << MaturityCapColNamesTable[Spread] << "[i]+";
	INFINEInterestDesc << MaturityCapColNamesTable[Index] << "[i]*";
	INFINEInterestDesc << MaturityCapColNamesTable[Volatility] << "[i]*";
	INFINEInterestDesc << "Sqrt(" << MaturityCapColNamesTable[Time] << "[i])*";
	INFINEInterestDesc << MaturityCapColNamesTable[Shock] << "[i])";
	INFINEInterestDesc << "*" << MaturityCapColNamesTable[IT] << "[i]";

	rowDescVec[INFINEInterest] = INFINEInterestDesc.str();
	rowTypeVec[INFINEInterest] = ARM_STRING;

	/// Estim TRI Yield & Nominal description
	CC_Ostringstream EstimINFINEYieldDesc;
	CC_Ostringstream EstimINFINENominalDesc;

	if (isFirstEvent)
	{
		EstimINFINEYieldDesc << itsInitTRI;
		EstimINFINENominalDesc << itsInitNominal;
	}
	else
	{
		EstimINFINEYieldDesc << MaturityCapColNamesTable[INFINEYield] << "[i]" ;
		EstimINFINENominalDesc << MaturityCapColNamesTable[INFINENominal] << "[i]";
	}

	rowDescVec[EstimINFINEYield] = EstimINFINEYieldDesc.str();
	rowTypeVec[EstimINFINEYield] = ARM_STRING;

	rowDescVec[EstimINFINENominal] = EstimINFINENominalDesc.str();
	rowTypeVec[EstimINFINENominal] = ARM_STRING;

	/// INFINE Payoff description
	if (isLastEvent)
	{
		CC_Ostringstream INFINEMaturityCapDesc;

		if (itsLongShort == K_PAY)
		{
			INFINEMaturityCapDesc <<	"-";
		}

		INFINEMaturityCapDesc << "Max(" << MaturityCapColNamesTable[INFINENominal] << "[i],0)*";
		INFINEMaturityCapDesc << MaturityCapColNamesTable[DF] << "[i]";
		rowDescVec[INFINEMaturityCap] = INFINEMaturityCapDesc.str();
		rowTypeVec[INFINEMaturityCap] = ARM_STRING;
	}

	if (beforeEndDate)
	{
		/// TRI Yield description
		CC_Ostringstream YieldDesc;
		if (isFirstEvent)
		{
			YieldDesc << CC_NS(std,fixed) << itsInitTRI;

			rowDescVec[Yield] = YieldDesc.str();
			rowTypeVec[Yield] = ARM_DOUBLE;
		}
		else
		{
			YieldDesc << "PriceToYield(" << MaturityCapColNamesTable[StartDate] << "[i],";
			YieldDesc << MaturityCapColNamesTable[FirstCouponDate] << "[i],";
			YieldDesc << MaturityCapColNamesTable[UnderlyingEndDate] << "[i],";
			YieldDesc << MaturityCapColNamesTable[Annuity] << "[i]*";
			YieldDesc << itsPayFreq << ",0,";
			YieldDesc << MaturityCapColNamesTable[TRINominal] << "[i],";
			YieldDesc << payFreq << ",";
			YieldDesc << "30/360)";

			rowDescVec[Yield] = YieldDesc.str();
			rowTypeVec[Yield] = ARM_STRING;
		}

		/// Yield0 description
		CC_Ostringstream Yield0Desc;
		
		Yield0Desc << CC_NS(std,fixed) << itsInitTRI;

		rowDescVec[Yield0] = Yield0Desc.str();
		rowTypeVec[Yield0] = ARM_STRING;

		/// TRI Interest description
		CC_Ostringstream TRIInterestDesc;

		TRIInterestDesc << "Min((" << MaturityCapColNamesTable[Index] << "[i]+";
		TRIInterestDesc << MaturityCapColNamesTable[Spread] << "[i])*";
		TRIInterestDesc << MaturityCapColNamesTable[IT] << "[i],";
		TRIInterestDesc << MaturityCapColNamesTable[Yield] << "[i]*";
		TRIInterestDesc << MaturityCapColNamesTable[Coeff] <<"[i])";
		
		

		rowDescVec[TRIInterest] = TRIInterestDesc.str();
		rowTypeVec[TRIInterest] = ARM_STRING;

		/// TRI Nominal description
		CC_Ostringstream TRInominalDesc;

		TRInominalDesc << "Max(";

		if (isFirstEvent)
		{
			TRInominalDesc << CC_NS(std,fixed) << itsInitNominal;
		}
		else
		{
			TRInominalDesc << "MAX(";
			TRInominalDesc << MaturityCapColNamesTable[TRINominal] << "[i-1]";
			TRInominalDesc << "*(1+" << MaturityCapColNamesTable[TRIInterest] << "[i-1])-";
			TRInominalDesc << MaturityCapColNamesTable[Annuity] << "[i-1]";
			TRInominalDesc << ",";
			TRInominalDesc << "YieldToPrice(" << MaturityCapColNamesTable[StartDate] << "[i],";
			TRInominalDesc << MaturityCapColNamesTable[FirstCouponDate] << "[i],";
			TRInominalDesc << MaturityCapColNamesTable[UnderlyingEndDate] << "[i],";
			TRInominalDesc << MaturityCapColNamesTable[Annuity] << "[i]*";
			TRInominalDesc << itsPayFreq << ",0,";
			TRInominalDesc << CC_NS(std,fixed) << CC_NS(std,setprecision)(8) << MaxYield << ",";
			TRInominalDesc << payFreq << ",";
			TRInominalDesc << "30/360))";
		}

		TRInominalDesc << ",0)";

		rowDescVec[TRINominal] = TRInominalDesc.str();
		rowTypeVec[TRINominal] = ARM_STRING;

		/// TRI Payoff description
		CC_Ostringstream TRIMaturityCapDesc;

		if (itsLongShort == K_PAY)
		{
			TRIMaturityCapDesc <<	"-";
		}

		TRIMaturityCapDesc << "Max((" << MaturityCapColNamesTable[Index] << "[i]+";
		TRIMaturityCapDesc << MaturityCapColNamesTable[Spread] << "[i])*";
		TRIMaturityCapDesc << MaturityCapColNamesTable[IT] << "[i]-";
		TRIMaturityCapDesc << MaturityCapColNamesTable[Yield] << "[i]*";
		TRIMaturityCapDesc << MaturityCapColNamesTable[Coeff] << "[i]";
		TRIMaturityCapDesc <<  ",0)*";
		TRIMaturityCapDesc << MaturityCapColNamesTable[DF] << "[i]*";
		TRIMaturityCapDesc << MaturityCapColNamesTable[TRINominal] << "[i]*";
		TRIMaturityCapDesc << MaturityCapColNamesTable[Amortizing] << "[i]";

		rowDescVec[TRIMaturityCap] = TRIMaturityCapDesc.str();
		rowTypeVec[TRIMaturityCap] = ARM_STRING;

		/// RefStdCap description
		CC_Ostringstream RefStdCapDesc;

		if (itsLongShort == K_PAY)
		{
			RefStdCapDesc <<	"-";
		}

		RefStdCapDesc << "Max((" << MaturityCapColNamesTable[Index] << "[i]+";
		RefStdCapDesc << MaturityCapColNamesTable[Spread] << "[i])*";
		RefStdCapDesc << MaturityCapColNamesTable[IT] << "[i]-";
		RefStdCapDesc << MaturityCapColNamesTable[Yield0] << "[i]*";
		RefStdCapDesc << MaturityCapColNamesTable[Coeff] << "[i]";
		RefStdCapDesc << ",0)*";
		RefStdCapDesc << MaturityCapColNamesTable[DF] << "[i]*";
		RefStdCapDesc << MaturityCapColNamesTable[RefStdCapNominal] << "[i]*";
		RefStdCapDesc << MaturityCapColNamesTable[Amortizing] << "[i]";

		rowDescVec[RefStdCap] = RefStdCapDesc.str();
		rowTypeVec[RefStdCap] = ARM_STRING;

		/// Estim TRI Yield & Nominal description
		CC_Ostringstream EstimTRIYieldDesc;
		CC_Ostringstream EstimTRINominalDesc;

		if (isFirstEvent)
		{
			EstimTRIYieldDesc << itsInitTRI;
			EstimTRINominalDesc << itsInitNominal;
		}
		else
		{
			EstimTRIYieldDesc << MaturityCapColNamesTable[Yield] << "[i]" ;
			EstimTRINominalDesc << MaturityCapColNamesTable[TRINominal] << "[i-1]";
		}

		rowDescVec[EstimTRIYield] = EstimTRIYieldDesc.str();
		rowTypeVec[EstimTRIYield] = ARM_STRING;

		rowDescVec[EstimTRINominal] = EstimTRINominalDesc.str();
		rowTypeVec[EstimTRINominal] = ARM_STRING;
	}

    return ARM_RowInfo(rowDescVec,rowTypeVec);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_MaturityCapCalculator
///	Routine: InitPriceableColumns
///	Returns: nothing
///	Action : Feed the MDM with input market objects
/////////////////////////////////////////////////////////////////
void ARM_MaturityCapCalculator::InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const
{
    string zeroValue("0");

	rowDescVec[TRINominal] = zeroValue;
    rowTypeVec[TRINominal] = ARM_DOUBLE;

	rowDescVec[INFINENominal] = zeroValue;
    rowTypeVec[INFINENominal] = ARM_DOUBLE;

	rowDescVec[INFINENominalFinal] = zeroValue;
    rowTypeVec[INFINENominalFinal] = ARM_DOUBLE;

	rowDescVec[TRIMaturityCap] = zeroValue;
    rowTypeVec[TRIMaturityCap] = ARM_DOUBLE;

	rowDescVec[INFINEMaturityCap] = zeroValue;
    rowTypeVec[INFINEMaturityCap] = ARM_DOUBLE;

	rowDescVec[EstimTRIYield] = zeroValue;
    rowTypeVec[EstimTRIYield] = ARM_DOUBLE;

	rowDescVec[EstimTRINominal] = zeroValue;
    rowTypeVec[EstimTRINominal] = ARM_DOUBLE;

	rowDescVec[EstimINFINEYield] = zeroValue;
    rowTypeVec[EstimINFINEYield] = ARM_DOUBLE;

	rowDescVec[EstimINFINENominal] = zeroValue;
    rowTypeVec[EstimINFINENominal] = ARM_DOUBLE;

	rowDescVec[RefStdCap] = zeroValue;
    rowTypeVec[RefStdCap] = ARM_DOUBLE;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MaturityCapCalculator
///	Routine: CustomDatesStructure
///	Returns: ARM_DateStripCombiner
///	Action : create the list of all event dates of the Maturity Cap. The
///          DateStripCombiner merges event dates of each
///          legs
///			customized Schedule Not implemented yet.
/////////////////////////////////////////////////////////////////
ARM_DateStripCombiner ARM_MaturityCapCalculator::CustomDatesStructure(const ARM_DateStripVector& dateStrips) const
{
	return ARM_DateStripCombiner();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_MaturityCapCalculator
///	Routine: DatesStructure
///	Returns: ARM_DateStripCombiner
///	Action : create the list of all event dates of the Maturity Cap. The
///          DateStripCombiner merges event dates of each
///          legs
/////////////////////////////////////////////////////////////////
ARM_DateStripCombiner ARM_MaturityCapCalculator::DatesStructure() const
{

	/// The MaturityCap is mono-currency at the moment : the currency
    /// is then the payement one & it is saved at the ARM_Security level

    /// Get reset & payment calendars
    const char* resetCalendar   = itsResetCal.c_str();
    const char* payCalendar     = itsPayCal.c_str();

    int fwdRule		=   K_MOD_FOLLOWING;	// for forward dates
    int stubRule	=   K_SHORTSTART;

    int resetTiming =   K_ADVANCE;
	int payGap      =   GETDEFAULTVALUE;
    int payTiming   =   K_ARREARS;


    /// Build the Loan schedule
    ARM_DateStrip LoanResetLegSched(
			itsStartDate,
			itsUnderlyingEndDate,
			itsResetFreq,
			itsDayCount,
            resetCalendar,
            fwdRule,
			itsIntRule,
			stubRule,
            itsResetGap,
            itsResetFreq,
			payGap,
            payCalendar,
            resetTiming,
			payTiming);

	/// Build the Loan schedule
    ARM_DateStrip LoanPayLegSched(
			itsStartDate,
			itsUnderlyingEndDate,
			itsPayFreq,
			itsDayCount,
            resetCalendar,
            fwdRule,
			itsIntRule,
			stubRule,
            itsResetGap,
            itsPayFreq,
			payGap,
            payCalendar,
            resetTiming,
			payTiming);

    /// Merge schedules on "ResetDate"
    ARM_DateStripVector SchedVect(NB_MATCAP_SCHED,NULL);
    SchedVect[LOANRESET_SCHED]     = &LoanResetLegSched;
	SchedVect[LOANPAY_SCHED]     = &LoanPayLegSched;

    ARM_DateStripCombiner EventSchedule(SchedVect,"ResetDate");

    /// Initialise the memorisation of the first call line
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	ARM_VectorPtr eventDates = EventSchedule.GetMergeData();

    return EventSchedule;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_MaturityCapCalculator
///	Routine: GetIndexType
///	Returns: the index type
///	Action : get the ARM index type of the coupon index
/////////////////////////////////////////////////////////////////
ARM_INDEX_TYPE ARM_MaturityCapCalculator::GetIndexType()
{
    string liborTypeName(string(GetCurrencyUnit()->GetCcyName()) == string("EUR") ? "EURIBOR" : "LIBOR");
    string indexTerm(itsIndexTerm);
    if(indexTerm=="12M")
        indexTerm="1Y"; // because we cant convert EURIBOR1Y and not EURIBOR12M
    liborTypeName += indexTerm;
    return static_cast< ARM_INDEX_TYPE > (ARM_ArgConv_IndexType.GetNumber(liborTypeName));
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_MaturityCapCalculator
///	Routine: CreateAndSetModel
///	Returns: void
///	Action : creates the model
/////////////////////////////////////////////////////////////////
void ARM_MaturityCapCalculator::CreateAndSetModel()
{
	/// Get yield curve
	ARM_ZeroCurve* curve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));

	const ARM_DateStripCombiner& dateStructure = DatesStructure();

	ARM_GP_Vector* resetDates = dateStructure.GetDateStrip(0)->GetResetDates();

	/// Creates default values for volatility & mean reversion
    /// (calibration will be called at pricing time) and set them
    /// in the reference model
	ARM_GP_Vector defaultTimes;
	ARM_GP_Vector defaultSigmas;

	/// Get asOfDate
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

	ARM_CurveModelParam volParam( ARM_ModelParamType::Volatility, &defaultSigmas, &defaultTimes );
	ARM_CurveModelParam* mrsParam = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(GetKeys()[MrsKey]));
    if(!mrsParam || mrsParam->GetType() != ARM_ModelParamType::MeanReversion)
    {
		CC_Ostringstream os;
		os << ARM_USERNAME << " : MRS Param for key=" << GetKeys()[MrsKey] << " is expected in the Market Data Manager";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
    }
	ARM_ModelParam* betaParam = NULL;	
	betaParam = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[BetaKey]));
	if(!betaParam || betaParam->GetType() != ARM_ModelParamType::Beta)
	{
		CC_Ostringstream os;
		os << ARM_USERNAME << " : Beta Param for key=" << GetKeys()[BetaKey] << " is expected in the Market Data Manager";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}
	ARM_ModelParam* correlParam = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[CorrelKey]));
    if(!correlParam || correlParam->GetType() != ARM_ModelParamType::Correlation)
    {
		CC_Ostringstream os;
		os << ARM_USERNAME << " : Correl Param for key=" << GetKeys()[CorrelKey] << " is expected in the Market Data Manager";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
    }

	ARM_ModelParamVector paramVector(4);
	paramVector[0] = &volParam;
	paramVector[1] = mrsParam;
	paramVector[2] = betaParam;
	paramVector[3] = correlParam;

	ARM_INDEX_TYPE liborType = GetIndexType();
	ARM_IRIndex IRIndex(liborType,K_DEF_FREQ,K_DEF_FREQ,GetCurrencyUnit());

	/// use the factory class to get the model params!
	ARM_ModelParamsSFRM* pSFRMModelParams = 
		ARM_ModelParamsSFRMFactory.Instance()->CreateModelParamsSFRM(paramVector,&IRIndex,SFRM_NB_FACTORS,SFRM_VOL_TYPE);

	/// Build the default stochastic model of the calculator : SFRM 2F
    ARM_PricingModelPtr model( ARM_PricingModelPtr(static_cast< ARM_PricingModel* >(new ARM_SFRM( CreateClonedPtr( curve ), *pSFRMModelParams))) );

	// Delte the model params because it is cloned in the model
	delete pSFRMModelParams;

	// MRGK5 Random Generator
	ARM_RandomGeneratorPtr  pBaseRandomGen( ARM_RandGenFactory.Instance()->CreateRandGen( 
				ARM_RandGenFactoryImp::MerseneStd,
				ARM_RandGenFactoryImp::UnknownTransformAlgo ) );

	/// antithetic box muller
	ARM_RandomGeneratorPtr normRandGen( ARM_RandGenFactory.Instance()->CreateRandGen( 
		ARM_RandGenFactoryImp::UnknownBaseGenAlgorithm,
		ARM_RandGenFactoryImp::BoxMuller,
		pBaseRandomGen ) );

	/// antithetic variates!
	ARM_RandomGeneratorPtr numRandGen = ARM_RandomGeneratorPtr( ARM_RandGenFactory.Instance()->CreateRandGen( 
		ARM_RandGenFactoryImp::UnknownBaseGenAlgorithm,
		ARM_RandGenFactoryImp::AntitheticOne,
		normRandGen ) );

	ARM_TimeStepPerYearScheduler scheduler(MC_NB_INTER_STEPS);
	ARM_NormalCentredSamplerND sampler(&scheduler);

    ARM_MCMethod* mcMethod = new ARM_MCMethod(itsNbIterations,numRandGen,&sampler);

    /// Create a Numeraire and set it
    ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::RollingEvent ) );
    model->SetNumeraire(numeraire);
	model->SetNumMethod(ARM_NumMethodPtr( mcMethod ) );

	/// Set the model
	SetPricingModel(model);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_MaturityCapCalculator
///	Routine: CreateAndSetCalibration
///	Returns: void
///	Action : create the calibration
/////////////////////////////////////////////////////////////////
void ARM_MaturityCapCalculator::CreateAndSetCalibration()
{
    bool isFreezeWeights = false;					// reinit non null weights in the portfolio
    bool isInitSigma = true;						// init sigma param for calibration (bounds & init guess)

	/// Create cap floor portfolio (sigma calibration)
    ARM_StdPortfolioPtr cfPortfolio = CreateCFPortfolio();

	ARM_CalibMethod* cfCalibMethod = new ARM_CalibMethod(
		cfPortfolio,
		ARM_ModelParamVector(),
		ARM_CalibMethodType::Bootstrap1D,
		ARM_MAX_ITER,
		ARM_CalibrationTarget::PriceTarget,
		NULL,
		NULL);
	
	SetCalibMethod( ARM_CalibMethodPtr( cfCalibMethod ) );

	if (itsCalibrationMode == EX_BOUNDARY)
	{
		ComputeCFStrikesEX_BOUNDARY();
	}
	else
	{
		ComputeCFStrikesATMOrFLAT();
	}

	ComputeCFPrices(isFreezeWeights,isInitSigma);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_MaturityCapCalculator
///	Routine: CreateCFPortfolio
///	Returns: nothing
///	Action : create the list of caplets
/////////////////////////////////////////////////////////////////
ARM_StdPortfolioPtr ARM_MaturityCapCalculator::CreateCFPortfolio()
{
	/// Get the B&S model and vol curves from the market data manager
    ARM_BSModel* cfBSModel = dynamic_cast< ARM_BSModel* >(GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );

    const ARM_DealDescription dealDesc = GetGenSecurity()->GetDealDescription();
    size_t nbEvents=dealDesc.GetRowsNb(); // + 1 because of column names

	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();

	int indexFreq=ARM_ArgConv_MatFrequency.GetNumber(itsIndexTerm);

	int step = 1;

	if (indexFreq < itsPayFreq)
	{
		step = itsPayFreq/indexFreq;
	}

    /// Build a caplet/floorlet security using the coupon currency
    ARM_INDEX_TYPE liborType = GetIndexType();

	vector<ARM_CapFloor*> assets;
    for(size_t eventIdx=1;eventIdx<nbEvents;eventIdx+=step)
    {
		ARM_Date endDate( atof(dealDesc.GetElem(eventIdx,EndDate).c_str()) );

		if( (itsProductMode==TRI) && (endDate <= itsEndDate.GetJulian() + DAY_TOLERANCE) || (itsProductMode==INFINE))
		{
			ARM_Date startDate( atof(dealDesc.GetElem(eventIdx,StartDate).c_str()) );

			 /// Compute not adjusted index end date to avoid problem in caplet building
			ARM_Date endDateIndex(startDate);
			endDateIndex.AddPeriod(itsIndexTerm,itsPayCal.c_str());

			long dayCount = GetCurrencyUnit()->GetLiborIndexDayCount();

			double fwdRate;
			fwdRate = cfBSModel->ExpectedFwdYield(startDate.GetJulian(),endDateIndex.GetJulian(),endDate.GetJulian(),dayCount);

			/// Build the cap (ARM_Security & VanillaArg versions)
			double strike = fwdRate;	

			ARM_CapFloor capFloor(
				startDate,
				endDateIndex,
				K_CAP,
				strike,
				liborType,
				0.0,
				K_DEF_FREQ,
				K_DEF_FREQ,
				GetCurrencyUnit());

			double expiry = (capFloor.GetExpiryDate().GetJulian()-asOfDate.GetJulian());

			if (fabs(expiry) > K_DOUBLE_TOL)
			{
				assets.push_back(static_cast<ARM_CapFloor*>(capFloor.Clone()));
			}
		}
    }

	ARM_StdPortfolioPtr capFloorPortfolio(new ARM_StdPortfolio(assets.size()));

	size_t i;
	for (i = 0; i < assets.size(); ++i)
	{
		capFloorPortfolio->SetAsset(assets[i], i);
		capFloorPortfolio->SetWeight(CF_DEFAULT_WEIGHT,i);
		capFloorPortfolio->SetPrice(eventIdx*CF_DEFAULT_PRICE,i);
	}

	ARM_CurveModelParam modelparamVol(ARM_ModelParamType::Volatility, SIGMA_DEFAULT_VALUE,"SIGMA",SIGMA_LOWER_BOUND,SIGMA_UPPER_BOUND );
	GetPricingModel()->GetModelParams()->SetModelParam(&modelparamVol);

	return capFloorPortfolio;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MaturityCapCalculator
///	Routine: Calibrate
///	Returns: void
///	Action : Calibrate the model of the calculator.
/////////////////////////////////////////////////////////////////
void ARM_MaturityCapCalculator::Calibrate()
{
	GetCalibMethod()->Calibrate(&*GetPricingModel());
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_MaturityCapCalculator
///	Routine: Price
///	Returns: a double
///	Action : price the Maturity Cap deal.
/////////////////////////////////////////////////////////////////
double ARM_MaturityCapCalculator::Price()
{
	CalibrateAndTimeIt();

    /// Price the implicit product according to internal flag
	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();

	CC_Ostringstream evalDateDesc;
    evalDateDesc << CC_NS(std,fixed) << asOfDate.GetJulian();
    string evalDateStr(evalDateDesc.str());

	ARM_ZeroCurve* zeroCurve = (ARM_ZeroCurve*) GetMktDataManager()->GetData(GetKeys()[YcKey]);

    for (size_t i = 0; i < NbProductsToPrice; ++i)
		if (itsProductsToPrice[i])
		{
			/// Select the right column that describes the product
			MaturityCapColAlias prodIdx;

			if(i == MaturityCapPrice)
			{
				if (itsProductMode == TRI)
					prodIdx = TRIMaturityCap;
				else if (itsProductMode == INFINE)
					prodIdx = INFINEMaturityCap;
			}
			else if(i == RefStdCapPrice)
			{
				prodIdx = RefStdCap;
			}
			else if(i == EstimatedTRI)
			{
				if (itsProductMode == TRI)
					prodIdx = EstimTRIYield;
				else if (itsProductMode == INFINE)
					prodIdx = EstimINFINEYield;
			}
			else if(i == EstimatedNominal)
			{
				if (itsProductMode == TRI)
					prodIdx = EstimTRINominal;
				else if (itsProductMode == INFINE)
					prodIdx = EstimINFINENominal;
			}


			ARM_DealDescription dealDesc = GetGenSecurity()->GetDealDescription();

			/// Build the sub-generic security
			size_t startRowIdx = 1; /// start with the first actual row (=> skip column names)
			size_t endRowIdx = dealDesc.GetRowsNb() - 1;
			ARM_DealDescriptionPtr subDealDesc = dealDesc.GetSubDescription(startRowIdx,endRowIdx,prodIdx+1);
			vector<string> names(1,SHOCKCST);
			vector<double> values(1,0.0);
			ARM_CstManagerPtr cstManager(new ARM_CstManager(names,values));
			ARM_GenSecurityPtr subGenSec = ARM_GenSecurityPtr(new ARM_GenSecurity(subDealDesc,"",cstManager));

			ARM_GenPricer* genPricer = new ARM_GenPricer( &*subGenSec,&*GetPricingModel());

			ARM_AutoCleaner<ARM_GenPricer> HoldGP(genPricer );

			double price = genPricer->Price();
			double stdDev = 0.0;
            if (genPricer->GetPricerInfo()->GetContents().IsDataExist("StdDev"))
                stdDev = genPricer->GetPricerInfo()->GetContents().GetData("StdDev").GetDouble();
			double duration = genPricer->GetPricerInfo()->GetDuration();

			if(i == MaturityCapPrice)
			{
				itsMaturityCapPrice = price;
				itsMaturityCapStdDev = stdDev;
				itsMaturityCapPricingDuration = duration;
			}
			else if(i == RefStdCapPrice)
			{
				itsRefStdCapPrice = price;
				itsRefStdCapStdDev = stdDev;
			}
			else if((i == EstimatedTRI) || (i == EstimatedNominal))
			{
				ARM_GramFctorArg gramFunctorArg = genPricer->GetPricerInfo()->GetContents().GetData("IntermediatePrices");
				ARM_VectorPtr itermPrices = gramFunctorArg.GetVector();

				for (size_t j=0; j <itermPrices->size(); ++j)
				{
					double df = 1.0;
					double eventDate = atof( dealDesc.GetElem(j+1,EventDate).c_str() );
					df = zeroCurve->DiscountPrice((eventDate-asOfDate.GetJulian())/K_YEAR_LEN);
					(*itermPrices)[j] = (*itermPrices)[j]/df;
				}
				
				if (i == EstimatedTRI)
				{
					itsEstimatedTRI = itermPrices;
				}
				else if (i == EstimatedNominal)
				{
					itsEstimatedNominal = itermPrices;
				}
			}
		}

	itsHasBeenPriced = true;
    
    return itsMaturityCapPrice;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_MaturityCapCalculator
///	Routine: ComputePricingData
///	Returns: a ARM_MultiTypeDict
///	Action : get pricing data of the calculaltor
/////////////////////////////////////////////////////////////////

void ARM_MaturityCapCalculator::ComputePricingData() const
{
	if (!itsHasBeenPriced)
		const_cast<ARM_MaturityCapCalculator*>(this)->PriceAndTimeIt();

	if (itsProductsToPrice[MaturityCapPrice])
	{
		GetPricingData()[ "MaturityCapPrice"		] = itsMaturityCapPrice;
		GetPricingData()[ "MaturityCapStdDev"		] = itsMaturityCapStdDev;
		GetPricingData()[ "MaturityCapPricingDuration"	] = itsMaturityCapPricingDuration;
	}
	if (itsProductsToPrice[RefStdCapPrice])
	{
		GetPricingData()[ "RefStdCapPrice"		] = itsRefStdCapPrice;
		GetPricingData()[ "RefStdCapStdDev"		] = itsRefStdCapStdDev;
	}
	if (itsProductsToPrice[EstimatedTRI])
	{
		GetPricingData()[ "EstimatedTRI"		] = itsEstimatedTRI;
	}
	if (itsProductsToPrice[EstimatedNominal])
	{
		GetPricingData()[ "EstimatedNominal"		] = itsEstimatedNominal;
	}
	if (itsCalibrationMode == EX_BOUNDARY)
	{
		GetPricingData()[ "Shock"		] = itsBoundaryShock;
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_MaturityCapCalculator
///	Routine: ComputeCFPrices
///	Returns: nothing
///	Action : compute market target prices of the C/F portfolio
/////////////////////////////////////////////////////////////////
void ARM_MaturityCapCalculator::ComputeCFPrices(
	bool isFreezeWeights, 
	bool isInitParam)
{
	/// Get the B&S model and vol curves from the market data manager
    ARM_BSModel* cfBSModel = dynamic_cast< ARM_BSModel* >(GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );

	/// Get asOfDate
	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();

	ARM_CalibMethod* cfCalibMethod = GetCFCalibMethod();
	ARM_StdPortfolioPtr port = GetCFPortfolio();

	ARM_GP_Vector initTimes;
	ARM_GP_Vector initSigmas;
	ARM_GP_Vector sigmaLowerBound;
	ARM_GP_Vector sigmaUpperBound;

	for (size_t i = 0; i < port->GetSize(); ++i)
	{
		ARM_CapFloor* capFloor = static_cast<ARM_CapFloor*>(port->GetAsset(i));
		capFloor->SetModel(cfBSModel);
		double price = capFloor->ComputePrice();
		double vega = 0.0;
		if (!isFreezeWeights)
			vega= capFloor->ComputeSensitivity(K_VEGA);

		double fwdRate;
		if (capFloor->GetSwapLeg()->GetFwdRates() 
				&& capFloor->GetSwapLeg()->GetFwdRates()->GetSize() >= 1)
		{
			fwdRate = capFloor->GetSwapLeg()->GetFwdRates()->Elt(0);
		}
		else
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : cap floor fwd rates are not available !" );
		}

		double expiry = (capFloor->GetExpiryDate().GetJulian()-asOfDate.GetJulian());
		double tenor = capFloor->GetSwapLeg()->GetIRIndex()->GetYearTerm();

		double vol = cfBSModel->ComputeVol(expiry/K_YEAR_LEN,tenor,fwdRate,fwdRate)/100.0;

		if (isInitParam)
		{
			initTimes.push_back(expiry);
			initSigmas.push_back(vol);
			sigmaLowerBound.push_back(SIGMA_LOWER_BOUND);
			sigmaUpperBound.push_back(SIGMA_UPPER_BOUND);
		}

		port->SetPrice(price,i);
		if (!isFreezeWeights)
			port->SetPrecision(vega,i);
	}

	if (isInitParam)
	{
		/// Test if volatility must be initialise. In this case,
		/// its schedule contents the expiry date of each swaption product
		size_t sigmaIdx,paramSize=cfCalibMethod->GetCalibParams().size();
		for(sigmaIdx=0;sigmaIdx<paramSize;++sigmaIdx)
			if(cfCalibMethod->GetCalibParam(sigmaIdx) &&
			(cfCalibMethod->GetCalibParams())[sigmaIdx]->GetType() == ARM_ModelParamType::Volatility)
				break;

		ARM_CurveModelParam* sigma = new ARM_CurveModelParam(
				ARM_ModelParamType::Volatility,
				&initSigmas,
				&initTimes,
				"SIGMA",
				"STEPUPRIGHT",
				&sigmaLowerBound,
				&sigmaUpperBound);

		if(sigmaIdx >= paramSize || paramSize == 0)
			cfCalibMethod->GetCalibParams().push_back(sigma);
		else
		{
			delete cfCalibMethod->GetCalibParam(sigmaIdx);
			(cfCalibMethod->GetCalibParams())[sigmaIdx] = sigma;
		}
	}
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MaturityCapCalculator
///	Routine: ComputeCFStrikes
///	Returns: nothing
///	Action : compute C/F Portfolios strikes for ATM and FLAT modes
/////////////////////////////////////////////////////////////////
void ARM_MaturityCapCalculator::ComputeCFStrikesATMOrFLAT()
{
	/// Get the B&S model and vol curves from the market data manager
    ARM_BSModel* cfBSModel = dynamic_cast< ARM_BSModel* >(GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );

	ARM_StdPortfolioPtr port = GetCFPortfolio();

	for (size_t i = 0; i < port->GetSize(); ++i)
	{
		ARM_CapFloor* capFloor = static_cast<ARM_CapFloor*>(port->GetAsset(i));
		capFloor->SetModel(cfBSModel);

		if (itsCalibrationMode == ATM)
		{
			double fwdRate;
			if (capFloor->GetSwapLeg()->GetFwdRates() 
					&& capFloor->GetSwapLeg()->GetFwdRates()->GetSize() >= 1)
			{
				fwdRate = capFloor->GetSwapLeg()->GetFwdRates()->Elt(0);
			}
			else
			{
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : cap floor fwd rates are not available !" );
			}
			capFloor->SetStrike(fwdRate);
		}
		else if (itsCalibrationMode == FLAT)
		{
			double strikeAdjust = 1.0;

			if (itsDayCount == KACTUAL_360)
				strikeAdjust = 360.0/365;

			capFloor->SetStrike((itsInitTRI*strikeAdjust-itsSpread)*100);
		}
	}
}

class InfineMaturityCapToInverse : public ARM_GP::UnaryFunc<double,double> 
{
public: 
		InfineMaturityCapToInverse(
			ARM_GenPricer* matCapGenpricer,
			ARM_CstManager::iterator shockIt) :
		itsMatCapGenPricer(matCapGenpricer),
		itsShockIt(shockIt)
		{
		};

		virtual double operator() (double shock) const 
		{
			itsShockIt->second->SetDouble (shock) ;
			return itsMatCapGenPricer->Price();
		}
private:
	ARM_GenPricer* itsMatCapGenPricer;
	ARM_CstManager::iterator itsShockIt;
};

/////////////////////////////////////////////////////////////////
///	Class  : ARM_MaturityCapCalculator
///	Routine: ComputeCFStrikes
///	Returns: nothing
///	Action : compute C/F Portfolios strikes for EX_BOUNDARY mode
/////////////////////////////////////////////////////////////////
void ARM_MaturityCapCalculator::ComputeCFStrikesEX_BOUNDARY()
{
	/// Get the B&S model and vol curves from the market data manager
    ARM_BSModel* cfBSModel = dynamic_cast< ARM_BSModel* >(GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );

	// Zc Curve;
	ARM_ZeroCurve* zcCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));

	/// Get asOfDate
	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();

	/// Make a copy of the initial description to update strike spreads
    ARM_DealDescriptionPtr dealDesc( (ARM_DealDescription*) const_cast< ARM_DealDescription& >(GetGenSecurity()->GetDealDescription()).Clone() );

	size_t nbEvents=dealDesc->GetRowsNb();

	ARM_StdPortfolioPtr port = GetCFPortfolio();

	for (size_t eventIdx = 1; eventIdx < nbEvents; ++eventIdx)
	{
		ARM_Date resetDate( atof(dealDesc->GetElem(eventIdx,EventDate).c_str()) );
		ARM_Date startDate( atof(dealDesc->GetElem(eventIdx,StartDate).c_str()) );

         /// Compute not adjusted index end date to avoid problem in caplet building
        ARM_Date endDate(startDate);
        endDate.AddPeriod(itsIndexTerm,itsPayCal.c_str());

		long dayCount = GetCurrencyUnit()->GetLiborIndexDayCount();

		double fwdRate;
		fwdRate = cfBSModel->ExpectedFwdYield(startDate.GetJulian(),endDate.GetJulian(),endDate.GetJulian(),dayCount);

		double expiry = (resetDate.GetJulian()-asOfDate.GetJulian())/K_YEAR_LEN;
		double tenor = (endDate.GetJulian()-startDate.GetJulian())/K_YEAR_LEN;

		double vol = cfBSModel->ComputeVol(expiry,tenor,fwdRate,fwdRate)/100.0;

		CC_Ostringstream volStr;
		volStr << vol;

		dealDesc->SetElem(eventIdx,Volatility,volStr.str(),ARM_DOUBLE);
	}

	double minShock = MINSHOCK, maxShock = MAXSHOCK, shockInit;
	double minInFine, maxInFine;

	size_t startRowIdx = 1; /// start with the first actual row (=> skip column names)
	size_t endRowIdx = dealDesc->GetRowsNb() - 1;
	size_t prodIdx = INFINENominalFinal;

	ARM_DealDescriptionPtr subDealDesc = dealDesc->GetSubDescription(startRowIdx,endRowIdx,prodIdx+1);

	ARM_IrFwdMod irFwdMod( CreateClonedPtr( zcCurve) );
	ARM_FINumMethod* fiNumMethod = new ARM_FINumMethod;
	irFwdMod.SetNumMethod(ARM_NumMethodPtr( fiNumMethod ) );

	vector<string> names(1,SHOCKCST);
	vector<double> values(1,0.0);
	ARM_CstManagerPtr cstManager(new ARM_CstManager(names,values));
	ARM_CstManager::iterator shockIt= cstManager->find(SHOCKCST);
	ARM_GenSecurityPtr subGenSec = ARM_GenSecurityPtr(new ARM_GenSecurity(subDealDesc,"",cstManager));
	ARM_GenPricer* genPricer = new ARM_GenPricer( &*subGenSec,&irFwdMod);
	ARM_AutoCleaner<ARM_GenPricer> HoldGP(genPricer );

	do {
		shockIt->second->SetDouble(minShock);
		minInFine = genPricer->Price();

		minShock -= SHOCKINC;
	}
	while (minInFine > 0);

	do {
		shockIt->second->SetDouble(maxShock);
		maxInFine = genPricer->Price();

		maxShock += SHOCKINC;
	}
	while (maxInFine < 0);

	shockInit = (minShock+maxShock)/2;

	InfineMaturityCapToInverse func(genPricer, shockIt);
	UnaryFuncWithNumDerivative<double> funcWithDev(func);

	T_NewtonRaphsonSolver< UnaryFuncWithNumDerivative<double> > solver(funcWithDev,0,DEFAULT_PREC,DEFAULT_PREC);
	
	///To initialize departure point
	solver.setInitialGuess(shockInit);

	itsBoundaryShock = solver.Solve();

	prodIdx = EstimINFINEYield;
	subDealDesc = dealDesc->GetSubDescription(startRowIdx,endRowIdx,prodIdx+1);
	subGenSec = ARM_GenSecurityPtr(new ARM_GenSecurity(subDealDesc,"",cstManager));
	ARM_GenPricer* triGenPricer = new ARM_GenPricer( &*subGenSec,&irFwdMod);
	ARM_AutoCleaner<ARM_GenPricer> HoldTRIGP(triGenPricer );
	triGenPricer->Price();
	ARM_GramFctorArg gramFunctorArg = triGenPricer->GetPricerInfo()->GetContents().GetData("IntermediatePrices");
	ARM_VectorPtr itermPrices = gramFunctorArg.GetVector();

	double strikeAdjust = 1.0;

	if (itsDayCount == KACTUAL_360)
		strikeAdjust = 360.0/365;

	for (size_t i = 0; i < port->GetSize(); ++i)
	{
		ARM_CapFloor* capFloor = static_cast<ARM_CapFloor*>(port->GetAsset(i));

		double df = 1.0;
		double eventDate = atof( dealDesc->GetElem(i+1,EventDate).c_str() );

		df = zcCurve->DiscountPrice((eventDate-asOfDate.GetJulian())/K_YEAR_LEN);

		capFloor->SetStrike(((*itermPrices)[i]/df*strikeAdjust-itsSpread)*100);
	}	
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_MaturityCapCalculator
///	Routine: UpdateModel
///	Returns: void
///	Action : update the model datas w.r.t. MarketDataManager
/////////////////////////////////////////////////////////////////
void ARM_MaturityCapCalculator::UpdateModel()
{
	const ARM_DateStripCombiner& dateStructure = DatesStructure();

	ARM_GP_Vector* resetDates = dateStructure.GetDateStrip(0)->GetResetDates();

	ARM_ModelParam* mrsParam = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[MrsKey]));
    if(!mrsParam || mrsParam->GetType() != ARM_ModelParamType::MeanReversion)
    {
		CC_Ostringstream os;
		os << ARM_USERNAME << " : MRS Param for key=" << GetKeys()[MrsKey] << " is expected in the Market Data Manager";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
    }
	ARM_ModelParam* betaParam = NULL;
	betaParam = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[BetaKey]));
	if(!betaParam || betaParam->GetType() != ARM_ModelParamType::Beta)
	{
		CC_Ostringstream os;
		os << ARM_USERNAME << " : Beta Param for key=" << GetKeys()[BetaKey] << " is expected in the Market Data Manager";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}
	ARM_ModelParam* correlParam = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[CorrelKey]));
    if(!correlParam || correlParam->GetType() != ARM_ModelParamType::Correlation)
    {
		CC_Ostringstream os;
		os << ARM_USERNAME << " : Correl Param for key=" << GetKeys()[CorrelKey] << " is expected in the Market Data Manager";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
    }

	ARM_SFRM* model = static_cast<ARM_SFRM*>( &*GetPricingModel() );

	ARM_ModelParamVector paramVector(4);

	paramVector[0] = static_cast<ARM_ModelParam*>(model->GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).Clone());

    /// Update the MRS
	paramVector[1] = mrsParam;
	/// Update the Beta
	paramVector[2] = betaParam;
	/// Update the Correl
	paramVector[3] = correlParam;

	ARM_INDEX_TYPE liborType = GetIndexType();

	ARM_IRIndex IRIndex(liborType,K_DEF_FREQ,K_DEF_FREQ,GetCurrencyUnit());

	/// use the factory class to get the model params!
	ARM_ModelParamsSFRM* pSFRMModelParams = 
		ARM_ModelParamsSFRMFactory.Instance()->CreateModelParamsSFRM(paramVector,&IRIndex,SFRM_NB_FACTORS,SFRM_VOL_TYPE);

	model->SetModelParams(*pSFRMModelParams);
	delete pSFRMModelParams;

	ARM_ZeroCurve* zeroCurve     = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));

	model->SetZeroCurve( CreateClonedPtr( zeroCurve) );

	ARM_Portfolio* port=NULL;
	model->ConvertToShiftorBetaParam(*port);

}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_MaturityCapCalculator
///	Routine: UpdateCalibration
///	Returns: void
///	Action : update the calibration datas w.r.t. MarketDataManager
///          The context is an hedge ratios computation
/////////////////////////////////////////////////////////////////
void ARM_MaturityCapCalculator::UpdateCalibration(bool isUpdateStrike)
{
    /// Update cap floor prices
    bool isFreezeWeights = true; // keep the portfolio size to avoid hedge jumps
    bool isInitSigma = true;       // keep current sigma param to init calibration

	if (itsCalibrationMode == EX_BOUNDARY)
	{
		ComputeCFStrikesEX_BOUNDARY();
	}
	else
	{
		ComputeCFStrikesATMOrFLAT();
	}

	ComputeCFPrices(isFreezeWeights,isInitSigma);
}

////////////////////////////////////////////////////
///	Class   : ARM_MaturityCapCalculator
///	Routines: Clone,View
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_MaturityCapCalculator::Clone() const
{
	return new ARM_MaturityCapCalculator(*this);
}

void ARM_MaturityCapCalculator::View(char* id, FILE* ficOut) const
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

    /// Maturity Cap Calculator specific datas viewing
    fprintf(fOut,"\n\n =======> MATURITY CAP <====== \n");

	CC_Ostringstream mcapData;
    mcapData << "\nStartDate = " <<  itsStartDate.toString() << "\n";
    mcapData << "EndDate = " << itsEndDate.toString() << "\n";
	mcapData << "UnderlyingEndDate = " << itsUnderlyingEndDate.toString() << "\n";
    mcapData << "Long/Short = " << ARM_ParamView::GetMappingName(S_LONG_SHORT_TYPES, itsLongShort) << "\n";
	mcapData << "Cap/Floor = " << ARM_ParamView::GetMappingName(S_CAPLIKE_TYPE, itsCapFloor) << "\n";
	mcapData << "Reset Frequency = " << ARM_ParamView::GetMappingName(S_FREQUENCY, itsResetFreq ) << "\n";
	mcapData << "Payment Frequency = " << ARM_ParamView::GetMappingName(S_FREQUENCY, itsPayFreq ) << "\n";
	mcapData << "Index Term = " << itsIndexTerm << "\n";
    mcapData << "Day Count = " << ARM_ParamView::GetMappingName( S_DAYCOUNT,  itsDayCount ) << "\n";
	mcapData << "Int Rule = " << ARM_ParamView::GetMappingName( S_INTEREST_RULES,  itsIntRule ) << "\n";
	mcapData << "Init Nominal = " << itsInitNominal << "\n";
	mcapData << "Init TRI = " << itsInitTRI << "\n";
	mcapData << "Annuity = " << itsAnnuity << "\n";
	mcapData << "Product Mode = " << ARM_ParamView::GetMappingName( S_MCAP_MODE_TYPES,  itsProductMode ) << "\n";
	mcapData << "Coeff = " << itsCoeff << "\n";
	mcapData << "Reset Gap = " << itsResetGap << "\n";
	mcapData << "Reset Calendar = " << itsPayCal << "\n";
    mcapData << "Pay Calendar = " << itsPayCal << "\n";
	mcapData << "Calibration Mode = " << ARM_ParamView::GetMappingName( S_MCAP_CALIB_MODE_TYPES,  itsCalibrationMode ) << "\n";
	mcapData << "Nb Iterations = " << itsNbIterations << "\n";
   
	/// part common to gencalculator
	mcapData << "\n\nCommon Calculator Part\n" << ARM_GenCalculator::toString() << "\n\n";

    fprintf(fOut,"%s",mcapData.str().c_str());

	ARM_StdPortfolioPtr capFloorPF = const_cast<ARM_MaturityCapCalculator*>(this)->GetCFPortfolio();

	if (capFloorPF != ARM_StdPortfolioPtr(NULL))
		capFloorPF->View(id,fOut);
	
    string prodTxt("Product to be priced = ");
    if(IsMaturityCapToPrice())				prodTxt += "Maturity Cap";
	else if(IsRefStdCapToPrice())			prodTxt += "Ref Standard Cap";
	else if(IsEstimatedTRIToPrice())		prodTxt += "Estimated TRI";
	else if(IsEstimatedNominalToPrice())	prodTxt += "Estimated Nominal";
    else prodTxt += "Unknown !";
    fprintf(fOut,"\n\n %s\n",prodTxt.c_str());

    /// Common viewing
    ARM_GenCalculator::View(id,fOut);


    if ( ficOut == NULL )
       fclose(fOut);
}

CC_END_NAMESPACE()
