/*!
 /** Copyright (c) CDC IXIS CM July 2003 Paris
 *	\file bermudaswaptioncalculator.cpp
 *  \brief file for the Bermuda Swaption Calculator
 *	\author  H. BAKHTRI & A. TRIKI
 *	\version 1.0
 *	\date June 2005
 */

/// this header comes first as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalculators/bermudaswaptioncalculator.h"
#include "gpcalculators/meanrevfinder.h"
#include "gpcalculators/argconvdefault.h"

/// gpbase
#include "gpbase/env.h"
#include "gpbase/datestripcombiner.h"
#include "gpbase/autocleaner.h"
#include "gpbase/singleton.h"
#include "gpbase/datestrip.h"
#include "gpbase/stringconvert.h"
#include "gpbase/gptrigomatrix.h"
#include "gpbase/gpvector.h"
#include "gpbase/utilityport.h"

#include "gpbase/ostringstream.h"
#include "gpbase/curve.h"
#include "gpbase/gplinalgconvert.h"
#include "gpbase/numericconstant.h"
#include "gpbase/surface.h"
#include "gpbase/surfacetypedef.h"
#include "gpbase/curveconvert.h"

/// gpinfra
#include "gpinfra/argconvdefault.h"
#include "gpinfra/gensecurity.h"
#include "gpinfra/dealdescription.h"
#include "gpinfra/gensecmanipulator.h"
#include "gpinfra/mktdatamanagerrep.h"
#include "gpinfra/pricingadviser.h"
#include "gpinfra/cstmanager.h"
#include "gpinfra/pricingmodel.h"
#include "gpinfra/gramfunctorargdict.h"
#include "gpinfra/genpricer.h"
#include "gpinfra/pricerinfo.h"
#include "gpinfra/numeraire.h"
#include "gpinfra/numerairefactory.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/correlmatparam.h"

/// gpcalib
#include "gpcalib/modelparamsfactory.h"
#include "gpcalib/calibmethod.h"
#include "gpcalib/vanillacap.h"
#include "gpcalib/kerneltogp.h"
#include "gpcalib/modelfitter.h"
#include "gpcalib/modelfitterdes.h"
#include "gpcalib/vanillapricer.h"
#include "gpcalib/vanillaarg.h"
#include "gpcalib/basket.h"

/// gpmodels
#include "gpmodels/modelparamssfrm.h"
#include "gpmodels/modelparamssfrmdiag.h"
#include "gpmodels/modelparamssfrmfactory.h"
#include "gpmodels/sfrm.h"
#include "gpmodels/hw1f.h"
#include "gpmodels/modelparamshw1f.h"
#include "gpmodels/QGM1F.h"
#include "gpmodels/ModelParamsQGM1F.h"

/// gpnumlib
#include "gpnumlib/randomgenfactory.h"
#include "gpnumlib/compositegen.h"
#include "gpnumlib/antitheticgen.h"
#include "gpnumlib/argconvdefault.h"

#include "gpnumlib/solver.h"
#include "gpnumlib/newtonraphson.h"
#include "gpnumlib/numfunction.h"

/// gpnummethods
/// Tree
#include "gpnummethods/treemethod.h"
#include "gpnummethods/treebase.h"
#include "gpnummethods/treefactory.h"
#include "gpnummethods/pdemethod.h"
#include "gpnummethods/sampler.h"
#include "gpnummethods/normalcentredsampler.h"
#include "gpnummethods/scheduler.h"
#include "gpnummethods/truncator.h"
#include "gpnummethods/reconnector.h"
#include "gpnummethods/smoother.h"
#include "gpnummethods/argconvdefault.h"
#include "gpnummethods/pathscheme.h"
#include "gpnummethods/pathschemefactory.h"
/// AMC
#include "gpnummethods/amcmethod.h"
#include "gpnummethods/amc_andersen.h"
#include "gpnummethods/amc_ls.h"
#include "gpnummethods/cfmethod.h"


/// kernel
//#include <inst/swaption.h>
//#include <inst/fixleg.h>
//#include <inst/portfolio.h>
#include <glob/paramview.h>
//#include <mod/bssmiled.h>
//////#include <crv/volcube.h>

/// STL
#include <iomanip> /// for setprecision()
#include <memory>

CC_BEGIN_NAMESPACE( ARM )

/// Default MDM key names
const string YC_KEY_NAME            = "YC_";
const string OSWMODEL_KEY_NAME      = "OSWMOD_";
const string CFMODEL_KEY_NAME       = "CFMOD_";
const string BETA_KEY_NAME          = "BETA_";
const string NMOD_KEY_NAME			= "NMOD_";

/// Call schedule for Bermuda Swaption date structure
const unsigned int CALL_BS_SCHED			= 0;

/// SFRM sigma range [10bp,1000bp] with a 20% default value
const double SIGMA_LOWER_BOUND				= 0.001;
const double SIGMA_UPPER_BOUND				= 1.0;

///SFRM2F Correlation range [-100%,100%]
const double CORREL_LOWER_BOUND				= -1.0;
const double CORREL_UPPER_BOUND				= 1.0;

/// SFRM BETA range [1%,150%] with a 100% default value
const double BETA_LOWER_BOUND				= 0.01;
const double BETA_UPPER_BOUND				= 1.5;

/// Anual frequency for correlation
const unsigned int CORREL_FREQ				= 1;

/// QGM Calibration
const unsigned int QGM_VOL_CALIB_ITER		= 5;

const unsigned int NB_PROBA_COL				= 44;
const unsigned int NB_CALCULATED_PROBA      = 21;

//Newton raphson precision
const double DEF_PRECISION = 1e-4;

const double NON_CALL_FEE			= 1e15;
const double DEFAULT_PRICE			= 1.0e+100;
const double DEFAULT_PRECISION		= 1.0;
const double DEFAULT_WEIGHT			= 1.0;
const double VEGA_MIN				= 1.0e-6; 

const string ARM_BermudaSwaptionCalculator::ControlVariateColNamesOutputs [] =
{
	"ControlVariate1Price",
	"ControlVariate2Price",
	"ControlVariate3Price",
	"ControlVariate4Price",
	"ControlVariate5Price",
	"ControlVariate6Price",
	"ControlVariate7Price",
	"ControlVariate8Price",
	"ControlVariate9Price",
	"ControlVariate10Price",
};


//GenSec 1
const string ARM_BermudaSwaptionCalculator::BermudaSwaptionColNamesTable1 [] =
{
     "ResetDate",
	 "StartDate",
	 "PayDate",
	 "SwapStartDate",
	 "NextSwapStartDate",
	 "LastSwapStartDate",
	 "Strike",
	 "Fees",
	 "Notional",
	 "StdFixSwaplet",
	 "StdVarSwaplet",
	 "StdSwaplet",
	 "StdSwap",
	 "Option",
	 "Bermuda",
};

const string ARM_BermudaSwaptionCalculator::ControlVariateColNamesTable1 [] = 
{
	"ControlVariate1",
	"ControlVariate2",
	"ControlVariate3",
	"ControlVariate4",
	"ControlVariate5",
	"ControlVariate6",
	"ControlVariate7",
	"ControlVariate8",
	"ControlVariate9",
	"ControlVariate10",
};

//GenSec 2
const string ARM_BermudaSwaptionCalculator::BermudaSwaptionColNamesTable2 [] =
{
     "ResetDate",
	 "StartDate",
	 "EndDate",
	 "SwapStartDate",
	 "NextSwapStartDate",
	 "LastSwapStartDate",
	 "Strike",
	 "Fees",
	 "Notional",
	 "StdFixSwaplet",
	 "StdVarSwaplet",
	 "StdSwaplet",
	 "StdSwap",
	 "Option",
	 "Bermuda",
	 "SwapRate",
	 "Frontier",
	 "DealEndDate",
	 "ExerciseCondition",
};

const string ARM_BermudaSwaptionCalculator::BermudaSwaptionProbaColNamesTable2 [] =
{
	"ResetPlus1",
	"FinalDate",
	"ExerciseIndicator",
	"Prob10",
	"Prob2",
	"Prob20",
	"Prob3",
	"Prob30",
	"Prob4",
	"Prob40",
	"Prob5",
	"Prob50",
	"Prob6",
	"Prob60",
	"Prob7",
	"Prob70",
	"Prob8",
	"Prob80",
	"Prob9",
	"Prob90",
	"Proba10",
	"Proba100",
	"Proba11",
	"Proba110",
	"Proba12",
	"Proba120",
	"Proba13",
	"Proba130",
	"Proba14",
	"Proba140",
	"Proba15",
	"Proba150",
	"Proba16",
	"Proba160",
	"Proba17",
	"Proba170",
	"Proba18",
	"Proba180",
	"Proba19",
	"Proba190",
	"Proba20",
	"Proba200",
	"ProbFinal",
	"ProbFinal0",
};
const string ARM_BermudaSwaptionCalculator::BermudaSwaptionProbaPricedColNamesTable2 [] =
{
	"Prob10",
	"Prob20",
	"Prob30",
	"Prob40",
	"Prob50",
	"Prob60",
	"Prob70",
	"Prob80",
	"Prob90",
	"Proba100",
	"Proba110",
	"Proba120",
	"Proba130",
	"Proba140",
	"Proba150",
	"Proba160",
	"Proba170",
	"Proba180",
	"Proba190",
	"Proba200",
	"ProbFinal0",
};

//GenSec 3
const string ARM_BermudaSwaptionCalculator::BermudaSwaptionColNamesTable3 [] =
{
     "ResetDate",
	 "StartDate",
	 "EndDate",
	 "SwapStartDate",
	 "Strike",
	 "Fees",
	 "StdSwaplet",
	 "ReverseStdSwaplet",
	 "OptionAnnulation",
	 "Bermuda",
};

const string ARM_BermudaSwaptionCalculator::ControlVariateColNamesTable3 [] = 
{
	"CVOptionAnnulation1",
	"ControlVariate1",
	"CVOptionAnnulation2",
	"ControlVariate2",	
	"CVOptionAnnulation3",
	"ControlVariate3",	
	"CVOptionAnnulation4",
	"ControlVariate4",	
	"CVOptionAnnulation5",
	"ControlVariate5",	
	"CVOptionAnnulation6",
	"ControlVariate6",
	"CVOptionAnnulation7",
	"ControlVariate7",
	"CVOptionAnnulation8",
	"ControlVariate8",
	"CVOptionAnnulation9",
	"ControlVariate9",
	"CVOptionAnnulation10",
	"ControlVariate10",
};

//GenSecType = 4
const string ARM_BermudaSwaptionCalculator::BermudaSwaptionColNamesTable4 [] =
{
	 "ResetDate",
	 "StartDate",
	 "PastStartDate",
	 "SwapStartDate",
	 "EndDate",
	 "Strike",
	 "Fees",
	 "DFStart",
	 "CvgStart",
	 "DFEnd",
	 "CvgEnd",
	 "Annuity",
	 "Flows",
	 "CashFlow",
	 "Bermuda",
 	 "SwapRate",
	 "Frontier",
	 "DealEndDate",
};

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator //////////////////////
///	Routine: Constructor ////////////////////////////////////////
///	Returns: void ///////////////////////////////////////////////
///	Action : builds the object //////////////////////////////////
/////////////////////////////////////////////////////////////////
ARM_BermudaSwaptionCalculator::ARM_BermudaSwaptionCalculator(
								  ARM_Currency& ccy,
								  ARM_Date&	startDate,
								  ARM_Date&	endDate,
								  ARM_ReferenceValue&  nominal,
								  ARM_ReferenceValue&  strike,
								  int payRec,
								  int stubRule,
								  int callFreq,
								  int callNotice,
								  string callCal,
								  ARM_Date firstCallDate,
								  ARM_Date lastCallDate,
								  ARM_ReferenceValue& fees,
								  int fixFreq,
								  int fixBasis,
								  string fixPayCal,
								  int fixAdjRule,
								  int fixRule,
								  int fixPayGap,
								  bool isZc,
								  int varFreq,
								  int varBasis,
								  string varResetCal,
								  string varPayCal,
								  string varIndexTerm,
								  int varAdjRule,
								  int varRule,
								  int varResetGap,
								  int varPayGap,
								  ARM_ReferenceValue& varSpread,
								  int genSecType,
								  vector<int>* controlVariate,
								  vector<double>* controlPrices,
								  ARM_StringVector& mdmKeys,
								  const ARM_MarketData_ManagerRep& mktDataManager,
								  int modelType,
								  vector<ARM_ReferenceValue*>* modelParams,
								  bool mrsCalibFlag,
								  bool atmDiagonalFlag,
								  int numMethodType,
								  int amcIter,
								  int mcIter,
								  int maxBucketSize,
								  string genType1,
								  string genType2,
								  string pathOrder,
								  string pathScheme,
								  int firstNbTimes,
								  int treeSteps,
								  vector<int>& portfolioMode,
								  bool fixBoundaryFlag,
								  bool approxMarginFlag,
								  bool freezeBetasFlag,
								  bool calculateProbaFlag)   							  
: ARM_GenCalculator(mktDataManager),
	itsCcy(ccy),
	itsPastZC(false),
	itsPastStartDate(startDate),
	itsStartDate(startDate),
	itsEndDate(endDate),
	itsNominal(nominal),
	itsStrike(strike),
	itsPayRec(payRec),
	itsStubRule(stubRule),
	itsCallFreq(callFreq),
	itsCallNotice(callNotice),
	itsCallCal(callCal),
	itsFirstCallDate(firstCallDate),
	itsLastCallDate(lastCallDate),
	itsFees(fees),
	itsFixFreq(fixFreq),
	itsFixBasis(fixBasis),
	itsFixPayCal(fixPayCal),
	itsFixAdjRule(fixAdjRule),
	itsFixRule(fixRule),
	itsFixPayGap(fixPayGap),
	itsZC(isZc),
	itsVarFreq(varFreq),
	itsVarBasis(varBasis),
	itsVarResetCal(varResetCal),
	itsVarPayCal(varPayCal),
	itsVarIndexTerm(varIndexTerm),
	itsVarAdjRule(varAdjRule),
	itsVarRule(varRule),
	itsVarResetGap(varResetGap),
	itsVarPayGap(varPayGap),
	itsVarSpread(varSpread),
	itsGenSecType(genSecType),
	itsControlVariate(controlVariate),
	itsControlVariatePrices(controlPrices),
	itsExoSwaption(NULL),
	itsUnderSwap(NULL),

	//Model Params
	itsModelParams(NULL),
	itsMrsCalibFlag(mrsCalibFlag),
	itsATMDiagonalFlag(atmDiagonalFlag),
	itsModelType(modelType),
	itsNumMethodType(numMethodType),
	itsAmcIter(amcIter),
	itsMcIter(mcIter),
	itsMaxBucketSize(maxBucketSize),
	itsGenType1(genType1),
	itsGenType2(genType2),
	itsPathOrder(pathOrder),
	itsPathScheme(pathScheme),
	itsFirstNbTimes(firstNbTimes),
	itsTreeSteps(treeSteps),
	itsPortfolioMode(portfolioMode),
	itsHasBeenPriced(false),
	itsBermudaSwaptionPrice(0.0),
	itsBermudaSwaptionStdDev(0.0),
	itsFixBoundaryFlag(fixBoundaryFlag),
	itsApproxMarginFlag(approxMarginFlag),
	itsFreezeBetasFlag(freezeBetasFlag),
	itsCalculateProbaFlag(calculateProbaFlag),
	itsFundDateStrip(0),
	itsStructDateStrip(0),
	itsBetasCoeff(0),
	itsCallDateStrip(0),
	itsPreCalibMethod(0)
{
	SetModelParams(modelParams);

	/// Creation of the Exercise Style for the Exo-Swaption Creation
	ARM_ExerciseStyle* exoDates  = CreateExerciseStyle();

	itsNbCalculatedProbabilities = (NB_CALCULATED_PROBA<exoDates->GetExerciseStartDates()->size() )?NB_CALCULATED_PROBA:exoDates->GetExerciseStartDates()->size();
	itsProbaMatrix				 = ARM_GP_Matrix(itsNbCalculatedProbabilities,itsNbCalculatedProbabilities,0.0);
	itsCumProba					 = ARM_VectorPtr(new std::vector<double>(itsNbCalculatedProbabilities,0.0));

	SetCurrencyUnit(&ccy); 

	SetKeys(mdmKeys);

	//Update the GenSecType.
	UpdateGenSecType(GetGenSecType());

	//Pertinent checks before starting...
	CheckData();
	CheckSwaptionInputs();

	//Propagate the right diffusion.
	SetEquivalentVarIndexTerm();	
	CreateEquivalentSwaption(exoDates);

	//Create the const manager
	ARM_CstManagerPtr cstManagerPtr = CreateCstManager();

	itsAtmMrsCalibType  = itsPortfolioMode[2];
	itsDefaultMrs		= (*itsModelParams)[1]->Interpolate(mktDataManager.GetAsOfDate().GetJulian());

	//Priced Columns.
	ARM_StringVector pricedColumns = PricedColumnNames();

	/// Flag for other payoffs in the creation of the generic security
	bool CVflag = false;
	if (itsControlVariate->size() > 0)
	{
		CVflag = true;
	}

	bool otherPayoffs = ((itsCalculateProbaFlag && (GetGenSecType() == 2)) || CVflag || (GetOSWPortfolioMode()==6));

	//Create the generic security with the const manager.
	CreateAndSetDealDescriptionAndTimeIt("",pricedColumns, cstManagerPtr,false,otherPayoffs);

	//Create the model.
	CreateAndSetModel();

	//Caplet/floorlet strike spread adjustment
    CreateAndSetCalibration();

	//Free memory
	delete exoDates;
	exoDates = NULL;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: Summit constructor
///	Returns: a built object
///	Action : builds the object
/////////////////////////////////////////////////////////////////
ARM_BermudaSwaptionCalculator::ARM_BermudaSwaptionCalculator(
								  const ARM_Date&	asOfDate,
								  ARM_Currency& ccy,
								  ARM_Date&	startDate,
								  ARM_Date&	endDate,
								  ARM_ReferenceValue&  nominal,
								  ARM_ReferenceValue&  strike,
								  int payRec,
								  int stubRule,
								  int callFreq,
								  int callNotice,
								  string callCal,
								  ARM_Date firstCallDate,
								  ARM_Date lastCallDate,
								  ARM_ReferenceValue& fees,
								  int fixFreq,
								  int fixBasis,
								  string fixPayCal,
								  int fixAdjRule,
								  int fixRule,
								  int fixPayGap,
								  bool isZc,
								  int varFreq,
								  int varBasis,
								  string varResetCal,
								  string varPayCal,
								  string varIndexTerm,
								  int varAdjRule,
								  int varRule,
								  int varResetGap,
								  int varPayGap,
								  ARM_ReferenceValue& varSpread,
								  int genSecType,
								  bool isPastZc,
								  ARM_Date& pastStartDate)
: ARM_GenCalculator(asOfDate),
	itsCcy(ccy),
	itsPastZC(isPastZc),
	itsPastStartDate(pastStartDate),
	itsStartDate(startDate),
	itsEndDate(endDate),
	itsNominal(nominal),
	itsStrike(strike),
	itsPayRec(payRec),
	itsStubRule(stubRule),
	itsCallFreq(callFreq),
	itsCallNotice(callNotice),
	itsCallCal(callCal),
	itsFirstCallDate(firstCallDate),
	itsLastCallDate(lastCallDate),
	itsFees(fees),
	itsFixFreq(fixFreq),
	itsFixBasis(fixBasis),
	itsFixPayCal(fixPayCal),
	itsFixAdjRule(fixAdjRule),
	itsFixRule(fixRule),
	itsFixPayGap(fixPayGap),
	itsZC(isZc),
	itsVarFreq(varFreq),
	itsVarBasis(varBasis),
	itsVarResetCal(varResetCal),
	itsVarPayCal(varPayCal),
	itsVarIndexTerm(varIndexTerm),
	itsVarAdjRule(varAdjRule),
	itsVarRule(varRule),
	itsVarResetGap(varResetGap),
	itsVarPayGap(varPayGap),
	itsVarSpread(varSpread),
	itsGenSecType(genSecType),
	itsControlVariate(new vector<int>(0)),
	itsControlVariatePrices(new vector<double>(0)),
	itsModelParams(NULL),
	itsHasBeenPriced(false),
	itsBermudaSwaptionPrice(0.0),
	itsBermudaSwaptionStdDev(0.0),
	itsExoSwaption(NULL),
	itsUnderSwap(NULL),

	//Default calib params
	itsMrsCalibFlag(false),
	itsATMDiagonalFlag(false),
	itsModelType(0),
	itsNumMethodType(0),
	itsAmcIter(0),
	itsMcIter(0),
	itsMaxBucketSize(0),
	itsGenType1("MRGK5"),
	itsGenType2("Sobol"),
	itsPathOrder("BucketOrder"),
	itsPathScheme("Incremental"),
	itsFirstNbTimes(0),
	itsTreeSteps(0),
	itsPortfolioMode(vector<int>(0)),
	itsAtmMrsCalibType(0),
	itsFixBoundaryFlag(false),
	itsApproxMarginFlag(false),
	itsFreezeBetasFlag(false),
	itsDefaultMrs(0.0),
	itsFundDateStrip(0),
	itsStructDateStrip(0),
	itsBetasCoeff(new std::vector<double>(NULL))
{
	/// Creation of the Exercise Style for the Exo-Swaption Creation
	ARM_ExerciseStyle* exoDates  = CreateExerciseStyle();
	itsPreCalibMethod			 = ARM_CalibMethodPtr(NULL);
	itsCallDateStrip			 = ARM_DateStripPtr(NULL);
	itsNbCalculatedProbabilities = (NB_CALCULATED_PROBA<exoDates->GetExerciseStartDates()->size() )?NB_CALCULATED_PROBA:exoDates->GetExerciseStartDates()->size();
	itsProbaMatrix				 = ARM_GP_Matrix(itsNbCalculatedProbabilities,itsNbCalculatedProbabilities,0.0);
	itsCumProba					 = ARM_VectorPtr(new std::vector<double>(itsNbCalculatedProbabilities,0.0));

	//Update the GenSecType.?????????????
	SetCurrencyUnit(&ccy); 

	//Propagate the right diffusion.
	SetEquivalentVarIndexTerm();
	
	//Equivalent swaption.
	CreateEquivalentSwaption(exoDates);

	//Free memory
	delete exoDates;
	exoDates = NULL;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: Summit constructor
///	Returns: a built object
///	Action : builds the object
/////////////////////////////////////////////////////////////////
ARM_BermudaSwaptionCalculator::ARM_BermudaSwaptionCalculator(
								  const ARM_Date&	asOfDate,
								  ARM_Swaption* swaption,
								  int fixFreqIfZC)
:	ARM_GenCalculator(asOfDate),
	itsControlVariate(new vector<int>(0)),
	itsControlVariatePrices(new vector<double>(0)),
	itsModelParams(NULL),
	itsHasBeenPriced(false),
	itsBermudaSwaptionPrice(0.0),
	itsBermudaSwaptionStdDev(0.0),
	itsGenSecType(-1),
	itsFixFreq(fixFreqIfZC),
	itsPastZC(false),
	itsUnderSwap(NULL),
	itsFees(*(swaption->GetFee())),

	//Default calib params
	itsMrsCalibFlag(false),
	itsATMDiagonalFlag(false),
	itsModelType(0),
	itsNumMethodType(0),
	itsAmcIter(0),
	itsMcIter(0),
	itsMaxBucketSize(0),
	itsGenType1("MRGK5"),
	itsGenType2("Sobol"),
	itsPathOrder("BucketOrder"),
	itsPathScheme("Incremental"),
	itsFirstNbTimes(0),
	itsTreeSteps(0),
	itsPortfolioMode(vector<int>(0)),
	itsAtmMrsCalibType(0),
	itsFixBoundaryFlag(false),
	itsApproxMarginFlag(false),
	itsFreezeBetasFlag(false),
	itsDefaultMrs(0.0),
	itsBetasCoeff(new std::vector<double>(NULL))
{
	SetCurrencyUnit(swaption->GetCurrencyUnit()); 

	itsExoSwaption = CreateForwardSwaption(swaption); //also set itsUnderSwap
	GenerateProductDescription();

	/// Creation of the Exercise Style for the Exo-Swaption Creation
	itsPreCalibMethod			 = ARM_CalibMethodPtr(NULL);
	itsCallDateStrip			 = ARM_DateStripPtr(NULL);
	itsNbCalculatedProbabilities = MAX(NB_CALCULATED_PROBA, itsExoSwaption->GetExerciseStyle()->GetExerciseStartDates()->size());
	itsProbaMatrix				 = ARM_GP_Matrix(itsNbCalculatedProbabilities, itsNbCalculatedProbabilities, 0.0);
	itsCumProba					 = ARM_VectorPtr(new std::vector<double>(itsNbCalculatedProbabilities,0.0));
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: CreateForwardSwaption
///	Returns: swaption
///	Action : keep forward flows of the swaption
/////////////////////////////////////////////////////////////////
ARM_Swaption* ARM_BermudaSwaptionCalculator::CreateForwardSwaption(ARM_Swaption* swaption)
{
	try
	{
		ARM_Swaption* newSwaption = NULL;

		ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();
		// keep only forward flows !!!
		ARM_SwapLeg* floatLeg = swaption->GetFloatLeg();
		ARM_SwapLeg* fixedLeg = swaption->GetFixedLeg();
		
		ARM_Date newStartDate;
		int i=0, j=0;
		// WARNING : if floatLeg frequency <> fixedLeg frequency, we may create a stub !
		if (floatLeg->GetIRIndex()->GetResetFrequency() >= fixedLeg->GetIRIndex()->GetPayFrequency())
		{
			while (asOfDate.GetJulian() > floatLeg->GetResetDates()->Elt(i))
			{
				i++;
				if (i >= floatLeg->GetResetDates()->size())
					ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : last floating leg period already fixed. Swaption should have at least one forward period.");
			}
			newStartDate = ARM_Date(floatLeg->GetFlowStartDates()->Elt(i));
		}
		else
		{
			while (asOfDate.GetJulian() > fixedLeg->GetResetDates()->Elt(j))
			{
				j++;
				if (j >= fixedLeg->GetResetDates()->size())
					ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : last fixed leg period already fixed. Swaption should have at least one forward period.");
			}
			newStartDate = ARM_Date(fixedLeg->GetFlowStartDates()->Elt(j));
		}

		//Exercise Dates
		ARM_ExerciseStyle* newExerStyle = NULL;
		ARM_Vector* newExerDates = new ARM_Vector();
		ARM_Vector* newEndDates = new ARM_Vector();
		for (i=0; i<swaption->GetExerciseStyle()->GetExerciseStartDates()->size(); i++)
		{
			if (asOfDate < ARM_Date(swaption->GetExerciseStyle()->GetExerciseStartDates()->Elt(i)))
			{
				newExerDates->push_back(swaption->GetExerciseStyle()->GetExerciseStartDates()->Elt(i));
				newEndDates->push_back(swaption->GetExerciseStyle()->GetExerciseEndDates()->Elt(i));
			}
		}
		newExerStyle = new ARM_ExerciseStyle(newExerDates,newEndDates);
		newExerStyle->SetExerciseType(K_BERMUDAN);
		delete newExerDates;
		delete newEndDates;

		newStartDate = ARM_Date(newExerStyle->GetExerciseEndDates()->Elt(0)); //TMP

		ARM_FixLeg* fixLeg = new ARM_FixLeg( newStartDate, 
											 fixedLeg->GetEndDateNA(),
											 swaption->GetStrikes(), //fixedLeg->GetFixedRate(), 
											 swaption->GetOptionType(),
											 fixedLeg->GetPaymentFreq(), 
											 fixedLeg->GetDayCount(), 
											 fixedLeg->GetDecompFreq(),
											 fixedLeg->GetIRIndex()->GetPayTiming(),
											 fixedLeg->GetIRIndex()->GetIntRule(),
											 fixedLeg->GetStubMeth(),
											 fixedLeg->GetCurrencyUnit(),
											 fixedLeg->GetPayCalName(),
											 fixedLeg->GetNxFlag(),
											 NULL, //ref date
											 1);   //adjust start

		ARM_SwapLeg* varLeg = new ARM_SwapLeg(newStartDate,
											 floatLeg->GetEndDateNA(),
											 floatLeg->GetIRIndex()->GetIndexType(),
											 -swaption->GetOptionType(),
											 floatLeg->GetSpread(),
											 floatLeg->GetIRIndex()->GetResetFrequency(),
											 floatLeg->GetIRIndex()->GetPayFrequency(),
											 floatLeg->GetIRIndex()->GetResetTiming(),
											 floatLeg->GetIRIndex()->GetPayTiming(),
											 floatLeg->GetCurrencyUnit(),
											 floatLeg->GetIRIndex()->GetIntRule(),
											 floatLeg->GetIRIndex()->GetResetGap(),
											 floatLeg->GetResetCalName(),
											 floatLeg->GetPayCalName(),
											 floatLeg->GetDecompPricingFlag(),
											 floatLeg->GetNxFlag(),
											 floatLeg->GetStubMeth(),
											 NULL, //ref date
											 1,    //adjust start
											 floatLeg->GetDayCount(),
											 floatLeg->GetIRIndex()->GetFwdRule(),
											 floatLeg->GetIRIndex()->GetPayGap());

		if (floatLeg->GetVariableSpread())
		{
			ARM_ReferenceValue* newSpreads = NULL;
			if (floatLeg->GetVariableSpread()->GetSize()>1)
			{
				newSpreads = floatLeg->GetVariableSpread()->CptReferenceValues(floatLeg->GetResetDates());
				varLeg->SetVariableSpread(newSpreads);
			}
			else
			{
				varLeg->SetVariableSpread(floatLeg->GetVariableSpread());
			}
			delete newSpreads;
		}

		// rebuild variable strikes to match fixleg's schedule
		ARM_ReferenceValue* newStrikes = NULL;
		if (swaption->GetStrikes()->GetSize() > 1)
		{
			newStrikes = swaption->GetStrikes()->CptReferenceValues(fixLeg->GetResetDates());
			fixLeg->SetVarCoupons(newStrikes);
		}
		else
			newStrikes = (ARM_ReferenceValue*)swaption->GetStrikes()->Clone();

		ARM_Swap* swap = new ARM_Swap(fixLeg, varLeg);
		SetUnderSwap(swap);
		delete fixLeg;
		delete varLeg;
	
		newSwaption = new ARM_Swaption(swap, 
									   swaption->GetOptionType(),
									   newExerStyle,
									   newStrikes,
									   0);
		delete newExerStyle;
		delete newStrikes;

		newSwaption->SetAmount(floatLeg->GetAmount());
		newSwaption->SetPorS(swaption->GetPorS());
		newSwaption->SetFee(swaption->GetFee());

		return newSwaption;
	}
	catch (Exception& x)
	{
	    x.DebugPrint();
		throw x;
	}
	catch (...)
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in ARM_BermudaSwaptionCalculator::CreateForwardSwaption()" );
	}
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: GenerateProductDescription
///	Returns: void
///	Action : set calculator's attributes from swaption
/////////////////////////////////////////////////////////////////
void ARM_BermudaSwaptionCalculator::GenerateProductDescription()
{
	try
	{
		ARM_SwapLeg* floatLeg = itsExoSwaption->GetFloatLeg();
		ARM_SwapLeg* fixedLeg = itsExoSwaption->GetFixedLeg();

		itsStrike = *((ARM_ReferenceValue*)(itsExoSwaption->GetStrikes()->Clone()));
		itsStrike /= 100.0;
		itsPayRec = itsExoSwaption->GetOptionType();
		itsStubRule = floatLeg->GetStubMeth();

		int callSize = itsExoSwaption->GetExerciseStyle()->GetExerciseStartDates()->size();
		itsFirstCallDate = itsExoSwaption->GetExerciseStyle()->GetExerciseEndDates()->Elt(0);
		itsLastCallDate = itsExoSwaption->GetExerciseStyle()->GetExerciseEndDates()->Elt(callSize-1);
		// callFreq = moyenne sur l'ensemble du deal
		itsCallFreq = (callSize > 1) ? (int)(K_YEAR_LEN / ((itsLastCallDate - itsFirstCallDate) / (callSize-1)) + 0.5) : K_ANNUAL;
		itsCallNotice = 9999;
		itsCallCal = string(floatLeg->GetResetCalName());

		itsCcy = *floatLeg->GetCurrencyUnit();
		itsStartDate = floatLeg->GetStartDateNA();
		itsEndDate = floatLeg->GetEndDateNA();

		int floatSize = floatLeg->GetAmount()->size();
		if (floatSize == 1)
		{
			ARM_Vector* vecValues = (ARM_Vector*)(floatLeg->GetAmount()->GetDiscreteValues())->Clone();
			ARM_Vector* vecDates = new ARM_Vector(1, itsStartDate.GetJulian());
			itsNominal = ARM_ReferenceValue(vecDates, vecValues);
		}
		else
		{
			itsNominal = *((ARM_ReferenceValue*)(floatLeg->GetAmount())->Clone());
		}

		itsZC = (fixedLeg->GetPaymentFreq() == K_ZEROCOUPON);
		// WARNING : if fix freq = ZC, take "COMP_compFreq" already stored in itsFixFreq,
		// otherwise, replace by the real freq of the leg
		if (!itsZC)
			itsFixFreq = fixedLeg->GetPaymentFreq();

		itsPastStartDate = fixedLeg->GetStartDate();
		if (itsZC && (itsPastStartDate.GetJulian() != itsExoSwaption->GetExerciseStyle()->GetExerciseEndDates()->Elt(0)))
			itsPastZC = true; //expiry date doesn't match start date

		itsFixBasis = fixedLeg->GetDayCount();
		itsFixPayCal = string(fixedLeg->GetPayCalName());
		itsFixAdjRule = fixedLeg->GetIRIndex()->GetFwdRule();
		itsFixRule = itsZC ? K_UNADJUSTED : K_ADJUSTED; //TMP
		itsFixPayGap = fixedLeg->GetIRIndex()->GetPayGap();
		itsVarFreq = floatLeg->GetIRIndex()->GetResetFrequency();
		itsVarBasis = floatLeg->GetIndexDayCount();
		itsVarResetCal = string(floatLeg->GetResetCalName());
		itsVarPayCal = string(floatLeg->GetPayCalName());

		int period = 12/itsVarFreq;
		if ( period == 12 )
		{
			itsVarIndexTerm = string("1Y");
		}
		else
		{
			char str[2];
			sprintf(str, "%d", period);
			itsVarIndexTerm = string(str) + string("M");
		}

		itsVarAdjRule = floatLeg->GetIRIndex()->GetFwdRule();
		itsVarRule = floatLeg->GetIRIndex()->GetIntRule();
		itsVarResetGap = floatLeg->GetIRIndex()->GetResetGap();
		itsVarPayGap = floatLeg->GetIRIndex()->GetPayGap();

		if (floatLeg->GetVariableSpread())
		{
			ARM_ReferenceValue* tempVarFundSpread = (ARM_ReferenceValue*) floatLeg->GetVariableSpread()->Clone();
			*(tempVarFundSpread->GetDiscreteValues()) /= 100.0; 
			itsVarSpread = (*tempVarFundSpread);
			delete tempVarFundSpread;
		}
		else
		{
			double spread = (floatLeg->GetSpread())/100;
			ARM_Vector* vecDates = new ARM_Vector(1, itsStartDate.GetJulian());
			ARM_Vector* vecFundSpread = new ARM_Vector(1, spread);
			ARM_ReferenceValue fundSpread = ARM_ReferenceValue(vecDates, vecFundSpread);
			itsVarSpread = fundSpread;
		}	
	}
	catch (Exception& x)
	{
	    x.DebugPrint();
		throw x;
	}
	catch (...)
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in ARM_BermudaSwaptionCalculator::GenerateProductDescription()" );
	}
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: InitBermudaSwaptionForSummit
///	Returns: a built object
///	Action : builds the object
/////////////////////////////////////////////////////////////////
void ARM_BermudaSwaptionCalculator::InitBermudaSwaptionForSummit( ARM_StringVector&	mdmKeys,
																  ARM_MarketData_ManagerRep* mktDataManager,
																  vector<int>* controlVariates,
																  vector<double>* controlPrices,
																  vector<ARM_ReferenceValue*>* modelParams,
																  bool mrsCalibFlag,
																  bool atmDiagonalFlag,
																  int numMethodType,
																  int amcIter,
																  int mcIter,
																  int maxBucketSize,
																  string genType,
																  int treeSteps,
																  vector<int> portfolioMode,
																  bool boundaryFlag,
																  bool approxMarginFlag,
																  bool freezeBetasFlag,
																  int modelType,
																  bool calculateProbaFlag)
{
	//Fullfill calibParams
	SetModelParams(modelParams);
	SetMrsCalibFlag(mrsCalibFlag);
	SetATMDiagonalFlag(atmDiagonalFlag);
	SetModelType(modelType);
	SetNumMethodType(numMethodType);
	SetAmcIter(amcIter);
	SetMcIter(mcIter);
	SetMaxBucketSize(maxBucketSize);
	SetGenType(genType);
	SetTreeSteps(treeSteps);
	SetPortfolioMode(portfolioMode);
	SetFixBoundaryFlag(boundaryFlag);
	SetApproxMarginFlag(approxMarginFlag);
	SetFreezeBetasFlag(freezeBetasFlag);
	SetCalculateProbaFlag(calculateProbaFlag);

	itsAtmMrsCalibType  = itsPortfolioMode[2] ;
	itsDefaultMrs		= (*itsModelParams)[1]->Interpolate((&*GetMktDataManager())->GetAsOfDate().GetJulian());

	//Fullfill control variates
	SetControlVariate(controlVariates);
	SetControlVariatePrices(controlPrices);

	//Update the genSecType
	UpdateGenSecType(GetGenSecType());

	//Update mdmKeys
	SetKeys(mdmKeys);
		
	//Update mktDataManager 
	InitializeMktDataManagerOnly(*mktDataManager);

	//Equivalent swaption.       SCOTCH!!!!!!!!!!!!!
	if (GetZCFlag() == true)
	{
		SetGenSecType(4);
		//TMP
		//ARM_ExerciseStyle* exerciseDates = (GetExoSwaption())->GetExerciseStyle();
		//delete itsExoSwaption;
		//CreateEquivalentSwaption(exerciseDates);
	}

	//Pertinent checks before starting...
	CheckData();
	CheckSwaptionInputs();

	//Priced Columns.
	ARM_StringVector pricedColumns  = PricedColumnNames();

	//Create the const manager
	ARM_CstManagerPtr cstManagerPtr = CreateCstManager();

	// Flag for other payoffs in the creation of the generic security
	bool CVflag = false;
	if (itsControlVariate->size()>0)
	{
		CVflag = true;
	}
	bool otherPayoffs = ((itsCalculateProbaFlag && (GetGenSecType() == 2)) || CVflag || (GetOSWPortfolioMode()==6));

	//Create the generic security with the const manager.
    CreateAndSetDealDescriptionAndTimeIt("",pricedColumns, cstManagerPtr,false, otherPayoffs);
	
	//Create the model.
	CreateAndSetModel();

	//Caplet/floorlet strike spread adjustment
    CreateAndSetCalibration();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: InitBermudaSwaptionForSummit
///	Returns: a built object
///	Action : builds the object
/////////////////////////////////////////////////////////////////
void ARM_BermudaSwaptionCalculator::InitBermudaSwaptionForSummit(  vector<int>* controlVariates,
																   vector<double>* controlPrices,
																   vector<ARM_ReferenceValue*>* modelParams,
																   bool mrsCalibFlag,
																   bool atmDiagonalFlag,
																   int numMethodType,
																   int amcIter,
																   int mcIter,
																   int maxBucketSize,
																   string genType,
																   int treeSteps,
																   vector<int> portfolioMode,
																   bool boundaryFlag,
																   bool approxMarginFlag,
																   bool freezeBetasFlag,
																   int modelType,
																   bool calculateProbaFlag,
																   ARM_ZeroCurve* zcCpn, 
																   ARM_VolCurve* swoptVC, 
																   ARM_VolCurve* capVC,
																   ARM_VolLInterpol* capRo, 
																   ARM_VolLInterpol* capNu,
																   ARM_VolLInterpol* capBeta,
																   ARM_VolLInterpol* swoptRo, 
																   ARM_VolLInterpol* swoptNu,
																   ARM_VolLInterpol* swoptBeta,
																   ARM_MarketIRModel* normalModel,
																   int hedgeUpdate,
																   int SABRSigmaOrAlpha)
{
	//Fullfill calib parameters ------------------------------------------------------>
	SetModelParams(modelParams);
	SetMrsCalibFlag(mrsCalibFlag);
	SetATMDiagonalFlag(atmDiagonalFlag);
	SetNumMethodType(numMethodType);
	SetAmcIter(amcIter);
	SetMcIter(mcIter);
	SetMaxBucketSize(maxBucketSize);
	SetGenType(genType);
	SetTreeSteps(treeSteps);
	SetPortfolioMode(portfolioMode);
	SetFixBoundaryFlag(boundaryFlag);
	SetApproxMarginFlag(approxMarginFlag);
	SetFreezeBetasFlag(freezeBetasFlag);
	SetModelType(modelType);
	SetCalculateProbaFlag(calculateProbaFlag);

	itsAtmMrsCalibType  = itsPortfolioMode[2] ;
	itsDefaultMrs		= (*itsModelParams)[1]->Interpolate((&*GetMktDataManager())->GetAsOfDate().GetJulian());
	
	//Fullfill control variates
	SetControlVariate(controlVariates);
	SetControlVariatePrices(controlPrices);

	//Update the genSecType
	UpdateGenSecType(GetGenSecType());

	//Model Keys -------------------------------------------------------------------->
	int nbKeys = 0;
	if (normalModel)
	{
		nbKeys = 4;
	}
	else
	{
		nbKeys = 3;
	}

	ARM_StringVector keys(nbKeys);
	string ccyName(zcCpn->GetCurrencyUnit()->GetCcyName());

	keys[YcKey]       = YC_KEY_NAME			+ ccyName;
    keys[OswModelKey] = OSWMODEL_KEY_NAME	+ ccyName;
    keys[CfModelKey]  = CFMODEL_KEY_NAME	+ ccyName;
  
	//The normal model is used if we calibrate the mean reversion speed.
	if (normalModel)
	{
		keys[NormalModel] = NMOD_KEY_NAME	+ ccyName;
	}

    SetKeys(keys);
	
	//Market Data Manager ----------------------------------------------------------->
	ARM_Date asof = GetMktDataManager()->GetAsOfDate();

	///Swopt Volatility Model
	ARM_BSModel* swoptBSmod = NULL;	
	ARM_VolCurve* ATMSwoptVol = (swoptVC->GetName() == ARM_VOL_CUBE) ? ((ARM_VolCube*)swoptVC)->GetATMVol() : swoptVC;

	//Sabr with Beta = 1: 
	if (swoptRo && swoptNu && !swoptBeta)
	{
		swoptBSmod = new ARM_BSSmiledModel( asof, 
											0.0, 
											zcCpn, 
											zcCpn, 
											ATMSwoptVol, 
											K_YIELD, 
											swoptRo, 
											swoptNu, 
											K_SABR_ARITH,
											NULL,
											0.5, // SABR Weight: Irrelevant
											SABRSigmaOrAlpha);
	}
	//Complete Sabr:
	else if (swoptRo && swoptNu && swoptBeta)
	{
		int methodType;

		if (swoptBeta->IsEqualToOne())
			methodType = K_SABR_ARITH;
		else
			methodType = K_SABR_IMPLNVOL;

		swoptBSmod = new ARM_BSSmiledModel( asof, 
											0.0, 
											zcCpn, 
											zcCpn, 
											ATMSwoptVol, 
											K_YIELD, 
											swoptRo, 
											swoptNu, 
											methodType, 
											swoptBeta,
											0.5, // SABR Weight: Irrelevant
											SABRSigmaOrAlpha);
	}
	else if (!swoptRo && !swoptNu && !swoptBeta)
	{
		swoptBSmod = new ARM_BSModel(asof, 0.0, zcCpn, zcCpn, swoptVC, K_YIELD);
	}

	//IRG Volatility Model 
	ARM_BSModel* capBSmod = NULL;
	ARM_VolCurve* ATMCapVol = (capVC->GetName() == ARM_VOL_CUBE) ? ((ARM_VolCube*)capVC)->GetATMVol() : capVC;

	if (capRo && capNu && !capBeta)
	{
		capBSmod = new ARM_BSSmiledModel(asof, 
										 0.0, 
										 zcCpn, 
										 zcCpn, 
										 ATMCapVol, 
										 K_YIELD, 
										 capRo, 
										 capNu, 
										 K_SABR_ARITH,
										 NULL,
										 0.5, // SABR Weight: Irrelevant
										 SABRSigmaOrAlpha);
	}
	//Complete Sabr :
	else if (capRo && capNu && capBeta)
	{
		int methodType;

		if (capBeta->IsEqualToOne())
			methodType = K_SABR_ARITH;
		else
			methodType = K_SABR_IMPLNVOL;

		capBSmod = new ARM_BSSmiledModel(asof, 
										 0.0, 
										 zcCpn, 
										 zcCpn, 
										 ATMCapVol, 
										 K_YIELD, 
										 capRo, 
										 capNu, 
										 methodType, 
										 capBeta,
										 0.5, // SABR Weight: Irrelevant
										 SABRSigmaOrAlpha);
	}
	//No Sabr --> Volatility Cube
	else if (!capRo && !capNu && !capBeta)
	{
		capBSmod = new ARM_BSModel( asof, 0.0, zcCpn, zcCpn, capVC, K_YIELD);
	}

	//Create mktData Manager & fill it with curves
	vector <ARM_Object*> marketDatas;

	marketDatas.push_back(zcCpn);
	marketDatas.push_back(swoptBSmod);
	marketDatas.push_back(capBSmod);
	
	if (normalModel)	
	{
		marketDatas.push_back(normalModel);
	}

	 Init(marketDatas,false);

	//Equivalent swaption.       SCOTCH!!!!!!!!!!!!!
	if (GetZCFlag() == true)
	{
		SetGenSecType(4);
		// TMP
		//ARM_ExerciseStyle* exerciseDates = (GetExoSwaption())->GetExerciseStyle();
		//delete itsExoSwaption;
		//CreateEquivalentSwaption(exerciseDates);
	}
	
	//Flag for other payoffs in the creation of the generic security
	bool CVflag = false;
	if (itsControlVariate->size()>0)
	{
		CVflag = true;
	}
	bool otherPayoffs = ((itsCalculateProbaFlag && (GetGenSecType() == 2)) || CVflag || (GetOSWPortfolioMode()==6));

	//Pertinent checks before starting ----------------------------------------------->
	CheckSwaptionInputs();

	//Priced Columns.
	ARM_StringVector pricedColumns = PricedColumnNames();

	//Create the const manager ------------------------------------------------------->
	ARM_CstManagerPtr cstManagerPtr = CreateCstManager();

	//Create the generic security with the const manager. ---------------------------->
    CreateAndSetDealDescriptionAndTimeIt("", pricedColumns, cstManagerPtr, false, otherPayoffs);
	
	//Mkt Data Manager + Model + Calibration ----------------------------------------->
    if (hedgeUpdate) 
    {
       Update(marketDatas);
    }
    else
    {
       Init(marketDatas); 
    }
	CheckData();

	//Free memory
	if (swoptBSmod)
	{
	    delete swoptBSmod;
		swoptBSmod = NULL;
	}
	if (capBSmod)
	{
		delete capBSmod;
		capBSmod = NULL;
	}
/*	if (normalModel)
	{
		delete normalModel;
		normalModel = NULL;
	}*/
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: Copy constructor
///	Returns: a built object
///	Action : builds the object
/////////////////////////////////////////////////////////////////
ARM_BermudaSwaptionCalculator::ARM_BermudaSwaptionCalculator(const ARM_BermudaSwaptionCalculator& rhs)
: ARM_GenCalculator(rhs),
itsCcy(rhs.itsCcy),
itsPastZC(rhs.itsPastZC),
itsPastStartDate(rhs.itsPastStartDate),
itsStartDate(rhs.itsStartDate),
itsEndDate(rhs.itsEndDate),
itsNominal(rhs.itsNominal),
itsStrike(rhs.itsStrike),
itsPayRec(rhs.itsPayRec),
itsStubRule(rhs.itsStubRule),
itsCallFreq(rhs.itsCallFreq),
itsCallNotice(rhs.itsCallNotice),
itsCallCal(rhs.itsCallCal),
itsFirstCallDate(rhs.itsFirstCallDate),
itsLastCallDate(rhs.itsLastCallDate),
itsFees(rhs.itsFees),
itsFixFreq(rhs.itsFixFreq),
itsFixBasis(rhs.itsFixBasis),
itsFixPayCal(rhs.itsFixPayCal),
itsFixAdjRule(rhs.itsFixAdjRule),
itsFixRule(rhs.itsFixRule),
itsFixPayGap(rhs.itsFixPayGap),
itsZC(rhs.itsZC),
itsVarFreq(rhs.itsVarFreq),
itsVarBasis(rhs.itsVarBasis),
itsVarResetCal(rhs.itsVarResetCal),
itsVarPayCal(rhs.itsVarPayCal),
itsVarIndexTerm(rhs.itsVarIndexTerm),
itsVarAdjRule(rhs.itsVarAdjRule),
itsVarRule(rhs.itsVarRule),
itsVarResetGap(rhs.itsVarResetGap),
itsVarPayGap(rhs.itsVarPayGap),
itsVarSpread(rhs.itsVarSpread),
itsGenSecType(rhs.itsGenSecType),
itsControlVariate(rhs.itsControlVariate),
itsControlVariatePrices(rhs.itsControlVariatePrices),
itsMrsCalibFlag(rhs.itsMrsCalibFlag),
itsATMDiagonalFlag(rhs.itsATMDiagonalFlag),
itsNumMethodType(rhs.itsNumMethodType),
itsModelType(rhs.itsModelType),
itsAmcIter(rhs.itsAmcIter),
itsMcIter(rhs.itsMcIter),
itsMaxBucketSize(rhs.itsMaxBucketSize),
itsGenType1(rhs.itsGenType1),
itsGenType2(rhs.itsGenType2),
itsPathOrder(rhs.itsPathOrder),
itsPathScheme(rhs.itsPathScheme),
itsFirstNbTimes(rhs.itsFirstNbTimes),
itsTreeSteps(rhs.itsTreeSteps),
itsPortfolioMode(rhs.itsPortfolioMode),
itsHasBeenPriced(rhs.itsHasBeenPriced),
itsBermudaSwaptionPrice(rhs.itsBermudaSwaptionPrice),
itsBermudaSwaptionStdDev(rhs.itsBermudaSwaptionStdDev),
itsFixBoundaryFlag(rhs.itsFixBoundaryFlag),
itsAtmMrsCalibType(rhs.itsAtmMrsCalibType),
itsApproxMarginFlag(rhs.itsApproxMarginFlag),
itsFreezeBetasFlag(rhs.itsFreezeBetasFlag),
itsBetasCoeff(rhs.itsBetasCoeff),
itsPreCalibMethod(rhs.itsPreCalibMethod),
itsCalculateProbaFlag(rhs.itsCalculateProbaFlag),
itsProbaMatrix(rhs.itsProbaMatrix),
itsCumProba(rhs.itsCumProba),
itsNbCalculatedProbabilities(rhs.itsNbCalculatedProbabilities),
itsDefaultMrs(rhs.itsDefaultMrs),
itsCallDateStrip(rhs.itsCallDateStrip),
itsFundDateStrip(rhs.itsFundDateStrip),
itsStructDateStrip(rhs.itsStructDateStrip)
{
	itsExoSwaption			= (ARM_Swaption*) (rhs.itsExoSwaption)->Clone();
	itsUnderSwap			= (ARM_Swap*) (rhs.itsUnderSwap)->Clone();
	itsControlVariatePrices = new vector<double>(*rhs.itsControlVariatePrices);
	itsControlVariate		= new vector<int>(*rhs.itsControlVariate);
	itsModelParams			= NULL;
	if (rhs.itsModelParams != NULL)
		SetModelParams(rhs.itsModelParams);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: Destructor
///	Returns: void
///	Action : destroys the object
/////////////////////////////////////////////////////////////////
ARM_BermudaSwaptionCalculator::~ARM_BermudaSwaptionCalculator()
{
	delete itsExoSwaption;
	itsExoSwaption = NULL;

	delete itsUnderSwap;
	itsUnderSwap = NULL;

	delete itsControlVariate;
	itsControlVariate = NULL;

	delete itsControlVariatePrices;
	itsControlVariatePrices = NULL;

	if (itsModelParams)
	{
		for (int i=0; i<itsModelParams->size(); i++)
			delete (*itsModelParams)[i];

		delete itsModelParams;
		itsModelParams = NULL;
	}
}

vector<ARM_ReferenceValue*>* ARM_BermudaSwaptionCalculator::GetModelParams() const
{
	return itsModelParams;
}

void ARM_BermudaSwaptionCalculator::SetModelParams(vector<ARM_ReferenceValue*>* modelParams)
{
	if (itsModelParams)
	{
		for (int i=0; i<itsModelParams->size(); i++)
			delete (*itsModelParams)[i];

		delete itsModelParams;
		itsModelParams = NULL;
	}
	itsModelParams = new vector<ARM_ReferenceValue*>(modelParams->size());
	for (int i=0; i<modelParams->size(); i++)
	{
		(*itsModelParams)[i] = (ARM_ReferenceValue*) (*modelParams)[i]->Clone();
	}
}

double ARM_BermudaSwaptionCalculator::GetDefaultMrs() const
{
	return itsDefaultMrs;
}

void ARM_BermudaSwaptionCalculator::SetInitMrs(double initMrs)
{
	delete (*itsModelParams)[1];
	(*itsModelParams)[1] = new ARM_ReferenceValue(initMrs);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: ARM_ExerciseStyle
///	Returns: ARM_ExerciseStyle
///	Action : 
/////////////////////////////////////////////////////////////////
ARM_ExerciseStyle* ARM_BermudaSwaptionCalculator::CreateExerciseStyle() const
{
	ARM_DateStripCombiner calculatedExerciseDates = DatesStructure();	
	////// The Calib Portfolio Construction takes StartDates
	std::vector<double>& resetDates = calculatedExerciseDates.GetDateStrip(CALL_BS_SCHED)->GetFlowStartDates();
	std::vector<double>& startDates = calculatedExerciseDates.GetDateStrip(CALL_BS_SCHED)->GetFlowStartDates();
	//// Conversion to ARM_Vector
	size_t sizeVector = resetDates->size();
	ARM_Vector resetDates_bis(sizeVector);
	ARM_Vector startDates_bis(sizeVector);
	for (size_t j = 0; j<sizeVector; ++j)
	{
		resetDates_bis[j] = (*resetDates)[j];
		startDates_bis[j] = (*startDates)[j];
	}
	ARM_ExerciseStyle* exoDates = new ARM_ExerciseStyle(&resetDates_bis,&startDates_bis);
	exoDates->SetExerciseType(K_BERMUDAN);
	return exoDates;
}
/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: ComputeMeanNotionalCurve
///	Returns: double
///	Action : Convert the Notional to the standard Frequency starting from the first exercise date
/// This notional will be used with the std portfolio when we calibrate the mean reversion
/////////////////////////////////////////////////////////////////
ARM_ReferenceValue*	ARM_BermudaSwaptionCalculator::ComputeMeanNotionalCurve(ARM_ReferenceValue& notionalCurve,int newFreq )
{	
	ARM_DateStrip FixSched(GetFirstCallDate(), GetEndDate(),GetFixFreq(),GetFixBasis(),
							GetVarResetCal().c_str(), GetFixAdjRule(), GetFixRule(), GetStubRule(), GetVarResetGap(),
							GetFixFreq(), GETDEFAULTVALUE, GetFixPayCal().c_str(), K_ADVANCE, K_ARREARS, true);	
	ARM_DateStrip newSched(GetFirstCallDate(), GetEndDate(),newFreq,GetFixBasis(),
							GetVarResetCal().c_str(), GetFixAdjRule(), GetFixRule(), GetStubRule(), GetVarResetGap(),
							newFreq, GETDEFAULTVALUE, GetFixPayCal().c_str(), K_ADVANCE, K_ARREARS, true);

	std::vector<double>& endDates = FixSched.GetFlowEndDates();
	std::vector<double>& newDates = newSched.GetFlowEndDates();

	int times = (int) (GetFixFreq()/newFreq) ;
	if (times<=1)
		return (ARM_ReferenceValue*) (notionalCurve.Clone());

	int sizeVector = newDates->size();
	ARM_Vector*	values = new ARM_Vector(sizeVector,0.0);
	ARM_Vector*	dates  = new ARM_Vector(sizeVector);

	int index = 0;
	for(int j=0;j<sizeVector;++j)
	{
		dates->Elt(j) = newDates->Elt(j);
		for(int i=0; (i<times && index<endDates->size());++i)
		{
			values->Elt(j)+=notionalCurve.Interpolate(endDates->Elt(index));
			index++;
		}
		values->Elt(j)/=times;
	}
	ARM_ReferenceValue*	result = new ARM_ReferenceValue(dates,values);
	result->SetCalcMethod(notionalCurve.GetCalcMethod());
	return result;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: ConvertVariableMarginTo Strike
///	Returns: double
///	Action : Convert the Variable Margen RefValue to a strike (for the creation of the calib portfolio)
/// WARNING :Function to be called only for Bermuda ZC options
///			 The margin curve have to be converted to the fix freq first
/////////////////////////////////////////////////////////////////
ARM_ReferenceValue* ARM_BermudaSwaptionCalculator::ConvertVariableMarginToStrike(ARM_ReferenceValue& marginCurve )
{
	/// Only Constant Strikes are Valid for ZC Bermudan Option
	double strikeCst = GetStrike().GetDiscreteValues()->Elt(0);

	//// Transform the margin Curve to the new frequency
	ARM_DateStrip FixSched(GetStartDate(), GetEndDate(),GetFixFreq(),GetFixBasis(),
							GetVarResetCal().c_str(), GetFixAdjRule(), GetFixRule(), GetStubRule(), GetVarResetGap(),
							GetFixFreq(), GETDEFAULTVALUE, GetFixPayCal().c_str(), K_ADVANCE, K_ARREARS, true);
	

	std::vector<double>& resetDates = FixSched.GetResetDates();
	std::vector<double>& startDates = FixSched.GetFlowStartDates();
	std::vector<double>& endDates = FixSched.GetFlowEndDates();

	size_t FixScheduleSize = FixSched.GetResetDates()->size();
	ARM_Vector FactorVector(FixScheduleSize,1.0);
	for(size_t i = 1 ; i<FixScheduleSize;i++)
	{
		FactorVector.Elt(i) = FactorVector.Elt(i-1)*1/(1+strikeCst/GetFixFreq());
	}

	ARM_ZeroCurve* pCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));

	ARM_Vector	marginFlows(FixScheduleSize,1.0);
	ARM_Vector*	addStrikes = new ARM_Vector(FixScheduleSize,1.0);
	ARM_Vector*	dates	   = new ARM_Vector(FixScheduleSize,1.0);
	for(i = 0 ; i<FixScheduleSize;i++)
	{
		addStrikes->Elt(i) = -marginCurve.Interpolate(resetDates->Elt(i))*FactorVector.Elt(i)+strikeCst;
		dates->Elt(i) = resetDates->Elt(i);
	}
	ARM_ReferenceValue*	result = new ARM_ReferenceValue(dates,addStrikes);
	result->SetCalcMethod(marginCurve.GetCalcMethod());
	return result;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: ConvertVariableMargin
///	Returns: ARM_ReferenceValue
///	Action : Convert the Variable Margen RefValue to the new frequency
/////////////////////////////////////////////////////////////////
//// A revoir
ARM_ReferenceValue* ARM_BermudaSwaptionCalculator::ConvertVariableMargin(ARM_ReferenceValue& marginCurve,int initFreq, int newFreq, int varBasis , int fixBasis )
{
	//Float Leg schedule
    ARM_DateStrip VarSched(GetStartDate(), GetEndDate(),initFreq,varBasis,
							GetVarResetCal().c_str(), GetVarAdjRule(), GetVarAdjRule(), GetStubRule(), GetVarResetGap(),
							initFreq, GETDEFAULTVALUE, GetVarPayCal().c_str(), K_ADVANCE, K_ARREARS, true);

	
	//// Transform the margin Curve to the new frequency
	ARM_DateStrip FixSched(GetStartDate(), GetEndDate(),newFreq,fixBasis ,
							GetVarResetCal().c_str(), GetFixAdjRule(), GetFixRule(), GetStubRule(), GetVarResetGap(),
							newFreq, GETDEFAULTVALUE, GetFixPayCal().c_str(), K_ADVANCE, K_ARREARS, true);


	size_t FixScheduleSize = FixSched.GetResetDates()->size();
	ARM_Vector* newDates = new ARM_Vector(FixScheduleSize,0.0);
	ARM_Vector* values = new ARM_Vector(FixScheduleSize,0.0);
	
	//values->Elt(0) = marginCurve.Interpolate(FixSched.GetResetDates()->Elt(0));
	//newDates->Elt(0) =  FixSched.GetResetDates()->Elt(0);
	/////	Interpolate the margin values
	int varIndex = 0;
	ARM_ZeroCurve* pCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));

	double maturityTime, DF, m, term, newTerm;
	
	for (size_t i=0; i<FixScheduleSize;i++)
	{
		double marginFlow = 0;
		/// 7 days for adjustments
		while((varIndex<VarSched.GetResetDates()->size()) && (((ARM_Date)(VarSched.GetFlowEndDates()->Elt(varIndex))).GetJulian()<=((ARM_Date) (FixSched.GetFlowEndDates()->Elt(i))).GetJulian()+7) )
		{
			/// Margin are paid each end period
			maturityTime = ((ARM_Date) (VarSched.GetFlowEndDates()->Elt(varIndex))).GetJulian() - GetMktDataManager()->GetAsOfDate().GetJulian();
			DF = pCurve->DiscountPrice(maturityTime/K_YEAR_LEN);
			term = VarSched.GetInterestTerms()->Elt(varIndex);
			m = marginCurve.Interpolate(VarSched.GetResetDates()->Elt(varIndex));
			marginFlow += m*DF*term;
			varIndex++;
		}
		newTerm = FixSched.GetInterestTerms()->Elt(i);
		values->Elt(i) = marginFlow/DF/newTerm;
		newDates->Elt(i) =  FixSched.GetResetDates()->Elt(i);
	}
	ARM_ReferenceValue* result =  new ARM_ReferenceValue(newDates,values);
	/// newDates and values are not cloned 
	result->SetCalcMethod(marginCurve.GetCalcMethod());
	return result;
}


////////////////////////////////////////////////////
///	Class   : ARM_BermudaSwaptionCalculator
///	Routines: Clone,View
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_BermudaSwaptionCalculator::Clone() const
{
	return new ARM_BermudaSwaptionCalculator(*this);
}

const int ARM_BermudaSwaptionCalculator::GetGenSecType() const
{
	return itsGenSecType;
}

void ARM_BermudaSwaptionCalculator::SetGenSecType(int genSecType)
{
	itsGenSecType = genSecType;
}

void ARM_BermudaSwaptionCalculator::UpdateGenSecType(int genSecType)
{
	//Zero-Coupon bermuda
	if (genSecType == 4)
		SetZc(true);
	
	//Default bermuda
	if (genSecType == -1)
	{
		if ((GetNumMethodType() == Tree) || (GetNumMethodType() == PDE))
		{
			SetGenSecType(2);
		}
		else if (GetNumMethodType() == MonteCarlo)
		{
			SetGenSecType(1);
		}
		if (GetZCFlag() == true)
		{
			SetGenSecType(4);
		}
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: GetBeta & SetBeta
///	Returns: 
///	Action : get & set the Beta param for SUMMIT interface
/////////////////////////////////////////////////////////////////
const ARM_ModelParam& ARM_BermudaSwaptionCalculator::GetBeta() const
{
    return *(static_cast< ARM_ModelParam* >(GetMktDataManager()->GetData(GetKeys()[BetaKey])));
}

void ARM_BermudaSwaptionCalculator::SetBeta(ARM_ModelParam* betaParam)
{
    if(!betaParam || betaParam->GetType() != ARM_ModelParamType::Beta)
		    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : an Beta Param is expected !" );

    GetMktDataManager()->RegisterData(GetKeys()[BetaKey],static_cast< ARM_Object* >(betaParam));
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: GetIndexType
///	Returns: the index type
///	Action : get the ARM index type of the coupon index
/////////////////////////////////////////////////////////////////
ARM_INDEX_TYPE ARM_BermudaSwaptionCalculator::GetIndexType()
{
    string liborTypeName(string(GetCcy().GetCcyName()) == string("EUR") ? "EURIBOR" : "LIBOR");
    string varIndexTerm(itsVarIndexTerm);
    if(varIndexTerm=="12M")
        varIndexTerm="1Y"; // because we cant convert EURIBOR1Y and not EURIBOR12M
    liborTypeName += varIndexTerm;
    return static_cast< ARM_INDEX_TYPE > (ARM_ArgConv_IndexType.GetNumber(liborTypeName));
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: CheckSwaptionInputs
///	Returns: void
///	Action : check if Bermuda Swaption inputs are consistent
/////////////////////////////////////////////////////////////////
void ARM_BermudaSwaptionCalculator::CheckSwaptionInputs()
{	
	if(GetZCFlag() && ((GetGenSecType()!=4)&&(GetGenSecType()!=-1)))
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : ZCFlag is activated. Please choose GenSecType 4 or -1");
	}

	if(GetZCFlag() && (GetFixRule()!=K_UNADJUSTED))
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Fix Rule Must be Unadjusted");
	}

	if ((GetFixRule() == K_ADJUSTED)&&(GetGenSecType() == 4))
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The FixRule must be UNADJ for ZC Bermudan Swaptions.");
	}
	//Dates.
	if ( GetEndDate() < GetStartDate())
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The Start Date of the deal has to be before the end date.");
	}

	if ( GetLastCallDate() < GetFirstCallDate())
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The first call Date of the deal has to be before the last call date.");
	}

	if ( GetLastCallDate() > GetEndDate())
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The last call Date of the deal has to be before the end date.");
	}

	if (GetCallFreq() > GetVarFreq())
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The call frequency cannot be greater than variable frequency.");
	}

	if (GetZCFlag() && !GetApproxMarginFlag())
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The ZC bermuda should be priced with the margin approximation.");
	}


	//Gen/NumMethod.
	if (GetGenSecType() == 2)
	{
		if (!(GetNumMethodType() == Tree || GetNumMethodType() == PDE))
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : You must price using a tree method.");
	}
	if ((GetNbControlVariate() !=0) && (GetNumMethodType() == Tree || GetNumMethodType() == PDE))
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : You must use a MC method when using control variates.");

	//to fullfill when necessary.

}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: CheckMktData & CheckData
///	Returns: void
///	Action : check if Bermuda Swaption datas are consistent
/////////////////////////////////////////////////////////////////
void ARM_BermudaSwaptionCalculator::CheckMktData()
{
    /// MdM datas checking
	ARM_ZeroCurve* ycCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));
    if(!ycCurve)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Yield Curve for key=" + GetKeys()[YcKey] + " is expected in the Market Data Manager");
    string Ccy(ycCurve->GetCurrencyUnit()->GetCcyName());

    ARM_BSModel* oswBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswModelKey]) );
    if(!oswBSModel)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : swaption B&S model for key=" + GetKeys()[OswModelKey] + " is expected in the Market Data Manager");

	ARM_BSModel* cfBSModel = dynamic_cast< ARM_BSModel* >(GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );
    if(!cfBSModel)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : cap/floor B&S model for key=" + GetKeys()[CfModelKey] + " is expected in the Market Data Manager");

	/*ARM_ModelParam* mrsParam = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[MrsKey]));
    if(!mrsParam || mrsParam->GetType() != ARM_ModelParamType::MeanReversion)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : MRS Param for key=" + GetKeys()[MrsKey] + " is expected in the Market Data Manager");*/
}


void ARM_BermudaSwaptionCalculator::CheckData()
{
	CheckMktData();
}

//////////////////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: PricedColumnNames
///	Returns: ARM_StringVector
///	Action : create the control variates column names of the deal description
//////////////////////////////////////////////////////////////////////////////
ARM_StringVector ARM_BermudaSwaptionCalculator::PricedColumnNames() const
{
	int genSecType = GetGenSecType();

	ARM_StringVector pricedColumns;

	if((genSecType == 1) || (genSecType == 3))
	{
		size_t nbCv = (*itsControlVariate).size();

		pricedColumns.resize(nbCv+1);
		ARM_StringVector cv = ControlVariateColumnNames();
		pricedColumns[0] = "Bermuda";
		for (size_t i=0; i<nbCv;i++)
		{
			pricedColumns[i+1] = cv[i];
		}
	}
	else if ((genSecType == 2) || (genSecType == 4))
	{
		/*int addColumns =0;
		if(itsCalculateProbaFlag && (genSecType == 2) )
			addColumns = itsNbCalculatedProbabilities;
		pricedColumns.resize(1+addColumns);
		pricedColumns[0] = "Bermuda";
		size_t j=3;
		for(size_t p=1;p<=addColumns;p++)
		{
			pricedColumns[p] = (BermudaSwaptionProbaColNamesTable2[j].c_str());
			j+=2;
		}*/
		int addColumns =0;
		if(itsCalculateProbaFlag && (genSecType == 2) )
			addColumns = itsNbCalculatedProbabilities;
		pricedColumns.resize(2+addColumns);
		pricedColumns[0] = "Bermuda";
		pricedColumns[1] = "Frontier";
		size_t j=3;
		for(size_t p=0;p<addColumns;p++)
		{
			pricedColumns[p+2] = (BermudaSwaptionProbaColNamesTable2[j].c_str());
			j+=2;
		}
	}

	return pricedColumns;
}

//////////////////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: ControlVariatesColumnNames
///	Returns: ARM_StringVector
///	Action : create the control variates column names of the deal description
//////////////////////////////////////////////////////////////////////////////
ARM_StringVector ARM_BermudaSwaptionCalculator::ControlVariateColumnNames() const
{
	//Generic Security Type:
	int genSecType = GetGenSecType();

	ARM_StringVector columnNames(0);
	size_t ctrlVarSize = (*itsControlVariate).size(); 
	size_t i;

	if ((genSecType == 1) || (genSecType == 3))
	{
		columnNames.resize(ctrlVarSize); //control variate
		for(i = 0; i< ctrlVarSize; i++)
		{
			columnNames[i] = ControlVariateColNamesTable1[i] ;
		}
	}
	return columnNames;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: ColumnNames
///	Returns: ARM_RowInfo
///	Action : create the column names of the deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_BermudaSwaptionCalculator::ColumnNames() const
{
	//Generic Security Type:
	int genSecType = GetGenSecType();

	//Number of Columns
	size_t colNamesInitSize;
	size_t colNamesSize;
	size_t ctrlVarSize = (*itsControlVariate).size();

	if (genSecType == 1)
	{
		colNamesInitSize = sizeof(BermudaSwaptionColNamesTable1)/sizeof(BermudaSwaptionColNamesTable1[0]);
		colNamesSize	 = colNamesInitSize + ctrlVarSize;
	}
	else if (genSecType == 2)
	{
		colNamesInitSize = sizeof(BermudaSwaptionColNamesTable2)/sizeof(BermudaSwaptionColNamesTable2[0]);
		colNamesSize     = colNamesInitSize;
		if (itsCalculateProbaFlag )
		{
			/// To be changed to take into account the number of lines
			colNamesSize += 4 + 2*(itsNbCalculatedProbabilities-1); 
		}
	}
	else if (genSecType == 3)
	{
		colNamesInitSize = sizeof(BermudaSwaptionColNamesTable3)/sizeof(BermudaSwaptionColNamesTable3[0]);
		colNamesSize = colNamesInitSize + 2*ctrlVarSize;
	}
	else if (genSecType == 4)
	{
		colNamesInitSize = sizeof(BermudaSwaptionColNamesTable4)/sizeof(BermudaSwaptionColNamesTable4[0]);
		colNamesSize     = colNamesInitSize;
	}

    vector< string > colNamesVec(colNamesSize);
    vector< ARM_GP_VALUE_TYPE > colTypeVec(colNamesSize, ARM_STRING); 

	//Standard columns
    for(size_t i=0;i<colNamesInitSize; ++i)
	{
		if (genSecType == 1)
		{
			colNamesVec[i] = BermudaSwaptionColNamesTable1[i];
		}
		else if (genSecType == 2)
		{
			colNamesVec[i] = BermudaSwaptionColNamesTable2[i];
		}
		else if (genSecType == 3)
		{
			colNamesVec[i] = BermudaSwaptionColNamesTable3[i];
		}
		else if (genSecType == 4)
		{
			colNamesVec[i] = BermudaSwaptionColNamesTable4[i];
		}
	}
	//Control Variate columns or Calculate probabilities
	for(i=colNamesInitSize; i<colNamesSize; ++i)
	{
		if (genSecType == 1)
		{
			colNamesVec[i] = ControlVariateColNamesTable1[i-colNamesInitSize] ;
		}
		else if (genSecType == 3)
		{
			colNamesVec[i]	 = ControlVariateColNamesTable3[i-colNamesInitSize];
			colNamesVec[i+1] = ControlVariateColNamesTable3[i+1-colNamesInitSize];
			i++;
		}
		else if (genSecType == 2)
		{
			colNamesVec[i]	 = BermudaSwaptionProbaColNamesTable2[i-colNamesInitSize];
		}
	}

    ARM_RowInfo rowInfo(colNamesVec,colTypeVec);

    return rowInfo;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: MiddleRows
///	Returns: ARM_RowInfo
///	Action : create a row of a deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_BermudaSwaptionCalculator::MiddleRows( size_t eventIdx, const ARM_DateStripCombiner& datesStructure ) const
{
	//Generic Security Type
	int genSecType = GetGenSecType();
	
	//Number of columns
	size_t descSize; 
	size_t ctrlVarSize = (*itsControlVariate).size();
	if (genSecType == 1) 
	{
		descSize = sizeof(BermudaSwaptionColNamesTable1)/sizeof(BermudaSwaptionColNamesTable1[0]);
		descSize += ctrlVarSize;
	}
	else if (genSecType == 2)
	{
		descSize = sizeof(BermudaSwaptionColNamesTable2)/sizeof(BermudaSwaptionColNamesTable2[0]);
		if(itsCalculateProbaFlag)
			descSize += 4 + 2*(itsNbCalculatedProbabilities-1);
	
	}
	else if (genSecType == 3)
	{
		descSize = sizeof(BermudaSwaptionColNamesTable3)/sizeof(BermudaSwaptionColNamesTable3[0]);
		descSize = descSize +2*ctrlVarSize;
	}
	else if (genSecType == 4)
	{
		descSize = sizeof(BermudaSwaptionColNamesTable4)/sizeof(BermudaSwaptionColNamesTable4[0]);
	}

	size_t callEventSize = datesStructure.GetDateStrip(CALL_BS_SCHED)->GetResetDates()->size();
    vector< string > rowDescVec(descSize);
    vector< ARM_GP_VALUE_TYPE > rowTypeVec(descSize, ARM_MISSING_TYPE); 
	
	//Set default 0 value for each column to be able to sum it
    InitPriceableColumns(rowDescVec,rowTypeVec);
	
	//Useful datas
	int intPayRec = itsPayRec;
	string payRec; 
	string invPayRec;
	if(intPayRec == -1)
	{
		payRec = "P";
		invPayRec = "R";
	}
	else
	{
		payRec = "R";
		invPayRec = "P";
	}
	char* ccy = GetCcy().GetCcyName();
	int intFixFreq = itsFixFreq;
	string fixFreq = ARM_ArgConvReverse_MatFrequency.GetString(intFixFreq);
	int intFixBasis = itsFixBasis;
	string fixBasis = ARM_ArgConvReverse_DayCount.GetString(intFixBasis);
	int intVarFreq = itsVarFreq;
	string varFreq = ARM_ArgConvReverse_MatFrequency.GetString(intVarFreq);
	int intVarBasis = itsVarBasis;
	string varBasis = ARM_ArgConvReverse_DayCount.GetString(intVarBasis);

	int callFreq = itsCallFreq;
	int callAdjRule = itsVarAdjRule;
	int fixRule = itsFixRule;
	int fixAdjRule = itsFixAdjRule;

	if (genSecType == 1) // AMC
	{
		//1)RESETDATE description (Notice date)
		double noticeDate = (*(datesStructure.GetDateStrip(CALL_BS_SCHED)->GetResetDates()))[eventIdx] ;
		CC_Ostringstream noticeDateDesc;
		noticeDateDesc << CC_NS(std,fixed) << noticeDate;
		rowDescVec[ResetDate1] = noticeDateDesc.str();
		rowTypeVec[ResetDate1] = ARM_DATE_TYPE;

		//2)STARTDATE description (Call date)
		double callDate = (*(datesStructure.GetDateStrip(CALL_BS_SCHED)->GetFlowStartDates()))[eventIdx] ;
		CC_Ostringstream callDateDesc;
		callDateDesc << CC_NS(std,fixed) << callDate;
		rowDescVec[StartDate1] = callDateDesc.str();
		rowTypeVec[StartDate1] = ARM_DATE_TYPE;

		//3)PAYDATE description (Underlying End date)
		double payDate = itsEndDate.GetJulian(); 
		CC_Ostringstream payDateDesc;
		payDateDesc << CC_NS(std,fixed) << payDate;
		rowDescVec[PayDate1] = payDateDesc.str();
		rowTypeVec[PayDate1] = ARM_DATE_TYPE;
	
		//4)SWAPSTARTDATE 
		double startDate = itsStartDate.GetJulian(); 
		double swapStartDate=0;
		if(startDate < callDate)
		{
			swapStartDate = callDate;
		}
		else
		{
			swapStartDate = startDate;
		}

		CC_Ostringstream swapStartDateDesc;
		swapStartDateDesc << CC_NS(std,fixed) << swapStartDate;
		rowDescVec[SwapStartDate1] = swapStartDateDesc.str();
		rowTypeVec[SwapStartDate1] = ARM_DATE_TYPE;

		// To handle case CallFreq > FixFreq
		int nbCallPeriods = CC_Round(((itsEndDate.GetJulian()-swapStartDate)/K_YEAR_LEN)*callFreq);
		int nbFixPeriods = ceil((double)nbCallPeriods/callFreq*intFixFreq);
		if ((intFixFreq < callFreq) && (nbCallPeriods%callFreq != 0))
		{
			ARM_Date nextSwapStartDate(itsEndDate.GetJulian());
			nextSwapStartDate.AddPeriodMult(intFixFreq, -nbFixPeriods/intFixFreq+1);

			if (fixRule == K_ADJUSTED)
			{
				nextSwapStartDate.AdjustToBusDate(ccy,fixAdjRule);
			}

			//5)NEXTSWAPSTARTDATE 
			CC_Ostringstream nextSwapStartDateDesc;
			nextSwapStartDateDesc << CC_NS(std,fixed) << nextSwapStartDate.GetJulian();
			rowDescVec[NextSwapStartDate2] = nextSwapStartDateDesc.str();
			rowTypeVec[NextSwapStartDate2] = ARM_DATE_TYPE;

			ARM_Date lastSwapStartDate(itsEndDate);
			lastSwapStartDate.AddPeriodMult(intFixFreq, -nbFixPeriods/intFixFreq);

			if (fixRule == K_ADJUSTED)
			{
				lastSwapStartDate.AdjustToBusDate(ccy,fixAdjRule);
			}
			
			//6)LASTSWAPSTARTDATE 
			CC_Ostringstream lastSwapStartDateDesc;
			lastSwapStartDateDesc << CC_NS(std,fixed) << lastSwapStartDate.GetJulian();
			rowDescVec[LastSwapStartDate2] = lastSwapStartDateDesc.str();
			rowTypeVec[LastSwapStartDate2] = ARM_DATE_TYPE;

			ARM_Date resetFromCallDate(lastSwapStartDate);
			int defaultResetGap = GetCurrencyUnit()->GetSpotDays();
			resetFromCallDate.GapBusinessDay(-defaultResetGap, const_cast<char *>(GetVarResetCal().c_str()));

			//7)STRIKE
			double strike = const_cast< ARM_ReferenceValue& >(itsStrike).Interpolate(resetFromCallDate.GetJulian());
			CC_Ostringstream strikeDesc;
			strikeDesc << CC_NS(std,fixed) << strike;
			rowDescVec[Strike1] = strikeDesc.str();
			rowTypeVec[Strike1] = ARM_DOUBLE;

			// NOTIONAL
			CC_Ostringstream notionalDesc;
			double notional = const_cast< ARM_ReferenceValue& >(itsNominal).Interpolate(nextSwapStartDate.GetJulian());
			notionalDesc << CC_NS(std,fixed) << notional;
			rowDescVec[Notional1] = notionalDesc.str();
			rowTypeVec[Notional1] = ARM_DOUBLE;
			
			//9)STDFIXSWAPLET
			CC_Ostringstream stdFixSwapletDesc;
			stdFixSwapletDesc << "(";
			if (intPayRec == K_RCV)
			{
				stdFixSwapletDesc << "-";
			}
				 
			stdFixSwapletDesc << "DCF(" << BermudaSwaptionColNamesTable1[LastSwapStartDate1] << "[i],"; 
			stdFixSwapletDesc << BermudaSwaptionColNamesTable1[SwapStartDate1] << "[i],";
			stdFixSwapletDesc << fixBasis << ")*";
			stdFixSwapletDesc << BermudaSwaptionColNamesTable1[Strike1] << "[i]*";
			stdFixSwapletDesc << "DF(" << ccy << "," << BermudaSwaptionColNamesTable1[SwapStartDate1] << "[i])";

			if (intPayRec == K_PAY)
			{
				stdFixSwapletDesc << "-";
			}
			else
			{
				stdFixSwapletDesc << "+";
			}

			stdFixSwapletDesc << "DCF(" << BermudaSwaptionColNamesTable1[LastSwapStartDate1] << "[i],"; 
			stdFixSwapletDesc << BermudaSwaptionColNamesTable1[NextSwapStartDate1] << "[i],";
			stdFixSwapletDesc << fixBasis << ")*";
			stdFixSwapletDesc << BermudaSwaptionColNamesTable1[Strike1] << "[i]*";
			stdFixSwapletDesc << "DF(" << ccy << "," << BermudaSwaptionColNamesTable1[NextSwapStartDate1] << "[i]))*";
			stdFixSwapletDesc << BermudaSwaptionColNamesTable1[Notional1] << "[i]";

			rowDescVec[StdFixSwaplet1] = stdFixSwapletDesc.str();
			rowTypeVec[StdFixSwaplet1] = ARM_STRING;

			CC_Ostringstream stdVarSwapletDesc;

			//10)STDVARSWAPLET
			stdVarSwapletDesc << "Swap(" << ccy << ",";
			stdVarSwapletDesc << BermudaSwaptionColNamesTable1[SwapStartDate1] << "[i],";
			stdVarSwapletDesc << BermudaSwaptionColNamesTable1[NextSwapStartDate1] << "[i], 0,";
			stdVarSwapletDesc << payRec << ",";
			stdVarSwapletDesc << fixFreq << ", " << fixBasis << ", " << varFreq << ", " << varBasis << ", ";
			stdVarSwapletDesc << "MarginCurve, NotionalCurve" << ")";

			rowDescVec[StdVarSwaplet1] = stdVarSwapletDesc.str();
			rowTypeVec[StdVarSwaplet1] = ARM_STRING;

			//11)STDSWAPLET
			CC_Ostringstream stdSwapletDesc;
			stdSwapletDesc << BermudaSwaptionColNamesTable1[StdFixSwaplet1] << "[i]+";
			stdSwapletDesc << BermudaSwaptionColNamesTable1[StdVarSwaplet1] << "[i]";

			rowDescVec[StdSwaplet1] = stdSwapletDesc.str();
			rowTypeVec[StdSwaplet1] = ARM_STRING;

			//12)STDSWAP
			CC_Ostringstream stdSwapDesc;
			stdSwapDesc << BermudaSwaptionColNamesTable1[StdSwaplet1] << "[i]";
			if (nbCallPeriods > callFreq)
			{
				stdSwapDesc << "+Swap(" << ccy << ",";
				stdSwapDesc << BermudaSwaptionColNamesTable1[NextSwapStartDate1] << "[i],";
				stdSwapDesc << BermudaSwaptionColNamesTable1[PayDate1] << "[i], StrikesCurve,";
				stdSwapDesc << payRec << ",";
				stdSwapDesc << fixFreq << ", " << fixBasis << ", " << varFreq << ", " << varBasis << ", ";
				stdSwapDesc << "MarginCurve, NotionalCurve" << ")";
			}

			rowDescVec[StdSwap1] = stdSwapDesc.str();
			rowTypeVec[StdSwap1] = ARM_STRING;
		}
		else
		{
			ARM_Date resetFromCallDate(callDate);
			int defaultResetGap = GetCurrencyUnit()->GetSpotDays();
			resetFromCallDate.GapBusinessDay(-defaultResetGap, const_cast<char *>(GetVarResetCal().c_str()));

			//7)STRIKE
			double strike = const_cast< ARM_ReferenceValue& >(itsStrike).Interpolate(resetFromCallDate.GetJulian());
			CC_Ostringstream strikeDesc;
			strikeDesc << CC_NS(std,fixed) << strike;
			rowDescVec[Strike1] = strikeDesc.str();
			rowTypeVec[Strike1] = ARM_DOUBLE;

			//1)STDSwap
			CC_Ostringstream stdSwapDesc;
			stdSwapDesc << "Swap(" << ccy << ",";
			stdSwapDesc << BermudaSwaptionColNamesTable1[SwapStartDate1] << "[i],";
			stdSwapDesc << BermudaSwaptionColNamesTable1[PayDate1] << "[i], StrikesCurve,";
			stdSwapDesc << payRec << ",";
			stdSwapDesc << fixFreq << ", " << fixBasis << ", " << varFreq << ", " << varBasis << ", ";
			stdSwapDesc << "MarginCurve, NotionalCurve" << ")";

			rowDescVec[StdSwap1] = stdSwapDesc.str();
			rowTypeVec[StdSwap1] = ARM_STRING;
		}

		//8)FEES
		CC_Ostringstream feesDesc;
		double fees = const_cast< ARM_ReferenceValue& >(itsFees).Interpolate(noticeDate);
		feesDesc << CC_NS(std,fixed) << fees;
		rowDescVec[Fees1] = feesDesc.str();
		rowTypeVec[Fees1] = ARM_DOUBLE;


		//13)OPTION
		CC_Ostringstream optionDesc;
		if(eventIdx == callEventSize-1)
		{
			optionDesc << "MAX(";
			optionDesc << BermudaSwaptionColNamesTable1[StdSwap1] << "[i]-";
			optionDesc << BermudaSwaptionColNamesTable1[Fees1] << "[i],0)";
		}
		else
		{
			optionDesc << "Exercise(0,";
			optionDesc << BermudaSwaptionColNamesTable1[StdSwap1] << "[i]-";
			optionDesc << BermudaSwaptionColNamesTable1[Fees1] << "[i]," ;
			optionDesc << BermudaSwaptionColNamesTable1[Option1] << "[i+1])";
		}
		rowDescVec[Option1] = optionDesc.str();
		rowTypeVec[Option1] = ARM_STRING;

		//14)BERMUDA
		CC_Ostringstream bermudaDesc;
		if(eventIdx == 0)
		{
			if (itsExoSwaption->GetPorS() == K_PAY)
				bermudaDesc << "-";

			bermudaDesc << BermudaSwaptionColNamesTable1[Option1] << "[i]";
			rowDescVec[Bermuda1] = bermudaDesc.str();
			rowTypeVec[Bermuda1] = ARM_STRING;
		}
		else
		{
			bermudaDesc << "0";
			rowDescVec[Bermuda1] = bermudaDesc.str();
			rowTypeVec[Bermuda1] = ARM_DOUBLE;
		}
		
		//15)CONTROLVARIATES
		size_t currCtrlVar;
		size_t i;
		size_t nbCol =  sizeof(BermudaSwaptionColNamesTable1)/sizeof(BermudaSwaptionColNamesTable1[0]);
		size_t currCol;
		if (ctrlVarSize != 0) 
		{
			for (i=0; i<ctrlVarSize;i++)
			{
				CC_Ostringstream controlVariateDesc;
				currCtrlVar = (*itsControlVariate)[i];
				
				currCol = Bermuda1 + i + 1;
				if(eventIdx == currCtrlVar)
				{
				
					controlVariateDesc << "MAX(";
					controlVariateDesc << BermudaSwaptionColNamesTable1[StdSwap1] << "[i],0)";
					rowDescVec[currCol] = controlVariateDesc.str();
					rowTypeVec[currCol] = ARM_STRING;
				}
				else
				{
					controlVariateDesc << "0";
					rowDescVec[currCol] = controlVariateDesc.str();
					rowTypeVec[currCol] = ARM_DOUBLE;
				}
			}

		}
	}
	else if (genSecType == 2) // TREE or PDE
	{
		//1)RESETDATE description (Notice date)
		double noticeDate = (*(datesStructure.GetDateStrip(CALL_BS_SCHED)->GetResetDates()))[eventIdx] ;
		CC_Ostringstream noticeDateDesc;
		noticeDateDesc << CC_NS(std,fixed) << noticeDate;
		rowDescVec[ResetDate2] = noticeDateDesc.str();
		rowTypeVec[ResetDate2] = ARM_DATE_TYPE;

		 //2)STARTDATE description (Call date)
		double callDate = (*(datesStructure.GetDateStrip(CALL_BS_SCHED)->GetFlowStartDates()))[eventIdx];
		ARM_Date armCallDate = ARM_Date(callDate);
		CC_Ostringstream callDateDesc;
		callDateDesc << CC_NS(std,fixed) << callDate;
		rowDescVec[StartDate2] = callDateDesc.str();
		rowTypeVec[StartDate2] = ARM_DATE_TYPE;

		//3)ENDDATE description 
		CC_Ostringstream endDateDesc;
		double tempEndDate = (*(datesStructure.GetDateStrip(CALL_BS_SCHED)->GetFlowEndDates()))[eventIdx] ;
		double endDate;
		if (eventIdx == callEventSize-1)
		{
			endDate = itsEndDate.GetJulian();
		}
		else
		{
			endDate = tempEndDate; 
		}
		endDateDesc << CC_NS(std,fixed) << endDate;
		rowDescVec[EndDate2] = endDateDesc.str();
		rowTypeVec[EndDate2] = ARM_DATE_TYPE;

		//4)SWAPSTARTDATE description
		double startDate = itsStartDate.GetJulian(); 
		double swapStartDate=0;
		if(startDate < callDate)
		{
			swapStartDate = callDate;
		}
		else
		{
			swapStartDate = startDate;
		}
		CC_Ostringstream swapStartDateDesc;
		swapStartDateDesc << CC_NS(std,fixed) << swapStartDate;
		rowDescVec[SwapStartDate2] = swapStartDateDesc.str();
		rowTypeVec[SwapStartDate2] = ARM_DATE_TYPE;

		//5bis)DEAL ENDDATE
		ARM_Date DealEndDate(itsEndDate);
		CC_Ostringstream dealEndDateDesc;
		dealEndDateDesc << CC_NS(std,fixed) << DealEndDate.GetJulian();
		rowDescVec[DealEndDate2] = dealEndDateDesc.str();
		rowTypeVec[DealEndDate2] = ARM_DATE_TYPE;

		//6)FEES
		CC_Ostringstream feesDesc;
		double fees = const_cast< ARM_ReferenceValue& >(itsFees).Interpolate(noticeDate);
		feesDesc << CC_NS(std,fixed) << fees;
		rowDescVec[Fees2] = feesDesc.str();
		rowTypeVec[Fees2] = ARM_DOUBLE;

		// To handle case CallFreq > FixFreq
		int nbCallPeriods = CC_Round(((itsEndDate.GetJulian()-swapStartDate)/K_YEAR_LEN)*callFreq);
		int nbFixPeriods = ceil((double)nbCallPeriods/callFreq*intFixFreq);
		if (intFixFreq < callFreq)
		{
			//5)NEXTSWAPSTARTDATE
			ARM_Date nextSwapStartDate(itsEndDate.GetJulian());
			nextSwapStartDate.AddPeriodMult(intFixFreq, -nbFixPeriods/intFixFreq+1);

			if (fixRule == K_ADJUSTED)
			{
				nextSwapStartDate.AdjustToBusDate(ccy,fixAdjRule);
			}

			CC_Ostringstream nextSwapStartDateDesc;
			nextSwapStartDateDesc << CC_NS(std,fixed) << nextSwapStartDate.GetJulian();
			rowDescVec[NextSwapStartDate2] = nextSwapStartDateDesc.str();
			rowTypeVec[NextSwapStartDate2] = ARM_DATE_TYPE;

			//6)LASTSWAPSTARTDATE
			ARM_Date lastSwapStartDate(itsEndDate);
			lastSwapStartDate.AddPeriodMult(intFixFreq, -nbFixPeriods/intFixFreq);

			if (fixRule == K_ADJUSTED)
			{
				lastSwapStartDate.AdjustToBusDate(ccy,fixAdjRule);
			}

			CC_Ostringstream lastSwapStartDateDesc;
			lastSwapStartDateDesc << CC_NS(std,fixed) << lastSwapStartDate.GetJulian();
			rowDescVec[LastSwapStartDate2] = lastSwapStartDateDesc.str();
			rowTypeVec[LastSwapStartDate2] = ARM_DATE_TYPE;

			ARM_Date resetFromCallDate(lastSwapStartDate);
			int defaultResetGap = GetCurrencyUnit()->GetSpotDays();
			resetFromCallDate.GapBusinessDay(-defaultResetGap, const_cast<char *>(GetVarResetCal().c_str()));

			//7)STRIKE
			double strike = const_cast< ARM_ReferenceValue& >(itsStrike).Interpolate(resetFromCallDate.GetJulian());
			CC_Ostringstream strikeDesc;
			strikeDesc << CC_NS(std,fixed) << strike;
			rowDescVec[Strike2] = strikeDesc.str();
			rowTypeVec[Strike2] = ARM_DOUBLE;

			// NOTIONAL
			CC_Ostringstream notionalDesc;
			double notional = const_cast< ARM_ReferenceValue& >(itsNominal).Interpolate(nextSwapStartDate.GetJulian());
			notionalDesc << CC_NS(std,fixed) << notional;
			rowDescVec[Notional2] = notionalDesc.str();
			rowTypeVec[Notional2] = ARM_DOUBLE;

			
			//9)STDFIXSWAPLET
			CC_Ostringstream stdFixSwapletDesc;
			stdFixSwapletDesc << "(";
			if (intPayRec == K_RCV)
			{
				stdFixSwapletDesc << "-";
			}
				
			stdFixSwapletDesc << "DCF(" << BermudaSwaptionColNamesTable2[LastSwapStartDate2] << "[i],"; 
			stdFixSwapletDesc << BermudaSwaptionColNamesTable2[SwapStartDate2] << "[i],";
			stdFixSwapletDesc << fixBasis << ")*";
			stdFixSwapletDesc << BermudaSwaptionColNamesTable2[Strike1] << "[i]*";
			stdFixSwapletDesc << "DF(" << ccy << "," << BermudaSwaptionColNamesTable1[SwapStartDate2] << "[i])";

			if (intPayRec == K_PAY)
			{
				stdFixSwapletDesc << "-";
			}
			else
			{
				stdFixSwapletDesc << "+";
			}


			stdFixSwapletDesc << "DCF(" << BermudaSwaptionColNamesTable2[LastSwapStartDate2] << "[i],";
			stdFixSwapletDesc << BermudaSwaptionColNamesTable2[NextSwapStartDate2] << "[i],";
			stdFixSwapletDesc << fixBasis << ")*";
			stdFixSwapletDesc << BermudaSwaptionColNamesTable2[Strike2] << "[i]*";
			stdFixSwapletDesc << "DF(" << ccy << "," << BermudaSwaptionColNamesTable2[NextSwapStartDate2] << "[i]))*";
			stdFixSwapletDesc << BermudaSwaptionColNamesTable2[Notional2] << "[i]";

			rowDescVec[StdFixSwaplet2] = stdFixSwapletDesc.str();
			rowTypeVec[StdFixSwaplet2] = ARM_STRING;
			
			//10)STDVARSWAPLET
			CC_Ostringstream stdVarSwapletDesc;

			stdVarSwapletDesc << "Swap(" << ccy << ",";
			stdVarSwapletDesc << BermudaSwaptionColNamesTable2[SwapStartDate2] << "[i],";
			stdVarSwapletDesc << BermudaSwaptionColNamesTable2[NextSwapStartDate2] << "[i], 0,";
			stdVarSwapletDesc << payRec << ",";
			stdVarSwapletDesc << fixFreq << "," << fixBasis << "," << varFreq << "," << varBasis << ",";
			stdVarSwapletDesc << "MarginCurve,NotionalCurve,";

			// With this code we want to prevent STUBs
			// Compute the number of periods between the startDate and the endDate
			int nbPeriods = CC_Round((nextSwapStartDate.GetJulian()-swapStartDate)/K_YEAR_LEN*callFreq);
			ARM_Date unAdjEndDate = swapStartDate;
			unAdjEndDate.AddPeriodMult(callFreq,nbPeriods);
			if (endDate <= unAdjEndDate.GetJulian())
				stdVarSwapletDesc << "SE";
			else
				stdVarSwapletDesc << "LE";

			stdVarSwapletDesc << ")";

			rowDescVec[StdVarSwaplet2] = stdVarSwapletDesc.str();
			rowTypeVec[StdVarSwaplet2] = ARM_STRING;

			//11)STDSWAPLET
			CC_Ostringstream stdSwapletDesc;
			stdSwapletDesc << BermudaSwaptionColNamesTable1[StdFixSwaplet2] << "[i]+";
			stdSwapletDesc << BermudaSwaptionColNamesTable1[StdVarSwaplet2] << "[i]";

			rowDescVec[StdSwaplet2] = stdSwapletDesc.str();
			rowTypeVec[StdSwaplet2] = ARM_STRING;

			//12)STDSWAP
			CC_Ostringstream stdSwapDesc;
			stdSwapDesc << BermudaSwaptionColNamesTable2[StdSwaplet2] << "[i]";

			if (nbCallPeriods > callFreq)
			{
				int dec = nbPeriods%callFreq;
				if (dec == 0)
				{
					dec = callFreq;
				}
				stdSwapDesc << "+ PV(" << "StdSwap[i+" << dec << "])";
			}
			rowDescVec[StdSwap2] = stdSwapDesc.str();
			rowTypeVec[StdSwap2] = ARM_STRING;
		}
		else
		{
			ARM_Date resetFromCallDate(callDate);
			int defaultResetGap = GetCurrencyUnit()->GetSpotDays();
			resetFromCallDate.GapBusinessDay(-defaultResetGap, const_cast<char *>(GetVarResetCal().c_str()));

			//7)STRIKE
			double strike = const_cast< ARM_ReferenceValue& >(itsStrike).Interpolate(resetFromCallDate.GetJulian());
			CC_Ostringstream strikeDesc;
			strikeDesc << CC_NS(std,fixed) << strike;
			rowDescVec[Strike2] = strikeDesc.str();
			rowTypeVec[Strike2] = ARM_DOUBLE;

			//11)STDSWAPLET
			CC_Ostringstream stdSwapletDesc;
			stdSwapletDesc << "Swap(" << ccy << "," << "SwapStartDate[i]," << "EndDate[i]," << "StrikesCurve" <<",";
			stdSwapletDesc << payRec << "," << fixFreq << "," << fixBasis << "," << varFreq << "," << varBasis << ",";
			stdSwapletDesc << "MarginCurve, NotionalCurve,";

			// With this code we want to prevent STUBs
			// Compute the number of periods between the startDate and the endDate
			double nbPeriodsDbl = (endDate-swapStartDate)/K_YEAR_LEN*intFixFreq;
			int nbPeriods = CC_Round(nbPeriodsDbl);
			ARM_Date unAdjEndDate = swapStartDate;
			unAdjEndDate.AddPeriodMult(intFixFreq,nbPeriods);
			if (endDate <= unAdjEndDate.GetJulian())
				stdSwapletDesc << "SE";
			else
				stdSwapletDesc << "LE";

			stdSwapletDesc << ")";

			rowDescVec[StdSwaplet2] = stdSwapletDesc.str();
			rowTypeVec[StdSwaplet2] = ARM_STRING;

			//12)STDSWAP
			CC_Ostringstream stdSwapDesc;
			if(eventIdx == callEventSize-1)
			{
				stdSwapDesc << BermudaSwaptionColNamesTable2[StdSwaplet2] << "[i]";
			}
			else
			{
				stdSwapDesc << BermudaSwaptionColNamesTable2[StdSwaplet2] << "[i] + PV(";
				stdSwapDesc << BermudaSwaptionColNamesTable2[StdSwap2] << "[i+1])";
			}
			rowDescVec[StdSwap2] = stdSwapDesc.str();
			rowTypeVec[StdSwap2] = ARM_STRING;   
		}
		
		
		//13)OPTION
		CC_Ostringstream optionDesc;
		if(eventIdx == callEventSize-1)
		{
			optionDesc << "MAX(" << BermudaSwaptionColNamesTable2[StdSwap2] << "[i]";
			optionDesc << "-" << BermudaSwaptionColNamesTable2[Fees2] << "[i]," << "0" << ")";
		}
		else
		{
			optionDesc << "MAX(" << BermudaSwaptionColNamesTable2[StdSwap2] << "[i]";
			optionDesc << "-" << BermudaSwaptionColNamesTable2[Fees2] << "[i],";
			optionDesc << "PV("<< BermudaSwaptionColNamesTable2[Option2] << "[i+1])" << ")";
		}

		rowDescVec[Option2] = optionDesc.str();
		rowTypeVec[Option2] = ARM_STRING;

		//14)BERMUDA
		CC_Ostringstream bermudaDesc;
		if (eventIdx == 0) 
		{
			if (itsExoSwaption->GetPorS() == K_PAY)
				bermudaDesc << "-";

			bermudaDesc << BermudaSwaptionColNamesTable2[Option2] << "[i]";
			rowDescVec[Bermuda2] = bermudaDesc.str();
			rowTypeVec[Bermuda2] = ARM_STRING;
		}
		else 
		{
			bermudaDesc << "0";	
			rowDescVec[Bermuda2] = bermudaDesc.str();
			rowTypeVec[Bermuda2] = ARM_DOUBLE;
		}

		//15)SWAPRATE description
		CC_Ostringstream swaprateDesc;
		swaprateDesc << "SwapRate(" << ccy << ",";
		swaprateDesc << BermudaSwaptionColNamesTable2[StartDate2] << "[i],";
		swaprateDesc << BermudaSwaptionColNamesTable2[DealEndDate2] << "[i])";
		
		rowDescVec[SwapRate2] = swaprateDesc.str();
		rowTypeVec[SwapRate2] = ARM_STRING;

		//15bis)FRONTIER DESCRIPTION
		CC_Ostringstream frontierDesc;
		if(eventIdx == callEventSize-1)
		{
			frontierDesc << "Frontier(" << BermudaSwaptionColNamesTable2[StdSwap2] << "[i]";
			frontierDesc << "-" << BermudaSwaptionColNamesTable2[Fees2] << "[i],";
			frontierDesc << "0," << BermudaSwaptionColNamesTable2[SwapRate2] << "[i])";
		}
		else
		{
			frontierDesc << "Frontier(" << BermudaSwaptionColNamesTable2[StdSwap2] << "[i]";
			frontierDesc << "-" << BermudaSwaptionColNamesTable2[Fees2] << "[i],";
			frontierDesc << BermudaSwaptionColNamesTable2[Option2] << "[i+1]," << BermudaSwaptionColNamesTable2[SwapRate2] << "[i])";
		}
		rowDescVec[Frontier2] = frontierDesc.str();
		rowTypeVec[Frontier2] = ARM_STRING;


		//16)EXERCISECONDITION
		CC_Ostringstream exerciseconditionDesc;
		if(eventIdx == callEventSize-1)
		{
			exerciseconditionDesc << BermudaSwaptionColNamesTable2[SwapStartDate2] << "[i]-";
			exerciseconditionDesc << BermudaSwaptionColNamesTable2[SwapStartDate2] << "[i]";
		}
		else
			exerciseconditionDesc << BermudaSwaptionColNamesTable2[SwapStartDate2] << "[i]-";
			exerciseconditionDesc << BermudaSwaptionColNamesTable2[SwapStartDate2] << "[i]-";
			exerciseconditionDesc << "PV(" << BermudaSwaptionColNamesTable2[Option2] << "[i+1]" << ")";
		rowDescVec[ExerciseCondition2] = exerciseconditionDesc.str();
		rowTypeVec[ExerciseCondition2] = ARM_STRING;

		//17) Probabilities Description
		if (itsCalculateProbaFlag)
		{
			// 13) Reset Plus 1
			CC_Ostringstream reset_Plus1_Desc;
			double noticeDate_Plus1;
			if(eventIdx<datesStructure.GetDateStrip(CALL_BS_SCHED)->GetResetDates()->size()-1)
				noticeDate_Plus1 = (*(datesStructure.GetDateStrip(CALL_BS_SCHED)->GetResetDates()))[eventIdx+1] ;
			else
				noticeDate_Plus1 = itsEndDate.GetJulian();
			reset_Plus1_Desc << CC_NS(std,fixed) << noticeDate_Plus1;
			rowDescVec[ExerciseCondition2 + 1] = reset_Plus1_Desc.str();
			rowTypeVec[ExerciseCondition2 + 1] = ARM_DATE_TYPE;
			// 14) Final Date
			CC_Ostringstream FinalDate_Desc;
			FinalDate_Desc<< CC_NS(std,fixed) << itsEndDate.GetJulian();
			rowDescVec[ExerciseCondition2  + 2] = FinalDate_Desc.str();
			rowTypeVec[ExerciseCondition2  + 2] = ARM_DATE_TYPE;

			// 15) Exercise Indicator 
			CC_Ostringstream Indicator_Desc;
			if(eventIdx == callEventSize-1)
			{
				Indicator_Desc<< "if(" << BermudaSwaptionColNamesTable2[StdSwap2] << "[i]-";
				Indicator_Desc <<  BermudaSwaptionColNamesTable2[Fees2] <<"[i]>0,1,0)";
			}

			else
			{
				Indicator_Desc << "if(" << BermudaSwaptionColNamesTable2[StdSwap2] << "[i]-"; 
				Indicator_Desc << BermudaSwaptionColNamesTable2[Fees2] << "[i]-";
				Indicator_Desc << "PV(" << BermudaSwaptionColNamesTable2[Bermuda2] << "[i+1])>0,1,0)";
			}

			rowDescVec[ExerciseCondition2  + 3] = Indicator_Desc.str();
			rowTypeVec[ExerciseCondition2  + 3] = ARM_STRING;

			// 17) Proba 1
			CC_Ostringstream Prob1_Desc;
			if (eventIdx == 0)
			{
				Prob1_Desc<< "ExerciseIndicator[i]*DF("<<ccy<<",FinalDate[i])";
				rowTypeVec[ExerciseCondition2  + 4] = ARM_STRING;
			}
			else
			{
				Prob1_Desc<< " 0 ";
				rowTypeVec[ExerciseCondition2  + 4] = ARM_DOUBLE;
			}
			rowDescVec[ExerciseCondition2  + 4] = Prob1_Desc.str();
			

			// 16) Probabilities
			size_t currProb = 5;
			int nbCols = 4 + 2*(itsNbCalculatedProbabilities-1);
			for(int i = 5; i<=nbCols;++i)
			{
				CC_Ostringstream Prob0_Desc;
				size_t maxIndex = (i == nbCols-1)?(callEventSize-1) : currProb-4;
				/////// Probi0
				if(eventIdx<maxIndex)
				{
					Prob0_Desc<<"MAX(ExerciseIndicator[i], PV("<< (BermudaSwaptionProbaColNamesTable2[i-1].c_str())<<"[i+1])/DF("<<ccy<<",ResetPlus1[i]))";
					rowTypeVec[ExerciseCondition2  + i] = ARM_STRING;
				}
				else if (eventIdx==maxIndex)
				{
					Prob0_Desc<<"ExerciseIndicator[i]";
					rowTypeVec[ExerciseCondition2  + i] = ARM_STRING;
				}
				else
				{
					Prob0_Desc<<"0";
					rowTypeVec[ExerciseCondition2  + i] = ARM_DOUBLE;
				}
				rowDescVec[ExerciseCondition2  + i] = Prob0_Desc.str();

				CC_Ostringstream Prob1_Desc;
				if(eventIdx<=maxIndex)
				{
					Prob1_Desc<<(BermudaSwaptionProbaColNamesTable2[i-1].c_str())<<"[i]*DF("<<ccy<<",FinalDate[i])";
					i++;
					rowTypeVec[ExerciseCondition2  + i] = ARM_STRING;
				}
				else
				{
					Prob1_Desc<<"0";
					i++;
					rowTypeVec[ExerciseCondition2  + i] = ARM_DOUBLE;
				}
				rowDescVec[ExerciseCondition2  + i] = Prob1_Desc.str();
				currProb++;
			}
		}
	}
	else if (genSecType == 3) // TREE&MC (not used default !)
	{
		//1)RESETDATE description (Notice date)
		double noticeDate = (*(datesStructure.GetDateStrip(CALL_BS_SCHED)->GetResetDates()))[eventIdx] ;
		CC_Ostringstream noticeDateDesc;
		noticeDateDesc << CC_NS(std,fixed) << noticeDate;
		rowDescVec[ResetDate3] = noticeDateDesc.str();
		rowTypeVec[ResetDate3] = ARM_DATE_TYPE;

		 //2)STARTDATE description (Call date)
		double callDate = (*(datesStructure.GetDateStrip(CALL_BS_SCHED)->GetFlowStartDates()))[eventIdx] ;
		ARM_Date armCallDate;
		armCallDate = ARM_Date(callDate);
		CC_Ostringstream callDateDesc;
		callDateDesc << CC_NS(std,fixed) << callDate;
		rowDescVec[StartDate3] = callDateDesc.str();
		rowTypeVec[StartDate3] = ARM_DATE_TYPE;

		//3)ENDDATE description 
		CC_Ostringstream endDateDesc;
		if(callFreq == K_ANNUAL)
		{
			armCallDate.AddMonths(12);
			armCallDate.AdjustToBusDate(ccy,callAdjRule);	
		}
		if(callFreq == K_SEMIANNUAL)
		{
			armCallDate.AddMonths(6);
			armCallDate.AdjustToBusDate(ccy,callAdjRule);	
		}
		if(callFreq == K_QUARTERLY)
		{
			armCallDate.AddMonths(3);
			armCallDate.AdjustToBusDate(ccy,callAdjRule);	
		}
		if(callFreq == K_MONTHLY)
		{
			armCallDate.AddMonths(1);
			armCallDate.AdjustToBusDate(ccy,callAdjRule);		
		}
		double endDate;
		if (eventIdx == callEventSize-1)
		{
			endDate = itsEndDate.GetJulian();
		}
		else
		{
			endDate = armCallDate.GetJulian(); 
		}
		endDateDesc << CC_NS(std,fixed) << endDate;
		rowDescVec[EndDate3] = endDateDesc.str();
		rowTypeVec[EndDate3] = ARM_DATE_TYPE;

		//4)SWAPSTARTDATE 
		double startDate = itsStartDate.GetJulian(); 
		double swapStartDate=0;
		if(startDate < callDate)
		{
			swapStartDate = callDate;
		}
		else
		{
			swapStartDate = startDate;
		}
		CC_Ostringstream swapStartDateDesc;
		swapStartDateDesc << CC_NS(std,fixed) << swapStartDate;
		rowDescVec[SwapStartDate3] = swapStartDateDesc.str();
		rowTypeVec[SwapStartDate3] = ARM_DATE_TYPE;

		//5)STRIKE
		double strike = const_cast< ARM_ReferenceValue& >(itsStrike).Interpolate(callDate);
		CC_Ostringstream strikeDesc;
		strikeDesc << CC_NS(std,fixed) << strike;
		rowDescVec[Strike3] = strikeDesc.str();
		rowTypeVec[Strike3] = ARM_DOUBLE;

		//6)FEES
		CC_Ostringstream feesDesc;
		double fees = const_cast< ARM_ReferenceValue& >(itsFees).Interpolate(noticeDate);
		feesDesc << CC_NS(std,fixed) << fees;
		rowDescVec[Fees3] = feesDesc.str();
		rowTypeVec[Fees3] = ARM_DOUBLE;

		//7)STDSWAPLET
		CC_Ostringstream stdSwapletDesc;
		stdSwapletDesc << "Swap(" << ccy << "," << "SwapStartDate[i]," << "EndDate[i]," << "StrikesCurve" <<",";
		stdSwapletDesc << payRec << "," << fixFreq << "," << fixBasis << "," << varFreq << "," << varBasis << ",";
		stdSwapletDesc << "MarginCurve, NotionalCurve)";
		rowDescVec[StdSwaplet3] = stdSwapletDesc.str();
		rowTypeVec[StdSwaplet3] = ARM_STRING;

		//8)REVERSESTDSWAPLET
		CC_Ostringstream reversestdswapletDesc;
		reversestdswapletDesc << "Swap(" << ccy << "," << "SwapStartDate[i]" << "," << "EndDate[i]" << "," << "StrikesCurve" << ",";
		reversestdswapletDesc << invPayRec << "," << fixFreq << "," << fixBasis << "," << varFreq << "," << varBasis << ",";
		reversestdswapletDesc << "MarginCurve, NotionalCurve)";
		rowDescVec[ReverseStdSwaplet3] = reversestdswapletDesc.str();
		rowTypeVec[ReverseStdSwaplet3] = ARM_STRING;

		//9)OPTIONANNULATION
		CC_Ostringstream optionannulationDesc;
		if (eventIdx < callEventSize -1)
		{
			optionannulationDesc << "ReverseStdSwaplet[i]" << "+" << "Exercise(-" << "ReverseStdSwaplet[i] - Fees[i] , OptionAnnulation[i+1])";
		}
		else
		{
			optionannulationDesc << "ReverseStdSwaplet[i]" << "+" << "MAX(-" <<  "ReverseStdSwaplet[i] - Fees[i], 0)";
		}
		rowDescVec[OptionAnnulation3] = optionannulationDesc.str();
		rowTypeVec[OptionAnnulation3] = ARM_STRING;

		//10)BERMUDA
		CC_Ostringstream bermudaDesc;
		if (eventIdx == 0) 
		{
			if (itsExoSwaption->GetPorS() == K_PAY)
				bermudaDesc << "-";

			bermudaDesc << "OptionAnnulation[i] + StdSwaplet[i]";
			rowDescVec[Bermuda3] = bermudaDesc.str();
			rowTypeVec[Bermuda3] = ARM_STRING;
		}
		else 
		{
			bermudaDesc << "StdSwaplet[i]";
			rowDescVec[Bermuda3] = bermudaDesc.str();
			rowTypeVec[Bermuda3] = ARM_STRING;
		}

		//10)CONTROLVARIATES
		size_t ctrlVarSize = (*itsControlVariate).size();
		size_t currCtrlVar;
		size_t i;
		size_t j=0;
		size_t currCol1;
		size_t currCol2;
		currCol1 = Bermuda3 + j + 1;
		if (ctrlVarSize != 0) 
		{
			for (i=0; i<2*ctrlVarSize;i++)
			{
				CC_Ostringstream cvoptionAnnulationIDesc;
				CC_Ostringstream controlVariateIDesc;

				currCtrlVar = (*itsControlVariate)[j];
				//currCol1 = Bermuda3 + j + 1;
				currCol2 = currCol1 + 1; 
				if(eventIdx == 0) //First row
				{
					if(eventIdx != currCtrlVar-1) //Diagonal
					{
						cvoptionAnnulationIDesc << "ReverseStdSwaplet[i]" << "+" << "Exercise(0,";
						cvoptionAnnulationIDesc << "-ReverseStdSwaplet[i]-1000000000," << "CVOptionAnnulation" << j+1 << "[i+1])";			
					}
					else
					{
						cvoptionAnnulationIDesc << "ReverseStdSwaplet[i]" << "+" << "Exercise(0,";
						cvoptionAnnulationIDesc << "-ReverseStdSwaplet[i]," << "CVOptionAnnulation" << j+1 << "[i+1])";
					
					}
					rowDescVec[currCol1] = cvoptionAnnulationIDesc.str();
					rowTypeVec[currCol1] = ARM_STRING;

					controlVariateIDesc << "StdSwaplet[i]" << "+" << "CVOptionAnnulation" << j+1 << "[i]";
					rowDescVec[currCol2] = controlVariateIDesc.str();
					rowTypeVec[currCol2] = ARM_STRING;
				}
				else if((eventIdx != 0) && (eventIdx == currCtrlVar-1)) //Diagonal
				{
					cvoptionAnnulationIDesc << "ReverseStdSwaplet[i]" << "+" << "Exercise(0,";
					cvoptionAnnulationIDesc << "-ReverseStdSwaplet[i]," << "CVOptionAnnulation" << j+1 << "[i+1])";
					rowDescVec[currCol1] = cvoptionAnnulationIDesc.str();
					rowTypeVec[currCol1] = ARM_STRING;

					controlVariateIDesc << "StdSwaplet[i]"; 
					rowDescVec[currCol2] = controlVariateIDesc.str();
					rowTypeVec[currCol2] = ARM_STRING;

				}
				else if ((eventIdx != 0) && (eventIdx != currCtrlVar-1))
				{
					if (eventIdx == callEventSize-1) //Last row description
					{
						cvoptionAnnulationIDesc << "ReverseStdSwaplet[i]" << "+" << "MAX(";
						cvoptionAnnulationIDesc << "-ReverseStdSwaplet[i]-100000000,0)";

						rowDescVec[currCol1] = cvoptionAnnulationIDesc.str();
						rowTypeVec[currCol1] = ARM_STRING;

						controlVariateIDesc << "StdSwaplet[i]";
						rowDescVec[currCol2] = controlVariateIDesc.str();
						rowTypeVec[currCol2] = ARM_STRING;

					}
					else //Not in diagonal
					{
						cvoptionAnnulationIDesc << "ReverseStdSwaplet[i]" << "+" << "Exercise(0,";
						cvoptionAnnulationIDesc << "-ReverseStdSwaplet[i]-100000000," << "CVOptionAnnulation" << j+1 << "[i+1])";
						rowDescVec[currCol1] = cvoptionAnnulationIDesc.str();
						rowTypeVec[currCol1] = ARM_STRING;

						controlVariateIDesc << "StdSwaplet[i]";
						rowDescVec[currCol2] = controlVariateIDesc.str();
						rowTypeVec[currCol2] = ARM_STRING;
					}
				}
				i++;
				j++;
				currCol1 = currCol1 +2;	
			}
		}
	}
	else if (genSecType == 4) // (Zero Coupon)
	{
		//1)RESETDATE description (Notice date)
		double noticeDate = (*(datesStructure.GetDateStrip(CALL_BS_SCHED)->GetResetDates()))[eventIdx] ;
		CC_Ostringstream noticeDateDesc;
		noticeDateDesc << CC_NS(std,fixed) << noticeDate;
		rowDescVec[ResetDate4] = noticeDateDesc.str();
		rowTypeVec[ResetDate4] = ARM_DATE_TYPE;
		
		//2)STARTDATE description (Call date)
		//Fix schedule
		ARM_Date swapStartDate	 = itsStartDate;
		ARM_Date swapEndDate	 = itsEndDate;
		const char* fixCal		 = itsFixPayCal.c_str();
		int fixFreq				 = itsFixFreq;
		int fixAdjRule			 = itsFixAdjRule;
		int fixRule				 = itsFixRule;
		int stubRule			 = itsStubRule;
		int fixBas				 = itsFixBasis;
		ARM_DateStrip fixSched(swapStartDate, swapEndDate,fixFreq,fixBas,
			fixCal, fixAdjRule, fixRule, stubRule, GETDEFAULTVALUE,
			fixFreq, GETDEFAULTVALUE, fixCal, K_ADVANCE, K_ARREARS, true);
		
		std::vector<double> fixSchedVector = *(fixSched.GetFlowStartDates());
		double fixStartDate;
		for (int f=0; f<fixSchedVector.size();f++)
		{
			fixStartDate = fixSchedVector[f];
			if (fixStartDate > noticeDate)
				break;
		}
		CC_Ostringstream startDateDesc;
		startDateDesc << CC_NS(std,fixed) << fixStartDate;
		rowDescVec[StartDate4] = startDateDesc.str();
		rowTypeVec[StartDate4] = ARM_DATE_TYPE;
		
		//3) PASTSTARTDATE
		double pastStartDate = itsPastStartDate.GetJulian();
		CC_Ostringstream pastStartDateDesc;
		pastStartDateDesc << CC_NS(std,fixed) << pastStartDate;
		rowDescVec[PastStartDate4] = pastStartDateDesc.str();
		rowTypeVec[PastStartDate4] = ARM_DATE_TYPE;

		//4)SWAPSTARTDATE description
		double swapStartDateDble = itsStartDate.GetJulian();
		CC_Ostringstream swapstartDateDesc;
		swapstartDateDesc << CC_NS(std,fixed) << swapStartDateDble;
		rowDescVec[SwapStartDate4] = swapstartDateDesc.str();
		rowTypeVec[SwapStartDate4] = ARM_DATE_TYPE;
		
		//5)ENDDATE description (Underlying End date)
		double endDate = itsEndDate.GetJulian(); 
		CC_Ostringstream endDateDesc;
		endDateDesc << CC_NS(std,fixed) << endDate;
		rowDescVec[EndDate4] = endDateDesc.str();
		rowTypeVec[EndDate4] = ARM_DATE_TYPE;
		
		//6)STRIKE
		double callDate = (*(datesStructure.GetDateStrip(CALL_BS_SCHED)->GetFlowStartDates()))[eventIdx] ;
		double strike = const_cast< ARM_ReferenceValue& >(itsStrike).Interpolate(callDate);
		CC_Ostringstream strikeDesc;
		strikeDesc << CC_NS(std,fixed) << strike;
		rowDescVec[Strike4] = strikeDesc.str();
		rowTypeVec[Strike4] = ARM_DOUBLE;
		
		//7)FEES
		CC_Ostringstream feesDesc;
		double fees = const_cast< ARM_ReferenceValue& >(itsFees).Interpolate(noticeDate);
		feesDesc << CC_NS(std,fixed) << fees;
		rowDescVec[Fees4] = feesDesc.str();
		rowTypeVec[Fees4] = ARM_DOUBLE;
		
		//8)DFSTART
		CC_Ostringstream dfstartDesc;
		dfstartDesc << "DF(" << ccy << ", StartDate[i])";
		rowDescVec[DfStart4] = dfstartDesc.str();
		rowTypeVec[DfStart4] = ARM_STRING;
		
		//9)CVGSTART
		double cvgStart;
		CC_Ostringstream cvgstartDesc;
		if (itsPastZC)
		{	
			cvgStart= CountYears(intFixBasis, itsPastStartDate.GetJulian(), fixStartDate);
		}
		else
		{
			cvgStart= CountYears(intFixBasis, swapStartDateDble, fixStartDate);
		}
		cvgstartDesc << cvgStart;
		rowDescVec[CvgStart4] = cvgstartDesc.str();
		rowTypeVec[CvgStart4] = ARM_STRING;
		
		//10)DFEND
		CC_Ostringstream dfendDesc;
		dfendDesc << "DF(" << ccy << ", EndDate[i])";
		rowDescVec[DfEnd4] = dfendDesc.str();
		rowTypeVec[DfEnd4] = ARM_STRING;
		
		//11)CVGEND
		double cvgEnd;
		if (itsPastZC)
			cvgEnd= CountYears(intFixBasis, itsPastStartDate.GetJulian(), endDate);
		else
			cvgEnd= CountYears(intFixBasis, swapStartDateDble, endDate);
		CC_Ostringstream cvgendDesc;
		cvgendDesc << cvgEnd;
		rowDescVec[CvgEnd4] = cvgendDesc.str();
		rowTypeVec[CvgEnd4] = ARM_STRING;
		
		//12)ANNUITY
		CC_Ostringstream annuityDesc;
		annuityDesc << "ANNUITY(" << ccy << ", StartDate[i], EndDate[i]," << varFreq << "," << varBasis << ",";
		annuityDesc << "MarginCurve)";
		rowDescVec[Annuity4] = annuityDesc.str();
		rowTypeVec[Annuity4] = ARM_STRING;
		
		//13)FLOWS
		CC_Ostringstream flowsDesc;
		if (payRec == "R")
		{
			flowsDesc << "NotionalCurve*(-Pow(1+Strike[i]/" << intFixFreq << ",CvgStart[i]*" << intFixFreq << ")*";
			flowsDesc << "DFStart[i]";
			flowsDesc << "+";
			flowsDesc << "Pow(1+Strike[i]/" << intFixFreq << ",CvgEnd[i]*" << intFixFreq << ")*";
			flowsDesc << "DFEnd[i]";
			flowsDesc << "-";
			flowsDesc << "Annuity[i])-Fees[i]" ;
		}
		else if (payRec == "P")
		{	
			flowsDesc << "NotionalCurve*(POW(1+Strike[i]/" << intFixFreq << ",CvgStart[i]*" << intFixFreq << ")*";
			flowsDesc << "DFStart[i]";
			flowsDesc << "-";
			flowsDesc << "POW(1+Strike[i]/" << intFixFreq << ",CvgEnd[i]*" << intFixFreq << ")*";
			flowsDesc << "DFEnd[i]";
			flowsDesc << "+";
			flowsDesc << "Annuity[i])-Fees[i]";
		}
		rowDescVec[Flows4] = flowsDesc.str();
		rowTypeVec[Flows4] = ARM_STRING;
		
		//14)CASHFLOW
		CC_Ostringstream cashflowDesc;
		if (eventIdx == callEventSize - 1)
		{
			cashflowDesc << "MAX(Flows[i],0)";
		}
		else
		{
			cashflowDesc << "Exercise(0,Flows[i] ,CashFlow[i+1])";
		}
		
		rowDescVec[CashFlow4] = cashflowDesc.str();
		rowTypeVec[CashFlow4] = ARM_STRING;
		
		//15)BERMUDA
		CC_Ostringstream bermudaDesc;
		if (eventIdx == 0)
		{
			if (itsExoSwaption->GetPorS() == K_PAY)
				bermudaDesc << "-";

			bermudaDesc << "CashFlow[i]";
			rowDescVec[Bermuda4] = bermudaDesc.str();
			rowTypeVec[Bermuda4] = ARM_STRING;
		}
		else
		{
			bermudaDesc << "0";
			rowDescVec[Bermuda4] = bermudaDesc.str();
			rowTypeVec[Bermuda4] = ARM_DOUBLE;
		}

		//16)DEAL ENDDATE
		ARM_Date DealEndDate(itsEndDate);
		CC_Ostringstream dealEndDateDesc;
		dealEndDateDesc << CC_NS(std,fixed) << DealEndDate.GetJulian();
		rowDescVec[DealEndDate4] = dealEndDateDesc.str();
		rowTypeVec[DealEndDate4] = ARM_DATE_TYPE;

		//17)SWAPRATE description
		CC_Ostringstream swaprateDesc;
		swaprateDesc << "SwapRate(" << ccy << ",";
		swaprateDesc << BermudaSwaptionColNamesTable4[StartDate4] << "[i],";
		swaprateDesc << BermudaSwaptionColNamesTable4[DealEndDate4] << "[i])";
		
		rowDescVec[SwapRate4] = swaprateDesc.str();
		rowTypeVec[SwapRate4] = ARM_STRING;

		//18)FRONTIER DESCRIPTION
		CC_Ostringstream frontierDesc;
		if(eventIdx == callEventSize-1)
		{
			frontierDesc << "Frontier(Flows[i],0,SwapRate[i])";
		}
		else
		{
			frontierDesc << "Frontier(Flows[i],CashFlow[i+1],SwapRate[i])";
		}
		rowDescVec[Frontier4] = frontierDesc.str();
		rowTypeVec[Frontier4] = ARM_STRING;
	}

    return ARM_RowInfo(rowDescVec,rowTypeVec);
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: InitPriceableColumns
///	Returns: nothing
///	Action : Feed the MDM with input market objects
/////////////////////////////////////////////////////////////////
void ARM_BermudaSwaptionCalculator::InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const
{
    string zeroValue("0");

	int genSecType = GetGenSecType();
	if (genSecType == 1)
	{
		rowDescVec[Strike1] = zeroValue;
		rowTypeVec[Strike1] = ARM_DOUBLE;

		rowDescVec[Fees1] = zeroValue;
		rowTypeVec[Fees1] = ARM_DOUBLE;

		rowDescVec[StdFixSwaplet1] = zeroValue;
		rowTypeVec[StdFixSwaplet1] = ARM_DOUBLE;

		rowDescVec[StdVarSwaplet1] = zeroValue;
		rowTypeVec[StdVarSwaplet1] = ARM_DOUBLE;

		rowDescVec[StdSwaplet1] = zeroValue;
		rowTypeVec[StdSwaplet1] = ARM_DOUBLE;

		rowDescVec[StdSwap1] = zeroValue;
		rowTypeVec[StdSwap1] = ARM_DOUBLE;

		rowDescVec[Option1] = zeroValue;
		rowTypeVec[Option1] = ARM_DOUBLE;

		rowDescVec[Bermuda1] = zeroValue;
		rowTypeVec[Bermuda1] = ARM_DOUBLE;
	}
	else if (genSecType == 2)
	{

		rowDescVec[Strike2] = zeroValue;
		rowTypeVec[Strike2] = ARM_DOUBLE;

		rowDescVec[Fees2] = zeroValue;
		rowTypeVec[Fees2] = ARM_DOUBLE;

		rowDescVec[StdFixSwaplet2] = zeroValue;
		rowTypeVec[StdFixSwaplet2] = ARM_DOUBLE;

		rowDescVec[StdVarSwaplet2] = zeroValue;
		rowTypeVec[StdVarSwaplet2] = ARM_DOUBLE;

		rowDescVec[StdSwaplet2] = zeroValue;
		rowTypeVec[StdSwaplet2] = ARM_DOUBLE;

		rowDescVec[StdSwap2] = zeroValue;
		rowTypeVec[StdSwap2] = ARM_DOUBLE;

		rowDescVec[Option2] = zeroValue;
		rowTypeVec[Option2] = ARM_DOUBLE;

		rowDescVec[Bermuda2] = zeroValue;
		rowTypeVec[Bermuda2] = ARM_DOUBLE;

		rowDescVec[SwapRate2] = zeroValue;
		rowTypeVec[SwapRate2] = ARM_DOUBLE;
		
		rowDescVec[ExerciseCondition2] = zeroValue;
		rowTypeVec[ExerciseCondition2] = ARM_DOUBLE;
	}
	else if (genSecType == 3)
	{
		rowDescVec[Strike3] = zeroValue;
		rowTypeVec[Strike3] = ARM_DOUBLE;

		rowDescVec[Fees3] = zeroValue;
		rowTypeVec[Fees3] = ARM_DOUBLE;

		rowDescVec[StdSwaplet3] = zeroValue;
		rowTypeVec[StdSwaplet3] = ARM_DOUBLE;

		rowDescVec[ReverseStdSwaplet3] = zeroValue;
		rowTypeVec[ReverseStdSwaplet3] = ARM_DOUBLE;

		rowDescVec[OptionAnnulation3] = zeroValue;
		rowTypeVec[OptionAnnulation3] = ARM_DOUBLE;

		rowDescVec[Bermuda3] = zeroValue;
		rowTypeVec[Bermuda3] = ARM_DOUBLE;
	}
	else if (genSecType == 4)
	{
		rowDescVec[Strike4] = zeroValue;
		rowTypeVec[Strike4] = ARM_DOUBLE;

		rowDescVec[Fees4] = zeroValue;
		rowTypeVec[Fees4] = ARM_DOUBLE;

		rowDescVec[Annuity4] = zeroValue;
		rowTypeVec[Annuity4] = ARM_DOUBLE;

		rowDescVec[DfStart4] = zeroValue;
		rowTypeVec[DfStart4] = ARM_DOUBLE;

		rowDescVec[DfEnd4] = zeroValue;
		rowTypeVec[DfEnd4] = ARM_DOUBLE;

		rowDescVec[CashFlow4] = zeroValue;
		rowTypeVec[CashFlow4] = ARM_DOUBLE;
		
		rowDescVec[Bermuda4] = zeroValue;
		rowTypeVec[Bermuda4] = ARM_DOUBLE;
	}
}

/////////////////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: AdjustStubRule
///	Returns: void
///	Action : rustine: we adjust the stub rule depending on the dayS. 
/////////////////////////////////////////////////////////////////
int ARM_BermudaSwaptionCalculator::AdjustStubRule(ARM_Date& firstDate, ARM_Date& lastDate)
{
	int firstDay = firstDate.GetDay();
	int lastDay	 = lastDate.GetDay();
	int stub;

	if (firstDay < lastDay)
	{
		stub = K_LONGSTART;
	}
	else
	{
		stub = K_SHORTSTART;
	}

	//SetStubRule(stub);
	return stub;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: CustomDatesStructure
///	Returns: ARM_DateStripCombiner
///	Action : create the list of all event/notice dates of the BermudaSwaption
///			customized Schedule Not implemented yet.
/////////////////////////////////////////////////////////////////
ARM_DateStripCombiner ARM_BermudaSwaptionCalculator::CustomDatesStructure(const ARM_DateStripVector& dateStrips) const
{
	return ARM_DateStripCombiner();
}

/////////////////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: DatesStructure
///	Returns: ARM_DateStripVector
///	Action : create the list of all event dates of the BermudaSwaption. 
///			 The DateStripCombiner merges event dates of each
///          legs
/////////////////////////////////////////////////////////////////
ARM_DateStripCombiner ARM_BermudaSwaptionCalculator::DatesStructure() const
{
	//General datas
	ARM_Currency ccy = GetCcy();
	
	if (itsExoSwaption && (itsCallNotice == 9999)) //case swaption loaded from Summit !!!
	{
		ARM_DateStripVector SchedVect(1,NULL);

		ARM_SwapLeg* fundingLeg			= (ARM_SwapLeg*) itsExoSwaption->GetFloatLeg();
		ARM_Vector* vFlowStartDates		= fundingLeg->GetFlowStartDates();
		ARM_Vector* vFlowEndDates		= fundingLeg->GetFlowEndDates();
		ARM_Vector* vResetDates			= fundingLeg->GetResetDates();
		ARM_Vector* vFwdRateStartDates	= fundingLeg->GetFwdRateStartDates();
		ARM_Vector* vFwdRateEndDates	= fundingLeg->GetFwdRateEndDates();
		ARM_Vector* vPaymentDates		= fundingLeg->GetPaymentDates();
		ARM_Vector* vInterestTerms		= fundingLeg->GetInterestTerms();

		ARM_Vector* NoticeDates = itsExoSwaption->GetExerciseStyle()->GetExerciseStartDates();
		ARM_Vector* ExpiryDates = itsExoSwaption->GetExerciseStyle()->GetExerciseEndDates();

		int noticeSize = NoticeDates->GetSize();
		int legSize    = vResetDates->GetSize();

		std::vector<double>& FlowStartDates		= new std::vector<double>(noticeSize);
		std::vector<double>& FlowEndDates			= new std::vector<double>(noticeSize);
		std::vector<double>& ResetDates			= new std::vector<double>(noticeSize);
		std::vector<double>& FwdRateStartDates	= new std::vector<double>(noticeSize);
		std::vector<double>& FwdRateEndDates	    = new std::vector<double>(noticeSize);
		std::vector<double>& PaymentDates			= new std::vector<double>(noticeSize);
		std::vector<double>& InterestTerms		= new std::vector<double>(noticeSize);
		std::vector<double>& InterestDays			= new std::vector<double>(noticeSize);

		ARM_Vector* fees = NULL;
		ARM_Vector* StartDatesFees = new ARM_Vector(noticeSize);

		for (int i = 0; i < noticeSize; i++)
		{
			for (int j = 0; j < legSize; j++)
			{
				if (( ARM_Date(ExpiryDates->Elt(i))-ARM_Date(vFlowStartDates->Elt(j)) )
                    <
                    7
                   )
				{
					int idxEnd = j;
      
					if ( i < noticeSize-1 )
					{
                       try
                       {
						  while (ARM_Date(ExpiryDates->Elt(i+1)) != ARM_Date(vFlowEndDates->Elt(idxEnd)))
                          {
							  idxEnd++;

							  if ( idxEnd == legSize )
                              {
								char msg[200];
								char strDate[11];
								ARM_Date(ExpiryDates->Elt(i+1)).JulianToStrDate(strDate);
								sprintf(msg, "ARM_BermudaSwaptionCalculator::DateStructure : ExpiryDates[%d] (%s) does not match any corridor start date", i+1, strDate);
								throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg);
                              }
                          }
                       }

                       catch(Exception& )
                       {
                           idxEnd = j;

                           while ( (ARM_Date(ExpiryDates->Elt(i+1))-ARM_Date(vFlowEndDates->Elt(idxEnd)))
                                   >=
                                   7
                                 )
                           {
							  idxEnd++;

							  if ( idxEnd == legSize )
                              {
								char msg[200];
								char strDate[11];
								ARM_Date(ExpiryDates->Elt(i+1)).JulianToStrDate(strDate);
								sprintf(msg, "ARM_BermudaSwaptionCalculator::DateStructure : ExpiryDates[%d] (%s) does not match any corridor start date", i+1, strDate);
								throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg);
                              }
                           }
                       }
					}
					else
					   idxEnd = legSize-1;

					ResetDates->Elt(i)		  = NoticeDates->Elt(i);
					FlowStartDates->Elt(i)    = vFlowStartDates->Elt(j);
					StartDatesFees->Elt(i)    = ResetDates->Elt(i); //FlowStartDates->Elt(i);
					FlowEndDates->Elt(i)	  = vFlowEndDates->Elt(idxEnd);
					FwdRateStartDates->Elt(i) = vFwdRateStartDates->Elt(j);
					FwdRateEndDates->Elt(i)   = vFwdRateEndDates->Elt(idxEnd);
					PaymentDates->Elt(i)      = vPaymentDates->Elt(idxEnd);
					InterestTerms->Elt(i)     = vInterestTerms->Elt(j);
					break;
				}
			}
			
            if ( j == legSize )
			{
				char msg[200];
				char strDate[11];

				ARM_Date(NoticeDates->Elt(i)).JulianToStrDate(strDate);
				sprintf(msg, "ARM_BermudaSwaptionCalculator::DateStructure : NoticeDates[%d] (%s) does not match any underlying start date", i, strDate);
				
                throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg);
			}
		}
		
		for (i=0; i<noticeSize; i++)
		{
			InterestDays->Elt(i) = DaysBetweenDates(itsVarBasis,
													FlowStartDates->Elt(i), 
													FlowEndDates->Elt(i));
		}

		ARM_DateStrip CallSched(FlowStartDates, FlowEndDates, FwdRateStartDates, FwdRateEndDates,
								ResetDates, PaymentDates, InterestDays, InterestTerms);

		SchedVect[CALL_BS_SCHED] = &CallSched;

		// process fees
		if (itsExoSwaption->GetFee() != NULL)
		{
			//Adjust summit fees
			int feeSize = itsExoSwaption->GetFee()->size();
			ARM_ReferenceValue newFees;
			ARM_Vector* vecNoticeDates	= itsExoSwaption->GetFee()->GetDiscreteDates();
			ARM_Vector* vecFees			= itsExoSwaption->GetFee()->GetDiscreteValues();
			ARM_Vector* vecNewFees		= NULL; 
			double sumFee				= vecFees->Elt(0);
			double firstFee				= 0.0;
			double nextFee				= 0.0;
			bool isZero					= false;
			bool isOne					= false;

			for (i = 1; i < feeSize; i++)
			{
				firstFee = vecFees->Elt(i-1);
				nextFee  = vecFees->Elt(i);
				sumFee += nextFee;

				if ((firstFee == 1) && (firstFee == nextFee))
				{
					isOne = true;
				}
				else
				{
					isOne = false;
				}
			}
			
			//Keep the current fees
			if (sumFee == 0.0)
			{
				isZero = true;
				itsFees = ARM_ReferenceValue((ARM_Vector*)vecNoticeDates->Clone(), (ARM_Vector*) vecFees->Clone());
			}
			//Modify the fees from ONE to ZERO
			else if (isOne == true)
			{
				for (i = 0; i < feeSize; i++)
				{
					vecNewFees = new ARM_Vector(feeSize, 0.0); 
					itsFees = ARM_ReferenceValue((ARM_Vector*)vecNoticeDates->Clone(), (ARM_Vector*) vecNewFees->Clone());
				}
			}
			//Custom Fee
			else
			{
				vecNewFees = new ARM_Vector(feeSize, 0.0); 
				double currFee = 0.0;
				double currNoticeDate = 0.0;
				double currNominal = 0.0;
				for (i = 0; i < feeSize; i++)
				{
					currFee = vecFees->Elt(i);
					currNoticeDate = vecNoticeDates->Elt(i);
					if ( currFee == 1.0)
					{
						vecNewFees->Elt(i) = 0.0;
					}
					else if (currFee == 0.0)
					{
						vecNewFees->Elt(i) = 0.0;
					}
					else
					{
						if (itsExoSwaption->GetAmount()->size() == 1)
							currNominal = itsExoSwaption->Get1stLeg()->GetAmount()->GetDiscreteValues()->Elt(0);
						else
							currNominal = itsExoSwaption->Get1stLeg()->GetAmount()->Interpolate(currNoticeDate);

						if (itsExoSwaption->GetFixedLeg()->GetRcvOrPay() == K_PAY)
						{
							if (currFee < 1.0)
								vecNewFees->Elt(i) = (1.0 - currFee)*currNominal;				
							else //(currFee > 1) 
								vecNewFees->Elt(i) = -(currFee - 1.0)*currNominal;
						}
						else
						{
							if (currFee < 1.0)
								vecNewFees->Elt(i) = -(1.0 - currFee)*currNominal;				
							else //(currFee > 1) 
								vecNewFees->Elt(i) = (currFee - 1.0)*currNominal;
						}
					}
				}

				itsFees = ARM_ReferenceValue((ARM_Vector*)vecNoticeDates->Clone(), (ARM_Vector*) vecNewFees->Clone());
			}
			
			if (vecNewFees)
			{
				delete vecNewFees;
				vecNewFees = NULL;
			}
		}
		else
		{
			ARM_Vector* fees = new ARM_Vector(noticeSize,0.0);
			itsFees.SetDiscreteDates(StartDatesFees);
			itsFees.SetDiscreteValues(fees);
		}
		itsFees.SetCalcMethod(K_STEPUP_LEFT);

		delete FlowStartDates;
		FlowStartDates = NULL;

		delete FlowEndDates;
		FlowEndDates = NULL;
		
        delete ResetDates;
		ResetDates = NULL;
		
        delete FwdRateStartDates;
		FwdRateStartDates = NULL;
		
        delete FwdRateEndDates;
		FwdRateEndDates = NULL;
		
        delete PaymentDates;
		PaymentDates = NULL;
		
        delete InterestDays;
		InterestDays = NULL;
		
        delete InterestTerms;
		InterestTerms = NULL;

		delete StartDatesFees;
		StartDatesFees = NULL;
		
		itsCallDateStrip = ARM_DateStripPtr(new ARM_DateStrip(CallSched));

		return SchedVect;
	}
	else
	{
 		//General datas
		ARM_Currency ccy = itsCcy;
		int stubRule = itsStubRule;

		//Call datas
		int callFreq		= itsCallFreq;
		int callNotice		= itsCallNotice;
		const char* callCal	= itsCallCal.c_str();
		ARM_Date firstCall	= itsFirstCallDate;
		int firstDay		= firstCall.GetDay();
		ARM_Date lastCall	= itsLastCallDate;
		int lastDay			= lastCall.GetDay();
		int stub;
		if (firstDay < lastDay)
			stub = K_LONGSTART;
		else
			stub = K_SHORTSTART;

		int fixBasis	= itsFixBasis;
		int varAdjRule	= itsVarAdjRule;

		//Adjust the last call date.
		ARM_Date adjLastCallDate = ARM_Date(lastCall);
		if(callFreq == K_ANNUAL)
		{
			adjLastCallDate.AddMonths(12);
		}
		if(callFreq == K_SEMIANNUAL)
		{
			adjLastCallDate.AddMonths(6);
		}
		if(callFreq == K_QUARTERLY)
		{
			adjLastCallDate.AddMonths(3);
		}
		if(callFreq == K_MONTHLY)
		{
			adjLastCallDate.AddMonths(1);
		}

		//Funding schedule
		ARM_DateStrip fundingSched( itsStartDate, itsEndDate, itsVarFreq, itsVarBasis, callCal, itsVarAdjRule, itsVarRule, stub, itsVarResetGap,
								 itsVarFreq,  GETDEFAULTVALUE, callCal, K_ADVANCE,  K_ARREARS, true );
		const_cast< ARM_BermudaSwaptionCalculator* >(this)->itsFundDateStrip = ARM_DateStripPtr(new ARM_DateStrip(fundingSched));

		//coupon schedule
		ARM_DateStrip cpnSched( itsStartDate, itsEndDate, itsFixFreq, itsFixBasis, callCal, itsFixAdjRule, itsFixRule, stub, itsFixPayGap,
								 itsFixFreq,  GETDEFAULTVALUE, callCal, K_ADVANCE,  K_ARREARS, true );
		const_cast< ARM_BermudaSwaptionCalculator* >(this)->itsStructDateStrip = ARM_DateStripPtr(new ARM_DateStrip(cpnSched));

		ARM_Date endDate		= itsEndDate;
		int trueDay				= endDate.GetDay();
		int trueMonth			= adjLastCallDate.GetMonth();
		int trueYear			= adjLastCallDate.GetYear();
		ARM_Date lastCallDate   = ARM_Date(trueDay, trueMonth, trueYear);

		//Call schedule
		ARM_DateStrip CallSched( firstCall, 
								 lastCallDate, 
								 callFreq, 
								 fixBasis,
								 callCal, 
								 GetFixAdjRule(), 
								 GetFixRule(), 
								 stub, 
								 callNotice,
								 callFreq, 
								 GETDEFAULTVALUE, 
								 callCal, 
								 K_ADVANCE, 
								 K_ARREARS, 
								 true );

		ARM_DateStripVector SchedVect(1,NULL);
		SchedVect[CALL_BS_SCHED] = &CallSched;
		
		itsCallDateStrip = ARM_DateStripPtr(new ARM_DateStrip(CallSched));

		return SchedVect;
	}
}

/*/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: ComputeDomesticBasis
///	Returns: void
///	Action : Convert le basis from funding Ccy to domestic Ccy
/////////////////////////////////////////////////////////////////
void ARM_BermudaSwaptionCalculator::ComputeDomesticBasis()
{
	ARM_ZeroCurve* ycDomCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));
	ARM_ZeroCurve* ycBasisDomCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcBasisDomKey]));
	ARM_ZeroCurve* ycfundCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[FundingKey]));
	ARM_ZeroCurve* basisFundCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcBasisFundKey]));
	ARM_Forex* forex = dynamic_cast<ARM_Forex*>(GetMktDataManager()->GetData(GetKeys()[ForexKey]));

	ARM_BasisConverter basisConveter(GetDomesticCcy(),
					GetFundingCcy(),
					itsStructDateStrip,
					itsFundDateStrip,
					itsFundDateStrip,
					itsFundDayCount,
					itsFundFreq,
					itsFundDayCount,
					itsFundFreq,
					CreateClonedPtr(ycDomCurve),	
					CreateClonedPtr(ycfundCurve),
					CreateClonedPtr(ycBasisDomCurve),
					CreateClonedPtr(basisFundCurve),
					*forex,
					itsvCpnNominal,
					itsvInitialFundNominal,
					itsvInitialFundSpread);

	std::vector<double>& vdomMargin = basisConveter.ComputeDomMargin();
	CC_NS(std,auto_ptr)<std::vector<double>> HOLDCVDM(vdomMargin);	

	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();	
	std::vector<double> timeLags(*itsStructDateStrip->GetResetDates());
	timeLags -= asOfDate;

	/// we create a local curve to resize.....
	ARM_CurveInterpolator*  interpolator = (ARM_CurveInterpolator*)itsInitFundMarginProfile.GetInterpolator()->Clone();
	itsFundMarginProfile = ARM_Curve(timeLags,*vdomMargin, interpolator);
}
*/

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: CreateCstManager
///	Returns: void
///	Action : create the const manager (static data).
/////////////////////////////////////////////////////////////////
ARM_CstManagerPtr ARM_BermudaSwaptionCalculator::CreateCstManager()
{
	ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();

	vector<string>	cstNames(3,"");
	cstNames[0] = "StrikesCurve";
	cstNames[1] = "NotionalCurve";
	cstNames[2] = "MarginCurve";

	vector<ARM_GramFctorArg> cstVector;
	
	//Strike conversion
	ARM_Curve* tmpStrike;
	tmpStrike = RefValueToARM_Curve(&itsStrike, asOfDate, new ARM::ARM_StepUpRightOpenCstExtrapolDble());

	cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*>(tmpStrike))));

	//Notional conversion
	if (GetZCFlag() == true)
	{
		double initNotional = itsNominal.GetDiscreteValues()->Elt(0);
		cstVector.push_back(ARM_GramFctorArg(initNotional));
	}
	else
	{
		ARM_Curve* tmpNominal;
		tmpNominal = RefValueToARM_Curve(&itsNominal, asOfDate, new ARM::ARM_StepUpRightOpenCstExtrapolDble());
		cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(tmpNominal)));
	}

	//Margin conversion
	ARM_Curve* tmpMargin;

	if (GetApproxMarginFlag()&&(GetFixFreq()<4)&&(GetVarFreq()==4))
	{
		ARM_ReferenceValue* newMargin = ConvertVariableMargin(itsVarSpread,4,GetFixFreq(), GetVarBasis(), GetVarBasis());
		tmpMargin = RefValueToARM_Curve(newMargin, asOfDate, new ARM::ARM_StepUpRightOpenCstExtrapolDble());
		SetVarFreq(GetFixFreq());
		SetVarSpread(*newMargin);
		delete newMargin;
	}
	else
		tmpMargin = RefValueToARM_Curve(&itsVarSpread, asOfDate, new ARM::ARM_StepUpRightOpenCstExtrapolDble());
	
	cstVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(tmpMargin)));

	ARM_CstManagerPtr cstManagerPtr = ARM_CstManagerPtr(new ARM_CstManager(cstNames, cstVector));

	return cstManagerPtr;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: RefValueToARM_Curve
///	Returns: void
///	Action : convert a RefValue object to an arm_curve one.
/////////////////////////////////////////////////////////////////
ARM_Curve* ARM_BermudaSwaptionCalculator::RefValueToARM_Curve(ARM_ReferenceValue* refVal, ARM_Date& date, ARM_Interpolator<double,double>* interpolator)
{
	int size = 0;

	ARM_Vector*    dates    = NULL;
	ARM_Vector*    values   = NULL;
	std::vector<double>& gpDates  = NULL;
	std::vector<double>& gpValues = NULL;
	ARM_Curve*     curve    = NULL;

	dates = refVal->GetDiscreteDates();
	if (dates != NULL)
		size = dates->GetSize();

	if (size == 0)
	{
		gpDates = new std::vector<double>(1, 0.0);
	}
	else
	{
		gpDates  = new std::vector<double>(dates->GetSize(),dates->GetElt());
		*gpDates -= date.GetJulian();
	}

	values   = refVal->GetDiscreteValues();
	gpValues = new std::vector<double>(values->GetSize(),values->GetElt());

	// steptUpLeft car les strikes sont donnés avec les start dates
	curve = new ARM_Curve(*gpDates,*gpValues, interpolator);

	if (gpDates)
	   delete gpDates;
	gpDates = NULL;

	if (gpValues)
	   delete gpValues;
	gpValues = NULL;

	return curve;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: SetEquivalentVarIndexTerm
///	Returns: void
///	Action : set the right rate to propagate
/////////////////////////////////////////////////////////////////
void ARM_BermudaSwaptionCalculator::SetEquivalentVarIndexTerm()
{
	string varIndexTerm = GetVarIndexTerm();
	int fixFreq = GetFixFreq();
	int varFreq = GetVarFreq();
	char* ccy = GetCcy().GetCcyName();
	bool approxFlag = GetApproxMarginFlag();

	//We diffuse the smaller index between fixLeg (imagine) and VarLeg.
	int equivFreq = ((fixFreq <= varFreq) ? varFreq : fixFreq);
	if (equivFreq == K_ANNUAL)
		SetVarIndexTerm("12M");
	else if (equivFreq == K_SEMIANNUAL)
		SetVarIndexTerm("6M");
	else if (equivFreq == K_QUARTERLY)
		SetVarIndexTerm("3M");
	else if (equivFreq == K_MONTHLY)
		SetVarIndexTerm("1M");

	//Some various adjustments. 
	if (varIndexTerm == "12M")
	{
		if (strcmp(ccy,"EUR"))
		{
			SetVarIndexTerm("6M");
		}
	}
	else if (varIndexTerm == "3M")
	{
		if (fixFreq != K_QUARTERLY) //6M-12M
		{
			if (approxFlag)
			{
				if ((fixFreq == K_ANNUAL)&& strcmp(ccy,"USD"))
					SetVarIndexTerm("12M");
				else
					SetVarIndexTerm("6M");
			}
		}
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: CreateEquivalentNotional
///	Returns: void
///	Action : create the equivalent Notional for Bermuda Zero Coupon
/////////////////////////////////////////////////////////////////
ARM_ReferenceValue* ARM_BermudaSwaptionCalculator::CreateEquivalentNotional(double strike, int FixFreq, ARM_ReferenceValue* initNotionalCurve)
{
	double initNotional = initNotionalCurve->GetDiscreteValues()->Elt(0) * 100;
	///// Case of past valuation
	std::vector<double>& fixPayTimes;
	ARM_DateStrip* FixSched;
	if (itsPastZC)
	{
		FixSched = new ARM_DateStrip( itsPastStartDate,  GetEndDate(), GetFixFreq(), GetFixBasis(), GetVarResetCal().c_str(),
			K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTSTART, GETDEFAULTVALUE, GetFixFreq(), GETDEFAULTVALUE,
			GetFixPayCal().c_str());

		fixPayTimes = FixSched->GetFlowEndDates();

	}
	else
	{
		FixSched = new ARM_DateStrip( GetStartDate(),  GetEndDate(), GetFixFreq(), GetFixBasis(), GetVarResetCal().c_str(),
			K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTSTART, GETDEFAULTVALUE, GetFixFreq(), GETDEFAULTVALUE,
			GetFixPayCal().c_str());

		fixPayTimes = FixSched->GetFlowEndDates();
	}

	ARM_Vector* values= new ARM_Vector(fixPayTimes->size(), initNotional);
	ARM_Vector* dates = new ARM_Vector(fixPayTimes->size());

	for (int i=0;i<fixPayTimes->size();++i)
	{
		dates->Elt(i) = fixPayTimes->Elt(i);
	}
	ARM_ReferenceValue* NotionalCurve = new ARM_ReferenceValue(dates,values);
	NotionalCurve->SetCalcMethod(K_STEPUP_RIGHT);
	delete FixSched;

	//// We use the init notional curve to have good dates
	NotionalCurve->GetDiscreteValues()->Elt(0) = initNotional;
	int sizeVector = NotionalCurve->GetDiscreteValues()->size(); 
	if (sizeVector >1)
	{
		for (int i=1; i<sizeVector;++i)
		{
			NotionalCurve->GetDiscreteValues()->Elt(i) = NotionalCurve->GetDiscreteValues()->Elt(i-1)* (1.0 + strike/FixFreq);
		}
	}
	return NotionalCurve;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: CreateEquivalentSwaption
///	Returns: void
///	Action : create the equivalent swaption if not parsed from Summit
/////////////////////////////////////////////////////////////////
void ARM_BermudaSwaptionCalculator::CreateEquivalentSwaption(ARM_ExerciseStyle* exerciseDates)
{
	ARM_ReferenceValue strike = GetStrike();
	ARM_ReferenceValue varSpread = GetSpread();
	ARM_ReferenceValue nominal = GetNominal();

	ARM_ReferenceValue* newStrikes = NULL; 
	ARM_ReferenceValue* newMargin = NULL; 
	// Convert GP strikes to ARM Strikes
	if (GetGenSecType() == 4)
	{
		/// Margin Conversion to the fixFrequency
		ARM_ReferenceValue* newMargin = ConvertVariableMargin(varSpread,GetVarFreq(),GetFixFreq(), GetVarBasis(), GetFixBasis());
		newStrikes = ConvertVariableMarginToStrike(*newMargin);
		for (size_t j=0; j<newStrikes->GetDiscreteValues()->size();++j)
			newStrikes->GetDiscreteValues()->Elt(j) *= 100;
		delete newMargin;
	}
	else
	{
		ARM_Vector* values= new ARM_Vector(strike.GetDiscreteValues());
		ARM_Vector* dates = new ARM_Vector(strike.GetDiscreteDates());
		for (size_t j=0; j<values->size();++j)
			values->Elt(j) *= 100;

		ARM_Vector* marginValues= new ARM_Vector(varSpread.GetDiscreteValues());
		ARM_Vector* marginDates = new ARM_Vector(varSpread.GetDiscreteDates());
			for (j=0; j<marginValues->size();++j)
				marginValues->Elt(j) *= 100;

		newStrikes = new ARM_ReferenceValue(dates,values);
		newMargin = new ARM_ReferenceValue(marginDates,marginValues);

		newStrikes->SetCalcMethod(strike.GetCalcMethod());
		newMargin->SetCalcMethod(varSpread.GetCalcMethod());
	}	
	// Create Equivalent ExoSwaption  for portfolio construction
	ARM_Date startDate = GetStartDate();
	ARM_Date endDate = GetEndDate();
	int payRec = GetPayRec();
	int fixFreq = GetFixFreq();
	int modelFreq = ARM_ArgConv_MatFrequency.GetNumber(itsVarIndexTerm);
	int varFreq;
	if (GetApproxMarginFlag())
		varFreq = modelFreq;
	else
		varFreq = GetVarFreq();
	int fixBasis = GetFixBasis();
	int varBasis = GetVarBasis(); 
	ARM_Currency ccy = GetCcy();
	int stubRule = GetStubRule();

	ARM_FixLeg fixedLeg( startDate, 
						 endDate, 
		                 newStrikes, 
						 payRec,
						 fixFreq, 
						 fixBasis, 
						 K_COMP_PROP, 
						 K_ARREARS, 
						 K_ADJUSTED, 
						 stubRule, 
						 &ccy );

	ARM_SwapLeg floatLeg ( startDate, 
						   endDate, 
						   GetIndexType() ,
						   -payRec,
						   0.0, 
						   varFreq,
						   varFreq,
						   K_ADVANCE, 
						   K_ARREARS, 
						   &ccy,
						   K_ADJUSTED,
						   10000,
						   NULL,
						   NULL,
						   1,
						   K_NX_NONE,
						   GetStubRule(),
						   "NULL",
						   1,
						   varBasis );	
	
	if (GetGenSecType() != 4)
		floatLeg.SetVariableSpread(newMargin);


	ARM_Swap* equivalentSwap = new ARM_Swap( &fixedLeg, &floatLeg);	
	ARM_Swaption* exoSwaption = NULL;

	if (GetGenSecType() != 4)
	{
		// Convert GP nominal to ARM nominal
		ARM_Vector* nominalValues= new ARM_Vector(nominal.GetDiscreteValues());
		ARM_Vector* nominalDates = new ARM_Vector(nominal.GetDiscreteDates());
		for (size_t j=0; j<nominalValues->size();++j)
			nominalValues->Elt(j) *= 100;
		ARM_ReferenceValue newNominal(nominalDates,nominalValues);
		newNominal.SetCalcMethod(nominal.GetCalcMethod());
		newNominal.SetValueType(nominal.GetValueType());
		equivalentSwap->SetAmount(&newNominal);

		exoSwaption = new ARM_Swaption(equivalentSwap, payRec, exerciseDates, newStrikes,0);
		exoSwaption->SetAmount(&newNominal);
	}
	else //// else Bermuda ZC
	{
		ARM_ReferenceValue* notional = CreateEquivalentNotional(strike.GetDiscreteValues()->Elt(0),fixFreq,&nominal);
		equivalentSwap->SetAmount(notional);
		exoSwaption = new ARM_Swaption(equivalentSwap, payRec, exerciseDates, newStrikes,0);
		exoSwaption->SetAmount(notional);
	}

	SetUnderSwap(equivalentSwap);
	SetExoSwaption(exoSwaption);

	if (newStrikes)
		delete newStrikes;

	if (newMargin)
		delete newMargin;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: CreateEmptyCalibration
///	Returns: void
///	Action : create an empty calibration structure
/////////////////////////////////////////////////////////////////
void ARM_BermudaSwaptionCalculator::CreateEmptyCalibration()
{

}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: GetOSWCalibMethod
///	Returns: ARM_CalibMethod
///	Action : get calibration method for diagonal swaptions
/////////////////////////////////////////////////////////////////
ARM_CalibMethod* ARM_BermudaSwaptionCalculator::GetOSWCalibMethod() const
{
  
#ifdef __GP_STRICT_VALIDATION
    if( GetCalibMethod() == ARM_CalibMethodPtr(NULL))
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Optimisation calib method not found");
#endif

    return &(*GetCalibMethod());
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: GetOSWPortfolio
///	Returns: ARM_Portfolio
///	Action : get the diagonal swaptions calibration portfolio
/////////////////////////////////////////////////////////////////
const ARM_StdPortfolioPtr ARM_BermudaSwaptionCalculator::GetOSWPortfolio() const
{
	if (GetMrsCalibFlag() == false || (GetMrsCalibFlag() && GetATMDiagonalFlag()))
    {
		if( GetCalibMethod() == ARM_CalibMethodPtr(NULL) ||
			GetCalibMethod()->GetPortfolio() == ARM_StdPortfolioPtr(NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Optimisation calib method or portfolio not found");

		return GetCalibMethod()->GetPortfolio();
	}
	else
	{
		if( GetCalibMethod() == ARM_CalibMethodPtr(NULL) || GetCalibMethod()->GetlinkedMethod() == NULL ||
        GetCalibMethod()->GetlinkedMethod()->GetPortfolio() == ARM_StdPortfolioPtr(NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Optimisation calib method or portfolio not found");

		if(GetCalibMethod()->GetNextMethod() == NULL)
			return GetCalibMethod()->GetlinkedMethod()->GetPortfolio();
		else
			//// return the portfolio of the final calib Method (Case os standard OSW Calibration for Mean Reversion
			return GetCalibMethod()->GetNextMethod()->GetPortfolio();
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: GetStdPortfolio
///	Returns: ARM_Portfolio
///	Action : get the diagonal swaptions calibration portfolio (case of Mean Reversion Calibration)
/////////////////////////////////////////////////////////////////
const ARM_StdPortfolioPtr ARM_BermudaSwaptionCalculator::GetStdPortfolio() const
{
	if (GetMrsCalibFlag() == false || (GetMrsCalibFlag() && !GetATMDiagonalFlag()))
	{
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : MRS Calibration or ATM Flag is Off! ");
	}
	else
	{
		if( itsPreCalibMethod == ARM_CalibMethodPtr(NULL) ||  itsPreCalibMethod->GetlinkedMethod()->GetPortfolio() == ARM_StdPortfolioPtr(NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Optimisation calib method or portfolio not found");

		//// return the portfolio of the final calib Method (Case os standard OSW Calibration for Mean Reversion
		return itsPreCalibMethod->GetlinkedMethod()->GetPortfolio();
	}
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: SetOSWPortfolio
///	Returns: void
///	Action : affectation of a calibration portfolio made of digonal 
///			 swaptions.
/////////////////////////////////////////////////////////////////
void ARM_BermudaSwaptionCalculator::SetOSWPortfolio(const ARM_StdPortfolio& port)
{
    int pfSize=port.GetSize();
    if(pfSize<1)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Not any asset in the given portfolio");

	if (!GetMrsCalibFlag())
	{
#ifdef __GP_STRICT_VALIDATION
		if( GetCalibMethod() == ARM_CalibMethodPtr(NULL))
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Optimisation calib method not found");
#endif

		/// Update portfolio
		GetCalibMethod()->SetPortfolio( ARM_StdPortfolioPtr((ARM_StdPortfolio*)(const_cast< ARM_StdPortfolio& >(port).Clone())) );
	}
	else
	{
#ifdef __GP_STRICT_VALIDATION
		if( GetCalibMethod() == ARM_CalibMethodPtr(NULL)|| GetCalibMethod()->GetlinkedMethod() == NULL)
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Optimisation calib method not found");
#endif

		/// Update portfolio
		GetCalibMethod()->GetlinkedMethod()->SetPortfolio( ARM_StdPortfolioPtr((ARM_StdPortfolio*)(const_cast< ARM_StdPortfolio& >(port).Clone())) );

	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: GetSTMPortfolio
///	Returns: ARM_Portfolio
///	Action : get short term swaptions calibration portfolio
/////////////////////////////////////////////////////////////////
const ARM_StdPortfolioPtr ARM_BermudaSwaptionCalculator::GetSTMPortfolio() const
{
	if(!GetMrsCalibFlag())
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Mean Reversion Calibration is Off");
    
	if( GetCalibMethod() == ARM_CalibMethodPtr(NULL) ||
        GetCalibMethod()->GetPortfolio() == ARM_StdPortfolioPtr(NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Optimisation calib method or portfolio not found");

	if (GetMrsCalibFlag() && GetATMDiagonalFlag())
		return itsPreCalibMethod->GetPortfolio();
    return GetCalibMethod()->GetPortfolio();
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: SetSTMPortfolio
///	Returns: void
///	Action : affectation of a calibration portfolio made of digonal 
///			 swaptions.
/////////////////////////////////////////////////////////////////
void ARM_BermudaSwaptionCalculator::SetSTMPortfolio(const ARM_StdPortfolio& port)
{
	if(!GetMrsCalibFlag())
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Mean Reversion Calibration is Off");
    
	int pfSize=port.GetSize();
    if(pfSize<1)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Not any asset in the given portfolio");

    if( GetCalibMethod() == ARM_CalibMethodPtr(NULL))
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Optimisation calib method not found");

    /// Update portfolio
    GetCalibMethod()->SetPortfolio( ARM_StdPortfolioPtr((ARM_StdPortfolio*)(const_cast< ARM_StdPortfolio& >(port).Clone())) );
}
/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: CreateDiagonalSwaption
///	Returns: a portfolio
///	Action : create the list of diagonal swaptions
/////////////////////////////////////////////////////////////////
ARM_StdPortfolioPtr ARM_BermudaSwaptionCalculator::CreateOSWPortfolio()
{
    /// Built OSW portfolio and set default weight
	/// The portfolio is buiit in ARM Kernel
    /// A default price is set to pass the validation test in CalibMethod constructor
	//
	ARM_ZeroCurve* pCurve		= dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));
	ARM_BSModel* oswBSModel		= dynamic_cast< ARM_BSModel* >(GetMktDataManager()->GetData(GetKeys()[OswModelKey]));
	
	//ARM_VolCurve* vol			= oswBSModel->GetVolatility();
	//ARM_Portfolio* port		= itsExoSwaption->CalibrationPF(pCurve, vol, GetOSWPortfolioMode());
	ARM_Portfolio* port			= itsExoSwaption->CalibrationPF(pCurve, oswBSModel, GetOSWPortfolioMode());

    ARM_StdPortfolio* stdPort	= new ARM_StdPortfolio (*port);
	
	delete port; 
	port = NULL;
    
	return ARM_StdPortfolioPtr(stdPort);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: CreateDiagonalSwaption 
///	Returns: a portfolio
///	Action : create the list of std swaptions for MRS calibration
/// We determine the level of mean reversion using standard swaptions 
/// and then calibrate the model using the portfolio of diagonal swaptions  
/////////////////////////////////////////////////////////////////
ARM_StdPortfolioPtr ARM_BermudaSwaptionCalculator::CreateStdOSWPortfolio()
{
    /// Built OSW portfolio and set default weight
	/// The portfolio is buiit in ARM Kernel
    /// A default price is set to pass the validation test in CalibMethod constructor
	//
	ARM_ZeroCurve* pCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));
	ARM_BSModel* oswBSModel = dynamic_cast< ARM_BSModel* >(GetMktDataManager()->GetData(GetKeys()[OswModelKey]));
	ARM_VolCurve* vol = oswBSModel->GetVolatility();

	///// Test : We keep constant nominal to have good values of Mean Reversion
//	ARM_Vector* values= new ARM_Vector(1, 100.0);
//	ARM_Vector* dates = new ARM_Vector(GetNominal().GetDiscreteDates()->Elt(0));
//	ARM_ReferenceValue*	notional = new ARM_ReferenceValue(dates,values);
//	notional->SetCalcMethod(K_STEPUP_LEFT);

	ARM_Portfolio* port = itsExoSwaption->StdCalibrationPF(pCurve, oswBSModel, &GetNominal());//ComputeMeanNotionalCurve(GetNominal(),2));
	//ARM_Portfolio* port = itsExoSwaption->StdCalibrationPF(pCurve, oswBSModel, ComputeMeanNotionalCurve((*notional),2));
    ARM_StdPortfolio* stdPort = new ARM_StdPortfolio (*port);
	/// Conversion from Portfolio to StdPortfolio
	

	/////  try to price forward swap options with an analytic model
	ARM_MarketIRModel* marketModel = dynamic_cast<ARM_MarketIRModel*>(GetMktDataManager()->GetData(GetKeys()[NormalModel]));
	size_t sizePF = stdPort->GetSize();
	for (int i=0; i<sizePF; ++i)
	{
		ARM_Swaption* swaption = dynamic_cast<ARM_Swaption*> (stdPort->GetAsset(i));
		ARM_VanillaArg* swaptionArg = ARM_ConverterFromKernel::ConvertSecuritytoArgObject( swaption, GetMktDataManager()->GetAsOfDate().GetJulian(), "TOTO");
		double price = swaptionArg->Price(marketModel);
		stdPort->SetPrice(price, i);
		stdPort->SetWeight(1, i);
		delete swaptionArg;
	}
	
	delete port; 
	port = NULL;
    return ARM_StdPortfolioPtr(stdPort);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: CreateSTMPortfolio
///	Returns: a portfolio
///	Action : create STM Portfolio for Mean Reversion Calibration 
/////////////////////////////////////////////////////////////////
ARM_StdPortfolioPtr ARM_BermudaSwaptionCalculator::CreateSTMPorftolio(ARM_StdPortfolioPtr OSWPortfolio)
{
    /// Built STM portfolio and set default weight
	/// The portfolio is buiit in ARM Kernel
    /// A default price is set to pass the validation test in CalibMethod constructor

	//Update the function with Calib Flags
	ARM_ZeroCurve* pCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));;
	ARM_BSModel* bsModel = NULL;
	if ((GetMrsCalibFreq() < 2) || (GetMrsCalibMode()!= Column))
	{
		bsModel= dynamic_cast< ARM_BSModel* >(GetMktDataManager()->GetData(GetKeys()[OswModelKey]));
	}
	else
	{
		/// Cap BS Model
		bsModel = dynamic_cast< ARM_BSModel* >(GetMktDataManager()->GetData(GetKeys()[CfModelKey]));
	}

	/// To be changed with new inputs
	ARM_Portfolio* port = itsExoSwaption->MRSCalibrationPF(pCurve, bsModel, &*OSWPortfolio, GetMrsCalibFreq(), GetAtmMrsCalibType(), GetMrsCalibMode(), GetVolatilityLag());
    ARM_StdPortfolio* stdPort = new ARM_StdPortfolio (*port);
	/// Conversion from Portfolio to StdPortfolio

	delete port; 
	port = NULL;
	
	if (GetMrsCalibMode() == Correl)
	{
		/////  try to price forward swap options with an analytic model
		ARM_MarketIRModel* marketModel = dynamic_cast<ARM_MarketIRModel*>(GetMktDataManager()->GetData(GetKeys()[NormalModel]));
		size_t sizePF = stdPort->GetSize();
		for (int i=0; i<sizePF; ++i)
		{
			ARM_Swaption* swaption = dynamic_cast<ARM_Swaption*> (stdPort->GetAsset(i));
			ARM_VanillaArg* swaptionArg = ARM_ConverterFromKernel::ConvertSecuritytoArgObject( swaption, GetMktDataManager()->GetAsOfDate().GetJulian(), "TOTO");
			double price = swaptionArg->Price(marketModel);
			stdPort->SetPrice(price, i);
			stdPort->SetWeight(1, i);
			delete swaptionArg;
		}
	}
    return ARM_StdPortfolioPtr(stdPort);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: ComputeControlVariable
///	Returns: double
///	Action : computes the price of the control variable (before compute the estimated strikes)
/////////////////////////////////////////////////////////////////
void ARM_BermudaSwaptionCalculator::ComputeControlVariate(std::vector<double>& cvprices)
{

}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: UpdateCalibration
///	Returns: void
///	Action : update the calibration datas w.r.t. MarketDataManager
///          The context is an hedge ratios computation
/////////////////////////////////////////////////////////////////
void ARM_BermudaSwaptionCalculator::UpdateCalibration(bool isUpdateStrike)
{
	CreateAndSetCalibration();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: MeanReversionRoot
///	Returns: void
///	Action : calculate the mean reversion corrresponding to a given price
/////////////////////////////////////////////////////////////////
double ARM_BermudaSwaptionCalculator::MeanReversionRoot(double targetPrice, double fTolerance, int maxIter)
{
	double initialMrs = GetInitMrs();

	if (GetMrsCalibFlag() == true )
	{
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : MRS Calibration flag should be off! ");
	}

	ARM_MeanRevFinder func(const_cast<ARM_BermudaSwaptionCalculator*> (this), initialMrs);

	UnaryFuncWithNumDerivative<double> funcWithDeriv(func);

	T_NewtonRaphsonSolver<UnaryFuncWithNumDerivative<double> > solver(funcWithDeriv, targetPrice, fTolerance, SolverConstant::DefaultFxTolerance, maxIter);

	solver.setInitialGuess(initialMrs);
	
	double targetMrs = solver.Solve();
	
	return targetMrs;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: ComputePricingData
///	Returns: a ARM_MultiTypeDict
///	Action : get pricing data of the calculaltor
/////////////////////////////////////////////////////////////////
void ARM_BermudaSwaptionCalculator::ComputePricingData() const
{
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: CreateAndSetModel
///	Returns: void
///	Action : creates the model
/////////////////////////////////////////////////////////////////
void ARM_BermudaSwaptionCalculator::CreateAndSetModel()
{
	/// Build the default stochastic model of the calculator : SFRM 1F
	ARM_ZeroCurve* pCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));

	/// Creates default values for volatility & mean reversion
    /// (calibration will be called at pricing time) and set them
    /// in the reference model
	std::vector<double> defaultTimes(1,0.0);
	std::vector<double> defaultSigmas(1,GetInitVol());
	std::vector<double> defaultMeanReversion(1,GetInitMrs());

	ARM_PricingModelPtr refModel(NULL);

	///// Create the Model
    if (GetModelType()== ARM_PricingModelType::SFRM1F || GetModelType()==ARM_PricingModelType::SFRM2F)
	{
		std::vector<double> defaultBetas(1,GetInitBeta());
		ARM_ModelParamVector sfrmParamVect;
		ARM_CurveModelParam  paramvol(ARM_ModelParamType::Volatility, &defaultSigmas, &defaultTimes );
		sfrmParamVect.push_back(&paramvol);
		ARM_CurveModelParam  paramMRS(ARM_ModelParamType::MeanReversion, &defaultMeanReversion, &defaultTimes );
		sfrmParamVect.push_back(&paramMRS);
		ARM_CurveModelParam  paramBeta(ARM_ModelParamType::Beta, &defaultBetas, &defaultTimes );
		sfrmParamVect.push_back(&paramBeta);                
		ARM_INDEX_TYPE liborType = GetIndexType();
		ARM_IRIndex irIndex(liborType,K_DEF_FREQ,K_DEF_FREQ,GetCurrencyUnit());
		
        if (GetModelType()== ARM_PricingModelType::SFRM1F)
		{
			ARM_ModelParamsSFRM* sfrmParams = ARM_ModelParamsSFRMFactory.Instance()->CreateModelParamsSFRM(sfrmParamVect,&irIndex,1,K_ROW);
			CC_NS(std,auto_ptr)<ARM_ModelParamsSFRM> HoldSFRMParams(sfrmParams);
			
			/// Build the default stochastic model of the calculator : SFRM 1F
			
			ARM_SFRM* sfrm = new ARM_SFRM( CreateClonedPtr( pCurve ), *sfrmParams);
			
			//We adjust the sfrm first start date in a specific case   
			int callFreq					= GetCallFreq();
			ARM_DateStripPtr callDateStrip	= GetCallDateStrip(); 	
			double firstResetDate			= (*(callDateStrip->GetResetDates()))[0];
			double firstStartDate			= (*(callDateStrip->GetFlowStartDates()))[0];
			ARM_Date newStartDate			= ARM_Date(firstStartDate);
			double nbPeriods				= (firstStartDate - firstResetDate) / K_YEAR_LEN * callFreq;   

			if (nbPeriods > 1)
			{	
				int nbMonths = 0;
				int counter = (int) nbPeriods;

				if (callFreq ==	K_MONTHLY)
				{
					nbMonths = 1;
				}
				else if (callFreq == K_QUARTERLY)
				{
					nbMonths = 3;
				}
				else if (callFreq == K_SEMIANNUAL)
				{
					nbMonths = 6;
				}
				else if (callFreq == K_ANNUAL)
				{
					nbMonths = 12;		
				}

				while (counter != 0)
				{	
					newStartDate.AddMonths(-nbMonths);
					counter--;
				}
			}
			
			sfrm->SetFixStartDate(newStartDate);

			refModel = ARM_PricingModelPtr(ARM_PricingModelPtr(static_cast< ARM_PricingModel* >(sfrm)) );
		}
		else
		{
			std::vector<double> correlLowerBound(1,CORREL_LOWER_BOUND);
			std::vector<double> correlUpperBound(1,CORREL_UPPER_BOUND);
			ARM_DateStrip datestrip( GetStartDate(),
									 GetEndDate(),
									 CORREL_FREQ,
									 GetVarBasis(),
									 GETDEFAULTVALUESTR,
									 K_MOD_FOLLOWING,
									 K_ADJUSTED,
									 K_SHORTSTART,
									 GETDEFAULTVALUE,
									 CORREL_FREQ );

			ARM_CorrelTrigoMatParam correl = ARM_CorrelTrigoMatParam(GetInitTheta(),
														datestrip,
														GetMktDataManager()->GetAsOfDate().GetJulian(),
														"LINEAR",
														&correlLowerBound,
														&correlUpperBound);
			sfrmParamVect.push_back(((ARM_CurveModelParam*) &correl));
			ARM_ModelParamsSFRM* sfrmParams = ARM_ModelParamsSFRMFactory.Instance()->CreateModelParamsSFRM(sfrmParamVect,&irIndex,2,K_ROW);
			CC_NS(std,auto_ptr)<ARM_ModelParamsSFRM> HoldSFRMParams(sfrmParams);
			
			/// Build the default stochastic model of the calculator : SFRM 1F
			ARM_SFRM* sfrm = new ARM_SFRM( CreateClonedPtr( pCurve ), *sfrmParams);
			sfrm->SetFixStartDate(GetStartDate());
			refModel = ARM_PricingModelPtr(ARM_PricingModelPtr(static_cast< ARM_PricingModel* >(sfrm)) );

		}
	}
    else if(GetModelType()==ARM_PricingModelType::QGM1F)
	{
		ARM_CurveModelParam  paramvol(ARM_ModelParamType::Volatility, &defaultSigmas, &defaultTimes );
		ARM_CurveModelParam  paramMRS(ARM_ModelParamType::MeanReversion, &defaultMeanReversion, &defaultTimes );
		std::vector<double> defaultSkews(1,GetInitSkew());
		ARM_CurveModelParam paramSkew( ARM_ModelParamType::Skew, &defaultSkews, &defaultTimes );

		ARM_ModelParamVector paramVector(3);
		paramVector[0] = &paramvol;
		paramVector[1] = &paramMRS;
		paramVector[2] = &paramSkew;

		refModel = ARM_PricingModelPtr( ARM_PricingModelPtr (static_cast< ARM_PricingModel* >(new ARM_QGM1F( CreateClonedPtr( pCurve), paramVector ) ) ) );
	}
    else if(GetModelType()==ARM_PricingModelType::HWM1F)
	{
		ARM_CurveModelParam  paramvol(ARM_ModelParamType::Volatility, &defaultSigmas, &defaultTimes );
		ARM_CurveModelParam  paramMRS(ARM_ModelParamType::MeanReversion, &defaultMeanReversion, &defaultTimes );

		ARM_ModelParamVector paramVector(2);
		paramVector[0] = &paramvol;
		paramVector[1] = &paramMRS;

		refModel= ARM_PricingModelPtr( ARM_PricingModelPtr (static_cast< ARM_PricingModel* >(new ARM_HullWhite1F( CreateClonedPtr( pCurve )) ) ) );
		refModel->SetModelParams( ARM_ModelParamsHW1FStd(paramVector) );
	}
	////// Attach the Numerical Method
	if (GetNumMethodType() == Tree)
	{
		std::vector<double> schedulerDatas(3);
		schedulerDatas[0] = GetTreeSteps();
		schedulerDatas[1] = 1;
		schedulerDatas[2] = 1.0e-2;
		std::vector<double> samplerDatas(1,1.0e-3);
		std::vector<double> truncatorDatas(1,5.0);
		int reconnectorType=ARM_ReconnectorBase::Mean;
		int smootherType=ARM_SmootherBase::DoNothing;
		int truncatorType=ARM_TruncatorBase::StandardDeviation;
		int schedulerType,samplerType;
        if(GetModelType()== ARM_PricingModelType::SFRM1F)
		{
			schedulerType = ARM_SchedulerBase::ConstantVariance;			
			samplerType=ARM_SamplerBase::NormalCentred;
		}
		else //// Mean Reverting Tree
		{
			schedulerType=ARM_SchedulerBase::ConstantVarianceMeanReverting;
			samplerType=ARM_SamplerBase::MeanReverting;		
		}
		ARM_TreeBase* tree = ARM_TreeFactory.Instance()->CreateTreeND(GetFactorNb(),schedulerType,schedulerDatas,
				samplerType,samplerDatas,truncatorType,truncatorDatas,false,reconnectorType,smootherType);
		refModel->SetNumMethod( ARM_NumMethodPtr( tree ) );
	}
	else if (GetNumMethodType() == PDE)
	{

		int numSchemeType = ARM_ArgConv_PDENumSchemeType.GetNumber("CN1F");
		ARM_PDENumericalScheme* numScheme = ARM_PDENumericalScheme::getNumericalSchemeInstanceById(numSchemeType);

		int nbSteps   = 1000;
		int gridSize  = 501;
		double nbStdDevs = 10;

		ARM_PDEMethod* pPDEMethod = new ARM_PDEMethod(numScheme,nbSteps, gridSize, 2*nbStdDevs);
		
		refModel->SetNumMethod( ARM_NumMethodPtr( pPDEMethod ) );
	}
	///// else AMC Method
	else
	{	
		ARM_RandomGeneratorPtr numRandGen( ARM_RandGenFactory.Instance()->CreateSimpleRandGen( 
				itsGenType1,
				itsGenType2,
				"BoxMuller",
				"InvNormCum",
				itsPathOrder,
				itsFirstNbTimes,
				1,
				true));

		ARM_TimeStepPerYearScheduler scheduler(0);
		ARM_NormalCentredSamplerND sampler(&scheduler);

		// Path Scheme
		int pathSchemeType = ARM_ArgConv_PathSchemeType.GetNumber(itsPathScheme);
		ARM_PathSchemePtr pathScheme(ARM_PathSchemeFactory.Instance()->CreatePathScheme(pathSchemeType));

		ARM_ExerciseBoundaryCalc* exerBoundCal = new ARM_AMCAndersen(GetAmcIter(),1);
		ARM_AMCMethod* amcMethod = new ARM_AMCMethod(GetMcIter(),numRandGen,&sampler,exerBoundCal,GetMaxBucketSize(),ARM_ImpSamplerPtr(NULL),pathScheme);
		refModel->SetNumMethod( ARM_NumMethodPtr( amcMethod ) );

		delete exerBoundCal;
	} 
	
    /// Create a Numeraire and set it
    ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::TerminalZc ) );
    refModel->SetNumeraire(numeraire);

	/// Set the model
	SetPricingModel(refModel);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: CreateSwaptionPortfolio_Frontier
///	Returns: ARM_StdPortfolioPtr
///	Action : create swaption portfolio
/////////////////////////////////////////////////////////////////
ARM_StdPortfolioPtr ARM_BermudaSwaptionCalculator::SwaptionPortfolio_Frontier()
{
	ARM_DateStripCombiner datesStructure = DatesStructure();	
	std::vector<double>& resetDates	= datesStructure.GetDateStrip(CALL_BS_SCHED)->GetResetDates();
	std::vector<double>& startDates	= datesStructure.GetDateStrip(CALL_BS_SCHED)->GetFlowStartDates();
	
	ARM_BSModel* BSModel		= dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswModelKey]) );
	double asOfDate				= GetMktDataManager()->GetAsOfDate().GetJulian();

	int i,size = startDates->size();
	list< ARM_Security* > swaptionList;

	double fees;
	for (i = 0; i < size; ++i)
	{
		ARM_Date startDate ((*startDates)[i]);

		fees = itsFees.Interpolate((*resetDates)[i]);

		if (((*resetDates)[i] > asOfDate) && (fees < NON_CALL_FEE))
		{
			ARM_Swap stdSwap(startDate,
							itsEndDate,
							GetIndexType(),
							0.0,
							1.0,
							GetPayRec(),
							K_DEF_FREQ,
							K_DEF_FREQ,
							GetCurrencyUnit());

			ARM_Date expiryDate((*(stdSwap.GetFloatLeg()->GetResetDates()))[0]);
			stdSwap.SetModel(BSModel);

			double atmStrike = stdSwap.PriceToRate((ARM_Date) asOfDate, 0.0);
					
			ARM_Swaption swaption((&stdSwap),K_RCV,K_EUROPEAN,atmStrike,expiryDate);
					
			swaptionList.push_back(static_cast<ARM_Swaption*>(swaption.Clone()));
		}
	}

	ARM_StdPortfolioPtr swaptionPF(new ARM_StdPortfolio(swaptionList));
    	
	for(i=0;i<swaptionPF->size();++i)
	{
		ARM_Swaption* swopt = static_cast< ARM_Swaption* >(swaptionPF->GetAsset(i));
		swopt->SetModel(BSModel);
		double price = swopt->ComputePrice();
		swaptionPF->SetPrecision(DEFAULT_PRECISION,i);
        swaptionPF->SetWeight(DEFAULT_WEIGHT,i);
		swaptionPF->SetPrice(price,i);
	}
	return swaptionPF;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: CreateAndSetCalibration_Frontier
///	Returns: void
///	Action : create the calibration
/////////////////////////////////////////////////////////////////
void ARM_BermudaSwaptionCalculator::CreateAndSetCalibration_Frontier()
{
	ARM_BSModel* BSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswModelKey]) );

	std::vector<double> defaultTimes(1,0.0);
	std::vector<double> defaultSigmas(1,GetInitVol());
	std::vector<double> volLowerBound(1,SIGMA_LOWER_BOUND);
	std::vector<double> volUpperBound(1,SIGMA_UPPER_BOUND);

	ARM_CurveModelParam volatility = ARM_CurveModelParam(ARM_ModelParamType::Volatility,&defaultSigmas,&defaultTimes,
			"VOLATILITY","STEPUPRIGHT",&volLowerBound,&volUpperBound);
	ARM_ModelParamVector volatilityVector(1,&volatility);


// PreCalib ATM
	ARM_StdPortfolioPtr swaptionPF = SwaptionPortfolio_Frontier();

	ARM_CalibMethod volCalibAtm(swaptionPF,volatilityVector,ARM_CalibMethodType::Bootstrap1D,20,ARM_CalibrationTarget::PriceTarget,NULL,NULL,true);
	SetCalibMethod( ARM_CalibMethodPtr(static_cast< ARM_CalibMethod* >(volCalibAtm.Clone()) ) );

	GetCalibMethod()->Calibrate(&*GetPricingModel());

// Price and frontier
    ARM_GenSecurityPtr genSec = GetGenSecurity();
	ARM_GenPricer* genPricer = NULL;
	genPricer = new ARM_GenPricer( &*genSec,&*GetPricingModel() );
	ARM_AutoCleaner<ARM_GenPricer> HoldGP(genPricer );
	genPricer->Price();

	double price = genPricer->GetPricerInfo()->GetContents("Bermuda").GetData("Price").GetDouble();
	ARM_VectorPtr frontier = genPricer->GetPricerInfo()->GetContents("Frontier").GetData("Intermediateprices").GetVector();

// Set Calib at frontier
	for(int i=0;i<swaptionPF->size();++i)
	{
		ARM_Swaption* swopt = static_cast< ARM_Swaption* >(swaptionPF->GetAsset(i));
		double old = swopt->GetStrike();
		swopt->UpdateStrike((*frontier)[i]*100.);
		double price = swopt->ComputePrice();
		double vega = swopt->ComputeSensitivity(K_VEGA)/100.;
		if (vega<VEGA_MIN)
		{
			swopt->UpdateStrike(old);
			price = swopt->ComputePrice();
		}
		swaptionPF->SetPrecision(DEFAULT_PRECISION,i);
        swaptionPF->SetWeight(DEFAULT_WEIGHT,i);
		swaptionPF->SetPrice(price,i);
	}

	ARM_CalibMethod volCalibFrontier(swaptionPF,volatilityVector,ARM_CalibMethodType::Bootstrap1D,20,ARM_CalibrationTarget::PriceTarget,NULL,NULL,true);
	SetCalibMethod( ARM_CalibMethodPtr(static_cast< ARM_CalibMethod* >(volCalibFrontier.Clone()) ) );

}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: CreateAndSetCalibration
///	Returns: void
///	Action : create the calibration
/////////////////////////////////////////////////////////////////
void ARM_BermudaSwaptionCalculator::CreateAndSetCalibration()
{
	if (GetOSWPortfolioMode()==4 || GetOSWPortfolioMode()==5)
	{
		CreateAndSetCalibration_Basket();
	}
	else if (GetOSWPortfolioMode()==6)
	{
		CreateAndSetCalibration_Frontier();
	}
	else
	{
		//// Case of Mean Reversion calibration using the Flag ATMDiagonal

		/// Create standard diagonal swaption portfolio (sigma calibration)
		ARM_StdPortfolioPtr diagonalSwaptionPF(NULL);
		
		if (GetATMDiagonalFlag() && GetMrsCalibFlag())
		{
			diagonalSwaptionPF = ARM_StdPortfolioPtr(CreateStdOSWPortfolio());
		}
		else
			diagonalSwaptionPF = ARM_StdPortfolioPtr(CreateOSWPortfolio());

		std::vector<double> defaultTimes(1,0.0);
		std::vector<double> defaultSigmas(1,GetInitVol());
		std::vector<double> volLowerBound(1,SIGMA_LOWER_BOUND);
		std::vector<double> volUpperBound(1,SIGMA_UPPER_BOUND);

		ARM_CurveModelParam volatility = ARM_CurveModelParam(ARM_ModelParamType::Volatility,&defaultSigmas,&defaultTimes,
				"VOLATILITY","STEPUPRIGHT",&volLowerBound,&volUpperBound);
		ARM_ModelParamVector volatilityVector(1,&volatility);

		/// Build a volatility bootstrap calibration on the latter portfolio
		/// The CalibMethod object will be cloned because by default it is not shared
		ARM_CalibMethod volCalib(diagonalSwaptionPF,volatilityVector,ARM_CalibMethodType::Bootstrap1D,20,ARM_CalibrationTarget::PriceTarget,NULL,NULL,true);
		
		if (dynamic_cast<ARM_QGM1F*> (&*GetPricingModel()))
		{
			volCalib.SetNbIteration(QGM_VOL_CALIB_ITER);
		}
		if(GetMrsCalibFlag())
		{
				//// Model Param for Mean Reversion Calibration
			std::vector<double> mrsLowerBound(1,GetInitMrsL());
			std::vector<double> mrslUpperBound(1,GetInitMrsU());
			std::vector<double> initMrs (1,GetInitMrs());		
			ARM_CurveModelParam mrsCurve = ARM_CurveModelParam(ARM_ModelParamType::MeanReversion,&initMrs,&defaultTimes,
			"MEANREVERSION","STEPUPRIGHT",&mrsLowerBound,&mrslUpperBound);
			ARM_ModelParamVector mrsVector(1,&mrsCurve);

			/// Create standard short term vanillas portfolio (MRS calibration)
			ARM_StdPortfolioPtr shortTermVanillaPF(CreateSTMPorftolio(diagonalSwaptionPF));

			/// Build a MRS optimisation and embed the volatility bootstrapping
			/// The CalibMethod object will be cloned because by default it is not shared

			ARM_CalibMethod* mrsCalib = new ARM_CalibMethod(shortTermVanillaPF,mrsVector,ARM_CalibMethodType::Optimize,
									 ARM_MAX_ITER,ARM_CalibrationTarget::PriceTarget,&volCalib);

			if (GetATMDiagonalFlag())
			{
				SetPreCalibMethod(ARM_CalibMethodPtr(mrsCalib));
				volCalib.SetPortfolio(CreateOSWPortfolio());
				SetCalibMethod( ARM_CalibMethodPtr(static_cast< ARM_CalibMethod* >(volCalib.Clone()) ) );   
			}
			else	
				SetCalibMethod( ARM_CalibMethodPtr(mrsCalib) );

			GetCalibMethod()->SetNextMethod(static_cast< ARM_CalibMethod* >(volCalib.Clone()));
			if (GetATMDiagonalFlag())
				GetCalibMethod()->GetNextMethod()->SetPortfolio(CreateOSWPortfolio());
			else
				GetCalibMethod()->GetNextMethod()->SetPortfolio(GetOSWPortfolio());
		}
		else
			SetCalibMethod( ARM_CalibMethodPtr(static_cast< ARM_CalibMethod* >(volCalib.Clone()) ) );    
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: CreateAndSetCalibration_Basket
///	Returns: void
///	Action : create the calibration
/////////////////////////////////////////////////////////////////
void ARM_BermudaSwaptionCalculator::CreateAndSetCalibration_Basket()
{
	ARM_BSModel* oswBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswModelKey]) );
	ARM_MarketIRModel* marketModel = dynamic_cast<ARM_MarketIRModel*>(GetMktDataManager()->GetData(GetKeys()[NormalModel]));

	ARM_BasketCalib::BasketStrike strike = (GetOSWPortfolioMode()==4?ARM_BasketCalib::ATM:ARM_BasketCalib::EQUIVALENT);
	
	vector<ARM_Security*> securities;
	securities.push_back(itsExoSwaption->GetFixedLeg());//GetUnderSwap()->GetFirstLeg());
	securities.push_back(itsExoSwaption->GetFloatLeg());//GetUnderSwap()->GetSecondLeg());

	vector<ARM_Model*> models;
	models.push_back(oswBSModel);
	models.push_back(oswBSModel);

	std::vector<double> weights;
	weights.push_back(-1.);
	weights.push_back(1.);

	double AsOf(GetMktDataManager()->GetAsOfDate().GetJulian());
	ARM_ReferenceValue fixNotional(itsExoSwaption->GetFixedLeg()->GetAmount()->Interpolate(AsOf));//GetUnderSwap()->GetFirstLeg()->GetAmount()->Interpolate(AsOf));

	ARM_BasketCalib basket(	&*GetCallDateStrip(), *itsExoSwaption->GetFixedLeg()->GetAmount(),//*GetUnderSwap()->GetFirstLeg()->GetAmount(),
							GetFees(), 1, (ARM_BasketCalib::FULL),strike);
	basket.Compute(securities,models,weights);
	basket.Price(marketModel);

	std::vector<double> defaultTimes(1,0.0);
	std::vector<double> defaultSigmas(1,GetInitVol());
	std::vector<double> volLowerBound(1,SIGMA_LOWER_BOUND);
	std::vector<double> volUpperBound(1,SIGMA_UPPER_BOUND);

	ARM_CurveModelParam volatility = ARM_CurveModelParam(ARM_ModelParamType::Volatility,&defaultSigmas,&defaultTimes,
			"VOLATILITY","STEPUPRIGHT",&volLowerBound,&volUpperBound);
	ARM_ModelParamVector volatilityVector(1,&volatility);

	ARM_CalibMethod volCalib(basket.GetPortfolio(),volatilityVector,ARM_CalibMethodType::Bootstrap1D,20,ARM_CalibrationTarget::PriceTarget,NULL,NULL,true);
	SetCalibMethod( ARM_CalibMethodPtr(static_cast< ARM_CalibMethod* >(volCalib.Clone()) ) );
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: UpdateModel
///	Returns: void
///	Action : update the model datas w.r.t. MarketDataManager
/////////////////////////////////////////////////////////////////
void ARM_BermudaSwaptionCalculator::UpdateModel()
{
	SetInitMrs(GetDefaultMrs());
	CreateAndSetModel();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: Calibrate
///	Returns: void
///	Action : calibrate the model 
/////////////////////////////////////////////////////////////////
void ARM_BermudaSwaptionCalculator::Calibrate()
{
	//// We do the calibration process in 2 steps to obtain the mean reversion level and clear the model params schedule
	if (GetATMDiagonalFlag() && GetMrsCalibFlag())
	{
		GetPreCalibMethod()->Calibrate(&*GetPricingModel());
		double initMrs = ((ARM_CurveModelParam&) GetPricingModel()->GetModelParams()->GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];
		SetInitMrs(initMrs);
		CreateAndSetModel();
	}
	GetCalibMethod()->Calibrate(&*GetPricingModel());   
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: PortfolioSwaptionPrice
///	Returns: a double
///	Action : price the Bermuda Swaption deal.
/////////////////////////////////////////////////////////////////
double ARM_BermudaSwaptionCalculator::GetCVPrice(int cvi)
{
	ARM_StdPortfolioPtr swaptionPortfolio = GetOSWPortfolio();

	double price = swaptionPortfolio->GetAsset(cvi-1)->GetPrice()/100;

    return price;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: CVPrices
///	Returns: a double
///	Action : Returns a vector of Control Variate Prices.
/////////////////////////////////////////////////////////////////
std::vector<double> ARM_BermudaSwaptionCalculator::GetCVPrices()
{
	ARM_StdPortfolioPtr swaptionPortfolio = GetOSWPortfolio();
	size_t nbCv = GetNbControlVariate();
	std::vector<double> CVPrices(nbCv,0.0);
	int varIndexFreq = ARM_ArgConv_MatFrequency.GetNumber(itsVarIndexTerm);
	for (size_t i=0; i<nbCv; ++i)
	{
		(CVPrices)[i] = swaptionPortfolio->GetAsset(((*GetControlVariate()))[i]/**varIndexFreq/itsCallFreq*/)->GetPrice()/100;
	}
    return CVPrices;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: Price
///	Returns: a double
///	Action : price the Bermuda Swaption deal.
/////////////////////////////////////////////////////////////////
double ARM_BermudaSwaptionCalculator::Price()
{
	/// Calibrate the model
	CalibrateAndTimeIt();

	/// Price the implicit product according to internal flag
    ARM_GenSecurityPtr genSec = GetGenSecurity();
	
	ARM_GenPricer* genPricer = NULL;
	double price;
	int nbProductsToPrice;	
	int genSecType = GetGenSecType();
	if (GetFixBoundaryFlag())
		GetGenSecurity()->SetExercBoundaryResetFlag(false);
	if ((genSecType == 2) || (genSecType == 4))
	{
		genPricer = new ARM_GenPricer( &*genSec,&*GetPricingModel() );

		nbProductsToPrice = 1;
		ARM_AutoCleaner<ARM_GenPricer> HoldGP(genPricer );
		genPricer->Price();
		price	= genPricer->GetPricerInfo()->GetContents( "Bermuda").GetData("Price").GetDouble();
		if (GetFixBoundaryFlag() && (GetNumMethodType() == MonteCarlo))
		{
			ARM_GenSecManipulator genSecManipulator;
			ARM_DealDescriptionPtr dealDesc( new ARM_DealDescription(GetGenSecurity()->GetDealDescription()));
			genSecManipulator.ChangeAmericanIntoTrigger( *GetGenSecurity(), dealDesc );
			ARM_GenSecurityPtr newGenSec = ARM_GenSecurityPtr(new ARM_GenSecurity(dealDesc,"",GetGenSecurity()->GetCstManager(),!GetFixBoundaryFlag(),GetGenSecurity()->GetOtherPayoffsFlag()));
			SetGenSecurity( newGenSec );
			itsFixBoundaryFlag = false;
		}
		GetPricingData()["BermudaPrice"] = price;
		SetBermudaSwaptionPrice(price);
		if (GetNumMethodType() == MonteCarlo)
		{
			double stdDev;
			stdDev	= genPricer->GetPricerInfo()->GetContents( "Bermuda").GetData("StdDev").GetDouble();
			GetPricingData()["BermudaStdDev"] = stdDev;
			SetBermudaSwaptionStdDev(stdDev);
			itsBermudaSwaptionRA = genPricer->GetPricerInfo()->GetContents( "Bermuda").GetData("RunningAvg").GetVector();
			GetPricingData()["BermudaRunningAvg"] = itsBermudaSwaptionRA;
			
		}
		if (itsCalculateProbaFlag && (genSecType == 2))
		{
			ARM_ZeroCurve* pCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));
			double DF = pCurve->DiscountPrice(((GetEndDate().GetJulian())-(GetMktDataManager()->GetAsOfDate().GetJulian()))/K_YEAR_LEN);
			for(size_t i = 0; i < itsCumProba->size(); i++)
			{
				itsCumProba->Elt(i) = genPricer->GetPricerInfo()->GetContents( (BermudaSwaptionProbaPricedColNamesTable2[i].c_str())).GetData("IntermediatePrices").GetVector()->Elt(0)/ DF;
			}
			GetPricingData()["Proba"] = itsCumProba;


			///// Compute weights for each line of the matrix P(No 1), P(No1 && No2),......
			std::vector<double> wRows(itsProbaMatrix.rows(),0.0);
			wRows[0] = (*itsCumProba)[0] ;
			double normR = wRows[0];
			for(int h = 1; h < itsCumProba->size(); ++h)
			{
				wRows[h] = ((*itsCumProba)[h] - (*itsCumProba)[h-1]);
			}
			for(h = 0; h < itsCumProba->size();++h)
			{
				wRows[h] /= normR;
			}

			//// Fill In Matrix Weights P(2/No1)*P(No 1) , P(No 2 && 3/ No1)*P(No 1).....
			//// Ans Transpose  the matrix
			int vProbaMatrixRow = itsProbaMatrix.rows() - 2;
			for(int j = 0 ; j < vProbaMatrixRow;j++)
			{
				double w0 = 1 - (genPricer->GetPricerInfo()->GetContents( (BermudaSwaptionProbaPricedColNamesTable2[j].c_str())).GetData("IntermediatePrices").GetVector()->Elt(j))/DF;
				for(int l=j;l<vProbaMatrixRow;l++)
				{
					double wl1= genPricer->GetPricerInfo()->GetContents( (BermudaSwaptionProbaPricedColNamesTable2[l].c_str())).GetData("IntermediatePrices").GetVector()->Elt(j);
					double wl2 = genPricer->GetPricerInfo()->GetContents( (BermudaSwaptionProbaPricedColNamesTable2[l+1].c_str())).GetData("IntermediatePrices").GetVector()->Elt(j);
					itsProbaMatrix.Elt(j,l-(j)) = (wl2-wl1); 
				}
				double norm = w0 / wRows[j];
				for (l=0;l<(vProbaMatrixRow-(j));l++)
				{
					itsProbaMatrix.Elt(j,l)/=norm;
				}
			}

			//// Renorm the Matrix
			double norm = itsProbaMatrix.Elt(0,0);
			for(j=0; j < vProbaMatrixRow; j++)
				for(int l=0; l < vProbaMatrixRow + 1; l++)
						itsProbaMatrix.Elt(j,l)/=norm;
		}
		
	}
	else if ((genSecType == 1) || (genSecType == 3))
	{
		const std::vector<double> cvPrices = GetCVPrices();
		ARM_StringVector cvNames  = ControlVariateColumnNames();
		/// Get the stored betas
		if (itsHasBeenPriced && GetFreezeBetasFlag() && cvPrices.size()>0)
		{			
#ifdef __GP_STRICT_VALIDATION
			if( itsBetasCoeff == ARM_VectorPtr(NULL) || itsBetasCoeff->size()!= cvPrices.size())
				 ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Betas Coeff are NULL or not of a good size");
#endif
			genPricer = new ARM_GenPricer( &*genSec,&*GetPricingModel(), cvNames, cvPrices, "Bermuda",*itsBetasCoeff);
			genPricer->Price();
		}
		else
		{
			genPricer = new ARM_GenPricer( &*genSec,&*GetPricingModel(), cvNames, cvPrices, "Bermuda");
			/// store betas
			genPricer->Price();
			if (GetFreezeBetasFlag() && cvPrices.size()>0)
			{
				itsBetasCoeff = ARM_VectorPtr(new std::vector<double>(genPricer->GetCVInfo()->GetBetas()));
			}
		}
		ARM_AutoCleaner<ARM_GenPricer> HoldGP(genPricer );
		nbProductsToPrice = cvNames.size() + 1;
		if (GetNbControlVariate()>0)
			price	= genPricer->GetPricerInfo()->GetContents( ARM_GenPricer::CVTag).GetData("Price").GetDouble();
		else
			price	= genPricer->GetPricerInfo()->GetContents( "Bermuda").GetData("Price").GetDouble();
		
		if (GetFixBoundaryFlag() && (GetNumMethodType() == MonteCarlo) )
		{
			ARM_GenSecManipulator genSecManipulator;
			ARM_DealDescriptionPtr dealDesc( new ARM_DealDescription(GetGenSecurity()->GetDealDescription()));
			genSecManipulator.ChangeAmericanIntoTrigger( *GetGenSecurity(), dealDesc );
			ARM_GenSecurityPtr newGenSec = ARM_GenSecurityPtr(new ARM_GenSecurity(dealDesc,"",GetGenSecurity()->GetCstManager(),!GetFixBoundaryFlag(),GetGenSecurity()->GetOtherPayoffsFlag()));
			SetGenSecurity( newGenSec );
			itsFixBoundaryFlag = false;
		}
		GetPricingData()["BermudaPrice"] = price;
		SetBermudaSwaptionPrice(price);
		if (GetNumMethodType() == MonteCarlo)
		{
			double stdDev;
			if (GetNbControlVariate()>0)
				stdDev	= genPricer->GetPricerInfo()->GetContents( ARM_GenPricer::CVTag).GetData("StdDev").GetDouble();
			else
				stdDev	= genPricer->GetPricerInfo()->GetContents( "Bermuda").GetData("StdDev").GetDouble();
			GetPricingData()["BermudaStdDev"] = stdDev;
			SetBermudaSwaptionStdDev(stdDev);
			itsBermudaSwaptionRA = genPricer->GetPricerInfo()->GetContents( "Bermuda").GetData("RunningAvg").GetVector();
			GetPricingData()["BermudaRunningAvg"] = itsBermudaSwaptionRA;
		}

		//Control Variates
		int nbCtrl = cvNames.size();
		double currCVPrice;	
		string currCVName;
		for (int i =0; i< nbCtrl;i++)
		{
			currCVName = cvNames[i];
			currCVPrice = genPricer->GetPricerInfo()->GetContents(currCVName).GetData("Price").GetDouble();
			GetPricingData()[ControlVariateColNamesOutputs[i]] = currCVPrice;
			SetControlVariatePrice(i,currCVPrice);
		}
	}
	itsHasBeenPriced = true;

    return price*GetPorS();
}

////////////////////////////////////////////////////
///	Class   : ARM_CRFCalculator
///	Routines: toString
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
string ARM_BermudaSwaptionCalculator::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream bermudaData;

    return bermudaData.str();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: View
///	Returns: 
///	Action : .
/////////////////////////////////////////////////////////////////
void ARM_BermudaSwaptionCalculator::View(char* id, FILE* ficOut) const
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

    /// Bermuda Swaption Calculator specific datas viewing
    fprintf(fOut,"\n\n =============================================== \n");
	fprintf(fOut,"\n\n =======> BERMUDA SWAPTION CALCULATOR <========= \n");
    fprintf(fOut,"\n\n =============================================== \n");

	//Product Description
	fprintf(fOut,"\n\n =======> Product Description <========= \n");
   
	CC_Ostringstream generalData;

	generalData << "Currency			= " << itsCcy.GetCcyName() << "\n";
	generalData << "StartDate			= " << itsStartDate.toString() << "\n";
    generalData << "EndDate				= " << itsEndDate.toString() << "\n";
    generalData << "Pay/Rec				= " << ARM_ParamView::GetMappingName(S_RECEIVE_PAY, itsPayRec) << "\n";
	generalData << "Stub Rule			= " << ARM_ParamView::GetMappingName(S_STUB_RULES, itsStubRule) << "\n";
	if (GetZCFlag() == true)
	{
		generalData << "Product Type:          Zero-Coupon Bermuda  \n";  
	}

	generalData << "\nCall Frequency	= " << ARM_ParamView::GetMappingName(S_FREQUENCY, itsCallFreq ) << "\n";
   	generalData << "Call Notice			= " << itsCallNotice << "\n";
	generalData << "Call Calendar		= " << itsCallCal << "\n";
	generalData << "FirstCallDate		= " << itsFirstCallDate.toString() << "\n";
	generalData << "LastCallDate		= " << itsLastCallDate.toString() << "\n";

    generalData << "\nFix Frequency		= " << ARM_ParamView::GetMappingName(S_FREQUENCY, itsFixFreq ) << "\n";
	generalData << "Fix Basis			= " << ARM_ParamView::GetMappingName( S_DAYCOUNT, itsFixBasis )  << "\n";
	generalData << "Fix Pay Calendar	= " << itsFixPayCal << "\n";
	generalData << "Fix Adj Rule		= " << ARM_ParamView::GetMappingName( S_FORWARD_RULES, itsFixAdjRule )  << "\n";
	generalData << "Fix Rule			= " << ARM_ParamView::GetMappingName( S_INTEREST_RULES, itsFixRule )  << "\n";
	generalData << "Fix Pay Gap			= " << itsFixPayGap << "\n";

    generalData << "\nVar Frequency		= " << ARM_ParamView::GetMappingName(S_FREQUENCY, itsVarFreq ) << "\n";
   	generalData << "Var Basis			= " << ARM_ParamView::GetMappingName( S_DAYCOUNT, itsVarBasis )  << "\n";
	generalData << "Var Reset Calendar	= " << itsVarResetCal << "\n";
	generalData << "Var Pay Calendar	= " << itsVarPayCal << "\n";
	generalData << "Var Index Term		= " << itsVarIndexTerm << "\n";
	generalData << "Var Adj Rule		= " << ARM_ParamView::GetMappingName( S_FORWARD_RULES, itsVarAdjRule )  << "\n";
	generalData << "Var Rule			= " << ARM_ParamView::GetMappingName( S_INTEREST_RULES, itsVarRule )  << "\n";
	generalData << "Var Reset Gap		= " << itsVarResetGap << "\n";
	generalData << "Var Pay Gap			= " << itsVarPayGap << "\n";
	
	generalData << "\n";
	int i, size;
	double currDate, currValue;
	//Strike profile
	generalData << "\n\n ===========Strike Profile============= \n";
	if( itsStrike.IsConstant() )
	{
		generalData << "ConstRefValue = " << (*(itsStrike.GetDiscreteValues()))[0] << " ; \n";
	}
	else
	{
		size = itsStrike.GetDiscreteValues()->size();
		for (i=0; i<size; i++)
		{
			currDate =(*(itsStrike.GetDiscreteDates()))[i];
			currValue =(*(itsStrike.GetDiscreteValues()))[i];
			generalData << "StartDate[" << i+1 << "]= " <<  ARM_Date(currDate).toString() << "  ;  ";
  			generalData << "Strike[ " << i+1 << "]= " <<  currValue << "\n";
		}
	}
	generalData << "\n";
	
	//Spread profile
	generalData << "\n\n ===========Spread Profile============= \n";
	if( itsVarSpread.IsConstant() )
	{
		generalData << "ConstRefValue = " << (*(itsVarSpread.GetDiscreteValues()))[0] << " ; \n";
	}
	else
	{
		size = itsVarSpread.GetDiscreteValues()->size();
		for (i=0; i<size; i++)
		{
			currDate =(*(itsVarSpread.GetDiscreteDates()))[i];
			currValue =(*(itsVarSpread.GetDiscreteValues()))[i];
			generalData << "StartDate[" << i+1 << "]= " <<  ARM_Date(currDate).toString() << "  ;  ";
  			generalData << "Spread[ " << i+1 << "]= " <<  currValue << "\n" ;
		}
	}
	generalData << "\n";
	
	//Notional profile
	generalData << "\n\n ===========Notional Profile============= \n";
	if( itsNominal.IsConstant() )
	{
		generalData << "ConstRefValue = " << (*(itsNominal.GetDiscreteValues()))[0] << " ; \n";
	}
	else
	{
		size = itsNominal.GetDiscreteValues()->size();
		for (i=0; i<size; i++)
		{
			currDate =(*(itsNominal.GetDiscreteDates()))[i];
			currValue =(*(itsNominal.GetDiscreteValues()))[i];
			generalData << "Date[" << i+1 << "]= " <<  ARM_Date(currDate).toString() << "  ;  ";
  			generalData << "Nominal[ " << i+1 << "]= " <<  currValue << "\n";
		}
	}
	generalData << "\n";

	//Fees profile
	generalData << "\n\n ===========Fees Profile============= \n";
	if( itsFees.IsConstant() )
	{
		generalData << "ConstRefValue = " << (*(itsFees.GetDiscreteValues()))[0] << " ; \n";
	}
	else
	{
		size = itsFees.GetDiscreteValues()->size();
		for (i=0; i<size; i++)
		{
			currDate =(*(itsFees.GetDiscreteDates()))[i];
			currValue =(*(itsFees.GetDiscreteValues()))[i];
			generalData << "NoticeDate[" << i+1 << "]= " <<  ARM_Date(currDate).toString() << "  ;  ";
  			generalData << "Fee[ " << i+1 << "]= " <<  currValue << "\n";
		}
	}
	generalData << "\n";

	fprintf(fOut,"%s",generalData.str().c_str());
    
	//Calibration parameters
	fprintf(fOut,"\n\n =======> Calibration Parameters <========= \n");
	CC_Ostringstream calibParameters;

	//Model Type
	double asOf = GetMktDataManager()->GetAsOfDate().GetJulian();
	calibParameters << "\n ------ Model Type	= " <<  GetModelType() << "\n";

	if (itsModelParams)
	{
		calibParameters << "\n ------ Init Volatility ------\n";
		ARM_GP_CurvePtr initVol(RefValueToCurve(*(*itsModelParams)[0], asOf));
		calibParameters << initVol->toString();

		calibParameters << "\n ------ Init MRS ------\n";
		ARM_GP_CurvePtr initMRS(RefValueToCurve(*(*itsModelParams)[1], asOf));
		calibParameters << initVol->toString();

		calibParameters << "\n ------ Default MRS	= " << itsDefaultMrs << " ------\n";

		calibParameters << "\n ------ Init Beta ------\n";
		ARM_GP_CurvePtr initBeta(RefValueToCurve(*(*itsModelParams)[2], asOf));
		calibParameters << initBeta->toString();

		calibParameters << "\n ------ Init Theta ------\n";
		ARM_GP_CurvePtr initTheta(RefValueToCurve(*(*itsModelParams)[3], asOf));
		calibParameters << initTheta->toString();

		calibParameters << "\n ------ Init Skew ------\n";
		ARM_GP_CurvePtr initSkew(RefValueToCurve(*(*itsModelParams)[4], asOf));
		calibParameters << initSkew->toString();

		calibParameters << "\n ------ Init MRS Low ------\n";
		ARM_GP_CurvePtr initMRSLow(RefValueToCurve(*(*itsModelParams)[5], asOf));
		calibParameters << initMRSLow->toString();

		calibParameters << "\n ------ Init MRS Up ------\n";
		ARM_GP_CurvePtr initMRSUp(RefValueToCurve(*(*itsModelParams)[6], asOf));
		calibParameters << initMRSUp->toString() << "\n\n";
	}

	//Calibrate MRS
	bool calibMRS = GetMrsCalibFlag();
	if (true == calibMRS)
		calibParameters << "Calibrate MRS ?		= " <<  "TRUE" << "\n";
	else if (false == calibMRS)
		calibParameters << "Calibrate MRS ?		= " <<  "FALSE" << "\n";
	//ATM Diagonal
	calibParameters << "ATM Diagonal		= " <<  GetATMDiagonalFlag() << "\n";
	//Method Type
	int numMethod = GetNumMethodType();
	if (MonteCarlo == numMethod)
		calibParameters << "Method Type			= " <<  "MONTECARLO" << "\n";
	else if (Tree == numMethod)
		calibParameters << "Method Type			= " <<  "TREE" << "\n";
	else if (PDE == numMethod)
		calibParameters << "Method Type			= " <<  "PDE" << "\n";
	//AMC Iter
	calibParameters << "AMC Iter			= " <<  GetAmcIter() << "\n";
	//MC Iter
	calibParameters << "MC Iter				= " <<  GetMcIter()	<< "\n";
	//Max Bucket Size
	calibParameters << "Max Bucket Size		= " <<  GetMaxBucketSize() << "\n";
	//Tree Steps
	calibParameters << "Tree Steps			= " <<  GetTreeSteps() << "\n";
	//MRS Portfolio Mode
	int mrsPortfolioMode = GetMrsPortfolioMode();
	if (Summit == mrsPortfolioMode)
		calibParameters << "MRS Portfolio Mode	= " <<  "SUMMIT" << "\n";
	else if (Manual == mrsPortfolioMode)
		calibParameters << "MRS Portfolio Mode	= " <<  "MANUAL" << "\n";
	else if (Faster == mrsPortfolioMode)
		calibParameters << "MRS Portfolio Mode	= " <<  "FASTER" << "\n";
	else if (NewMode == mrsPortfolioMode)
		calibParameters << "MRS Portfolio Mode	= " <<  "NEWMODE" << "\n";
	//MRS Calib Freq
	calibParameters << "MRS Calib Freq		= " <<  ARM_ParamView::GetMappingName(S_FREQUENCY, GetMrsCalibFreq()) << "\n";
	//MRS Calib Strikes
	calibParameters << "MRS Calib Strikes	= " <<  GetMrsCalibStrikes()	<< "\n";
	//Fix Boundary
	bool fixBoundary = GetFixBoundaryFlag();
	if (true == fixBoundary)
		calibParameters << "Fix Boundary ?		= " <<  "TRUE" << "\n";
	else if (false == fixBoundary) 
		calibParameters << "Fix Boundary ?		= " <<  "FALSE" << "\n";
	//Approx Margin Flag
	bool approxMargin = GetApproxMarginFlag();
	if (true == approxMargin)
		calibParameters << "Approx Margin ?		= " <<  "TRUE" << "\n";
	else if (false == approxMargin)
		calibParameters << "Approx Margin ?		= " <<  "FALSE" << "\n";
	//Freeze Betas
	bool freezeBetas = GetFreezeBetasFlag();
	if (true == freezeBetas)
		calibParameters << "Freeze Betas ?		= " <<  "TRUE" << "\n";
	else if (false == freezeBetas)
		calibParameters << "Freeze Betas ?		= " <<  "FALSE" << "\n";
	//Generator Type 1
	calibParameters << "Generator Type	1	= " <<  GetGenType1() << "\n";
	//Generator Type 2
	calibParameters << "Generator Type	2	= " <<  GetGenType2() << "\n";
	//Path Order
	calibParameters << "PathOrder			= " <<  GetPathOrder() << "\n";
	//Path Scheme
	calibParameters << "PathScheme			= " <<  GetPathScheme() << "\n";
	//First Nb Times
	calibParameters << "FirstNbTimes		= " <<  GetFirstNbTimes() << "\n";
	//MRS Calib Mode
	int mrsCalibMode = GetMrsCalibMode();
	if (Column == mrsCalibMode)
		calibParameters << "MRS Calib Mode		= " <<  "COLUMN" << "\n";
	else if (Correl == mrsCalibMode)
		calibParameters << "MRS Portfolio Mode	= " <<  "CORREL" << "\n";
	else if (Surdiag == mrsCalibMode)
		calibParameters << "MRS Portfolio Mode	= " <<  "SURDIAG" << "\n";
	else if (AntiDiag == mrsCalibMode)
		calibParameters << "MRS Portfolio Mode	= " <<  "ANTIDIAG" << "\n";
	//Calculate Proba Flag
	bool calcProba = itsCalculateProbaFlag;
	if (true == calcProba)
		calibParameters << "Calculate Proba ?	= " <<  "TRUE" << "\n";
	else if (false == calcProba)
		calibParameters << "Calculate Proba ?	= " <<  "FALSE" << "\n";
	//Vol Time Lag
	calibParameters << "Vol Time Lag		= " <<  GetVolatilityLag() << "\n";

	fprintf(fOut,"%s",calibParameters.str().c_str());

	//Calib Portfolio
	CC_Ostringstream bermudaData;
	
	fprintf(fOut,"\n\n ============>    CALIB PORTFOLIO   <=========== \n");
	if(GetCalibMethod()!= ARM_CalibMethodPtr(NULL))
		bermudaData << GetCalibMethod()->GetPortfolio()->toString()<< "\n\n";
	fprintf(fOut,"%s",bermudaData.str().c_str());
	fprintf(fOut,"\n\n ============> EOF CALIB PORTFOLIO  <=========== \n");

	if (itsCalculateProbaFlag && (itsGenSecType == 2) )
	{
		CC_Ostringstream probaData;
		probaData<<"\n\n ===========> Cumulative Probabilities (in %)<============ \n\n";
		double previousProba = 0.;
		for(size_t i =0; i<itsCumProba->size()-1;i++)
		{
			probaData<<"Exercise" << i+1 << "     ";
			probaData<<setprecision(0)<<setw(0)<<(itsCumProba->Elt(i) -  previousProba) * 100<<"\n";
			previousProba = itsCumProba->Elt(i);
		}
		probaData<<endl;
		probaData<<"Total Proba"<< "     ";
		probaData<<setprecision(0)<<setw(0)<<itsCumProba->Elt(itsCumProba->size()-1)*100<<"\n\n";
		probaData<<itsProbaMatrix.toString()<<"\n\n";
		fprintf(fOut,"%s",probaData.str().c_str());
		
	}

	//GenSecType
	CC_Ostringstream genSecType;
	genSecType << "Generic Security Type = " <<  itsGenSecType << "\n";
	fprintf(fOut,"%s",genSecType.str().c_str());


	//Gen Calculator
	fprintf(fOut,"\n\n ============>    GEN CALCULATOR    <=========== \n");
	bermudaData << "\n\nCommon Calculator Part\n" << ARM_GenCalculator::toString() << "\n\n";
	fprintf(fOut,"%s",bermudaData.str().c_str());
	fprintf(fOut,"\n\n ============>  EOF GEN CALCULATOR  <=========== \n");

	try
	{
	    ARM_MarketIRModel* marketModel = dynamic_cast<ARM_MarketIRModel*>(GetMktDataManager()->GetData(GetKeys()[NormalModel]));

		ARM_BSModel* bsModel = dynamic_cast<ARM_BSModel*>(marketModel->GetMktDataManager()->GetData(marketModel->GetBsModelKey()));

		ARM_VolCurve* correlCube  = bsModel->GetCorrelationMatrix();

		fprintf(fOut,"\n\n ============>    Used Correlation    <=========== \n");
		if (correlCube)
		   correlCube->View(id, fOut); 
	}
	catch(...)
	{
		fprintf(fOut,"\n\n ============>    Used Correlation Not Found   <=========== \n");
	}

	fprintf(fOut,"\n\n ======END OF VIEW=============================== \n");

	if ( ficOut == NULL )
       fclose(fOut);
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
