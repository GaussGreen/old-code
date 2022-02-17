/*!
 /** Copyright (c) CDC IXIS CM July 2003 Paris
 *	\file CallableQuantoCalculator.cpp
 *  \brief file for the Snow Range Calculator
 *	\author  J.M
 *	\version 1.0
 *	\date Apr 2007
 */


#include "gpcalculators/callablequantocalculator.h"

/// gpbase
#include "gpbase/env.h"
#include "gpbase/typedef.h"
#include "gpbase/argconvdefault.h"
#include "gpbase/datestripcombiner.h"
#include "gpbase/autocleaner.h"
#include "gpbase/singleton.h"
#include "gpbase/datestrip.h"

/// gpinfra
#include "gpinfra/argconvdefault.h"
#include "gpinfra/gensecurity.h"
#include "gpinfra/dealdescription.h"
#include "gpinfra/gensecmanipulator.h"
#include "gpinfra/mktdatamanagerrep.h"
#include "gpinfra/pricingadviser.h"
#include "gpinfra/genpricer.h"
#include "gpinfra/numeraire.h"
#include "gpinfra/numerairefactory.h"
#include "gpinfra/pricerinfo.h"
#include "gpinfra/gramfunctorargdict.h"
#include "gpinfra/cstmanager.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/correlmatparam.h"
#include "gpinfra/modelnamemap.h"

/// gpnumlib
#include "gpnumlib/argconvdefault.h"
#include "gpnumlib/randomgenfactory.h"

/// gpnummethods
#include "gpnummethods/treemethod.h"
#include "gpnummethods/treebase.h"
#include "gpnummethods/treefactory.h"
#include "gpnummethods/scheduler.h"
#include "gpnummethods/sampler.h"
#include "gpnummethods/truncator.h"
#include "gpnummethods/reconnector.h"
#include "gpnummethods/smoother.h"
#include "gpnummethods/mcmethod.h"
#include "gpnummethods/meanrevertingsampler.h"
#include "gpnummethods/amcmethod.h"
#include "gpnummethods/amc_andersen.h"
#include "gpnummethods/amc_ls.h"

/// gpmodels
#include "gpmodels/modelparamshw1f.h"
#include "gpmodels/MultiAssetsFactory.h"
#include "gpmodels/multiassets.h"
#include "gpmodels/hw1f.h"
#include "gpmodels/HWHWQtoModel.h"
#include "gpmodels/eqfx_modelfactory.h"
#include "gpmodels/ModelParams_EqFxBase.h"
#include "gpmodels/forwardmarginbasis.h"
#include "gpmodels/local_normal_model.h"
#include "gpmodels/local_normal_modelparams.h"




CC_BEGIN_NAMESPACE( ARM )

/// Default MDM key names
const string YC_KEY_NAME		= "YC_";
const string BASIS_KEY_NAME		= "YC_BASIS_";
const string OSWMODEL_KEY_NAME	= "OSWMOD_";
const string CFMODEL_KEY_NAME	= "CFMOD_";
const string MRS_KEY_NAME		= "MRS_";
const string CORREL_KEY_NAME	= "CORREL";
const string FOREX_KEY_NAME		= "FOREX";
const string FXVOL_KEY_NAME		= "FXVOL";

const string ARM_CallableQuantoCalculator::CallableQuantoColNamesTable [] =
{  
	"ResetDate",
	"StartDate",
	"EndDate",
	"Fees",
	"Bermuda",
	"Option"
};

const int ARM_CallableQuantoCalculator::CallableQuantoProductToPriceColumns [] =
{  
	Option
};

//////////////////////////////////////////////////////////////////////////////
/// Class  : ARM_CallableQuantoCalculator
/// Routines: 
/// Returns :
/// Action  : constructor
//////////////////////////////////////////////////////////////////////////////
ARM_CallableQuantoCalculator::ARM_CallableQuantoCalculator(		const ARM_Currency& ccyDom,
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
																const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice)
:	ARM_GenCalculator(),
	itsStartDate(startDate),
	itsEndDate(endDate),
	itsPayRec(payRec),
	itsNotional(notional),
	itsCallFreq(callFreq),
	itsCallNotice(callNotice),
	itsCallCal(callCal),
	itsCallFees(callFees),
	itsFundFreq(fundFreq),
	itsFundSpread(fundSpread),
	itsFundDayCount(fundDayCount),
	itsCpnPayFreq(cpnPayFreq),
	itsCpnDayCount(cpnDayCount),
	itsCpnResetCal(cpnResetCal),
	itsCpnPayCal(cpnPayCal),
	itsProductsToPrice(productsToPrice)
{
	SetDomesticCcy(ARM_Currency(ccyDom));
	SetForeignCcy(ARM_Currency(ccyFor));

	string domCcyName(GetDomesticCcy().GetCcyName());
    string forCcyName(GetForeignCcy().GetCcyName());

	// create internal keys
	ARM_StringVector keys(NbKeys);
	keys[YcDomKey]			= YC_KEY_NAME + domCcyName;
	keys[YcForKey]			= YC_KEY_NAME + forCcyName;
	keys[ForexKey]			= FOREX_KEY_NAME;
	keys[YcBasisDomKey]		= BASIS_KEY_NAME + domCcyName;
	keys[YcBasisForKey]		= BASIS_KEY_NAME + forCcyName;
	keys[OswDomModelKey]	= OSWMODEL_KEY_NAME + domCcyName;
	keys[OswForModelKey]	= OSWMODEL_KEY_NAME + forCcyName;
	keys[CfDomModelKey]		= CFMODEL_KEY_NAME + domCcyName;
	keys[CfForModelKey]		= CFMODEL_KEY_NAME + forCcyName;
	keys[CorrelMatrixKey]	= CORREL_KEY_NAME;
	keys[MrsDomKey]			= MRS_KEY_NAME + domCcyName;
	keys[MrsForKey]			= MRS_KEY_NAME + forCcyName;
	keys[FXVolKey]			= FXVOL_KEY_NAME;

	SetKeys(keys);


	itsSchedulerDatas.resize(3);
	itsSchedulerDatas[0]=20;	// Nb min steps before 1st notice
	itsSchedulerDatas[1]=2;		// Nb max steps before 1st notice
    itsSchedulerDatas[2]=1.e-3; // Nb steps per year before 1st notice

	itsTruncatorDatas.resize(1);
	itsTruncatorDatas[0]=5;     // Nb max std dev

	itsProductPrices.resize(itsProductsToPrice.size());
}


//////////////////////////////////////////////////////////////////////////////
/// Class  : ARM_CallableQuantoCalculator
/// Routines: 
/// Returns :
/// Action  : copy constructor
//////////////////////////////////////////////////////////////////////////////
ARM_CallableQuantoCalculator::ARM_CallableQuantoCalculator(const ARM_CallableQuantoCalculator& rhs)
: ARM_GenCalculator(rhs)
{
	CopyNoCleanUp(rhs);
}
	
//////////////////////////////////////////////////////////////////////////////
/// Class  : ARM_CallableQuantoCalculator
/// Routine: CopyNoCleanUp
/// Returns: void
/// Action : Arguments copy
//////////////////////////////////////////////////////////////////////////////
void ARM_CallableQuantoCalculator::CopyNoCleanUp(const ARM_CallableQuantoCalculator& rhs)
{
	itsStartDate = rhs.itsStartDate;
	itsEndDate = rhs.itsEndDate;
	itsPayRec = rhs.itsPayRec;
	itsNotional = rhs.itsNotional;
	itsCallFreq = rhs.itsCallFreq;
	itsCallNotice = rhs.itsCallNotice;
	itsCallCal = rhs.itsCallCal;
	itsCallFees = rhs.itsCallFees;
	itsFundFreq = rhs.itsFundFreq;
	itsFundSpread = rhs.itsFundSpread;
	itsFundDayCount = rhs.itsFundDayCount;
	itsCpnPayFreq = rhs.itsCpnPayFreq;
	itsCpnDayCount = rhs.itsCpnDayCount;
	itsCpnResetCal = rhs.itsCpnResetCal;
	itsCpnPayCal = rhs.itsCpnPayCal;
	itsProductsToPrice = rhs.itsProductsToPrice;
	itsProductPrices = rhs.itsProductPrices;
}

//////////////////////////////////////////////////////////////////////////////
/// Class  : ARM_CallableQuantoCalculator
/// Routines: operator=
/// Returns :
/// Action  : Call copy constructor
//////////////////////////////////////////////////////////////////////////////
ARM_CallableQuantoCalculator& ARM_CallableQuantoCalculator::operator=(const ARM_CallableQuantoCalculator& rhs)
{
	if( this != & rhs )
	{
		ARM_GenCalculator::operator=(rhs);
        CleanUp();
        CopyNoCleanUp(rhs);
	}
	return *this;
}

//////////////////////////////////////////////////////////////////////////////
/// Class  : ARM_CallableQuantoCalculator
/// Routine: CleanUp
/// Returns: void
/// Action : delete pointors
//////////////////////////////////////////////////////////////////////////////
void ARM_CallableQuantoCalculator::CleanUp()
{
}

//////////////////////////////////////////////////////////////////////////////
/// Class  : ARM_CallableQuantoCalculator
/// Routine: Destructor
/// Returns: void
/// Action : destroys the object
//////////////////////////////////////////////////////////////////////////////
ARM_CallableQuantoCalculator::~ARM_CallableQuantoCalculator()
{
	CleanUp();
}

//////////////////////////////////////////////////////////////////////////////
/// Class  : ARM_CallableQuantoCalculator
/// Routine: ColumnNames
/// Returns: ARM_RowInfo
/// Action : returns column names
//////////////////////////////////////////////////////////////////////////////
ARM_RowInfo	ARM_CallableQuantoCalculator::ColumnNames() const
{
	// Number of Columns
	size_t colNamesSizeCallableQuanto;
	size_t colNamesSize;

	colNamesSizeCallableQuanto = sizeof(CallableQuantoColNamesTable)/sizeof(CallableQuantoColNamesTable[0]);
    colNamesSize = colNamesSizeCallableQuanto + ExtraDescSize();

	vector< string > colNamesVec(colNamesSize);
    vector< ARM_GP_VALUE_TYPE > colTypeVec(colNamesSize, ARM_STRING); 

	for (size_t i = 0; i < colNamesSizeCallableQuanto; ++i)
		colNamesVec[i] = CallableQuantoColNamesTable[i];

	UpdateExtraColNames(colNamesVec);
    ARM_RowInfo rowInfo(colNamesVec,colTypeVec);

    return rowInfo;
}

//////////////////////////////////////////////////////////////////////////////
/// Class  : ARM_CallableQuantoCalculator
/// Routine: MiddleRows
/// Returns: ARM_RowInfo
/// Action : create a row of a deal description
//////////////////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_CallableQuantoCalculator::MiddleRows(size_t eventIdx, 
                                          const ARM_DateStripCombiner& datesStructure) const
{
	try
	{
		double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

		size_t EXER_SCHED = ScheduleCount::ExerSched;
		size_t eventSize = datesStructure.GetDateStrip(EXER_SCHED)->GetResetDates()->size();
		size_t descSize = sizeof(CallableQuantoColNamesTable)/sizeof(CallableQuantoColNamesTable[0]);

		descSize += ExtraDescSize();

		vector< string > rowDescVec(descSize);
		vector< ARM_GP_VALUE_TYPE > rowTypeVec(descSize, ARM_MISSING_TYPE); 

		/// Set default 0 value for each column to be able to sum it
		InitPriceableColumns(rowDescVec,rowTypeVec);

		std::vector<double>& exerciseDates = datesStructure.GetDateStrip(EXER_SCHED)->GetResetDates();
		size_t nbPastNoCall=0;
		while(nbPastNoCall < eventSize && (*exerciseDates)[nbPastNoCall] < asOfDate)
			++nbPastNoCall;

		bool isFirstEvent = (eventIdx==nbPastNoCall);
		bool isLastEvent  = (eventIdx==(eventSize-1));
		
		//RESET DATE
		double	noticeDate = (*(datesStructure.GetDateStrip(EXER_SCHED)->GetResetDates()))[eventIdx];
		CC_Ostringstream noticeDateDesc;
		noticeDateDesc << CC_NS(std,fixed) << noticeDate;
		rowDescVec[ResetDate] = noticeDateDesc.str();
		rowTypeVec[ResetDate] = ARM_DATE_TYPE;

		//START DATE
		CC_Ostringstream startDateDesc;
		double startDate = (*(datesStructure.GetDateStrip(EXER_SCHED)->GetFlowStartDates()))[eventIdx] ;
		startDateDesc << CC_NS(std,fixed) << startDate;
		rowDescVec[StartDate] = startDateDesc.str();
		rowTypeVec[StartDate] = ARM_DATE_TYPE;
		
		//END DATE
		CC_Ostringstream endDateDesc;
		double endDate = (*(datesStructure.GetDateStrip(EXER_SCHED)->GetFlowEndDates()))[eventIdx] ;
		endDateDesc << CC_NS(std,fixed) << endDate;
		rowDescVec[EndDate] = endDateDesc.str();
		rowTypeVec[EndDate] = ARM_DATE_TYPE;

		//FEES
		double fees = const_cast< ARM_ReferenceValue& >(itsCallFees).Interpolate(noticeDate);
		CC_Ostringstream feesDesc;
		feesDesc << CC_NS(std,fixed) << fees;
		rowDescVec[Fees] = feesDesc.str();
		rowTypeVec[Fees] = ARM_DOUBLE;

		//BERMUDA ALGO
		CC_Ostringstream bermudaDesc;
		if(isLastEvent)
			bermudaDesc << "Swap[i] + MAX(-Swap[i]-Fees[i],0)";
		else
			bermudaDesc << "Exercise(Swap[i],-Swap[i]-Fees[i],Bermuda[i+1])";
		rowDescVec[Bermuda] = bermudaDesc.str();
		rowTypeVec[Bermuda] = ARM_STRING;

		//OPTION
		CC_Ostringstream optionDesc;
		if(isFirstEvent)
			optionDesc << "Bermuda[i]";
		else
			optionDesc << "0";
		rowDescVec[Option] = optionDesc.str();
		rowTypeVec[Option] = ARM_STRING;


		UpdateExtraDesc(eventIdx,datesStructure,rowDescVec,rowTypeVec);

		return ARM_RowInfo(rowDescVec,rowTypeVec);
	}

	catch (Exception& x)
	{
	    x.DebugPrint();

		throw x;
	}
	catch (...)
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in CallableQuantoCalculator::MiddleRows" );
	}
}

//////////////////////////////////////////////////////////////////////////////
/// Class  : ARM_CallableQuantoCalculator
/// Routine: DatesStructure
/// Returns: ARM_DateStripVector
/// Action : create the list of all event dates of the calculator. 
//////////////////////////////////////////////////////////////////////////////
ARM_DateStripCombiner ARM_CallableQuantoCalculator::DatesStructure() const
{
	try
	{
		ARM_DateStripVector SchedVect(1,NULL);

		double startDate = itsStartDate.GetJulian();
		double endDate = itsEndDate.GetJulian();

		size_t EXER_SCHED = ScheduleCount::ExerSched;

		SchedVect[EXER_SCHED] = new ARM_DateStrip(	itsStartDate.GetJulian(), 
													itsEndDate.GetJulian(),
													itsCallFreq, 
													itsCpnDayCount,
													itsCallCal.c_str(), 
													GETDEFAULTVALUE, 
													K_UNADJUSTED, 
													K_SHORTSTART, 
													itsCallNotice,
													itsCallFreq, 
													GETDEFAULTVALUE, 
													itsCallCal.c_str(), 
													K_ADVANCE, 
													K_ARREARS, 
													true);

		return SchedVect;
	}
	catch (Exception& x)
	{
	    x.DebugPrint();

		throw x;
	}
	catch (...)
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in SnowRangeCalculator::DatesStructure" );
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableQuantoCalculator
///	Routine: CustomDatesStructure
///	Returns: ARM_DateStripCombiner
///	Action : create the list of all event dates
/////////////////////////////////////////////////////////////////
ARM_DateStripCombiner ARM_CallableQuantoCalculator::CustomDatesStructure(const ARM_DateStripVector& dateStrips) const
{
	return ARM_DateStripCombiner();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableQuantoCalculator
///	Routine: InitPriceableColumns
///	Returns: void
///	Action : init
/////////////////////////////////////////////////////////////////
void ARM_CallableQuantoCalculator::InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const
{
    string zeroValue("0");

	rowDescVec[0] = zeroValue;
	rowTypeVec[0] = ARM_DOUBLE;

}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableQuantoCalculator
///	Routine: ComputePricingData
///	Returns: a ARM_MultiTypeDict
///	Action : get pricing data of the calculaltor
/////////////////////////////////////////////////////////////////
void ARM_CallableQuantoCalculator::ComputePricingData() const
{
	if (!itsHasBeenPriced)
		const_cast<ARM_CallableQuantoCalculator*>(this)->PriceAndTimeIt();

	int colSize = MIN(NbProductsToPrice + ExtraNbProductsToPrice(), GetProductsToPrice().size());

	for (size_t i(0); i < colSize; ++i)
		GetPricingData()[ ColNames(i) ] = itsProductPrices[i];

}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableQuantoCalculator
///	Routine: Price
///	Returns: a double
///	Action : price the ARM_CallableQuantoCalculator.
/////////////////////////////////////////////////////////////////
double ARM_CallableQuantoCalculator::Price()
{
	try
	{
		CalibrateAndTimeIt();

		/// Price the implicit product according to internal flag
		ARM_GenSecurityPtr genSec = GetGenSecurity();

		size_t nbPricedColumns = genSec->GetDealDescription().GetPricedColumnNames().size();

		double price = 0.0;
		ARM_GenPricer* genPricer = NULL;
		genPricer = new ARM_GenPricer( &*genSec,&*GetPricingModel() );
		ARM_AutoCleaner<ARM_GenPricer> HoldGP(genPricer );
		genPricer->Price();

		// effective nb of columns to price :
		int colSize = MIN(NbProductsToPrice + ExtraNbProductsToPrice(), GetProductsToPrice().size());

		for (size_t i(0); i < colSize; ++i)
		{
			if (GetProductsToPrice()[i])
			{
				price	= genPricer->GetPricerInfo()->GetContents( ColNames(i) ).GetData("Price").GetDouble();
				itsProductPrices[i] = price;
			}
		}
		itsHasBeenPriced = true;
		return itsProductPrices[OptionPrice];

	}
	catch (Exception& x)
	{
	    x.DebugPrint();

		throw x;
	}
	catch (...)
	{
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in CraCalculator::Price" );
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableQuantoCalculator
///	Routine: UpdateModel
///	Returns: void
///	Action : update the model datas w.r.t. MarketDataManager
/////////////////////////////////////////////////////////////////
void ARM_CallableQuantoCalculator::UpdateModel()
{
	CreateAndSetModel();

	itsHasBeenPriced=false;
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableQuantoCalculator
///	Routine: CheckMktData & CheckData
///	Returns: void
///	Action : check if data are consistent
/////////////////////////////////////////////////////////////////
void ARM_CallableQuantoCalculator::CheckMktData()
{
//*************************************************************************************************************/
}


void ARM_CallableQuantoCalculator::CheckData()
{
	CheckMktData();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableQuantoCalculator
///	Routine: CreateAndSetModel
///	Returns: void
///	Action : create the HWHWQto model
/////////////////////////////////////////////////////////////////
void ARM_CallableQuantoCalculator::CreateAndSetModel()
{
	// MODELS
    int nbModels = ARM_HWHWQtoModel::NbModels;
	vector< ARM_PricingModelPtr > models(nbModels);
    ARM_StringVector names(nbModels);
	ARM_StringVectorVector depends(nbModels);
    names[ARM_HWHWQtoModel::DomModel]        = GetKeys()[YcDomKey];
    names[ARM_HWHWQtoModel::ForModel]        = GetKeys()[YcForKey];
    names[ARM_HWHWQtoModel::DomBasisModel]   = GetKeys()[YcBasisDomKey];
    names[ARM_HWHWQtoModel::ForBasisModel]   = GetKeys()[YcBasisForKey];
	names[ARM_HWHWQtoModel::FxModel]         = GetKeys()[ForexKey];
	names[ARM_HWHWQtoModel::ForLocalModel]	 = LOCAL_MODEL_NAME;
    
	depends[ARM_HWHWQtoModel::DomBasisModel] = ARM_StringVector(1,names[ARM_HWHWQtoModel::DomModel]);
	depends[ARM_HWHWQtoModel::ForBasisModel] = ARM_StringVector(1,names[ARM_HWHWQtoModel::ForModel]);
	depends[ARM_HWHWQtoModel::FxModel]       = ARM_StringVector(2);
	depends[ARM_HWHWQtoModel::FxModel][ARM_HWHWQtoModel::DomModel]	   = names[ARM_HWHWQtoModel::DomBasisModel];
	depends[ARM_HWHWQtoModel::FxModel][ARM_HWHWQtoModel::ForModel]	   = names[ARM_HWHWQtoModel::ForBasisModel];
	depends[ARM_HWHWQtoModel::ForLocalModel] = ARM_StringVector(1,names[ARM_HWHWQtoModel::ForModel]);

    ARM_ModelParamVector modelParams(2);
    ARM_CurveModelParam volParam( ARM_ModelParamType::Volatility,0.01,"Volatility");
    modelParams[0]  = &volParam;
	
	modelParams[1] = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(GetKeys()[MrsDomKey]));
	ARM_ZeroCurve* ycDomCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcDomKey]));
	models[ARM_HWHWQtoModel::DomModel] = ARM_PricingModelPtr( new ARM_HullWhite1F( CreateClonedPtr(ycDomCurve),new ARM_ModelParamsHW1FStd(modelParams)) );
	
    modelParams[1] = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(GetKeys()[MrsForKey]));
	ARM_ZeroCurve* ycForCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcForKey]));
    models[ARM_HWHWQtoModel::ForModel] = ARM_PricingModelPtr( new ARM_HullWhite1F( CreateClonedPtr(ycForCurve),new ARM_ModelParamsHW1FStd(modelParams)) );

	ARM_ZeroCurve* ycBasisDomCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcBasisDomKey]));
	models[ARM_HWHWQtoModel::DomBasisModel] = ARM_PricingModelPtr( new ARM_ForwardMarginBasis(CreateClonedPtr(ycBasisDomCurve)) );
    
	ARM_ZeroCurve* ycBasisForCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcBasisForKey]));
    models[ARM_HWHWQtoModel::ForBasisModel] = ARM_PricingModelPtr( new ARM_ForwardMarginBasis(CreateClonedPtr(ycBasisForCurve)) );

	ARM_CurveModelParam mrsParam( ARM_ModelParamType::MeanReversion,0.0,"mrs");
	modelParams[0] = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(GetKeys()[FXVolKey]));
	modelParams[1]  = &mrsParam;
	ARM_Forex* forex = dynamic_cast<ARM_Forex*>(GetMktDataManager()->GetData(GetKeys()[ForexKey]));
    ARM_GP_Matrix* correlMatrix = dynamic_cast<ARM_GP_Matrix*>(GetMktDataManager()->GetData(GetKeys()[CorrelMatrixKey]));
	models[ARM_HWHWQtoModel::FxModel] = ARM_PricingModelPtr( 
								ARM_EqFx_ModelFactory.Instance()->CreateModel(	
										CreateClonedPtr(ycBasisDomCurve), 
										modelParams,
										forex->GetMarketPrice(), 
										CreateClonedPtr(ycBasisForCurve),
										*correlMatrix,ARM_EqFx_ModelFactoryImp::BS_Model)	);

	ARM_Local_Normal_ModelParams noLocalParams(ARM_ModelParamVector(0));
	models[ARM_HWHWQtoModel::ForLocalModel] = ARM_PricingModelPtr( new ARM_Local_Normal_Model(CreateClonedPtr(ycForCurve),noLocalParams) );
    
	ARM_ModelNameMap modelMap( names, models, depends );
	ARM_PricingModelPtr hwhwqtoModel = ARM_PricingModelPtr( new ARM_HWHWQtoModel(modelMap,*correlMatrix,false) );

	// TREE

	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	double lastEventTime = GetEndDate().GetJulian(); 

	/*size_t nbStepBefore1st = 4;
	int nbSteps = static_cast<int>(floor(20*lastEventTime/K_YEAR_LEN));

	int smootherType=ARM_SmootherBase::Linear;
	int schedulerType=ARM_SchedulerBase::ConstantVarianceMeanReverting;
	std::vector<double> schedulerDatas(3);
	schedulerDatas[0] = 20;
	schedulerDatas[1] = nbStepBefore1st;
	schedulerDatas[2] = 1.0e-3;
	int samplerType=ARM_SamplerBase::MeanReverting;
	std::vector<double> samplerDatas(1,1.0e-3);
	int truncatorType=ARM_TruncatorBase::StandardDeviation;
	std::vector<double> truncatorDatas(1,5.0);
	int reconnectorType=ARM_ReconnectorBase::Mean;
	bool probaFlag = false; // to allow spot probabilities price computation (... and PV calibration under terminal ZC numeraire)
	ARM_TreeBase* tree = ARM_TreeFactory.Instance()->CreateTreeND( 2, schedulerType,schedulerDatas,
		samplerType,samplerDatas,truncatorType,truncatorDatas,probaFlag,reconnectorType,smootherType);

	ARM_NumMethodPtr numMethod = ARM_NumMethodPtr(tree);
	hwhwqtoModel->SetNumMethod(numMethod);*/
	
	ARM_RandomGeneratorPtr  randGen = ARM_RandomGeneratorPtr(
		ARM_RandGenFactory.Instance()->CreateSimpleRandGen(	"NR_Ran2",
															"Sobol",
															"BoxMuller",
															"InvNormCum",
															"PathOrder",
															0,
															0,
															true) );
	ARM_TimeStepPerYearScheduler scheduler(0);
	ARM_MeanRevertingSamplerND sampler(&scheduler);

	ARM_ExerciseBoundaryCalc* exerBoundCal = new ARM_AMCAndersen(10000, true);
	ARM_NumMethodPtr numMethod = ARM_NumMethodPtr(new ARM_AMCMethod(10000,randGen,&sampler,exerBoundCal));
	delete exerBoundCal;

	hwhwqtoModel->SetNumMethod(numMethod);

    // NUMERAIRE
    ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::Cash ) );
    hwhwqtoModel->SetNumeraire(numeraire);

	// SET
	SetPricingModel(hwhwqtoModel);
}


CC_END_NAMESPACE()
