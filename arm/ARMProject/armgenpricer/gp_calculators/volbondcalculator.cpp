/*!
 *
 * Copyright (c) NATIXIS July 2007 Paris
 *
 *	\file gencsocalculator.cpp
 *
 *  \brief file for the virtual calculator for Vol Bond Calculator
 *
 *	\author  Frédéric Degraeve
 *	\version 1.1
 *	\date July 2007
 */

#include "gpbase/removeidentifiedwarning.h"
#include <list>
#include <vector>
#include <sstream>
#include <iostream>
#include <iterator>

#include "gpcalculators/volbondcalculator.h"

#include "gpbase/datestrip.h"
#include "gpbase/countedptr.h"
#include "gpbase/autocleaner.h"
#include "gpbase/gpvector.h"
#include "gpbase/singleton.h"
#include "gpbase/datestripcombiner.h"

#include "gpinfra/argconvdefault.h"
#include "gpinfra/gramnode.h"
#include "gpinfra/genpricer.h"
#include "gpinfra/gensecurity.h"
#include "gpinfra/dealdescription.h"
#include "gpinfra/modelparam.h"
#include "gpinfra/numerairefactory.h"
#include "gpinfra/cstmanager.h"
#include "gpinfra/pricingadviser.h"
#include "gpinfra/pricingstates.h"

#include "gpmodels/hwsv1f.h"
#include "gpmodels/ModelParamsHWSV.h"

#include "gpnummethods/mcmethod.h"
#include "gpnummethods/normalcentredsampler.h"
#include "gpnummethods/meanrevertingsampler.h"
#include "gpnummethods/scheduler.h"
#include "gpnummethods/argconvdefault.h"
#include "gpnummethods/pdemethod.h"
#include "gpnummethods/pathscheme.h"
#include "gpnummethods/pathschemefactory.h"
#include "gpnummethods/normalcentredsampler.h"
#include "gpnummethods/scheduler.h"
#include "gpnummethods/pathscheme.h"
#include "gpnummethods/impsampler.h"
#include "gpnummethods/schedulerfactory.h"
#include "gpnummethods/samplerfactory.h"
#include "gpnummethods/typedef.h"
#include "gpnummethods/impsamplerfactory.h"

#include "gpnumlib/argconvdefault.h"
#include "gpnumlib/randomgenfactory.h"
#include "gpnumlib/compositegen.h"
#include "gpnumlib/antitheticgen.h"
#include "gpnumlib/odefunctions.h"
#include "gpnumlib/random.h"

#include "gpcalib/calibmethod.h"
#include "gpcalib/argconvdefault.h"
#include "gpcalib/modelfitterdes.h"

// for the calibration

//#include <inst/swap.h>
//#include <inst/swaption.h>
//#include <inst/portfolio.h>

CC_BEGIN_NAMESPACE( ARM )


const string ZCCURVE        = "ZCCURVE";
const string VOLCURVE       = "VOLCURVE";
const string MRS			= "MRS";
const string VOLOFVOL       = "VOLOFVOL";
const string CORRELATION    = "CORRELATION";
const string VOLMRS			= "VOLMRS";

// for calibration of vol mean reversion
const string BSMODEL		= "BSMODEL";
const string VOVPARAM		= "VOVPARAM";
const string CORRELPARAM	= "CORRELPARAM";
const string VOLCALIB		= "VOLCALIB";


const int NB_VOLBOND_SCHED	= 2;

const int VOLBOND_SCHED_ADV	= 0;
const int VOLBOND_SCHED_ARR	= 1;


class SwaptionPortfolios
{
public:
	SwaptionPortfolios(
		const ARM_MarketData_ManagerRep&	mktDataManager, 
		const ARM_StringVector&				mktDataManagerKeys,
		const std::vector<double>&				startdate, 
		const std::vector<double>&				resetdate,
		const std::vector<double>&				list_strike_modificators,
		const int							SwapLiborType = K_EURIBOR6M,
		const int							SwaptionLiborType = K_EURIBOR1Y,
		const int							TenorInMonths = 120 // 120Months == 10Y
		)
		: itsStartDate(startdate), 
		itsResetDate(resetdate),
		itsMktDataManager(mktDataManager),
		itsMktDataManagerKeys(mktDataManagerKeys),
		itsListStrikeModificators(list_strike_modificators),
		itsSwapLiborType(SwapLiborType),
		itsSwaptionLiborType(SwaptionLiborType),
		itsTenor(TenorInMonths)
	{		
		itsDone = false;
	}	

	~SwaptionPortfolios()
	{
	}

protected:
						
	void BuildSwapList()
	{		
		ARM_ZeroCurve* zc = dynamic_cast<ARM_ZeroCurve*>( itsMktDataManager.GetData(itsMktDataManagerKeys[ARM_VolBondCalculator::mdmKeysAlias::ZcCurveKey]) );
		if (zc == 0) ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : ZCCURVE not of good type ");

		ARM_BSModel* BSModel = dynamic_cast<ARM_BSModel*>( itsMktDataManager.GetData(itsMktDataManagerKeys[ARM_VolBondCalculator::mdmKeysCalibrationAlias::BSModelKey]) );
		if (BSModel == 0) ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : BSMODEL not of good type ");

		double	Spread = 0.;			
		double	FixedRate = 5.;
		int		ReceiveOrPay = K_RCV;

		ARM_Date startDate(itsMktDataManager.GetAsOfDate().GetJulian());
		ARM_Date endDate(itsMktDataManager.GetAsOfDate().GetJulian() + 20);
		ARM_Date AsOfDate = (ARM_Date) itsMktDataManager.GetAsOfDate().GetJulian();

		ARM_Swap newSwap(startDate, endDate,
						   (ARM_INDEX_TYPE) itsSwapLiborType,
							Spread,
							FixedRate,
							ReceiveOrPay);
		newSwap.SetModel(BSModel);

		for (int i = 0; i < itsStartDate.size(); i++)
		{
			startDate = static_cast<ARM_Date>(itsStartDate[i]);
			endDate = startDate.AddMonths(itsTenor);

			newSwap.SetStartDate(startDate);
			newSwap.SetEndDate(endDate);

			double priceswap = newSwap.PriceToRate(AsOfDate, 0.);

			itsSwapPrice.push_back(priceswap);
		}		
	}
	
	void BuildSwaptionList()
	{
		ARM_BSModel* BSModel	= dynamic_cast< ARM_BSModel* >		( itsMktDataManager.GetData(itsMktDataManagerKeys[ARM_VolBondCalculator::mdmKeysCalibrationAlias::BSModelKey]) );
		if (BSModel == 0) ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : BSMODEL not of good type ");
		
		double	Spread = 0.;			
		double	FixedRate = 5.;
		int		ReceiveOrPay = K_PAY;
		long	ExerciseType = K_EUROPEAN;
		
		double			strike;
		ARM_Date		startDate;
		ARM_Date		endDate;
		ARM_Date		resetDate;
		ARM_Swaption*	newSwaption;
		int				i;

		for (i = 0; i < itsStartDate.size(); i++)
		{
			strike = itsSwapPrice[i];
			startDate = static_cast<ARM_Date>(itsStartDate[i]);
			endDate = startDate;
			endDate = endDate.AddMonths(itsTenor);
			resetDate = static_cast<ARM_Date>(itsResetDate[i]);

			newSwaption = new ARM_Swaption(startDate,
											endDate,
											ReceiveOrPay,
											ExerciseType,
											strike,
											resetDate,
											(ARM_INDEX_TYPE) itsSwaptionLiborType);

			newSwaption->SetModel(BSModel);
			newSwaption->ComputePrice();
			itsBootstrapSwaptionList.push_back(newSwaption);
		}		

		// we want the last one
		strike = itsSwapPrice[itsStartDate.size()-1];
		startDate = static_cast<ARM_Date>(itsStartDate[itsStartDate.size()-1]);
		endDate = startDate;
		endDate = endDate.AddMonths(itsTenor);
		resetDate = static_cast<ARM_Date>(itsResetDate[itsResetDate.size()-1]);

		for (i = 0; i < itsListStrikeModificators.size(); ++i)
		{
			newSwaption = new ARM_Swaption(startDate,
											endDate,
											ReceiveOrPay,
											ExerciseType,
											strike - itsListStrikeModificators[i], // FIXME: possible probleme de convention entre pourcentage -0.5 ou -0.005
											resetDate,
											(ARM_INDEX_TYPE) itsSwaptionLiborType);
			newSwaption->SetModel(BSModel);
			newSwaption->ComputePrice();
			itsCalibrationSwaptionList.push_back(newSwaption);			
		}		
	}	

	void ListsToPortfolios()
	{
		ARM_BSModel* BSModel	= dynamic_cast< ARM_BSModel* >		( itsMktDataManager.GetData(itsMktDataManagerKeys[ARM_VolBondCalculator::mdmKeysCalibrationAlias::BSModelKey]) );
		if (BSModel == 0) ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : BSMODEL not of good type ");

		const double weight = 1.;
		const double precision = 1.e-3;

		int i;
		std::list<ARM_Security*>::const_iterator ite;

		double marketprice;

	    //ctor ARM_StdPortfolio(int sz, double* weights, double* mktprices, double* precision = NULL)

		itsBootstrapSwaptionPortfolio = ARM_StdPortfolioPtr(new ARM_StdPortfolio(itsBootstrapSwaptionList));   
		
		ite = itsBootstrapSwaptionList.begin();
		for(i=0; ite != itsBootstrapSwaptionList.end(); ++i, ++ite)
		{
			itsBootstrapSwaptionPortfolio->SetWeight(weight,i);
			itsBootstrapSwaptionPortfolio->SetPrecision(precision,i);
			itsBootstrapSwaptionPortfolio->SetPrice((*ite)->GetPrice(), i);
		}

		marketprice = itsBootstrapSwaptionPortfolio->ComputePrice();
		itsBootstrapSwaptionPortfolio->SetPrice(marketprice);

		itsCalibrationSwaptionPortfolio = ARM_StdPortfolioPtr(new ARM_StdPortfolio(itsCalibrationSwaptionList));

		ite = itsCalibrationSwaptionList.begin();
		for(i=0; ite != itsCalibrationSwaptionList.end(); ++i, ++ite)
		{
			itsCalibrationSwaptionPortfolio->SetWeight(weight,i);
			itsCalibrationSwaptionPortfolio->SetPrecision(precision,i);
			itsCalibrationSwaptionPortfolio->SetPrice((*ite)->GetPrice(),i);
		}

		marketprice = itsCalibrationSwaptionPortfolio->ComputePrice();
		itsCalibrationSwaptionPortfolio->SetPrice(marketprice);
	}

	void Action()
	{
		BuildSwapList();
		BuildSwaptionList();
		ListsToPortfolios();
	}

public:		
	ARM_StdPortfolioPtr			getBootstrapSwaptionPortfolio()
	{
		if (itsDone == false)
		{
			Action();
			itsDone = true;		
		}

		ARM_StdPortfolio* ptr = dynamic_cast<ARM_StdPortfolio*>(itsBootstrapSwaptionPortfolio->Clone());
		if (ptr == 0) ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : problem inside portfolio creation ");
		return ARM_StdPortfolioPtr(ptr);
	}

	ARM_StdPortfolioPtr			getCalibrationSwaptionPortfolio()
	{
		if (itsDone == false)
		{
			Action();
			itsDone = true;
		}

		ARM_StdPortfolio* ptr = dynamic_cast<ARM_StdPortfolio*>(itsCalibrationSwaptionPortfolio->Clone());
		if (ptr == 0) ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : problem inside portfolio creation ");
		return ARM_StdPortfolioPtr(ptr);
	}
private:
	const std::vector<double>&				itsResetDate;
	const std::vector<double>&				itsStartDate;
	const ARM_MarketData_ManagerRep&	itsMktDataManager;
	const ARM_StringVector&				itsMktDataManagerKeys;

	std::vector<double>						itsSwapPrice;
	std::list<ARM_Security*>			itsBootstrapSwaptionList;
	std::list<ARM_Security*>			itsCalibrationSwaptionList;

	ARM_StdPortfolioPtr					itsBootstrapSwaptionPortfolio;
	ARM_StdPortfolioPtr					itsCalibrationSwaptionPortfolio;	

	int									itsTenor;

	const std::vector<double>&				itsListStrikeModificators;

	int									itsSwapLiborType;
	int									itsSwaptionLiborType;
	
	bool								itsDone;
};




const string ColNamesTableTypeI = "ResetDate StartDate FirstResetDate FirstStartDate PayDate Nominal Leverage InterestTerms ZcPay CMSRate Coupon";
const string ColNamesTableTypeII = "ResetDate StartDate FirstResetDate FirstStartDate PayDate Nominal Leverage InterestTerms ZcPay CMSRate MinRate MaxRate Coupon";
const string ColNamesTableLookBack = "ResetDate StartDate FirstResetDate FirstStartDate PayDate Nominal Leverage InterestTerms ZcPay CMSRate MinRate MaxRate MaxMinOption AccruedMaxMinOption Coupon";

	// 3 steps: 
	//		1 retrieve attributes from arguments
	//		2 create the deal description
	//		3 fill model's attributes (chosen by the developer)

	ARM_VolBondCalculator::ARM_VolBondCalculator(
			// market data manager
			const ARM_MarketData_ManagerRep&	MktDataManager,
			const ARM_StringVector&				MdmKeys,

			// deal parameters
			const double&						Nominal,
			const double&						Strike,
			const double&						Leverage,

			const ARM_Date	&					StartDate,
			const ARM_Date	&					EndDate,

			const long &						PayFreq,
			const long &						ResetFreq,

			const long &						PayDayCount,
			const long &						ResetDayCount,

			const string &						Tenor,

			const long &						IntRule,
			const long &						StubRule,

			const double &						ResetGap,

			const string &						PayCalendar,
			const string &						ResetCalendar,

			// modelparams: Runge-Kutta
			const long &	 					SolverType, 
			const std::vector<double> &				SolverParams, 

			// modelparams: MonteCarlo
			const long &						NbSteps,
			const long &						StepsPerYear,
			const long &		 				BucketSize,
			const vector<string> &				RandGenerator,

			const int&							payofftype,
			const string&						underlying,

			const bool&							CalibrationFlag,

			// output
			const std::deque<bool>&				ProductsToPrice
		)
		: // step 1

			ARM_GenCalculator(MktDataManager),
			itsHasBeenPriced(false),
			itsEventScheduleForwardStart(0),
			
			itsCcy(PayCalendar.c_str()),

			itsNominal(Nominal),
			itsStrike(Strike),
			itsLeverage(Leverage),

			itsStartDate(StartDate),
			itsEndDate(EndDate),

			itsPayFreq(PayFreq),
			itsResetFreq(ResetFreq),

			itsPayDayCount(PayDayCount),
			itsResetDayCount(ResetDayCount),

			itsTenor(Tenor),

			itsIntRule(IntRule),
			itsStubRule(StubRule),

			itsResetGap(ResetGap),

			itsPayCalendar(PayCalendar),
			itsResetCalendar(ResetCalendar),


			itsSolverType(SolverType), 
			itsSolverParams(SolverParams), 			

			itsNbSteps(NbSteps),
			itsStepsPerYear(StepsPerYear),
			itsBucketSize(BucketSize),
			itsRandGenerator(RandGenerator),

			itsPayoffType(payofftype),
			itsUnderlying(underlying),

			itsCalibrationFlag(CalibrationFlag),

			itsProductsToPrice(ProductsToPrice),

			itsPrice(0.)
	{	
		SetName(ARM_VOLBOND_CALCULATOR);		

		 //step 2
		SetKeys(MdmKeys);

		CheckData();

		SetCurrencyUnit(&itsCcy);

		SetColNamesTable();
		
		if (itsCalibrationFlag == true)
			CreateAndSetCalibration();		

		CreateAndSetDealDescriptionAndTimeIt(""); // security is setted
		CreateAndSetModel();
		
	}

	void ARM_VolBondCalculator::SetColNamesTable()
	{
		if (itsPayoffType == payoffType::TypeI)
		{
			istringstream iss(ColNamesTableTypeI);
			copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter(ColNamesTable));
		}
		else if (itsPayoffType == payoffType::TypeII)					
		{
			istringstream iss(ColNamesTableTypeII);
			copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter(ColNamesTable));
		}
		else if (itsPayoffType == payoffType::LookBack)					
		{
			istringstream iss(ColNamesTableLookBack);
			copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter(ColNamesTable));
		} 
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : PayOffType not implemented yet ");
	}

	ARM_VolBondCalculator::ARM_VolBondCalculator(const ARM_VolBondCalculator&rhs)
	{
		this->operator =(rhs);
	}

	ARM_VolBondCalculator& ARM_VolBondCalculator::operator=(const ARM_VolBondCalculator& rhs)
	{
		if (this != &rhs)
		{
			ARM_GenCalculator::operator=(rhs);
	
			ColNamesTable = rhs.ColNamesTable;

			itsCcy = rhs.itsCcy;

			itsHasBeenPriced = rhs.itsHasBeenPriced;

			itsNominal = rhs.itsNominal;
			itsStrike = rhs.itsStrike;
			itsLeverage = rhs.itsLeverage;

			itsStartDate = rhs.itsStartDate;
			itsEndDate = rhs.itsEndDate;

			itsPayFreq = rhs.itsPayFreq;
			itsResetFreq = rhs.itsResetFreq;

			itsPayDayCount = rhs.itsPayDayCount;
			itsResetDayCount = rhs.itsResetDayCount;

			itsTenor = rhs.itsTenor;

			itsIntRule = rhs.itsIntRule;
			itsStubRule = rhs.itsStubRule;

			itsResetGap = rhs.itsResetGap;

			itsPayCalendar = rhs.itsPayCalendar;
			itsResetCalendar = rhs.itsResetCalendar;

			itsPayoffType = rhs.itsPayoffType;
			itsUnderlying = rhs.itsUnderlying;

			itsSolverType = rhs.itsSolverType; 
			itsSolverParams = rhs.itsSolverParams; 


			itsNbSteps = rhs.itsNbSteps;
			itsStepsPerYear = rhs.itsStepsPerYear;
			itsBucketSize = rhs.itsBucketSize;
			itsRandGenerator = rhs.itsRandGenerator;


			itsProductsToPrice = rhs.itsProductsToPrice;

			itsCalibrationFlag = rhs.itsCalibrationFlag;

			itsPrice = rhs.itsPrice;

			itsEventScheduleForwardStart.operator =(rhs.itsEventScheduleForwardStart);
		}

		return *this;
	}

	ARM_VolBondCalculator::~ARM_VolBondCalculator()
	{
	}

	ARM_Object*	ARM_VolBondCalculator::Clone() const
	{
		return new ARM_VolBondCalculator(*this);
	}



    void ARM_VolBondCalculator::CheckData() 
	{
		// Check internal data consistency

		CheckMktData();

		if (itsCalibrationFlag == true)
			CheckCalibrationData();
	} 
	
    void ARM_VolBondCalculator::CheckCalibrationData() 
	{
		ARM_BSModel* BSModel	= dynamic_cast< ARM_BSModel* >		( GetMktDataManager()->GetData(GetKeys()[mdmKeysCalibrationAlias::BSModelKey]) );
		if (BSModel == 0) ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : BSMODEL not of good type ");

		ARM_ModelParam* volParam = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[mdmKeysCalibrationAlias::VolCalibKey]));
		if(!volParam || volParam->GetType() != ARM_ModelParamType::Volatility)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Volatility Param for key=" + GetKeys()[mdmKeysCalibrationAlias::VolCalibKey] + " is expected in the Market Data Manager");
		
		ARM_ModelParam* vovParam = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[mdmKeysCalibrationAlias::VovCalibKey]));
		if(!vovParam || vovParam->GetType() != ARM_ModelParamType::VolOfVol)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : VolOfVol Param for key=" + GetKeys()[mdmKeysCalibrationAlias::VovCalibKey] + " is expected in the Market Data Manager");
		
		ARM_ModelParam* correlParam = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[mdmKeysCalibrationAlias::CorrelCalibKey]));
		if(!correlParam || correlParam->GetType() != ARM_ModelParamType::Correlation)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Correlation Param for key=" + GetKeys()[mdmKeysCalibrationAlias::CorrelCalibKey] + " is expected in the Market Data Manager");
	}


	void ARM_VolBondCalculator::CheckMktData() 
	{
		// Check mkt data consistency
		ARM_ZeroCurve* ycCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[ZcCurveKey]));
		if(!ycCurve)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Yield Curve for key=" + GetKeys()[ZcCurveKey] + " is expected in the Market Data Manager");

		ARM_ModelParam* volParam = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[VolCurveKey]));
		if(!volParam || volParam->GetType() != ARM_ModelParamType::Volatility)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Volatility Curve Param for key=" + GetKeys()[VolCurveKey] + " is expected in the Market Data Manager");

		ARM_ModelParam* volOfVolParam = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[VolOfVolKey]));
		if(!volOfVolParam || volOfVolParam->GetType() != ARM_ModelParamType::VolOfVol)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Volatility of Volatility Param for key=" + GetKeys()[VolOfVolKey] + " is expected in the Market Data Manager");

		ARM_ModelParam* mrsParam = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[MrsKey]));
		if(!mrsParam || mrsParam->GetType() != ARM_ModelParamType::MeanReversion)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : MRS Param for key=" + GetKeys()[MrsKey] + " is expected in the Market Data Manager");

		ARM_ModelParam* correlationParam = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[CorrelationKey]));
		if(!correlationParam || correlationParam->GetType() != ARM_ModelParamType::Correlation)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Correlation Param for key=" + GetKeys()[CorrelationKey] + " is expected in the Market Data Manager");

		ARM_ModelParam* volmrsParam = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[VolMrsKey]));
		if(!mrsParam || volmrsParam->GetType() != ARM_ModelParamType::VolMeanReversion)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Volatility of Mean Reversion Param for key=" + GetKeys()[VolMrsKey] + " is expected in the Market Data Manager");

	}

	void ARM_VolBondCalculator::CreateAndSetModel() 
	{
		ARM_CountedPtr<ARM_HWSV1F> myModel;

		try
		{
			ARM_ModelParamVector* modelParams = new ARM_ModelParamVector();

		// Check mkt data consistency
			ARM_ZeroCurve* ycCurve = static_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[ZcCurveKey]));

			ARM_ModelParam* volCurveParam = static_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[VolCurveKey]));
			ARM_ModelParam* mrsParam = static_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[MrsKey]));
			ARM_ModelParam* volOfVolParam = static_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[VolOfVolKey]));
			ARM_ModelParam* correlationParam = static_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[CorrelationKey]));
			ARM_ModelParam* volMeanReversionParam = static_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[VolMrsKey]));

			modelParams->push_back(volCurveParam);
			modelParams->push_back(mrsParam);
			
			modelParams->push_back(volOfVolParam);
			modelParams->push_back(correlationParam);	
			modelParams->push_back(volMeanReversionParam);					

			ARM_ModelParamsHWSV HWSVModelParams(*modelParams);

			std::vector<double>							solverDatas(itsSolverParams);

			int formulatypeSO = 0;
			vector<double> formulaParamsSO;
			std::vector<double> formulaDatasSO(formulaParamsSO);
			double maxdecay = 0.;

			myModel = ARM_CountedPtr<ARM_HWSV1F>(new ARM_HWSV1F( CreateClonedPtr( ycCurve ), &HWSVModelParams, 
				itsSolverType, solverDatas, formulatypeSO, formulaDatasSO, maxdecay));
			/// creates the default numeraire for the time being
			ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::Cash ) );
			myModel->SetNumeraire( numeraire );

			SetPricingModel(myModel);

			/// Create MC Method

			// random generator
			ARM_RandomGenerator* randGen1=NULL;
			ARM_RandomGenerator* randGen2=NULL;
			ARM_RandomGenerator* randGen=NULL;

			ARM_RandomGeneratorPtrVector randGenVector;

			ARM_RandomGeneratorPtr pRandGen1, pRandGen2;
				
			randGen1 = ARM_RandGenFactory.Instance()->CreateRandGen(
					(ARM_RandGenFactoryImp::BaseGenType) ARM_ArgConv_BaseGenAlgoType.GetNumber(itsRandGenerator[0]));
			pRandGen1 = ARM_RandomGeneratorPtr(randGen1);

			randGen2 = ARM_RandGenFactory.Instance()->CreateRandGen(
					ARM_RandGenFactoryImp::UnknownBaseGenAlgorithm,
					(ARM_RandGenFactoryImp::TransformAlgo) ARM_ArgConv_TransformAlgoType.GetNumber(itsRandGenerator[1]),
					pRandGen1);

			if (itsRandGenerator[2] != "Nothing")
			{
				pRandGen2 = ARM_RandomGeneratorPtr(randGen2);
				randGen = ARM_RandGenFactory.Instance()->CreateRandGen(
					ARM_RandGenFactoryImp::UnknownBaseGenAlgorithm,
					(ARM_RandGenFactoryImp::TransformAlgo) ARM_ArgConv_TransformAlgoType.GetNumber("AntitheticOne"),
					pRandGen2);
			}
			else
				randGen = randGen2;

			randGenVector.push_back(static_cast<ARM_RandomGeneratorPtr>(randGen));			

	ARM_SchedulerBase* scheduler = NULL;
	ARM_SamplerBase* sampler = NULL;
	ARM_ImpSamplerPtr impSampler(NULL);
	ARM_PathSchemePtr pathScheme(NULL);
	ARM_MCMethod* mcMethod= NULL;
	ARM_MultiCurve* alpha;
	bool isAlpha = false;

        std::vector<double> schedulerDatas;
		schedulerDatas.push_back(itsStepsPerYear);

		int nbSteps;
		int schedulerType = 3;

		scheduler = ARM_SchedulerFactory.Instance()->CreateScheduler(
			schedulerType,
			schedulerDatas,
			nbSteps);

		//nbSteps=25?  by default but it seems to be obligatory actually

		// Sampler
		// We need a multi dimension sampler
		int multiDim = 2;
		int samplerType = 0;
		std::vector<double> samplerDatas;

		sampler = ARM_SamplerFactory.Instance()->CreateSampler(
			multiDim,
			samplerType,
			samplerDatas,
			scheduler);

		// Imp Sampler
		int impSamplerType;

		std::vector<double> abscisses(1, 0.0);
		ARM_GP_T_Vector<std::vector<double> > ordinates(1, std::vector<double>(1,0.0));

		alpha = new ARM_MultiCurve(abscisses,ordinates,new ARM_LinInterpCstExtrapolVec);
		isAlpha = true;

		impSamplerType = 0;

		impSampler = ARM_ImpSamplerPtr(ARM_ImpSamplerFactory.Instance()->CreateImpSampler(
			impSamplerType,
			alpha));


			myModel->SetNumMethod(ARM_CountedPtr<ARM_NumMethod> (new ARM_MCMethod(
				itsNbSteps,
				randGenVector,
				sampler, 
				itsNbSteps/itsBucketSize,
				impSampler,
				pathScheme)));
		}
		catch(Exception& x)
		{
			x.DebugPrint();
			throw x;
		}

	}


	ARM_RowInfo ARM_VolBondCalculator::ColumnNames() const
	{		
		vector< ARM_GP_VALUE_TYPE >		colTypeVec(ColNamesTable.size(), ARM_STRING); 
		return ARM_RowInfo(ColNamesTable,colTypeVec);		
	}

/*****************************************************************************
/ Class  : ARM_VolBondCalculator
/ Routine: DatesStructure
/ Returns: ARM_DateStripVector
/	 return reset start date scheduler
******************************************************************************/
	 ARM_DateStripCombiner ARM_VolBondCalculator::DatesStructure() const 
	{
		try
		{
			ARM_DateStrip DateStripAdv( 
			(ARM_Date) itsStartDate,
			(ARM_Date) itsEndDate,
			itsResetFreq,
			itsPayDayCount,
			GETDEFAULTVALUESTR,
			K_MOD_FOLLOWING,
			itsIntRule,
			itsStubRule,
			itsResetGap,
			itsPayFreq,
			GETDEFAULTVALUE,
			itsResetCalendar.c_str(),
			K_ADVANCE);

			ARM_DateStrip DateStripArr(
			(ARM_Date) itsStartDate,
			(ARM_Date) itsEndDate,
			itsResetFreq,
			itsPayDayCount,
			GETDEFAULTVALUESTR,
			K_MOD_FOLLOWING,
			itsIntRule,
			itsStubRule,
			itsResetGap,
			itsPayFreq,
			GETDEFAULTVALUE,
			itsResetCalendar.c_str(),
			K_ARREARS);

			ARM_DateStripVector SchedVect(NB_VOLBOND_SCHED, 0); 
			SchedVect[VOLBOND_SCHED_ADV] = &DateStripAdv;
			SchedVect[VOLBOND_SCHED_ARR] = &DateStripArr;

			ARM_DateStripCombiner EventScheduleResetDate(SchedVect,"ResetDate");

			itsEventScheduleForwardStart = ARM_CountedPtr<ARM_DateStripCombiner>(new ARM_DateStripCombiner(SchedVect,"FwdStartDate"));

			/// Initialise the memorisation of the first call line
			double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
			ARM_VectorPtr eventDates = EventScheduleResetDate.GetMergeData();
			ARM_VectorPtr eventDatesfwdstdates = itsEventScheduleForwardStart->GetMergeData();

			return EventScheduleResetDate;

		}
		catch (...)
		{
			ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in ARM_VolBondCalculator::DatesStructure" );
		}
	 }

/////////////////////////////////////////////////////////////////
///	Class  : ARM_VolBondCalculator
///	Routine: CustomDatesStructure
///	Returns: ARM_DateStripCombiner
///	 return forward start date scheduler
/////////////////////////////////////////////////////////////////
	 ARM_DateStripCombiner ARM_VolBondCalculator::CustomDatesStructure(const ARM_DateStripVector& dateStrips) const 
	 {
		return ARM_DateStripCombiner();
	 }


	 void ARM_VolBondCalculator::ComputePricingData() const
	 {
		if (!itsHasBeenPriced)		
			const_cast<ARM_VolBondCalculator*>(this)->PriceAndTimeIt();
	 }




/////////////////////////////////////////////////////////////////
///	Class  : ARM_VolBondCalculator
///	Routine: UpdateModel
///	Returns: void
///	Action : update the model datas w.r.t. MarketDataManager
/////////////////////////////////////////////////////////////////
     void ARM_VolBondCalculator::UpdateModel() 
	 {
		itsHasBeenPriced = false;

		/// get zc curve
		ARM_ZeroCurve* zcCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[ZcCurveKey]));

		/// set it to pricing model
		GetPricingModel()->SetZeroCurve(CreateClonedPtr( zcCurve ));
		
		/// this takes into account the new mean reversion if its has been changed
		/// un peu bourrin mais bon on est plus à ca près
		CreateAndSetModel();
	 }


// EXECUTION

	void ARM_VolBondCalculator::InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const
	{
		string zeroValue("0");
				
		rowDescVec[Index("Coupon")] = zeroValue;
		rowTypeVec[Index("Coupon")] = ARM_DOUBLE;
	}

	size_t ARM_VolBondCalculator::Index(const char* str) const
	{
		size_t index = 0;
		for (; index != ColNamesTable.size() && ColNamesTable[index] != str; ++index)
			;
		if (ColNamesTable[index] != str)
			ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : no " + str + " column" );

		return index;
	}

	ARM_RowInfo ARM_VolBondCalculator::MiddleRows( size_t eventIdx, const ARM_DateStripCombiner& datesStructure ) const
	{
		try
		{

			// create a row of a deal description

			// parameters:
			//	eventIdx means the i th row
			//	contains all dates

			vector< string >				rowDescVec(ColNamesTable.size(), "");
			vector< ARM_GP_VALUE_TYPE >		rowTypeVec(ColNamesTable.size(), ARM_MISSING_TYPE);

			ARM_GP_VectorPtr& resetDates = datesStructure.GetMergeData();
			ARM_GP_VectorPtr& fwdstartDates = itsEventScheduleForwardStart->GetMergeData();
					
			InitPriceableColumns(rowDescVec, rowTypeVec); // initializes with "0"

			//RESET DATE
			CC_Ostringstream resetDateDesc;
			resetDateDesc << CC_NS(std,fixed) << (*resetDates)[eventIdx];			
			rowDescVec[Index("ResetDate")] = resetDateDesc.str();
			rowTypeVec[Index("ResetDate")] = ARM_DATE_TYPE;

			//FORWARD START DATE
			CC_Ostringstream startDateDesc;
			startDateDesc << CC_NS(std,fixed) << (*fwdstartDates)[eventIdx];
			rowDescVec[Index("StartDate")] = startDateDesc.str();
			rowTypeVec[Index("StartDate")] = ARM_DATE_TYPE;

			static double oldreset;
			static double oldstart;

			if (eventIdx >= 1)
			{
				//FIRST RESET DATE
				CC_Ostringstream firstresetDateDesc;
				firstresetDateDesc << CC_NS(std,fixed) << oldreset;
				rowDescVec[Index("FirstResetDate")] = firstresetDateDesc.str();
				rowTypeVec[Index("FirstResetDate")] = ARM_DATE_TYPE;

				//FIRST FORWARD START DATE
				CC_Ostringstream firststartDateDesc;
				firststartDateDesc << CC_NS(std,fixed) << oldstart;
				rowDescVec[Index("FirstStartDate")] = firststartDateDesc.str();
				rowTypeVec[Index("FirstStartDate")] = ARM_DATE_TYPE;
			}
			oldstart = (*fwdstartDates)[eventIdx];
			oldreset = (*resetDates)[eventIdx];

			//PAYDATE	
			double payDate;
			if (eventIdx == 0)		
				payDate = (*(datesStructure.GetDateStrip(VOLBOND_SCHED_ADV)->GetPaymentDates()))[0];
			else
				payDate = (*(datesStructure.GetDateStrip(VOLBOND_SCHED_ADV)->GetPaymentDates()))[eventIdx - 1];			
			CC_Ostringstream payDateDesc;
			payDateDesc << CC_NS(std,fixed) << payDate;
			rowDescVec[Index("PayDate")] = payDateDesc.str();
			rowTypeVec[Index("PayDate")] = ARM_DATE_TYPE;


			//NOMINAL
			CC_Ostringstream nominalDesc;
			nominalDesc << CC_NS(std,fixed) << itsNominal;
			rowDescVec[Index("Nominal")] = nominalDesc.str();
			rowTypeVec[Index("Nominal")] = ARM_STRING;
			
			//leverage
			CC_Ostringstream leverageDesc;
			leverageDesc << CC_NS(std,fixed) << itsLeverage;
			rowDescVec[Index("Leverage")] = leverageDesc.str();
			rowTypeVec[Index("Leverage")] = ARM_STRING;

			//INTERESTTERMS
			CC_Ostringstream itDesc;
			if (eventIdx == 0)						
				itDesc << CC_NS(std,fixed) << (*(datesStructure.GetDateStrip(VOLBOND_SCHED_ADV)->GetInterestTerms()))[0] ;
			else
				itDesc << CC_NS(std,fixed) << (*(datesStructure.GetDateStrip(VOLBOND_SCHED_ADV)->GetInterestTerms()))[eventIdx - 1] ;
			rowDescVec[Index("InterestTerms")] = itDesc.str();
			rowTypeVec[Index("InterestTerms")] = ARM_STRING;

			// DF
			CC_Ostringstream zcpayDesc;
			zcpayDesc << "DF(" << itsPayCalendar << ",PayDate[i])";
			rowDescVec[Index("ZcPay")] = zcpayDesc.str();
			rowTypeVec[Index("ZcPay")] = ARM_STRING;

			// CMSRATE
			CC_Ostringstream cmsrateDesc;
			cmsrateDesc << "SWAPRATE(" << itsPayCalendar << ",StartDate[i]," << itsTenor << "," << ARM_ArgConvReverse_StdFrequency.GetString(itsResetFreq) << "," << ARM_ArgConvReverse_DayCount.GetString(itsPayDayCount) << ")";
			rowDescVec[Index("CMSRate")] = cmsrateDesc.str();
			rowTypeVec[Index("CMSRate")] = ARM_STRING;

			if (itsPayoffType == payoffType::TypeII || itsPayoffType == payoffType::LookBack)
			{
				// MINRATE
				CC_Ostringstream minrateDesc;
				if (eventIdx >= 1)
					minrateDesc << "MINMAX(" << itsPayCalendar << ",StartDate[i]," << itsTenor << ",FirstResetDate[i],FirstStartDate[i],CMSRate[i-1],MIN," << ARM_ArgConvReverse_StdFrequency.GetString(itsResetFreq) << ")";
				else
					minrateDesc << "0";
				rowDescVec[Index("MinRate")] = minrateDesc.str();
				rowTypeVec[Index("MinRate")] = ARM_STRING;

				// MAXRATE
				CC_Ostringstream maxrateDesc;
				if (eventIdx >= 1)
					maxrateDesc << "MINMAX(" << itsPayCalendar << ",StartDate[i]," << itsTenor << ",FirstResetDate[i],FirstStartDate[i],CMSRate[i-1],MAX," << ARM_ArgConvReverse_StdFrequency.GetString(itsResetFreq) << ")";
				else
					maxrateDesc << "0";
				rowDescVec[Index("MaxRate")] = maxrateDesc.str();
				rowTypeVec[Index("MaxRate")] = ARM_STRING;
			}

			// MaxMinOption
			if (itsPayoffType == payoffType::LookBack)
			{
				CC_Ostringstream MaxMinOptionDesc;
				if (eventIdx >= 1)
					MaxMinOptionDesc << "MINMAX(" << itsPayCalendar << ",StartDate[i]," << itsTenor << ",FirstResetDate[i],FirstStartDate[i],CMSRate[i-1],MAXMIN," << ARM_ArgConvReverse_StdFrequency.GetString(itsPayFreq) << ","<< itsStrike<<",FLOOR)";
				else
					MaxMinOptionDesc << "0";
				rowDescVec[Index("MaxMinOption")] = MaxMinOptionDesc.str();
				rowTypeVec[Index("MaxMinOption")] = ARM_STRING;			

				CC_Ostringstream AccruedMaxMinOptionDesc;
				if (eventIdx >= 1)
					AccruedMaxMinOptionDesc << "MINMAX(" << itsPayCalendar << ",StartDate[i]," << itsTenor << ",FirstResetDate[i],FirstStartDate[i],CMSRate[i-1],MAXMIN," << ARM_ArgConvReverse_StdFrequency.GetString(itsPayFreq) << "," << itsStrike << ",FLOOR,0.54, AccruedMinMax)";
				else
					AccruedMaxMinOptionDesc << "0";
				rowDescVec[Index("AccruedMaxMinOption")] = AccruedMaxMinOptionDesc.str();
				rowTypeVec[Index("AccruedMaxMinOption")] = ARM_STRING;
			}

			//COUPONRATE
			CC_Ostringstream couponRateDesc;
			if (eventIdx >= 1)
			{
				if (itsPayoffType == payoffType::TypeI)
					couponRateDesc << "ABS(Leverage[i]*(CMSRate[i]-CMSRate[i-1]))*InterestTerms[i]*ZcPay[i]*Nominal[i]";
				else if (itsPayoffType == payoffType::TypeII)					
					couponRateDesc << "Leverage[i]*(MaxRate[i]-MinRate[i])*InterestTerms[i]*ZcPay[i]*Nominal[i]";
				else if (itsPayoffType == payoffType::LookBack)				
					couponRateDesc << "Leverage[i]*MaxMinOption[i]*InterestTerms[i]*ZcPay[i]*Nominal[i]";			
			}
			else
				couponRateDesc << "0";
			rowDescVec[Index("Coupon")] = couponRateDesc.str();
			rowTypeVec[Index("Coupon")] = ARM_STRING;

			return ARM_RowInfo(rowDescVec,rowTypeVec);
		
		}
		catch (Exception& x)
		{
			x.DebugPrint();

			throw x;
		}
		catch (...)
		{
			ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in VolBondCalculator::MiddleRows" );
		}

	}


    double	ARM_VolBondCalculator::Price() 
	{
		// pas encore: CalibrateAndTimeIt();

		// only Coupon for now : size_t NbProductsToPrice = itsProductsToPrice.size(); // FIXMEFRED: faire MIN(nb_output max total, itsProductsToPrice.size())
			
		if (itsHasBeenPriced == true)
			return itsPrice;

		if (itsCalibrationFlag == true)
			CalibrateAndTimeIt();

		ARM_GenPricer genPricer(&*GetGenSecurity(), &*GetPricingModel());
		itsHasBeenPriced = true;
		itsPrice = genPricer.Price();
		return itsPrice;
	}



/////////////////////////////////////////////////////////////////
///	Routine: Calibrate
///	Returns: void
/////////////////////////////////////////////////////////////////
	void ARM_VolBondCalculator::CreateAndSetCalibration() 
	{
		ARM_ModelParam* volParam = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[VolCalibKey]));
		if(!volParam || volParam->GetType() != ARM_ModelParamType::Volatility)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Volatility Param for key=" + GetKeys()[VolCalibKey] + " is expected in the Market Data Manager");
		
		ARM_ModelParam* vovParam = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[VovCalibKey]));
		if(!vovParam || vovParam->GetType() != ARM_ModelParamType::VolOfVol)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : VolOfVol Param for key=" + GetKeys()[VovCalibKey] + " is expected in the Market Data Manager");
		
		ARM_ModelParam* correlParam = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[CorrelCalibKey]));
		if(!correlParam || correlParam->GetType() != ARM_ModelParamType::Correlation)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Correlation Param for key=" + GetKeys()[CorrelCalibKey] + " is expected in the Market Data Manager");

		try
		{
			ARM_GP_VectorPtr& resetDates = DatesStructure().GetMergeData();
			ARM_GP_VectorPtr& fwdstartDates = itsEventScheduleForwardStart->GetMergeData();

			std::vector<double> list_strike_modificators(2);
			list_strike_modificators[0] = -0.5;
			list_strike_modificators[1] = -1.0;

			SwaptionPortfolios portfolios(*GetMktDataManager(), GetKeys(), *fwdstartDates, *resetDates, list_strike_modificators);	
			
			std::string bootstrapTypeSolver("BRENT");
			std::string bootstrapTypeMethod("BOOTSTRAP1D");
			int bootstrap_max_iter = 20;
			double bootstrap_fxTol			= 1.e-6;
			double bootstrap_xTol		= SolverConstant::DefaultXTolerance;

			ARM_ModelFitterDes* bootstrapModelFitterDes = new ARM_ModelFitterDes((ARM_SolverType) ARM_ArgConv_SolverTypeMethod.GetNumber(bootstrapTypeSolver),
										(size_t)bootstrap_max_iter,
										bootstrap_fxTol,
										bootstrap_xTol);

			ARM_ModelParamVector volParamVec(1,volParam);

			ARM_CalibMethod* bootstrapMethod = new ARM_CalibMethod(portfolios.getBootstrapSwaptionPortfolio(),
										volParamVec,
										(ARM_MethodType) ARM_ArgConv_CalibMethod.GetNumber(bootstrapTypeMethod),
										bootstrapModelFitterDes);

			std::string calibTypeSolver("NAG_SOLVER");
			std::string calibTypeMethod("OPTIMIZE1D");
			int calib_max_iter = 15;
			double calib_xTol = SolverConstant::DefaultXTolerance;
			double calib_fxTol = 0.001;

			ARM_ModelFitterDes* calibModelFitterDes = new ARM_ModelFitterDes((ARM_SolverType) ARM_ArgConv_SolverTypeMethod.GetNumber(calibTypeSolver),
										(size_t)calib_max_iter,
										calib_fxTol,
										calib_xTol); 

			ARM_ModelParamVector vovParamVec(1,vovParam);

			ARM_CalibMethod* vovCalibMethod = new ARM_CalibMethod(portfolios.getCalibrationSwaptionPortfolio(),
										vovParamVec,
										(ARM_MethodType) ARM_ArgConv_CalibMethod.GetNumber(calibTypeMethod),
										calibModelFitterDes,
										ARM_CalibrationTarget::PriceTarget,
										bootstrapMethod);

			ARM_ModelParamVector correlParamVec(1,correlParam);

			bool validate_correl = false;
			ARM_CalibMethod* correlCalibMethod = new ARM_CalibMethod(portfolios.getCalibrationSwaptionPortfolio(),
										correlParamVec,
										(ARM_MethodType) ARM_ArgConv_CalibMethod.GetNumber(calibTypeMethod),
										calibModelFitterDes,
										ARM_CalibrationTarget::PriceTarget,
										vovCalibMethod,
										NULL,
										true,
										0,
										1,
										validate_correl
										);

	
			SetCalibMethod(ARM_CalibMethodPtr(correlCalibMethod));		
		}
		catch(Exception& x)
		{
			x.DebugPrint();
			throw x;
		}
		
		if (0 /* vov<-correl<-mrs*/)
		{
			// manque portfolio kernel7Bis
		}

	}
	void ARM_VolBondCalculator::Calibrate()
	{
		try
		{
			GetCalibMethod()->Calibrate(&*GetPricingModel()); 
		}
		catch(Exception& x)
		{
			x.DebugPrint();
		}

	}

	void ARM_VolBondCalculator::UpdateCalibration(bool isUpdateStrike) 
	{
	}
/////////////////////////////////////////////////////////////////



	  ARM_DateStripPtr ARM_VolBondCalculator::GetOriginalFundDateStrip() const  
	  {
		  ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : not implemented" ); return NULL;
	  }

	  ARM_DateStripPtr ARM_VolBondCalculator::GetRefFundDateStrip() const 
	  {
		  ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : not implemented" ); return NULL;
	  }

	  ARM_DateStripPtr ARM_VolBondCalculator::GetRefDateStrip() const 
	  {
		  ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : not implemented" ); return NULL;
	  }

	  std::vector<double> ARM_VolBondCalculator::GetvCpnNominal() const  
	  {
		  ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : not implemented" ); return NULL;
	  }

	  std::vector<double> ARM_VolBondCalculator::GetvFundNominal() const 
	  {
		  ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : not implemented" ); return NULL;
	  }

	  std::vector<double> ARM_VolBondCalculator::GetvFundSpread() const 
	  {
		  ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : not implemented" ); return NULL;
	  }

	  ARM_ZeroCurve* ARM_VolBondCalculator::GetDomesticZeroCurve() const  
	  {
		  ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : not implemented" ); return NULL;
	  }

	  ARM_ZeroCurve* ARM_VolBondCalculator::GetForeignZeroCurve() const  
	  {
		  ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : not implemented" ); return NULL;
	  }

	  ARM_ZeroCurve* ARM_VolBondCalculator::GetDomesticDiscountZeroCurve() const  
	  {
		  ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : not implemented" ); return NULL;
	  }

	  ARM_ZeroCurve* ARM_VolBondCalculator::GetForeignDiscountZeroCurve() const  
	  {
		  ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : not implemented" ); return NULL;
	  }

	  ARM_Forex* ARM_VolBondCalculator::GetForex() const  
	  {
		  ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : not implemented" ); return NULL;
	  }


CC_END_NAMESPACE()