#include "gpbase/removeidentifiedwarning.h"
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

#include "gpnumlib/argconvdefault.h"
#include "gpnumlib/randomgenfactory.h"
#include "gpnumlib/compositegen.h"
#include "gpnumlib/antitheticgen.h"
#include "gpnumlib/odefunctions.h"
#include "gpnumlib/random.h"

#include "gpcalib/calibmethod.h"



CC_BEGIN_NAMESPACE( ARM )

const string ZCCURVE        = "ZCCURVE";
const string VOLCURVE       = "VOLCURVE";
const string MRS			= "MRS";
const string VOLOFVOL       = "VOLOFVOL";
const string CORRELATION    = "CORRELATION";
const string VOLMRS			= "VOLMRS";

const unsigned int NB_VOLBOND_SCHED = 2;
const int VOLBOND_SCHED_ADV	= 0;
const int VOLBOND_SCHED_ARR	= 1;

const string ColNamesTableTypeI = "ResetDate StartDate FirstResetDate FirstStartDate PayDate Nominal InterestTerms ZcPay CMSRate Coupon";
const string ColNamesTableTypeII = "ResetDate StartDate FirstResetDate FirstStartDate PayDate Nominal InterestTerms ZcPay CMSRate MinRate MaxRate Coupon";
const string ColNamesTableLookBack = "ResetDate StartDate FirstResetDate FirstStartDate PayDate Nominal InterestTerms ZcPay CMSRate MinRate MaxRate MaxMinOption Coupon";

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

			const ARM_Date	&					StartDate,
			const ARM_Date	&					EndDate,

			const long &						PayFreq,
			const long &						ResetFreq,

			const long &						DayCount,
			const string &						Tenor,

			const long &						IntRule,
			const long &						StubRule,

			const double &						ResetGap,

			const string &						PayCalendar,
			const string &						ResetCalendar,

			// modelparams: Runge-Kutta
			const long &	 					SolverType, 
			const ARM_GP_Vector &				SolverParams, 

			// modelparams: MonteCarlo
			const long &						NbSteps,
			const long &						StepsPerYear,
			const long &		 				BucketSize,
			const vector<string> &				RandGenerator,

			const int&							payofftype,

			// output
			const std::deque<bool>&				ProductsToPrice
		)
		: // step 1

			ARM_GenCalculator(MktDataManager),
			itsHasBeenPriced(false),
			itsEventScheduleForwardStart(0),
			
			itsCcy(PayCalendar.c_str()),

			itsNominal(Nominal),

			itsStartDate(StartDate),
			itsEndDate(EndDate),

			itsPayFreq(PayFreq),
			itsResetFreq(ResetFreq),

			itsDayCount(DayCount),
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

			itsProductsToPrice(ProductsToPrice),

			itsPrice(0.),

			itsFundingPrice(0.),
			itsGlobalCapPrice(0.),
			itsCouponPrice(0.),
			itsSwapPrice(0.),
			itsProductPrice(0.),
			itsSumCoupPrice(0.),
			itsDFPrice(0.),
			itsFundRatePrice(0.),
			itsIndexPrice(0.)
	{	
		SetName(ARM_VOLBOND_CALCULATOR);		

		 //step 2
		SetKeys(MdmKeys);

		CheckData();

		SetCurrencyUnit(&itsCcy);

		SetColNamesTable();

		CreateAndSetDealDescriptionAndTimeIt(""); // security is setted
		CreateAndSetModel();
		CreateAndSetCalibration();
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
		// copy constructor
	}

	ARM_VolBondCalculator& ARM_VolBondCalculator::operator=(const ARM_VolBondCalculator& rhs)
	{
		if (this != &rhs)
		{
			ARM_GenCalculator::operator=(rhs);

			//itsCalibND					= rhs.itsCalibND;
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

		CheckCalibrationData();
	} 
	
    void ARM_VolBondCalculator::CheckCalibrationData() 
	{
		
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

		ARM_ModelParam* mrsParam = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[MrsKey]));
		if(!mrsParam || mrsParam->GetType() != ARM_ModelParamType::MeanReversion)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : MRS Param for key=" + GetKeys()[MrsKey] + " is expected in the Market Data Manager");

		ARM_ModelParam* volOfVolParam = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[VolOfVolKey]));
		if(!volOfVolParam || volOfVolParam->GetType() != ARM_ModelParamType::VolOfVol)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Volatility of Volatility Param for key=" + GetKeys()[VolOfVolKey] + " is expected in the Market Data Manager");

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
		// Check mkt data consistency
			ARM_ZeroCurve* ycCurve = static_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[ZcCurveKey]));

			ARM_ModelParam* volCurveParam = static_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[VolCurveKey]));
			ARM_ModelParam* mrsParam = static_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[MrsKey]));
			ARM_ModelParam* volOfVolParam = static_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[VolOfVolKey]));
			ARM_ModelParam* correlationParam = static_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[CorrelationKey]));
			ARM_ModelParam* volMeanReversionParam = static_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[VolMrsKey]));
	
			ARM_ModelParamVector* modelParams = new ARM_ModelParamVector();
			modelParams->push_back(volCurveParam);
			modelParams->push_back(mrsParam);
			modelParams->push_back(volOfVolParam);
			modelParams->push_back(correlationParam);
			modelParams->push_back(volMeanReversionParam);

			//FIXMEFRED warning: clone dans ARM_ModelParamsHWSV
			ARM_ModelParamsHWSV* HWSVModelParams = new ARM_ModelParamsHWSV(*modelParams);

			/// to avoid memory leak
			/// put an autocleaner on the model
			ARM_AutoCleaner<ARM_ModelParamsHWSV>	Hold(HWSVModelParams);
			ARM_GP_Vector							solverDatas(itsSolverParams);

			myModel = ARM_CountedPtr<ARM_HWSV1F>(new ARM_HWSV1F( CreateClonedPtr( ycCurve ), HWSVModelParams, 
				itsSolverType, solverDatas));

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

			ARM_TimeStepPerYearScheduler* scheduler = new ARM_TimeStepPerYearScheduler(itsStepsPerYear);
			ARM_SamplerBase * ptrsampler = new ARM_NormalCentredSamplerND(scheduler);			
			ARM_ImpSamplerPtr ptrimpsampler(new ARM_DummyImpSampler());
			ARM_PathSchemePtr  ptrpathscheme(new  ARM_IncrementalPathScheme());

			myModel->SetNumMethod(ARM_CountedPtr<ARM_NumMethod> (new ARM_MCMethod(
									itsNbSteps,
									randGenVector,
									ptrsampler, 
									itsNbSteps/itsBucketSize,
									ptrimpsampler,
									ptrpathscheme)));
		}
		catch(Exception& x)
		{
			x.DebugPrint();			
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
			itsDayCount,
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
			itsDayCount,
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

//		if (itsProductsToPrice[CSOPrice])
//			GetPricingData()[ "VBD" ] = itsVBDPrice;


	 }




/////////////////////////////////////////////////////////////////
///	Class  : ARM_VolBondCalculator
///	Routine: UpdateModel
///	Returns: void
///	Action : update the model datas w.r.t. MarketDataManager
/////////////////////////////////////////////////////////////////
     void ARM_VolBondCalculator::UpdateModel() 
	 {
		/// get zc curve
		ARM_ZeroCurve* zcCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[ZcCurveKey]));

		/// set it to pricing model
		GetPricingModel()->SetZeroCurve(CreateClonedPtr( zcCurve ));
		
		/// this takes into account the new mean reversion if its has been changed
		/// un peu bourrin mais bon on est pas à ca près
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
			ARM_VectorPtr& fwdstartDates = itsEventScheduleForwardStart->GetMergeData();
					
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
			cmsrateDesc << "SWAPRATE(" << itsPayCalendar << ",StartDate[i]," << itsTenor << "," << ARM_ArgConvReverse_StdFrequency.GetString(itsResetFreq) << "," << ARM_ArgConvReverse_DayCount.GetString(itsDayCount) << ")";
			rowDescVec[Index("CMSRate")] = cmsrateDesc.str();
			rowTypeVec[Index("CMSRate")] = ARM_STRING;

			if (itsPayoffType == payoffType::TypeII || itsPayoffType == payoffType::LookBack)
			{
				// MINRATE
				CC_Ostringstream minrateDesc;
				if (eventIdx >= 1)
					minrateDesc << "MINMAX(EUR,StartDate[i]," << itsTenor << ",FirstResetDate[i],FirstStartDate[i],CMSRate[i-1],MIN)";
				else
					minrateDesc << "0";
				rowDescVec[Index("MinRate")] = minrateDesc.str();
				rowTypeVec[Index("MinRate")] = ARM_STRING;

				// MAXRATE
				CC_Ostringstream maxrateDesc;
				if (eventIdx >= 1)
					maxrateDesc << "MINMAX(EUR,StartDate[i]," << itsTenor << ",FirstResetDate[i],FirstStartDate[i],CMSRate[i-1],MAX)";
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
					MaxMinOptionDesc << "MINMAX(EUR,StartDate[i]," << itsTenor << ",FirstResetDate[i],FirstStartDate[i],CMSRate[i-1],MAXMIN,D,0.008,FLOOR)";
				else
					MaxMinOptionDesc << "0";
				rowDescVec[Index("MaxMinOption")] = MaxMinOptionDesc.str();
				rowTypeVec[Index("MaxMinOption")] = ARM_STRING;			
			}

			//COUPONRATE
			CC_Ostringstream couponRateDesc;
			if (eventIdx >= 1)
			{
				if (itsPayoffType == payoffType::TypeI)
					couponRateDesc << "ABS(CMSRate[i]-CMSRate[i-1])*InterestTerms[i]*ZcPay[i]";
				else if (itsPayoffType == payoffType::TypeII)
					couponRateDesc << "(MaxRate[i]-MinRate[i])*InterestTerms[i]*ZcPay[i]";
				else if (itsPayoffType == payoffType::LookBack)				
					couponRateDesc << "MaxMinOption[i]*InterestTerms[i]*ZcPay[i]";			
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
			
		ARM_GenPricer genPricer(&*GetGenSecurity(), &*GetPricingModel());
		itsHasBeenPriced = true;
		return genPricer.Price();
	}




/////////////////////////////////////////////////////////////////
///	Routine: Calibrate
///	Returns: void
/////////////////////////////////////////////////////////////////
	void ARM_VolBondCalculator::CreateAndSetCalibration() 
	{
		try
		{
			// not used today but will be used for mean reversion parameter calibration	
		}
		catch(Exception& x)
		{
			x.DebugPrint();
		}


	}
	void ARM_VolBondCalculator::Calibrate()
	{
		try
		{
			// not used today but will be used for mean reversion parameter calibration	
		}
		catch(Exception& x)
		{
			x.DebugPrint();
		}

	}

	void ARM_VolBondCalculator::UpdateCalibration(bool isUpdateStrike) 
	{
		// not used today but will be used for mean reversion parameter calibration	
	}
/////////////////////////////////////////////////////////////////



	  ARM_DateStripPtr ARM_VolBondCalculator::GetOriginalFundDateStrip() const  
	  {
		  return NULL;
	  }

	  ARM_DateStripPtr ARM_VolBondCalculator::GetRefFundDateStrip() const 
	  {
		  return NULL;
	  }

	  ARM_DateStripPtr ARM_VolBondCalculator::GetRefDateStrip() const 
	  {
		  return NULL;
	  }

	  ARM_GP_Vector ARM_VolBondCalculator::GetvCpnNominal() const  
	  {
		  return NULL;
	  }

	  ARM_GP_Vector ARM_VolBondCalculator::GetvFundNominal() const 
	  {
		  return NULL;
	  }

	  ARM_GP_Vector ARM_VolBondCalculator::GetvFundSpread() const 
	  {
		  return NULL;
	  }

	  ARM_ZeroCurve* ARM_VolBondCalculator::GetDomesticZeroCurve() const  
	  {
		  return NULL;
	  }

	  ARM_ZeroCurve* ARM_VolBondCalculator::GetForeignZeroCurve() const  
	  {
		  return NULL;
	  }

	  ARM_ZeroCurve* ARM_VolBondCalculator::GetDomesticDiscountZeroCurve() const  
	  {
		  return NULL;
	  }

	  ARM_ZeroCurve* ARM_VolBondCalculator::GetForeignDiscountZeroCurve() const  
	  {
		  return NULL;
	  }

	  ARM_Forex* ARM_VolBondCalculator::GetForex() const  
	  {
		  return NULL;
	  }


CC_END_NAMESPACE()