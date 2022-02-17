//----------------------------------------------------------------------------
//
//   Group       : Credit Derivatives Research
//
//   Filename    : SCIDparameters.cpp
//
//   Date        : July, 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SCIDparameters.hpp"
#include "edginc/ICDSParSpreads.hpp"
#include "edginc/CDSParSpreads.hpp"
#include "edginc/EffectiveCurve.hpp"
#include "edginc/SwapTool.hpp"
#include <numeric>
#include <cmath>
#include <fstream>


DRLIB_BEGIN_NAMESPACE

/* Make sure the class links */
bool  SCIDparametersLoad() {
    return (SCIDparameters::TYPE != 0);
   }

SCIDparameters::SCIDparameters(const CClassConstSP& clazz) : MarketObject(clazz)
{
}

SCIDparameters::SCIDparameters()
	: MarketObject(SCIDparameters::TYPE),
	  name("sCID parameters Default Name"),
	  m_initialWeights(1,1.0)
{
}
/*
string SCIDparameters::discountYieldCurveName() const {
    return discount.getName();
}
*/

SCIDparameters::~SCIDparameters()
{
}

string SCIDparameters::getName() const
{
    return name;
}


/** populate from market cache - default implementation provided */
void SCIDparameters::getMarket(const IModel* model, const MarketData* market) {
    static const string method = "SCIDparameters::getMarket";

	try
	{
        market->GetReferenceDate(today);
	    /*=========================================================================
		 * GET THE DISCOUNT CURVE 
		 *=======================================================================
	    discount.getData(model, market);*/

		m_nbNames = parSpreads->size();
		m_nbWorlds = m_linearShift.size();
		survProb.resize(m_nbNames, 1.0);

		double sumNotional = 0;
		for (int m=0; m<m_nbNames; m++) sumNotional += m_notional[m];
		for (int m=0; m<m_nbNames; m++) m_notional[m] /= sumNotional ;

		DateTimeArray sgleNameDT;
		getSingleNameCalibrationDates(sgleNameDT);
		DoubleArray baseRates(m_nbNames);
		m_recovery.resize(m_nbNames);
		for (int m=0; m<parSpreads->size(); m++)
		{
		    (*parSpreads)[m].getData(model, market);
		    m_recovery[m]= (*parSpreads)[m].getSP()->getRecovery();
		}

		if ( (singleNameExpiries->size()!=baseWorldCalibSurvProba.numCols()) || (m_nbNames!=baseWorldCalibSurvProba.numRows()) )
		{
			baseWorldCalibSurvProba.resize(singleNameExpiries->size(), m_nbNames);
			for (int m=0; m<m_nbNames; m++)
			{
			    ICDSParSpreadsSP CDSm = (*parSpreads)[m].getSP();
			    for (int k=0; k<sgleNameDT.size(); k++)
				    baseWorldCalibSurvProba[k][m] = CDSm->survivalProb(sgleNameDT[k]);
			}
		}
		double firstDateAsDouble = today.yearFrac(sgleNameDT[0]);
		for (int m=0; m<parSpreads->size(); m++)
  	   baseRates[m] = - log( baseWorldCalibSurvProba[0][m] ) / firstDateAsDouble;
		

		IntArray realJumpProcesses;
		for (int i=0; i<jumpParameters->numCols(); i++)
			if ((*jumpParameters)[i][0] * (*jumpParameters)[i][1] * (*jumpParameters)[i][2] > 1e-15) realJumpProcesses.push_back(i);
		m_nbJumpProcesses = realJumpProcesses.size();
	    m_world.m_jumps.resize(m_nbJumpProcesses);
		for (int j=0; j<m_nbJumpProcesses; j++)
		{
			int i = realJumpProcesses[j];
			m_world.m_jumps[j].setParameters((*jumpParameters)[i][0], (*jumpParameters)[i][1], baseRates, (*jumpParameters)[i][2], (*jumpParameters)[i][3]);
		}

		DateTime cfDate = cfExpiry->toDate(today);
		m_world.CalibrateDiffusion(m_nbNames,
						   today,
						   sgleNameDT,
						   baseWorldCalibSurvProba,
						   cfDate,
						   (*diffParameters)[0],
						   (*diffParameters)[1],
						   (*diffParameters)[2],
						   (*diffParameters)[3]);
    }
    catch (exception& e)
	{
		throw ModelException(e, method);
	}
}

/**Check stuff here */
void SCIDparameters::validatePop2Object() 
{
    static const string routine("SCIDparameters::validatePop2Object");

	try 
    {
		if (m_parallelMultShift.size()!=m_linearShift.size())
		    throw ModelException("not the same number of parallel shifts and linear shifts", routine);
		if (m_initialWeights.size()!=m_linearShift.size())
			m_initialWeights.resize(m_linearShift.size(), 1.0/m_linearShift.size());
        //TODO: add this new check at the future release
        /*double weightSum = std::accumulate(m_initialWeights.begin(),m_initialWeights.end(),0);
        if (abs(weightSum-1.0) > .0001) 
		    throw ModelException("the sum of world weights is not equal to one", routine);
        */
		if (m_parallelMultShift.size()==0)
		    throw ModelException("worlds not defined", routine);
		if (m_notional.size()!=parSpreads->size())
		    throw ModelException("not the same number of notionals as number of spread curves", routine);
		if (m_notional.size()==0)
		    throw ModelException("No names defined", routine);
		if (diffParameters->size()!=4)
		    throw ModelException("Diffusion part badly defined", routine);
		if (jumpParameters->numRows()!=4)
		    throw ModelException("Jump part badly defined", routine);
	}
	catch (exception& e)
	{
		throw ModelException(e, routine);
	}

}


void SCIDparameters::getSingleNameCalibrationDates(DateTimeArray &dates)
{
	dates.resize(singleNameExpiries->size());
	for (int k=0; k<dates.size(); k++) dates[k] = (*singleNameExpiries)[k]->toDate(today);
}

void SCIDparameters::reInitializeSurvProba(DoubleMatrix &survProbaAtCalibrationSingleNameDates)
{
	static const string method = "SCIDparameters::reInitialize";
    try
	{
		DateTimeArray sgleNameDT;
		getSingleNameCalibrationDates(sgleNameDT);
		DateTime cfDate = cfExpiry->toDate(today);
		m_world.CalibrateDiffusion(
						   m_nbNames,
						   today,
						   sgleNameDT,
						   survProbaAtCalibrationSingleNameDates, 
						   cfDate,
						   (*diffParameters)[0],
						   (*diffParameters)[1],
						   (*diffParameters)[2],
						   (*diffParameters)[3]);
    }
    catch (exception& e)
	{
		throw ModelException(e, method);
	}
}



/**Check stuff here */



void SCIDparameters::load(CClassSP& clazz)
{
    clazz->setPublic(); 
    REGISTER(SCIDparameters, clazz);
    SUPERCLASS(MarketObject);
	EMPTY_SHELL_METHOD(defaultSCIDparameters);

//    FIELD(discount, "Discount curve");        
	FIELD(name, "Name of the scid parameter class");
    FIELD_MAKE_OPTIONAL(name); // default name is CIDParameters::DEFAULTNAME

	FIELD(parSpreads,  "Porfolio of par spreads");
	FIELD(m_notional,  "Notionals");
	FIELD(today,  "Today");
	FIELD_MAKE_TRANSIENT(today);
	FIELD(singleNameExpiries,  "Dates at which the base CID worlds calibrate to single name survival probabilities");

	FIELD(diffParameters,  "Diffusion parameters: Idio Vol, Common Factor Vol, Decay, Ratio");
	FIELD(cfExpiry,  "Diffusion parameters: Expiry used for Diffusion Survival Probability Ratio");
	FIELD(jumpParameters,  "Matrix of size nbMarkets times 4, 4 for frequency, jumpSize, impact, decay");

	FIELD(m_parallelMultShift,  "parallel shifts of the base World");
    FIELD(m_linearShift,  "linear shifts of the base World");
	FIELD(m_initialWeights,  "weights of the different worlds");
    FIELD_MAKE_OPTIONAL(m_initialWeights); // if not provided, size 1 and then set to 1/nbWorlds
	FIELD(baseWorldCalibSurvProba,  "Survival Probabilities for the base CID world. Coming from Calibration");
    FIELD_MAKE_OPTIONAL(baseWorldCalibSurvProba); // if not provided, size 1 and then set to 1/nbWorlds
	

    FIELD(m_nbNames, "");
    FIELD_MAKE_TRANSIENT(m_nbNames);
    FIELD(m_nbWorlds, "");
    FIELD_MAKE_TRANSIENT(m_nbWorlds);
    FIELD(m_nbJumpProcesses, "");
    FIELD_MAKE_TRANSIENT(m_nbJumpProcesses);
    FIELD(m_recovery, "");
    FIELD_MAKE_TRANSIENT(m_recovery);


}

IObject* SCIDparameters::defaultSCIDparameters(){
    return new SCIDparameters();
}

CClassConstSP const SCIDparameters::TYPE = CClass::registerClassLoadMethod(
   "SCIDparameters", typeid(SCIDparameters), load);



ICDSParSpreadsSP SCIDparameters::getCDSparSpread(int name)
{
	QLIB_VERIFY(name>=0 && name<m_nbNames, "name index out of range");
	return (*parSpreads)[name].getSP();
}

double SCIDparameters::getMarketPortfolioEL(DateTime date, bool notionalEL)
{
	double EL = 0;
	for (int m=0; m<m_nbNames; m++)
	{
		ICDSParSpreadsSP CDS = getCDSparSpread(m);
		if (!notionalEL) EL += (1-m_recovery[m])*m_notional[m]*(1-CDS->survivalProb(date));
		else EL += m_notional[m]*(1-CDS->survivalProb(date));
	}
	return EL;
}

double SCIDparameters::getMarketPortfolioEL(double date, bool notionalEL)
{
	DateTime _date = today.rollDate(int(date*365+0.5));
	return getMarketPortfolioEL(_date,notionalEL);
}

DoubleArraySP SCIDparameters::getMarketPortfolioEL(DateTimeArray dates, bool notionalEL)
{
    DoubleArraySP result(new DoubleArray(dates.size()));
    for (int i=0; i<dates.size(); i++) 
        (*result)[i]=getMarketPortfolioEL(dates[i],notionalEL);
    return result; 
}

void SCIDparameters::MarketSurvProb(DoubleArray &time, DoubleMatrix &survProbability)
{
	DateTimeArray dates(time.size());
	for (int i=0; i<time.size(); i++) dates[i] = today.rollDate(int(time[i]*365+0.5));
	MarketSurvProb(dates, survProbability);
}

void SCIDparameters::MarketSurvProb(DateTimeArray &dates, DoubleMatrix &survProbability)
{
	survProbability.resize(dates.size(), m_nbNames);
	for (int m=0; m<m_nbNames; m++)
	{
		ICDSParSpreadsSP CDS = getCDSparSpread(m);
		for (int k=0; k<dates.size(); k++)
			survProbability[k][m] = CDS->survivalProb(dates[k]);
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////RISKY ANNUITIES AND DEFAULT LEGS////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

void SCIDparameters::computeRiskyAnnuityLeg(DateTime& valueDate, YieldCurveSP discount, DateTimeArray& simDates,
                                DoubleArray& expectedLoss, double trancheDiv, DateTime& startDate, DateTimeArray& endDates,
                                int coupon, const DayCountConvention & dcc, DoubleArray & RAvalues)
{
    DateTimeArray curveDates(1, valueDate);
    curveDates.insert(curveDates.end(), simDates.begin(), simDates.end());
	DoubleArray riskyDiscount(1,1);
   
    //computes risky discounts
    for(int i=0; i < expectedLoss.size(); ++i) 
        riskyDiscount.push_back(Maths::max(1e-12, 1-expectedLoss[i]*trancheDiv));           
    EffectiveCurve trancheCurve(valueDate, discount, curveDates, riskyDiscount, EffectiveCurve::FLAT_FORWARD);
    long matSize = endDates.size();
    RAvalues.resize(matSize);
	double RAT, RAt;
    for (long i=0; i < matSize; ++i){	
        CashFlowArray cashflows = SwapTool::cashflows(valueDate, endDates[i], false, 1.0, coupon, "M", &dcc);
        cashflows[cashflows.size()-1].amount -= 1.0; // ugly
        RAT = trancheCurve.annuityPV(cashflows,valueDate,IDiscountCurveRisky::RECOVER_1);
        if (startDate == valueDate)
            RAvalues[i] = RAT;
        else{    
            cashflows = SwapTool::cashflows(valueDate, startDate, false, 1.0, coupon, "M", &dcc);
            cashflows[cashflows.size()-1].amount -= 1.0; // ugly
            RAt = trancheCurve.annuityPV(cashflows,valueDate,IDiscountCurveRisky::RECOVER_1);
            RAvalues[i] = RAT - RAt;        
        }
    }
}

void SCIDparameters::computeDefaultLeg(DateTime& valueDate, YieldCurveSP discount, DateTimeArray& simDates,
                                DoubleArray& expectedLoss, double trancheDiv, DateTime& startDate, DateTimeArray& endDates,
                                DoubleArray & DLvalues)
{
    DateTimeArray curveDates(1, valueDate);
    curveDates.insert(curveDates.end(), simDates.begin(), simDates.end());
	  DoubleArray riskyDiscount(1,1);
   
    //computes risky discounts
    for(int i=0; i < expectedLoss.size(); ++i) 
        riskyDiscount.push_back(Maths::max(1e-12, 1-expectedLoss[i]*trancheDiv));           
    EffectiveCurve trancheCurve(valueDate, discount, curveDates, riskyDiscount, EffectiveCurve::FLAT_FORWARD);
    long matSize = endDates.size();
	  double DLT, DLt;
    DLvalues.resize(matSize);
    for (long i=0; i < matSize; ++i){	
        if (startDate.equals(valueDate))
            DLvalues[i] = trancheCurve.protectionPV(valueDate,valueDate, endDates[i],IDiscountCurveRisky::RECOVER_1);
        else{
            DLT = trancheCurve.protectionPV(valueDate,valueDate, endDates[i],IDiscountCurveRisky::RECOVER_1);
            DLt = trancheCurve.protectionPV(valueDate,valueDate,startDate,IDiscountCurveRisky::RECOVER_1);
            DLvalues[i] = DLT-DLt;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////DATA STORAGE FOR FULL MONTE CARLO////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SCIDparameters::InitializeDataStorage(int timeSlices, int nrPaths, int nrStoredValues)
{
    dataStorage.resize(timeSlices);
    for (int i=0; i < dataStorage.size(); ++i)
        dataStorage[i] = CDoubleMatrixSP(new DoubleMatrix(timeSlices, nrPaths));
}
void SCIDparameters::StorePathValues(int timeSlice, int nrPath, DoubleArray & values)
{
    QLIB_VERIFY((*dataStorage[timeSlice]).numCols() == values.size(),"StorePathValues tries to store a vector data different than what it should be");
    QLIB_VERIFY((*dataStorage[timeSlice]).numRows() < nrPath,"StorePathValues tries to store a vector data for a path index bigger than initial number of paths");
    for(int i=0; i< values.size(); ++i)
        (*dataStorage[timeSlice])[i][nrPath] = values[i];
}

DoubleMatrixArray SCIDparameters::getStoredValues()
{
    return dataStorage;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////TIME LINE HELP /////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

void SCIDparameters::push_backTimeLine(DateTimeArray &timeLine,
		  				   DateTime firstDate, 
						   DateTime lastDate, 
						   MaturityPeriodSP period, 
						   bool includeFirstDate)
{
//	period.toYears();
	if (includeFirstDate) timeLine.push_back(firstDate);
	DateTime nextDate = period->toDate(firstDate);
	while (nextDate<lastDate)
	{
		timeLine.push_back(nextDate);
		nextDate = period->toDate(nextDate);
	}
	if ((timeLine.back()==firstDate)&&(firstDate<lastDate)) timeLine.push_back(lastDate);
	else
	{
		if (timeLine.back().yearFrac(lastDate)<0.2*period->toYears()) timeLine.back()=lastDate;
			else timeLine.push_back(lastDate);
	}
}


DoubleArray SCIDparameters::DateAsDouble(DateTimeArray &dates)
{
	DoubleArray t(dates.size());
	for (int i=0; i<t.size(); i++)
		t[i] = today.yearFrac(dates[i]);
	return t;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////CLOSED FORM/////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

void SCIDparameters::NamesSurvProbinGivenWorld(int world, DoubleArray &time, DoubleArray &survProbability)
			// survProbabilility(name,time) = survProbabilility(name*time.size()+index time)
{
	QLIB_VERIFY(world>=0 && world<m_nbWorlds, "world index out of range");

	survProbability.resize(time.size()* m_nbNames);
	m_world.survProba(true, true, true, m_linearShift[world], m_parallelMultShift[world], time, &survProbability[0]);
}

void SCIDparameters::NamesSurvProbinGivenWorld(int world, DateTimeArray &dates, DoubleArray &survProbability)
			// survProbabilility(name,time) = survProbabilility(name*time.size()+index time)
{
	array<double> t(dates.size());
	for (int i=0; i<t.size(); i++) t[i] = today.yearFrac(dates[i]);
	NamesSurvProbinGivenWorld(world,t,survProbability);
}

void SCIDparameters::NamesSurvProb(DoubleArray &time, DoubleMatrix &survProbability)
{
	survProbability.resize(time.size(),m_nbNames);
	survProbability.fill(0.0);
	DoubleArray survProba(time.size()*m_nbNames);
	for (int i=0; i<m_nbWorlds; i++)
	{
		m_world.survProba(true, true, true, m_linearShift[i], m_parallelMultShift[i], time, &survProba[0]);
		for (int j=0; j<time.size(); j++)
			for (int k=0; k<m_nbNames; k++)
				survProbability[j][k] += m_initialWeights[i]*survProba[k*time.size()+j];
	}
}

void SCIDparameters::NamesSurvProb(DateTimeArray &dates, DoubleMatrix &survProbability)
{
	array<double> t(dates.size());
	for (int i=0; i<t.size(); i++) t[i] = today.yearFrac(dates[i]);
	NamesSurvProb(t,survProbability);
}



void SCIDparameters::PortfolioEL(DoubleArray &times, double *pel)  // pel must be of of size times.size()
{
	DoubleArray survProba;
	for (int k=0; k<times.size(); k++) pel[k]=0;
	for (int i=0; i<m_nbWorlds; i++)
	{
		NamesSurvProbinGivenWorld(i, times, survProba);
		for (int k=0; k<times.size(); k++)
			for (int m=0; m<m_nbNames; m++)
				pel[k]+= (1-survProba[m*times.size()+k])*m_notional[m]*(1-m_recovery[m])*m_initialWeights[i];
	}
}

void SCIDparameters::PortfolioEL(DateTimeArray &dates, double *pel)  // pel must be of of size dates.size()
{
	array<double> t(dates.size());
	for (int i=0; i<t.size(); i++) t[i] = today.yearFrac(dates[i]);
	PortfolioEL(t,pel);
}


void SCIDparameters::PortfolioELinGivenWorld(int world, DoubleArray &times, double *pel)  // pel of size times.size().
{
	QLIB_VERIFY(world>=0 && world<m_nbWorlds, "world index out of range");

	DoubleArray survProba;
	for (int k=0; k<times.size(); k++) pel[k]=0;
	NamesSurvProbinGivenWorld(world, times, survProba);
	for (int k=0; k<times.size(); k++)
		for (int m=0; m<m_nbNames; m++)
			pel[k]+= (1-survProba[m*times.size()+k])*m_notional[m]*(1-m_recovery[m]);
}

void SCIDparameters::PortfolioELinGivenWorld(int world, DateTimeArray &dates, double *pel) // pel must be of of size dates.size()
{
	array<double> t(dates.size());
	for (int i=0; i<t.size(); i++) t[i] = today.yearFrac(dates[i]);
	PortfolioELinGivenWorld(world,t,pel);

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////FAST MC/////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

void SCIDparameters::setFastMC(int seed, double timeSteps, int nbSteps,	int noJumpPaths, int jumpsPaths)
{
	QLIB_VERIFY(timeSteps>0, "negative time steps");
	QLIB_VERIFY(nbSteps>=0,  "negative number of steps");

	m_lastTimeFast = nbSteps*timeSteps;
	if (m_world.HasCommonFactor)
	{
		DoubleArray timeSimulation;		  // discretization time for BM
		timeSimulation.resize(nbSteps+1);
		timeSimulation[0]=0;
		for (int k=0; k<nbSteps; k++)
			timeSimulation[k+1]=timeSimulation[k]+timeSteps;
		m_world.setTimes(timeSimulation);
		m_BMFast.setParameters(&timeSimulation[0], timeSimulation.size(), seed);
		m_BMincrementsFast.resize(nbSteps);
		nbPathsNoJump = noJumpPaths;
	}
	else nbPathsNoJump = 1;

	if (m_nbJumpProcesses>0)
	{
		m_nbJumpsFast.resize(m_nbJumpProcesses);
		if (jumpsPaths>0) nbPathsJumps = jumpsPaths;  
			else nbPathsJumps = 1;
		DoubleArray freq(m_nbJumpProcesses);
		for (int r=0; r<m_nbJumpProcesses; r++) freq[r] = m_world.m_jumps[r].freq;
		m_PoissonFast.setParameters(&freq[0], freq.size(), m_lastTimeFast, seed+1);
		probaZeroJump = m_PoissonFast.probaNoJump();
		m_jumpTimeFast.initializeObject(seed+2,1);
	}
	else probaZeroJump=1.0;
	resetIpathFastMC();
}

void SCIDparameters::resetIpathFastMC(int iPathBM, int iPathPoisson, int iPathJumpTime)
{
	if (m_world.HasCommonFactor) m_BMFast.setIpath(iPathBM);
	if (m_nbJumpProcesses>0)
	{
		m_PoissonFast.setIpath(iPathPoisson);
		m_jumpTimeFast.setIpath(iPathJumpTime);
	}
}

void SCIDparameters::getIpathFastMC(int &iPathBM, int &iPathPoisson, int &iPathJumpTime)
{
	if (m_world.HasCommonFactor) iPathBM = m_BMFast.getIpath();
	if (m_nbJumpProcesses>0)
	{
		iPathPoisson = m_PoissonFast.getIpath();
		iPathJumpTime = m_jumpTimeFast.getIpath();
	}
}

void SCIDparameters::setConvolution(DoubleArray &kmin,
									DoubleArray &kmax)
{
	convol.setParameters(m_nbNames,kmin,kmax);
	convol.setLGD(&m_recovery[0],&m_notional[0]);
}

// Compute, as seen as today, forward TEL between times[0] and times...... put times[0]=0 to get TEL
// attachment points have been set up in setConvolution
// Monte-Carlo parameters have been set up in setFastMC
// convolution parameters are 

void SCIDparameters::ComputeTELinAllWorlds( 
					    DoubleArray &times,  // times at which we compute these TEL,
						vector< DoubleMatrix > &TEL,  // size nbWorlds, DoubleMatrix of size times.size(), kmin.size().
						int convolutionNoJump, int convolutionJumps)
{
	// ALGORITHM TO BE CHANGED, NEED MORE CACHING..
	// RIGHT NOW, SIMULATING EXACTLY THE SAME RANDOM PATHS FOR EACH WORLDS.
	int i,k,l,m,mc;
	double firstJump;
	double mult1 = 1.0 / double(nbPathsNoJump), mult2=0;
	if (nbPathsJumps>0) mult2 = 1.0 / double(nbPathsJumps);

	DoubleMatrix work(times.size(), convol.getNbTranches());
	work.fill(0.0);  
	if (TEL.size()!=m_nbWorlds) TEL.resize(m_nbWorlds);

	for (i=0; i<m_nbWorlds; i++)
	{
		work.fill(0);
		resetIpathFastMC();
		m_world.InitializeCondSurvProb(m_linearShift[i], m_parallelMultShift[i], times);
		for (mc=0; mc<nbPathsNoJump; mc++)
		{
			if (m_world.HasCommonFactor)
			{
				m_BMFast.GetBMincrements(&m_BMincrementsFast[0]);
				m_world.getCondSurvProba(m_linearShift[i], m_parallelMultShift[i], times, &m_BMincrementsFast[0]);
			}
			else 
				m_world.getCondSurvProba(m_linearShift[i], m_parallelMultShift[i], times, 0);
			for (k=1; k<times.size(); k++)
			{
				for (m=0; m<m_nbNames; m++)
					survProb[m] = m_world.m_condSurvProb[m*times.size()+k]
							+ (1.0 - m_world.m_condSurvProb[m*times.size()]); // this line is 0 for classical tranches
				convol.computeTEL(convolutionNoJump, &survProb[0], &work[k][0]); 
			}
		}
		work.scale(mult1);
		TEL[i].resize(times.size(), convol.getNbTranches());
		TEL[i].fill(0);

		if (m_nbJumpProcesses>0)
		{
			for (mc=0; mc<nbPathsJumps; mc++)
			{
				m_PoissonFast.GetPoissonNbJumps(true,&m_nbJumpsFast[0],(0.5+mc)*mult2);
				if (m_world.HasCommonFactor)
				{
					m_BMFast.GetBMincrements(&m_BMincrementsFast[0]);
					firstJump = m_world.getCondSurvProba(m_linearShift[i], m_parallelMultShift[i], times, 
													&m_BMincrementsFast[0], &m_nbJumpsFast[0],&m_jumpTimeFast, m_lastTimeFast);
				}
				else 
				{
					firstJump = m_world.getCondSurvProba(m_linearShift[i], m_parallelMultShift[i], times, 
													0, &m_nbJumpsFast[0],&m_jumpTimeFast, m_lastTimeFast);
				}
				for (k=1; k<times.size(); k++) 
				{
					if (firstJump>=times[k])
						for (l=0; l<work.numRows(); l++) TEL[i][k][l] += work[k][l];
					else
					{
						for (m=0; m<m_nbNames; m++)
							survProb[m] = m_world.m_condSurvProb[m*times.size()+k]
									+ ( 1 - m_world.m_condSurvProb[m*times.size()]); // this line is 0 for classical tranches
						convol.computeTEL(convolutionJumps, &survProb[0], &TEL[i][k][0]); 
					}
				}
			}
			TEL[i].scale(mult2*(1-probaZeroJump));
		}
		TEL[i].add(work, probaZeroJump);
	}
}

void SCIDparameters::ComputeSurvProbinAllWorld_fastMCtest(DoubleArray &times,     // compute the survival probability of each names, using the fast MC methodology
											  vector< DoubleMatrix > &survProba)  // only used for testing
{
	int i,k,m, mc;
	double firstJump;
	double mult1 = 1.0 / double(nbPathsNoJump), mult2=0;
	if (m_nbJumpProcesses>0) mult2 = 1.0 / double(nbPathsJumps);

	DoubleMatrix work(times.size(), m_nbNames);
	survProba.resize(m_nbWorlds);
	work.fill(0);  

	for (i=0; i<m_nbWorlds; i++)
	{
		work.fill(0);
		resetIpathFastMC();
		m_world.InitializeCondSurvProb(m_linearShift[i], m_parallelMultShift[i], times);
		for (mc=0; mc<nbPathsNoJump; mc++)
		{
			if (m_world.HasCommonFactor)
			{
				m_BMFast.GetBMincrements(&m_BMincrementsFast[0]);
				m_world.getCondSurvProba(m_linearShift[i], m_parallelMultShift[i], times, &m_BMincrementsFast[0]);
			}
			else 
				m_world.getCondSurvProba(m_linearShift[i], m_parallelMultShift[i], times, 0);
			for (k=0; k<times.size(); k++)
				for (m=0; m<m_nbNames; m++)
					work[k][m] += m_world.m_condSurvProb[m*times.size()+k];
		}
		work.scale(mult1);
		survProba[i].resize(times.size(), m_nbNames);
		survProba[i].fill(0);

		if (m_nbJumpProcesses>0)
		{
			for (mc=0; mc<nbPathsJumps; mc++)
			{
				m_PoissonFast.GetPoissonNbJumps(true,&m_nbJumpsFast[0],(0.5+mc)*mult2);
				if (m_world.HasCommonFactor)
				{
					m_BMFast.GetBMincrements(&m_BMincrementsFast[0]);
					firstJump = m_world.getCondSurvProba(m_linearShift[i], m_parallelMultShift[i], times, 
													&m_BMincrementsFast[0], &m_nbJumpsFast[0],&m_jumpTimeFast, m_lastTimeFast);
				}
				else 
				{
					firstJump = m_world.getCondSurvProba(m_linearShift[i], m_parallelMultShift[i], times, 
													0, &m_nbJumpsFast[0],&m_jumpTimeFast, m_lastTimeFast);
				}
				for (k=0; k<times.size(); k++) 
				{
					if (firstJump>=times[k])
						for (m=0; m<m_nbNames; m++)
							survProba[i][k][m] += work[k][m];
					else
						for (m=0; m<m_nbNames; m++)
							survProba[i][k][m] += m_world.m_condSurvProb[m*times.size()+k];
				}
			}
			survProba[i].scale(mult2*(1-probaZeroJump));
		}
		survProba[i].add(work, probaZeroJump);
	}
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////FULL MC/////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

void SCIDparameters::setFullMC(int seed, DoubleArray &discTimes)
{
    QLIB_VERIFY(discTimes.size()>0, "no dates in the Monte-Carlo simulation");
	for (int k=0; k<m_discTimes.size()-1; k++)
		QLIB_VERIFY(discTimes[k+1]>discTimes[k], "Dates in Full MC not in increasing order.");

	m_discTimes = discTimes;
	m_lastDateFull = m_discTimes.back();
	if (m_world.HasCommonFactor)
	{
		m_expDecayCF.resize(m_discTimes.size()-1);
		for (int k=0; k<m_discTimes.size()-1; k++) 
			m_expDecayCF[k]=exp(-m_world.m_cfCIR.kappa*0.5*(m_discTimes[k+1]-m_discTimes[k]));
	}

	m_jumpIntensity.resize(m_nbJumpProcesses,m_nbNames*m_discTimes.size());
	m_nbJumpsFull.resize(m_nbJumpProcesses);
	m_weightsDynamic.resize(m_discTimes.size(), m_nbWorlds);
	m_defaultTime.resize(m_nbNames);
	m_loss.resize(m_discTimes.size());
	m_notionalLoss.resize(m_discTimes.size());
	m_nbDefaultedNames.resize(m_discTimes.size());

	m_lambdaFutureTimeJump.resize(m_nbJumpProcesses,m_nbNames);
	m_lambdaFutureTimeIdio.resize(m_nbNames);
	m_intensity.resize(m_nbWorlds);
	m_defaultTime.resize(m_nbNames);
	m_forUpdate.resize(m_nbWorlds);
	m_defaultFull.initializeObject(seed+4,m_nbNames);

	m_intensity.resize(m_nbWorlds); 
	for (int i=0; i<m_nbWorlds; i++) m_intensity[i].initialize(m_nbNames, m_nbJumpProcesses, m_discTimes.size());

	if (m_nbJumpProcesses>0)
	{
		DoubleArray freq(m_nbJumpProcesses);
		for (int r=0; r<freq.size(); r++) 
			freq[r] = m_world.m_jumps[r].freq;
		m_poissonFull.setParameters(&freq[0], m_nbJumpProcesses, m_lastDateFull, seed+1);
		m_jumpTimeFull.initializeObject(seed+2,1);
		m_jumpSizeFull.initializeObject(seed+3,m_nbNames);
	}
	m_nbJumpsFull.resize(m_nbJumpProcesses);
	if ((m_world.HasCommonFactor)||(m_world.HasIdioVol)) 
		m_bmFull.setParameters(&m_discTimes[0], m_discTimes.size(), seed);
	int n=0;
	if (m_world.HasIdioVol) n = m_nbNames;
	if (m_world.HasCommonFactor) n++;
	m_bmIncrementsFull.resize(m_discTimes.size()*n,0.0); // will be used to store the Brownian Motions.
	m_intIntensityTemp.resize(m_nbNames*m_discTimes.size());
	m_world.m_idioCIR.setThetaDiff(m_discTimes);
	resetSpreads();
}

void SCIDparameters::setFullMC(int seed, DateTimeArray &discDates)
{
	array<double> t(discDates.size());
	for (int i=0; i<t.size(); i++) t[i] = today.yearFrac(discDates[i]);
	setFullMC(seed,t);
}


void SCIDparameters::resetIpathFullMC(int iPathBM, int iPathPoisson, int iPathJumpTime, int iPathJumpSize, int iPathDeafult)
{
	m_bmFull.setIpath(iPathBM);
	m_poissonFull.setIpath(iPathPoisson);
	m_jumpTimeFull.setIpath(iPathJumpTime);
	m_jumpSizeFull.setIpath(iPathJumpSize);
	m_defaultFull.setIpath(iPathDeafult);
}

void SCIDparameters::resetSpreads()
{
	int i,k,m;
	for (i=0; i<m_nbWorlds; i++)
	{
		fill(m_intensity[i].idio.begin(),m_intensity[i].idio.end(),0);
		fill(m_intensity[i].cf.begin(),m_intensity[i].cf.end(),0);
		fill(m_intensity[i].baseTotal.begin(),m_intensity[i].baseTotal.end(),0);
		fill(m_intensity[i].baseIntegral.begin(),m_intensity[i].baseIntegral.end(),0);
		fill(m_intensity[i].total.begin(),m_intensity[i].total.end(),0);
		fill(m_intensity[i].integral.begin(),m_intensity[i].integral.end(),0);
		m_jumpIntensity.fill(0);
	}
	if (!m_world.HasIdioVol)
	{
		for (i=0; i<m_nbWorlds; i++)
		{
			m_world.m_idioCIR.diffuseNoVol(m_linearShift[i],m_parallelMultShift[i],&m_intensity[i].idio[0]);
			m_intensity[i].baseTotal = m_intensity[i].idio;
			double df;
			for (m=0; m<m_nbNames; m++)
			{
				m_intensity[i].baseIntegral[m*m_discTimes.size()]=0;
				for (k=1; k<m_discTimes.size(); k++)
				{
					df = 0.5*(m_intensity[i].idio[m*m_discTimes.size()+k-1] + m_intensity[i].idio[m*m_discTimes.size()+k]);
				    df *= (m_discTimes[k]-m_discTimes[k-1]);
					m_intensity[i].baseIntegral[m*m_discTimes.size()+k] = m_intensity[i].baseIntegral[m*m_discTimes.size()+k-1] + df;
				}
			}
		}
	} 
}

void SCIDparameters::SimulateSpread()
{
	int i,k,m,r;
	double df;
	if (m_world.HasIdioVol)
	{
		for (m=0; m<m_nbNames; m++) m_bmFull.GetBMincrements(&m_bmIncrementsFull[m*(m_discTimes.size()-1)]);
		for (i=0; i<m_nbWorlds; i++)
		{
			m_world.m_idioCIR.diffuse(m_linearShift[i],m_parallelMultShift[i],&m_bmIncrementsFull[0],&m_intensity[i].idio[0]);
			for (m=0; m<m_nbNames; m++)
			{
				m_intensity[i].integral[m*m_discTimes.size()]=0;
				for (k=1; k<m_discTimes.size(); k++)
				{
					df = 0.5*(m_intensity[i].idio[m*m_discTimes.size()+k-1] + m_intensity[i].idio[m*m_discTimes.size()+k]);
				    df *= (m_discTimes[k]-m_discTimes[k-1]);
					m_intensity[i].integral[m*m_discTimes.size()+k] = m_intensity[i].integral[m*m_discTimes.size()+k-1] + df;
				}
			}
			m_intensity[i].total = m_intensity[i].idio;
		}
	}
	else
	{
		for (i=0; i<m_nbWorlds; i++) 
		{
			m_intensity[i].total = m_intensity[i].baseTotal;
			m_intensity[i].integral = m_intensity[i].baseIntegral;
		}
	}

	// Common Factor
	if (m_world.HasCommonFactor)
	{
		double integrale;
		m_bmFull.GetBMincrements(&m_bmIncrementsFull[0]);
		for (i=0; i<m_nbWorlds; i++)
		{
			integrale = 0;
			m_world.m_cfCIR.diffuse(m_linearShift[i],m_parallelMultShift[i],
				&m_bmIncrementsFull[0],&m_expDecayCF[0],m_discTimes,&m_intensity[i].cf[0]);
			for (k=1; k<m_discTimes.size(); k++)
			{
				integrale += 0.5*(m_intensity[i].cf[k-1] + m_intensity[i].cf[k])*(m_discTimes[k]-m_discTimes[k-1]);
				for (m=0; m<m_nbNames; m++)
					m_intensity[i].integral[m*m_discTimes.size()+k] += integrale*m_world.m_ams[m];
			}
			for (k=0; k<m_discTimes.size(); k++)
				for (m=0; m<m_nbNames; m++)
					m_intensity[i].total[m*m_discTimes.size()+k] += m_intensity[i].cf[k]*m_world.m_ams[m];
		}
	};

	for (r=0; r<m_nbJumpProcesses; r++)
	{
		if (m_nbJumpsFull[r]!=0)
		{
			m_world.m_jumps[r].condSpreadAndHazard(
				m_discTimes,
			    m_nbJumpsFull[r],
				m_jumpTimeFull,
				m_lastDateFull,
				m_jumpSizeFull,       // gives uniform in R^m_nbNames;
				&m_jumpIntensity[r][0],    // (O) (name,index time) -> (index time)*m_nbNames+name
				&m_intIntensityTemp[0]);   // (O) (name,index time) -> (index time)*m_nbNames+name
			for (i=0; i<m_nbWorlds; i++)
				for (k=0; k<m_intIntensityTemp.size(); k++)
				{
					m_intensity[i].integral[k] += m_linearShift[i]*m_intIntensityTemp[k];
					m_intensity[i].total[k] += m_linearShift[i]*m_jumpIntensity[r][k];
				}
		}
		else m_jumpIntensity.fill(0);
	}
}

void SCIDparameters::SimulateNbJumps(bool conditional)
{
	if (m_nbJumpProcesses>0) 
		m_poissonFull.GetPoissonNbJumps(conditional,&m_nbJumpsFull[0]);
}

void SCIDparameters::SimulateDefaults()
{
	// we assume that the spread have already been simulated
	int i,k,m;
	double survProba, OneDivMeanForUpdate;
	for (i=0; i<m_nbWorlds; i++) m_weightsDynamic[0][i] = m_initialWeights[i];
	
	fill(m_defaultTime.begin(),m_defaultTime.end(),999); // no defaults yet
	m_loss[0]=0.0;			 // no defaults yet
	m_notionalLoss[0]=0.0;       // no defaults yet
	m_nbDefaultedNames[0]=0;

	for (k=0; k<m_discTimes.size()-1; k++)
	{
		fill(m_forUpdate.begin(),m_forUpdate.end(),1.0);
		for (m=0; m<m_nbNames; m++)
		{
			if (m_defaultTime[m]>=999)
			{
				survProb[m] = 0.0;
				for(i=0; i<m_nbWorlds; i++)
				{
					survProba = exp ( - m_intensity[i].integral[m*m_discTimes.size()+k+1]+m_intensity[i].integral[m*m_discTimes.size()+k]);
				    m_forUpdate[i] *= survProba;
					survProb[m] += m_weightsDynamic[k][i] * survProba;
				}
			}
			else survProb[m] = 1.0;
		}
		m_loss[k+1] = m_loss[k];
		m_notionalLoss[k+1] = m_notionalLoss[k];
		m_nbDefaultedNames[k+1] = m_nbDefaultedNames[k];

		double * unif = m_defaultFull.getVector();
		for (m=0; m<m_nbNames; m++)
			if (unif[m]>survProb[m])
			{
				m_defaultTime[m] = m_discTimes[k] + (unif[m]-survProb[m])/(1.0-survProb[m])*(m_discTimes[k+1]-m_discTimes[k]);
				m_notionalLoss[k+1] += m_notional[m];
				++m_nbDefaultedNames[k+1];
				m_loss[k+1] += m_notional[m]*(1-m_recovery[m]);
				for (i=0; i<m_nbWorlds; i++)
					m_forUpdate[i] *= exp ( m_intensity[i].integral[m*m_discTimes.size()+k+1] - m_intensity[i].integral[m*m_discTimes.size()+k] ) - 1;
			}

		// update weights
		OneDivMeanForUpdate = 0;
		for(i=0; i<m_nbWorlds; i++) OneDivMeanForUpdate += m_weightsDynamic[k][i]*m_forUpdate[i];
		OneDivMeanForUpdate = 1/OneDivMeanForUpdate;
		for(i=0; i<m_nbWorlds; i++) 
			m_weightsDynamic[k+1][i] = m_weightsDynamic[k][i]*m_forUpdate[i]*OneDivMeanForUpdate;
	}
}

void SCIDparameters::getSurvProba(int indexLittle_t,
				  vector<int> &indexBig_t,                  
				  DoubleMatrix &survProba,
				  bool add)	  // return the survival probability between discTimes[indexLittle_t] and discTimes[indexBig_t],
							  // as computed in the given simulation (no average is being taken
							  // survProba must be of size indexBig_t.size() * m_nbNames, otherwise will be resized
							  // if add, add to survProba the answer
{
	survProba.resize(indexBig_t.size(),m_nbNames); // does not do anything if already of the right sizes...
	if (!add) survProba.fill(0);
	double logSP_t, logSP_T, w;
	for (int m=0; m<m_nbNames; m++)
		for (int i=0; i<m_nbWorlds; i++)
		{
			logSP_t = m_intensity[i].integral[m*m_discTimes.size()+indexLittle_t];
			w = m_weightsDynamic[indexLittle_t][i];
			for (size_t j=0; j<indexBig_t.size(); j++)
			{
				logSP_T = m_intensity[i].integral[m*m_discTimes.size()+indexBig_t[j]];
				survProba[j][m] += w*exp(-logSP_T + logSP_t);
			}
		}
}

DoubleArray SCIDparameters::getNameIntensity(int nameIndex)
{
    DoubleArray nameIntensities(m_discTimes.size(),0.0);
	for (int k=0; k<m_discTimes.size()-1; k++)
		for (int i=0; i<m_nbWorlds; i++)
            nameIntensities[k]+= m_weightsDynamic[k][i] * m_intensity[i].total[nameIndex*m_discTimes.size()+k];
    return nameIntensities;
}

DoubleArrayArray SCIDparameters::getIntensities()
{
    DoubleArrayArray intensities(m_nbNames, DoubleArray(m_discTimes.size(),0.0));
	for (int m=0; m<m_nbNames; m++)
        for (int k=0; k<m_discTimes.size()-1; k++)
            for (int i=0; i<m_nbWorlds; i++)
            intensities[m][k]+= m_weightsDynamic[k][i] * m_intensity[i].total[m*m_discTimes.size()+k];
    return intensities;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// HYBRID MONTE CARLO /////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////


double SCIDparameters::getFuturePortfolioEL(int indexFutureTime,
		  					    DoubleArray &times,
								DoubleArray &EL)
{
	QLIB_VERIFY(indexFutureTime>=0 && indexFutureTime<m_discTimes.size(), "indexFutureTime out of range");

	int i,k,m,r;
	DoubleArray survProba1world(m_nbNames*times.size());
	for (r=0; r<m_nbJumpProcesses; ++r)
		for (m=0; m<m_nbNames; ++m)
			m_lambdaFutureTimeJump[r][m] = m_jumpIntensity[r][m*m_discTimes.size()+indexFutureTime];
	for (i=0; i<m_nbWorlds; ++i)
	{
		for (m=0; m<m_nbNames; m++) 
			m_lambdaFutureTimeIdio[m] = m_intensity[i].idio[m*m_discTimes.size()+indexFutureTime];
		m_world.survProba(true,true,true,m_linearShift[i],m_parallelMultShift[i],times,&survProba1world[0], 
				m_discTimes[indexFutureTime],&m_lambdaFutureTimeIdio[0],m_intensity[i].cf[indexFutureTime],&m_lambdaFutureTimeJump);
		for (m=0; m<m_nbNames; m++)
			if (m_defaultTime[m]>m_discTimes[indexFutureTime])
			{
				for (k=0; k<times.size(); k++)
					EL[k] += m_weightsDynamic[indexFutureTime][i]*(1-survProba1world[m*times.size()+k])*(1-m_recovery[m])*m_notional[m];
			}
	}
	return m_loss[indexFutureTime];
}

double SCIDparameters::getFuturePortfolioNotionalEL(int indexFutureTime,
		  					    DoubleArray &times,
								DoubleArray &notionalEL)
{
	QLIB_VERIFY(indexFutureTime>=0 && indexFutureTime<m_discTimes.size(), "indexFutureTime out of range");

	notionalEL.resize(times.size());
	int i,k,m,r;
	DoubleArray survProba1world(m_nbNames*times.size());
	for (r=0; r<m_nbJumpProcesses; ++r)
		for (m=0; m<m_nbNames; ++m)
			m_lambdaFutureTimeJump[r][m] = m_jumpIntensity[r][m*m_discTimes.size()+indexFutureTime];
	for (i=0; i<m_nbWorlds; ++i)
	{
		for (m=0; m<m_nbNames; m++) 
			m_lambdaFutureTimeIdio[m] = m_intensity[i].idio[m*m_discTimes.size()+indexFutureTime];
		m_world.survProba(true,true,true,m_linearShift[i],m_parallelMultShift[i],times,&survProba1world[0], 
				m_discTimes[indexFutureTime],&m_lambdaFutureTimeIdio[0],m_intensity[i].cf[indexFutureTime],&m_lambdaFutureTimeJump);
		for (m=0; m<m_nbNames; m++)
			if (m_defaultTime[m]>m_discTimes[indexFutureTime])
			{
				for (k=0; k<times.size(); k++)
					notionalEL[k] += m_weightsDynamic[indexFutureTime][i]*(1-survProba1world[m*times.size()+k])*m_notional[m];
			}
	}
	return m_notionalLoss[indexFutureTime];
}

void SCIDparameters::getFutureSurvProba(int indexFutureTime, DoubleArray &times, DoubleMatrix &survProba)
{
	QLIB_VERIFY(indexFutureTime>=0 && indexFutureTime<m_discTimes.size(), "indexFutureTime out of range");

	DoubleArray survProba1world(m_nbNames*times.size());
	survProba.resize(times.size(),m_nbNames); // does nothing if of the right size
	survProba.fill(0);
	int i,k,m,r;
	for (m=0; m<m_nbNames; m++)
		for (r=0; r<m_nbJumpProcesses; r++)
			m_lambdaFutureTimeJump[r][m] = m_jumpIntensity[r][m*m_discTimes.size()+indexFutureTime];
	for (i=0; i<m_nbWorlds; i++)
	{
		for (m=0; m<m_nbNames; m++) 
			m_lambdaFutureTimeIdio[m] = m_intensity[i].idio[m*m_discTimes.size()+indexFutureTime];
		m_world.survProba(true,true,true,
								m_linearShift[i], 
								m_parallelMultShift[i], 
								times, 
								&survProba1world[0],
								m_discTimes[indexFutureTime], 
								&m_lambdaFutureTimeIdio[0], 
								m_intensity[i].cf[indexFutureTime],
								&m_lambdaFutureTimeJump);

		for (k=0; k<times.size(); k++)
				for (m=0; m<m_nbNames; m++)
					survProba[k][m] += m_weightsDynamic[indexFutureTime][i]*survProba1world[m*times.size()+k];//*exp(-m_intensity[i].integral[m*m_discTimes.size()+indexFutureTime]);
	}
}


double SCIDparameters::getFutureETL(int indexFutureTime, 
								DoubleArray &times, int begin,
								DoubleMatrix &ETL,   // ETL of size begin+times.size() * kmin.size()
								bool countPastLosses,
								int convolutionNoJump, int convolutionJumps)
{
	QLIB_VERIFY(indexFutureTime>=0 && indexFutureTime<m_discTimes.size(), "indexFutureTime out of range");

	ETL.resize(begin+times.size(),convol.getNbTranches()); // does nothing if of the right size

	if (countPastLosses) convol.setPast(m_loss[indexFutureTime],m_notionalLoss[indexFutureTime]);
					else convol.setPast(0.0,0.0);

	int i, k, l, m, mc, r;
	double firstJump;
    
	telWork1.resize(times.size(),convol.getNbTranches());   // make sure these arrays are of the right size 
	telWork2 = telWork1;								    // it does not cost anything if of the right size

	m_PoissonFast.changeMaturity(times.back() - m_discTimes[indexFutureTime]);
	double mult1 = 1.0 / double(nbPathsNoJump), mult2=0.0;
	if (m_nbJumpProcesses>0) mult2= 1.0 / double(nbPathsJumps);

	for (r=0; r<m_nbJumpProcesses; r++)
		for (m=0; m<m_nbNames; m++)
			m_lambdaFutureTimeJump[r][m] = m_jumpIntensity[r][m*m_discTimes.size()+indexFutureTime];

	int bmIpath, poissonIpath, jumpTimeIpath;
	getIpathFastMC(bmIpath, poissonIpath, jumpTimeIpath);

	for (i=0; i<m_nbWorlds; i++)
	{
		telWork1.fill(0);
		resetIpathFastMC(bmIpath, poissonIpath, jumpTimeIpath);
		for (m=0; m<m_nbNames; m++) 
			m_lambdaFutureTimeIdio[m] = m_intensity[i].idio[m*m_discTimes.size()+indexFutureTime];
		m_world.InitializeCondSurvProb(m_linearShift[i], m_parallelMultShift[i], times, m_discTimes[indexFutureTime], &m_lambdaFutureTimeIdio[0]);

		for (mc=0; mc<nbPathsNoJump; mc++)
		{
			if (m_world.HasCommonFactor)
			{
				m_BMFast.GetBMincrements(&m_BMincrementsFast[0]);
				m_world.getCondSurvProba(m_linearShift[i],m_parallelMultShift[i],times,&m_BMincrementsFast[0],0,0,0,
								  m_discTimes[indexFutureTime],m_intensity[i].cf[indexFutureTime],&m_lambdaFutureTimeJump);
			}
			else
				m_world.getCondSurvProba(m_linearShift[i],m_parallelMultShift[i],times,0,0,0,0,
								  m_discTimes[indexFutureTime],0,&m_lambdaFutureTimeJump);

			for (k=0; k<times.size(); k++)
			{
				for (m=0; m<m_nbNames; m++) survProb[m] = m_world.m_condSurvProb[m*times.size()+k];
				convol.computeTELfuture(convolutionNoJump, &survProb[0], &m_defaultTime[0], m_discTimes[indexFutureTime], &telWork1[k][0]);  
			}
		}
		telWork1.scale(mult1);

		if (m_nbJumpProcesses>0)
		{
			telWork2.fill(0);
			for (mc=0; mc<nbPathsJumps; mc++)
			{
				m_PoissonFast.GetPoissonNbJumps(true,&m_nbJumpsFast[0], (mc+0.5)*mult2);
				if (m_world.HasCommonFactor)
				{
					m_BMFast.GetBMincrements(&m_BMincrementsFast[0]);
					firstJump = m_world.getCondSurvProba(m_linearShift[i], 
											      m_parallelMultShift[i], 
												  times, 
												  &m_BMincrementsFast[0], 
												  &m_nbJumpsFast[0],
												  &m_jumpTimeFast, 
												  times.back(),
												  m_discTimes[indexFutureTime], 
												  m_intensity[i].cf[indexFutureTime],
												  &m_lambdaFutureTimeJump);
				}
				else firstJump = m_world.getCondSurvProba(m_linearShift[i], 
											      m_parallelMultShift[i], 
												  times, 
												  0, 
												  &m_nbJumpsFast[0],
												  &m_jumpTimeFast, 
												  times.back(),
												  m_discTimes[indexFutureTime], 
												  m_intensity[i].cf[indexFutureTime],
												  &m_lambdaFutureTimeJump);
				for (k=0; k<times.size(); k++) 
				{
					if (firstJump>=times[k])
					{
						for (l=0; l<convol.getNbTranches(); l++) telWork2[k][l] += telWork1[k][l];
					}
					else
					{
						for (m=0; m<m_nbNames; m++) survProb[m] = m_world.m_condSurvProb[m*times.size()+k];
						convol.computeTELfuture(convolutionJumps, &survProb[0], &m_defaultTime[0], m_discTimes[indexFutureTime], &telWork2[k][0]);  
					}
				}
			}
		}	
		if (m_nbJumpProcesses>0)
			for (k=0; k<times.size(); k++) 
				for (l=0; l<convol.getNbTranches(); l++)
					ETL[begin+k][l] += m_weightsDynamic[indexFutureTime][i]*telWork2[k][l]*mult2*(1-probaZeroJump);

		for (k=0; k<times.size(); k++) 
			for (l=0; l<convol.getNbTranches(); l++)
				ETL[begin+k][l] += m_weightsDynamic[indexFutureTime][i]*telWork1[k][l]*probaZeroJump;
        
	}
	return m_loss[indexFutureTime];
}



///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////





DEFINE_TEMPLATE_TYPE(SCIDparametersWrapper);
DRLIB_END_NAMESPACE

