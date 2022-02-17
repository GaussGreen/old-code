//----------------------------------------------------------------------------
//
//   Group       : Credit Derivatives Research
//
//   Filename    : SCIDtreeParameters.cpp
//
//   Date        : July, 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SCIDtreeParameters.hpp"
#include "edginc/ICDSParSpreads.hpp"
#include "edginc/CDSParSpreads.hpp"
#include "edginc/SwapTool.hpp"

DRLIB_BEGIN_NAMESPACE

/* Make sure the class links */
bool  SCIDtreeParametersLoad() {
    return (SCIDtreeParameters::TYPE != 0);
   }

SCIDtreeParameters::SCIDtreeParameters(const CClassConstSP& clazz) : MarketObject(clazz)
{
}

SCIDtreeParameters::SCIDtreeParameters()
	: MarketObject(SCIDtreeParameters::TYPE),
	  name("sCID parameters Default Name")
{
}
/*
string SCIDtreeParameters::discountYieldCurveName() const {
    return discount.getName();
}
*/

SCIDtreeParameters::~SCIDtreeParameters()
{
}

string SCIDtreeParameters::getName() const
{
    return name;
}


/** populate from market cache - default implementation provided */
void SCIDtreeParameters::getMarket(const IModel* model, const MarketData* market) 
{
    static const string method = "SCIDtreeParameters::getMarket";
	try
	{
		m_nbNames = parSpreads->size();
		double sumNotional = 0;
		for (int m=0; m<m_nbNames; m++) sumNotional += m_notional[m];
		for (int m=0; m<m_nbNames; m++) m_notional[m] /= sumNotional ;

		DoubleMatrix survProbaSingleNames(m_nbNames, singleNameExpiries->size());
		m_recovery.resize(m_nbNames);
		DateTimeArray sgleNameDT;
		getSingleNameCalibrationDates(sgleNameDT);
		DoubleArray sgleNameTimes = DateAsDouble(sgleNameDT);

		for (int m=0; m<parSpreads->size(); m++)
		{
		   (*parSpreads)[m].getData(model, market);
		   ICDSParSpreadsSP CDSm = (*parSpreads)[m].getSP();
		   m_recovery[m]= CDSm->getRecovery();
		   for (int k=0; k<sgleNameDT.size(); k++)
			   survProbaSingleNames[m][k] = CDSm->survivalProb(sgleNameDT[k]);
		}
		m_world.setUnderlyingProcesses(m_nbNames, volDiff, decayDiff, jumpDecay, jumpImpact, jumpFreq, jumpRatios, jumpMeans, jumpWeights);
		m_world.CalibrateSingleNames(today, sgleNameDT, survProbaSingleNames, idioRatio, cfRatio, jumpRatios);
    }
    catch (exception& e)
	{
		throw ModelException(e, method);
	}
}

/**Check stuff here */
void SCIDtreeParameters::validatePop2Object() 
{
    static const string routine("SCIDtreeParameters::validatePop2Object");

	try 
    {
		if (m_notional.size()!=parSpreads->size())
		    throw ModelException("not the same number of notionals as number of spread curves", routine);
		if (m_notional.size()==0)
		    throw ModelException("No names defined", routine);
	}
	catch (exception& e)
	{
		throw ModelException(e, routine);
	}

}


void SCIDtreeParameters::getSingleNameCalibrationDates(DateTimeArray &dates)
{
	dates.resize(singleNameExpiries->size());
	for (int k=0; k<dates.size(); k++) dates[k] = (*singleNameExpiries)[k]->toDate(today);
}

void SCIDtreeParameters::CalibrateBaseWorld(DoubleMatrix &survProba)
{
	DateTimeArray sgleNameDT;
	getSingleNameCalibrationDates(sgleNameDT);
	m_world.CalibrateSingleNames(today, sgleNameDT, survProba, idioRatio, cfRatio, jumpRatios);
}


void SCIDtreeParameters::load(CClassSP& clazz)
{
    clazz->setPublic(); 
    REGISTER(SCIDtreeParameters, clazz);
    SUPERCLASS(MarketObject);
	EMPTY_SHELL_METHOD(defaultSCIDtreeParameters);

	FIELD(name, "Name of the scid parameter class");
    FIELD_MAKE_OPTIONAL(name); // default name is CIDParameters::DEFAULTNAME

	FIELD(parSpreads,  "Porfolio of par spreads");
	FIELD(m_notional,  "Notionals");
	FIELD(today,  "Today");
	FIELD(singleNameExpiries,  "Dates at which the base CID worlds calibrate to single name survival probabilities");

	FIELD(volDiff, "volatility for the diffusion part");
	FIELD(decayDiff,  "mean reversion speed for the diffusion part");
	FIELD(jumpDecay,  "jump mean reversion speed");
	FIELD(jumpImpact, "jump Impact");
	FIELD(jumpFreq,   "frequencies of jumps");
	FIELD(jumpWeights, "non parametric disitrution");
	FIELD(jumpMeans,   "possible mean jump sizes");

	FIELD(idioRatio, "ratio of the survival probability going into the idiosyncratic surv proba");
	FIELD(cfRatio, "ratio of the survival probability going into the common factor surv proba");
	FIELD(jumpRatios, "ratio of the survival probability going into the jump surv proba");

	FIELD(bounds, "minimum and maximum values of the scalar function defined by the tree");
	FIELD(stressPEL, "how much weight on the PEL");
	FIELD(positions, "pricing positions");
	FIELD(weights, "initial Guess for the weights");


	FIELD(m_nbNames, "");
    FIELD_MAKE_TRANSIENT(m_nbNames);
    FIELD(m_nbJumpProcesses, "");
    FIELD_MAKE_TRANSIENT(m_nbJumpProcesses);
	FIELD(m_recovery, "");
    FIELD_MAKE_TRANSIENT(m_recovery);
}

IObject* SCIDtreeParameters::defaultSCIDtreeParameters(){
    return new SCIDtreeParameters();
}

CClassConstSP const SCIDtreeParameters::TYPE = CClass::registerClassLoadMethod(
   "SCIDtreeParameters", typeid(SCIDtreeParameters), load);



ICDSParSpreadsSP SCIDtreeParameters::getCDSparSpread(int name)
{
    static const string method = "SCIDtreeParameters::getCDSparSpread";
    if ((name<0)||(name>=m_nbNames)) throw ModelException(method, "name index out of range");
	return (*parSpreads)[name].getSP();
}

double SCIDtreeParameters::getMarketPortfolioEL(DateTime date, bool notionalEL)
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

double SCIDtreeParameters::getMarketPortfolioEL(double date, bool notionalEL)
{
	DateTime _date = today.rollDate(int(date*365+0.5));
	return getMarketPortfolioEL(_date,notionalEL);
}

DoubleArraySP SCIDtreeParameters::getMarketPortfolioEL(DateTimeArray dates, bool notionalEL)
{
    DoubleArraySP result(new DoubleArray(dates.size()));
    for (int i=0; i<dates.size(); i++) 
        (*result)[i]=getMarketPortfolioEL(dates[i],notionalEL);
    return result; 
}

void SCIDtreeParameters::MarketSurvProb(DoubleArray &time, DoubleMatrix &survProbability)
{
	DateTimeArray dates(time.size());
	for (int i=0; i<time.size(); i++) dates[i] = today.rollDate(int(time[i]*365+0.5));
	MarketSurvProb(dates, survProbability);
}

void SCIDtreeParameters::MarketSurvProb(DateTimeArray &dates, DoubleMatrix &survProbability)
{
	survProbability.resize(m_nbNames,dates.size());
	for (int m=0; m<m_nbNames; m++)
	{
		ICDSParSpreadsSP CDS = getCDSparSpread(m);
		for (int k=0; k<dates.size(); k++)
			survProbability[m][k] = CDS->survivalProb(dates[k]);
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////TIME LINE HELP /////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

void SCIDtreeParameters::push_backTimeLine(DateTimeArray &timeLine,
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


DoubleArray SCIDtreeParameters::DateAsDouble(const DateTimeArray &dates)
{
	DoubleArray t(dates.size());
	for (int i=0; i<t.size(); i++)
		t[i] = today.yearFrac(dates[i]);
	return t;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////FAST MC/////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

void SCIDtreeParameters::setFastMC(int seed,               // seed used by the random number	
					  double timeSteps,		  // time steps to discretise the CIR process	
					  double lastTime,		  // until when do we simulate
					  int noJumpPaths,		  // how many paths conditional on no jumps
					  int jumpsPaths,		  // how many paths conditional on at least one jumps
					  int _convolutionNoJump, // convolution method conditional on no jumps
					  int _convolutionJumps,  // convolution method conditional on at least one jumps
					  DoubleArray &kmin,      // lower attachment points
					  DoubleArray &kmax,      // upper attachment points 
					  int		  coupon,	  // coupon paid every ? months
					  DateTimeArray maturities, // tree dates
					  YieldCurveSP &discount, // discount factor use
					  DayCountConventionSP &dcc)  // day count convention used
{
	QLIB_VERIFY(timeSteps>0, "negative time steps in setFastMC");
	QLIB_VERIFY(lastTime>0, "negative last Time in setFastMC");
	m_world.InitializeRNG(seed);
	if (m_world.hasCommonFactor())
		m_world.setDiscretizationTimesCF(lastTime, timeSteps);
	fast.nbPathsNoJump = noJumpPaths;
	fast.nbPathsJumps = jumpsPaths;
	fast.lastTime = lastTime;

	 // deal with jumps when available

	convol.setParameters(m_nbNames,kmin,kmax);
	convol.setLGD(&m_recovery[0],&m_notional[0]);
	fast.convolutionNoJump = _convolutionNoJump;
	fast.convolutionJumps = _convolutionJumps;
	m_tree.setTreeDates(maturities,today);
	m_tree.setTrancheData(kmin, kmax, coupon, discount, dcc);

}

// Compute, as seen as today, TEL between times[0] and times
// parameters have been set up in setFastMC

void  SCIDtreeParameters::computeTELinGivenWorlds(DoubleArray &times,
												  vector< LinearInterpolantSP > &f,
							  					  DoubleMatrixArray & TEL) 
{
	m_world.setWorlds(f);
	size_t nbWorlds = f.size();
	QLIB_VERIFY(nbWorlds == TEL.size(), "TEL of the wrong size");
	int nbTranches = convol.getNbTranches();
	DoubleMatrix survProba(times.size(),m_nbNames);
	vector< DoubleMatrix > minusLogSurvProbaIdio(nbWorlds,survProba);

	DoubleMatrixArray TELdiff(nbWorlds);

	for (size_t i=0; i<nbWorlds; i++)
	{
		TEL[i]->resize(times.size(), nbTranches);
		TEL[i]->fill(0.0);
		TELdiff[i] = CDoubleMatrixSP( new DoubleMatrix(times.size(), nbTranches));
		minusLogSurvProbaIdio[i].resize(m_nbNames,times.size());
		minusLogSurvProbaIdio[i].fill(0.0);
		for (int m=0; m<m_nbNames; m++) m_world.minusLogSurvProbaIdio(i, m, times, minusLogSurvProbaIdio[i][m]);
	}
	DoubleArray minusLogSurvProbaCond;
	double probaZeroJump=1.0;
	int totalNbPaths = 1;
	if ( (m_world.hasCommonFactor()) || (m_world.getNbMarkets()>0) )
	{
		m_world.setIpathFast(0);
		minusLogSurvProbaCond.resize(times.size());
	}
	if (m_world.hasCommonFactor()) totalNbPaths = fast.nbPathsNoJump;
	if (m_world.getNbMarkets()>0) 
	{
		probaZeroJump = exp(-fast.lastTime*m_world.getSumFreq());
		totalNbPaths += fast.nbPathsJumps;
	}

	if (!m_world.hasCommonFactor())
	{
		for (size_t i=0; i<nbWorlds; i++)
		{
			for (int m=0; m<m_nbNames; m++)
				for (int l=0; l<times.size(); l++)
					survProba[l][m] = exp(-minusLogSurvProbaIdio[i][m][l]);
			for (int l=0; l<times.size(); l++)
				convol.computeTEL(fast.convolutionNoJump, survProba[l], (*TELdiff[i])[l]);
		}
	}
	else
	{
		for (int mc=0; mc<fast.nbPathsNoJump; mc++)
		{
			m_world.SimulateCFFast();
			for (size_t i=0; i<nbWorlds; i++)
			{
				for (int m=0; m<m_nbNames; m++) 
				{
					m_world.minusLogCondSurvProbaCF(i,m,times,&minusLogSurvProbaCond[0],false);
					for (int l=0; l<times.size(); l++) 
						survProba[l][m] = exp(-minusLogSurvProbaIdio[i][m][l]-minusLogSurvProbaCond[l]);
				}
				for (int l=0; l<times.size(); l++) 
					convol.computeTEL(fast.convolutionNoJump, survProba[l], (*TELdiff[i])[l]); 
			}
		}
		for (size_t i=0; i<nbWorlds; i++)
			TELdiff[i]->scale(1.0/double(fast.nbPathsNoJump));
	}
	for (size_t i=0; i<nbWorlds; i++)
		DoubleMatrix::scale(probaZeroJump,*TELdiff[i],&(*TEL[i]));
	if (m_world.getNbMarkets()>0)
	{
		double weight = (1.0 - probaZeroJump) / double(fast.nbPathsJumps);
		for (int mc=0; mc<fast.nbPathsJumps; mc++)
		{
			m_world.SimulateCFFast();
			double firstJump = m_world.SimulateJumpsFast(fast.lastTime,true);
			int indexFirstJump = 0;
			while ((indexFirstJump<times.size())&&(times[indexFirstJump]<=firstJump)) 
				++indexFirstJump;
			if (indexFirstJump<times.size())
			{
				DoubleArray restrictedTimes(times.begin()+indexFirstJump,times.end());
				for (size_t i=0; i<nbWorlds; i++)
				{
					for (int m=0; m<m_nbNames; m++) 
					{
						fill(minusLogSurvProbaCond.begin()+indexFirstJump,minusLogSurvProbaCond.end(),0.0);
						m_world.minusLogCondSurvProbaCF  (i,m,restrictedTimes,&minusLogSurvProbaCond[indexFirstJump],true);
						m_world.minusLogCondSurvProbaJump(i,m,restrictedTimes,&minusLogSurvProbaCond[indexFirstJump],true);
						for (int l=indexFirstJump; l<times.size(); l++) 
							survProba[l][m] = exp(-minusLogSurvProbaIdio[i][m][l]-minusLogSurvProbaCond[l]);
					}
					for (int l=0; l<indexFirstJump; l++)
						for (int j=0; j<nbTranches; ++j)
							(*TEL[i])[l][j] += weight*(*TELdiff[i])[l][j];

					for (int l=indexFirstJump; l<times.size(); l++) 
						convol.computeTEL(fast.convolutionJumps, survProba[l], (*TEL[i])[l], weight); 
				}
			}
			else
				for (size_t i=0; i<nbWorlds; i++)
					(*TEL[i]).add( *TELdiff[i], weight);
		}
	}
}

// Compute, as seen as today, PEL between times[0] and times
void SCIDtreeParameters::computePELinGivenWorlds(double time,
									vector< LinearInterpolantSP > &f,
									vector<double *> & PEL)  // must be initialized
{
	size_t nbWorlds = f.size();
	DoubleArray minusLogSurvProba(1);
	DoubleArray times(1,time);
	m_world.setWorlds(f);
	for (size_t i=0; i<nbWorlds; i++)
	{
		*PEL[i]=0.0;
		for (int m=0; m<m_nbNames; m++) 
		{
			m_world.minusLogSurvProba(i,m,times,&minusLogSurvProba[0],false);
			*(PEL[i]) +=  ( 1.0 - exp(-minusLogSurvProba[0]) ) * m_notional[m]*(1-m_recovery[m]);
		}
	}
}

void SCIDtreeParameters::computeTELandPELOneGeneration(vector<SimpleTreeSP> &leaves)
{
	QLIB_VERIFY(leaves.size()>0, "nothing to do here in computeTELandPELOneGeneration");
	size_t generation = leaves[0]->generation;
	QLIB_VERIFY(generation>0 && generation<= m_tree.getMaxNbGeneration(), "generation outside range in computeTELinAllWorlds");
	vector<LinearInterpolantSP> f(leaves.size());
	DoubleArray timesTEL(m_tree.TELtimes.begin()+m_tree.TELindexes[generation-1]+1,
						 m_tree.TELtimes.begin()+m_tree.TELindexes[generation]+1);

	DoubleMatrixArray TEL(leaves.size());
	vector<double*> PEL(leaves.size());
	m_tree.getWorldsFromTree(f,leaves);
	for (size_t i=0; i<f.size(); i++) 
	{
		size_t gen = leaves[i]->generation;
		QLIB_VERIFY(gen==generation, "generation out of range in computeTELinAllWorlds");
		TEL[i] = leaves[i]->TEL;
		PEL[i] = &leaves[i]->PEL;
	}
	computeTELinGivenWorlds(timesTEL, f, TEL);
	computePELinGivenWorlds(m_tree.treeTimes[generation-1], f, PEL);
}


void SCIDtreeParameters::computeTELandPELNGeneration(vector<SimpleTreeSP> &leaves, int lastGen)
{
	int gen = leaves[0]->generation;
	if (gen<=lastGen) computeTELandPELOneGeneration(leaves);
	if (gen<lastGen) 
		for (size_t i=0; i<leaves.size(); ++i)
			computeTELandPELNGeneration(leaves[i]->children, lastGen);
}

void SCIDtreeParameters::computeAllTELandPEL(bool feedTheTree)
{
	vector< SimpleTreeSP > leaves;
	size_t maxNbGen = m_tree.getMaxNbGeneration();
	for (size_t gen=1; gen<=maxNbGen; ++gen) 
	{
		if (feedTheTree) m_tree.feed(gen-1,a_t[gen-1]);
		leaves.clear();
		m_tree.getLeaves(leaves,gen);
		computeTELandPELOneGeneration(leaves);
	}
}

void SCIDtreeParameters::computeAllLegsFromTEL()
{
	vector< SimpleTreeSP > leaves;
	for (size_t gen=1; gen<=m_tree.getMaxNbGeneration(); ++gen) 
	{
		leaves.clear();
		m_tree.getLeaves(leaves,gen);
		m_tree.fromTELtoLegsOneGeneration(leaves);
	}
}

void SCIDtreeParameters::computeAllTELeasy(DoubleMatrixArray *outputTEL)
{
	if (outputTEL!=0)
	{
		vector< SimpleTreeSP > leaves;
		m_tree.getLeaves(leaves);
		vector<LinearInterpolantSP> f(leaves.size());
		outputTEL->resize(leaves.size());
		m_tree.getWorldsFromTree(f,leaves);
		for (size_t i=0; i<f.size(); i++) 
			(*outputTEL)[i] = CDoubleMatrixSP(new DoubleMatrix());
		computeTELinGivenWorlds(m_tree.TELtimes, f, *outputTEL);	
	}
}

//void SCIDtreeParameters::ImprovePosition(const DoubleMatrix &parSpread,
//												 const DoubleArray &PEL,
//												 int lowGen,
//												 int highGen,
//												 bool reconcile,
//												 double epsilon,
//												 int maxNbLoops,
//												 double ftol)      // terminate when values move by less than ftol
//{
//	bool test = m_tree.check();
//	DoubleArray eps(2);
//	eps[0]=epsilon; eps[1]=-epsilon;
//	int nbLoops = 0;
//	double value = 0.0;
//	double newValue = value;
//	vector<SimpleTreeSP> leaves;
//	for (int gen=lowGen; gen<=highGen; ++gen)
//	{
//		leaves.clear();
//		m_tree.getLeaves(leaves,gen);
//		computeTELandPELOneGeneration(leaves);
//		m_tree.fromTELtoLegsOneGeneration(leaves);
//	}
//	test = m_tree.check();
//
//	while ( (nbLoops<2) || ( (nbLoops < maxNbLoops) && (newValue<value-ftol) ) )
//	{
//		if (lowGen==highGen)
//		{
//			leaves.clear();
//			m_tree.getLeaves(leaves,lowGen);
//			vector<SimpleTreeSP> newLeaves;
//			for (size_t i=0; i<leaves.size(); ++i)
//			{
//				test = m_tree.check();
//				leaves[i]->CreateNearbyLeaves(newLeaves, eps, bounds[0], bounds[1], true);
//				test = m_tree.check();
//			}
//			computeTELandPELOneGeneration(newLeaves);
//			m_tree.fromTELtoLegsOneGeneration(newLeaves);
//		}
//		else
//		{
//			for (int gen=lowGen; gen<=highGen; ++gen)
//			{
//				leaves.clear();
//				m_tree.getLeaves(leaves,gen);
//				vector<SimpleTreeSP> newLeaves;
//				for (size_t i=0; i<leaves.size(); ++i)
//					leaves[i]->CreateNearbyLeaves(newLeaves, eps, bounds[0], bounds[1], true);
//				computeTELandPELNGeneration(newLeaves,highGen);
//				m_tree.fromTELtoLegsNGeneration(newLeaves,highGen);
//			}
//		
//			//vector< vector<SimpleTreeSP> > newLeaves(highGen+1);
//			//for (size_t i=0; i<leaves.size(); ++i)
//			//	leaves[i]->CreateNearbyLeavesMultipleGens(newLeaves, highGen, eps, bounds[0], bounds[1], true);
//			//for (size_t gen=lowGen; gen<=highGen; ++gen)
//			//{
//			//	leaves.clear();
//			//	m_tree.getLeaves(leaves,gen);
//			//	computeTELandPELOneGeneration(leaves);
//			//	m_tree.fromTELtoLegs(leaves);
//			//	//computeTELandPELOneGeneration(newLeaves[gen]);
//			//	//m_tree.fromTELtoLegs(newLeaves[gen]);
//			//}
//		}
//		test = m_tree.check();
//		value = newValue;
//		//newValue = m_tree.findWeights(lowGen, highGen, parSpread, PEL, stressPEL, 0.0, true);
//		test = m_tree.check();
//		++nbLoops;
//	}
//	if (reconcile)
//	{
//		for (int gen=highGen; gen>=lowGen; --gen)
//		{
//			leaves.clear();
//			m_tree.getLeaves(leaves,gen-1);
//			for (size_t j=0; j<leaves.size(); ++j)
//				leaves[j]->reconcile(2.0*epsilon,true);
//		}
//		for (int gen=lowGen; gen<=highGen; ++gen)
//		{
//			leaves.clear();
//			m_tree.getLeaves(leaves,lowGen);
//			computeTELandPELOneGeneration(leaves);
//			m_tree.fromTELtoLegsOneGeneration(leaves);
//		}
//		m_tree.findWeights(lowGen, highGen, parSpread, PEL, stressPEL, 0.0, true);
//	}
//}
//



void SCIDtreeParameters::CalibrateTreeToTranches(const DoubleMatrix &parSpread,
												const DoubleArray &PEL,
												double smoothing,
												size_t lowGen,
												size_t highGen)
{
	QLIB_VERIFY(lowGen<=highGen, "highGen<lowGen!");
	QLIB_VERIFY(lowGen>0 && highGen<=m_tree.getMaxNbGeneration(), "generation out of range");
	vector<SimpleTreeSP> leaves; 
	DoubleArray invDist;
	for (size_t gen=lowGen; gen<=highGen; ++gen)
	{
	   m_tree.feed(gen-1,a_t[gen-1]);
	   leaves.clear();
	   m_tree.getLeaves(leaves,gen);
	   computeTELandPELOneGeneration(leaves);
	   m_tree.fromTELtoLegsOneGeneration(leaves);
	   DoubleArray invDistTemp;
	   m_tree.getInverseDistanceBetweenLeaves(leaves, invDistTemp);
	   invDist.insert(invDist.end(),invDistTemp.begin(),invDistTemp.end());
	}
	for (int i=0; i<invDist.size(); ++i) invDist[i] *= smoothing;
	m_tree.findWeights(lowGen, highGen, parSpread, PEL, stressPEL, invDist, true);
}


void SCIDtreeParameters::guessPositionAtAllMaturities(
							 const DoubleMatrix &parSpread,
							 const DoubleArray &PELmarket,
							 double Dx,
							 double smoothing,
							 DoubleArrayArray &selectedPos)
{
	int nbPoints = 0; 
	vector<LinearInterpolantSP> f;
	DoubleArray x(1,1.0), y(1, bounds[0]);
	QLIB_VERIFY(Dx>0 && (bounds[1]-bounds[0])/Dx<1000, "too many points selected");
	while (y[0]<bounds[1])
	{
		f.push_back(LinearInterpolantSP( new LinearInterpolant(x,y)));
		y[0] += Dx;
		++nbPoints;
	}
	DoubleMatrixArray TEL(nbPoints);
	vector<double> PELworlds(nbPoints);
	selectedPos.clear();
	selectedPos.resize(m_tree.treeDates.size());
	for (int i=0; i<nbPoints; ++i)
		TEL[i] = CDoubleMatrixSP( new DoubleMatrix());
	computeTELinGivenWorlds(m_tree.TELtimes, f, TEL);

	DoubleMatrix weightConstraints(1,nbPoints);
	DoubleArray sumWeightsConstraints(1, 1.0);
	for (int l=0; l<nbPoints; l++)
		weightConstraints[0][l] = 1.0;

	for (int i=0; i<m_tree.treeDates.size(); ++i)
	{
		computePELinGivenWorlds(m_tree.treeTimes[i], f, PELworlds);

		CashFlowArray cf;
		cf = SwapTool::cashflows(today, m_tree.treeDates[i], false, 1.0, m_tree.coupon, "M", &(*m_tree.dcc));
		cf.back().amount -= 1.0; // ugly

		DoubleArray DL(m_tree.kmin.size()),RA(m_tree.kmin.size());
		vector<DoubleArray> MTMs(nbPoints,DL), PELconstraints(nbPoints,DoubleArray(1));;

		for (int l=0; l<nbPoints; ++l)
		{
			SCIDtreeWorld::convertTELtoLegs(m_tree.kmin, m_tree.kmax, *TEL[l], today,
							     m_tree.TELdates, cf, m_tree.treeDates[i], m_tree.discount, &DL[0], &RA[0]); 
			for (int k=0; k<m_tree.kmin.size(); ++k) 
				MTMs[l][k] = DL[k] - parSpread[k][i] * RA[k]; // - upfront
			PELconstraints[l][0] = PELworlds[l] - PELmarket[i];
		}
		DoubleArraySP newWeights = DoubleArraySP (new DoubleArray(nbPoints,0.0));
		DoubleArray smoothingArray(nbPoints, smoothing/Dx);
		smoothingArray.back()=0.0;
		SCIDtreeWorld::calibrateQuadProg(newWeights, MTMs, PELconstraints, weightConstraints, sumWeightsConstraints, -1.0, smoothingArray);
		for (int l=0; l<nbPoints; ++l)
			if ((*newWeights)[l]>1e-8) selectedPos[i].push_back(bounds[0]+l*Dx);
	}
}


//void SCIDtreeParameters::CalibratorHelp(vector<SimpleTreeSP> &leaves, 
//										const DoubleArray &ratios, 
//										DoubleArray &MTMandPEL, // (O)
//										const DoubleMatrix &parSpreadMarket,
//										const DoubleArray PELmarket,
//										bool changeWeights,
//										double smoothing,
//										bool deleteZeroWorlds)
//{
//	QLIB_VERIFY(leaves.size()>0, "no leaves in CalibratorHelp");
//	int nbTranches = parSpreadMarket.numCols();
//	size_t gen = leaves[0]->generation;
//	QLIB_VERIFY(gen>0 && gen<=m_tree.getMaxNbGeneration(), "gen out of range");
//
//	MTMandPEL.resize(nbTranches+1);
//	fill(MTMandPEL.begin(),MTMandPEL.end(),0.0);
//
//	for (size_t i=0; i<leaves.size(); ++i)
//		leaves[i]->changePositions(ratios[i], bounds[0], bounds[1]);
//	leaves.clear();
//	m_tree.getLeaves(leaves,gen); // ??? DOES THIS MAKE SENSE???
//	computeTELandPELOneGeneration(leaves);
//	m_tree.fromTELtoLegsOneGeneration(leaves);
//
//	if (changeWeights) 
//		m_tree.findWeights(gen, gen, parSpreadMarket, PELmarket, stressPEL, smoothing, deleteZeroWorlds);
//
//	for (size_t i=0; i<leaves.size(); ++i)
//	{
//		QLIB_VERIFY(leaves[i]->DL->size()==nbTranches, "wrong number of tranches");
//		QLIB_VERIFY(leaves[i]->generation == gen, "leaves of different generation");
//		for (int j=0; j<nbTranches; ++j)
//		{
//			double DL = (*(leaves[i]->DL))[j];
//			double RA = (*(leaves[i]->RA))[j];
//			double ps = parSpreadMarket[j][gen-1];
//			MTMandPEL[j] += leaves[i]->weight * (DL - ps*RA);
//		}
//		MTMandPEL.back() +=  stressPEL * leaves[i]->weight*(leaves[i]->PEL - PELmarket[gen-1]);
//	}
//}


void scale(DoubleArrayArray &x, double alpha)
{
	for (int i=0; i<x.size(); ++i) for (int j=0; j<x[i].size(); ++j) x[i][j] *= alpha;
}
void add(DoubleArrayArray &x, const DoubleArrayArray &y, double alpha)
{
	for (int i=0; i<x.size(); ++i) for (int j=0; j<x[i].size(); ++j) x[i][j] += alpha*y[i][j];
}

void scale(DoubleArray &x, double alpha)
{
	for (int i=0; i<x.size(); ++i) x[i] *= alpha;
}
void add(DoubleArray &x, const DoubleArray &y, double alpha)
{
	for (int i=0; i<x.size(); ++i) x[i] += alpha*y[i];
}


void SCIDtreeParameters::readAverageLegs(vector < SimpleTreeSP > &leaves, DoubleMatrix &DL, DoubleMatrix &RA)
{
	if (leaves.size()==0) return;
	DoubleMatrix leafDL,leafRA;
	m_tree.readLegs(leaves[0],leafDL,leafRA);
	leafDL.scale(leaves[0]->weight);
	leafRA.scale(leaves[0]->weight);
	DL = leafDL; RA = leafRA;
	for (size_t i=1; i<leaves.size(); ++i)
	{
		m_tree.readLegs(leaves[i],leafDL,leafRA);
		DL.add(leafDL,leaves[i]->weight);
		RA.add(leafRA,leaves[i]->weight);
	}
}
void SCIDtreeParameters::readAveragePEL(vector < SimpleTreeSP > &leaves, DoubleArray &PEL)
{
	if (leaves.size()==0) return;
	DoubleArray leafPEL;
	m_tree.readPEL(leaves[0], leafPEL);
	scale(leafPEL,leaves[0]->weight);
	PEL = leafPEL;
	for (size_t i=1; i<leaves.size(); ++i)
	{
		m_tree.readPEL(leaves[i],leafPEL);
		add(PEL,leafPEL,leaves[i]->weight);
	}
}

void SCIDtreeParameters::readAverageTEL(vector < SimpleTreeSP > &leaves, DoubleMatrix &TEL)
{
	if (leaves.size()==0) return;
	DoubleMatrix leafTEL;
	m_tree.readTEL(leaves[0], leafTEL);
	leafTEL.scale(leaves[0]->weight);
	TEL = leafTEL;
	for (size_t i=1; i<leaves.size(); ++i)
	{
		m_tree.readTEL(leaves[i],leafTEL);
		TEL.add(leafTEL,leaves[i]->weight);
	}
}

void SCIDtreeParameters::readAllLegs(vector < SimpleTreeSP > &leaves, DoubleMatrixArray &DL, DoubleMatrixArray &RA)
{
	if (leaves.size()==0) return;
	DL.resize(leaves.size());
	RA.resize(leaves.size());
	for (size_t i=0; i<leaves.size(); ++i)
	{
		DL[i]= CDoubleMatrixSP(new DoubleMatrix());
		RA[i]= CDoubleMatrixSP(new DoubleMatrix());
		m_tree.readLegs(leaves[i],*DL[i],*RA[i]);
	}
}

void SCIDtreeParameters::readAllPEL(vector < SimpleTreeSP > &leaves, DoubleArrayArray &PEL)
{
	if (leaves.size()==0) return;
	PEL.resize(leaves.size());
	for (size_t i=0; i<leaves.size(); ++i)
		m_tree.readPEL(leaves[i],PEL[i]);
}

void SCIDtreeParameters::readAllTEL(vector < SimpleTreeSP > &leaves, DoubleMatrixArray &TEL)
{
	if (leaves.size()==0) return;
	TEL.resize(leaves.size());
	for (size_t i=0; i<leaves.size(); ++i)
	{
		TEL[i] = CDoubleMatrixSP ( new DoubleMatrix());
		m_tree.readTEL(leaves[i],*TEL[i]);
	}
}



void SCIDtreeParameters::TESTSurvProbinGivenWorldFastMC (DoubleArray &times,
														 LinearInterpolantSP f,
														 DoubleMatrix &survProbability)
{
	//m_world.changeXaxisToCommonTimesLinearInterpolant(f); // ensure world X axis is the right one
	//DoubleArray fy;
	//m_world.getYaxisFatTimes(f,fy); // ensure world X axis is the right one

	//survProbability.resize(m_nbNames, times.size()); 
	//survProbability.fill(0.0);
	//DoubleArray spMC(times.size());

	//if (m_world.hasCommonFactor())
	//{
	//	for (int mc=0; mc<fast.nbPathsNoJump; mc++)
	//	{
	//		m_world.SimulateBaseFast();
	//		for (int m=0; m<m_nbNames; m++)
	//		{
	//			m_world.minusLogCondSurvProbaCF(fy,m,times,&spMC[0],false);
	//			for (int i=0; i<times.size(); i++) survProbability[m][i] += exp(-spMC[i]);
	//		}
	//	}
	//	survProbability.scale(1.0/fast.nbPathsNoJump);
	//}
	//else survProbability.fill(1.0);
	//for (int m=0; m<m_nbNames; m++)
	//{
	//	m_world.minusLogSurvProbaIdio(f,m,times,&spMC[0],false);
	//	for (int i=0; i<times.size(); i++) 
	//		survProbability[m][i] *= exp(-spMC[i]);
	//}
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////FULL MC/////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////


void SCIDtreeParameters::setFullMC(long seed, DoubleArray &discTimes)
{

}

void SCIDtreeParameters::setFullMC(long seed, DateTimeArray &discDates)
{

}

/** simulate the spreads in all worlds
always call SimulateNbJumps first, as it assumes it knows how many jumps there are */
void SCIDtreeParameters::SimulateSpread()
{
	m_world.SimulateSpreadFull(full.spread);
}

/** simulate default times */
void SCIDtreeParameters::SimulateDefaults()
{
	// we assume that the spread have already been simulated

	double OneDivMeanForUpdate;
	for (int i=0; i<m_nbWorlds; i++) 
		full.weightsDynamic[0][i] = m_world.fweights[i];

	fill(full.defaultTime.begin(),full.defaultTime.end(),999); // no defaults yet
	full.loss[0]=0.0;			 // no defaults yet
	full.notionalLoss[0]=0.0;       // no defaults yet
	full.nbDefaultedNames[0]=0;

	DoubleArray forUpdate(m_nbWorlds), survProbaWorld(m_nbWorlds);
	int guessF = 0;
	for (int k=0; k<full.discTimesDefaults.size()-1; k++)
	{
		fill(forUpdate.begin(),forUpdate.end(),1.0);
		full.loss[k+1] = full.loss[k];
		full.notionalLoss[k+1] = full.notionalLoss[k];
		full.nbDefaultedNames[k+1] = full.nbDefaultedNames[k];

		for (int m=0; m<m_nbNames; m++)
		{
			double survProbName = 0.0;
			if (full.defaultTime[m]>=999)
			{
				for(int i=0; i<m_nbWorlds; i++)
				{
					double intensityInGivenWorld = 0.0;
					int j=full.spread.discTimesLossIndices[k];
					double spreadfj = m_world.f[i]->valueWithGuess(full.spread.discTimesSpread[j],guessF)
										*full.spread.intensity[m][j];
					double spreadfj1;
					for (; j<full.spread.discTimesLossIndices[k+1]; ++j)
					{
						spreadfj1 = m_world.f[i]->valueWithGuess(full.spread.discTimesSpread[j+1],guessF)
										*full.spread.intensity[m][j+1];
						intensityInGivenWorld += 0.5*(spreadfj+spreadfj1)*(full.spread.discTimesSpread[j+1]-full.spread.discTimesSpread[j]);
						spreadfj = spreadfj1;
					}
					survProbaWorld[i] = exp (-intensityInGivenWorld);
					forUpdate[i] *= survProbaWorld[i];
					survProbName += full.weightsDynamic[k][i] * survProbaWorld[i];
				}
			}
			else survProbName = 1.0;
			double unif = m_world.getRandomNb(false, true);
			if (unif>survProbName)
			{
				full.defaultTime[m] = full.spread.discTimesSpread[k] + (unif-survProbName)/(1.0-survProbName)
					* (full.spread.discTimesSpread[k+1] - full.spread.discTimesSpread[k]);
				full.notionalLoss[k+1] += m_notional[m];
				++full.nbDefaultedNames[k+1];
				full.loss[k+1] += m_notional[m]*(1-m_recovery[m]);
				OneDivMeanForUpdate = 0;
				for (int i=0; i<m_nbWorlds; i++)
					forUpdate[i] *= survProbaWorld[i] - 1.0;

			}
		}

			// update weights
		OneDivMeanForUpdate = 0;
		for(int i=0; i<m_nbWorlds; i++) OneDivMeanForUpdate += full.weightsDynamic[k][i]*forUpdate[i];
		OneDivMeanForUpdate = 1/OneDivMeanForUpdate;
		for(int i=0; i<m_nbWorlds; i++) 
			full.weightsDynamic[k+1][i] = full.weightsDynamic[k][i]*forUpdate[i]*OneDivMeanForUpdate;
	}

}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////






DEFINE_TEMPLATE_TYPE(SCIDtreeParametersWrapper);
DRLIB_END_NAMESPACE

