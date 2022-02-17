//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids DR
//
//   Filename    : MCPathConfigCCM.cpp
//
//   Description : Implementation of MCPathConfig class, which drives the
//				   simulation of default times for a CCMCopulaModel
//
//   Date        : Feb 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/MCPathConfigCCM.hpp"

#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/DependentCopulaModel.hpp"
#include "edginc/IndependentCopulaModel.hpp"
#include "edginc/CCMCopulaModel.hpp"
#include "edginc/BenchmarkDate.hpp"
#include "edginc/MCCache.hpp"
#include "edginc/DECLARE.hpp"
#include "edginc/DependenceGauss.hpp"
#include "edginc/GaussSampling.hpp"
#include "edginc/imsl.h"
#include "edginc/BaseSimulation.hpp"

DRLIB_BEGIN_NAMESPACE

class StratMCRandomNoCache : public MCRandomNoCache
{

	IMFSamplingSP sampling;
	int nPaths;

	DoubleArrayConstSP points;
	IntArrayConstSP  numSamples;

	int pointIndex;

	int sampleSize;

	DoubleMatrix buffer;

public:

	// sampling needs to have been configured for nPaths
	StratMCRandomNoCache(
		MCRandomNoCache mcRand,
		IMFSamplingSP sampling,
		int nPaths,
		int sampleSize) :
	MCRandomNoCache(mcRand),
	sampling(sampling),
	nPaths(nPaths),
	pointIndex(0),
	sampleSize(sampleSize),
	buffer(1, sampleSize)
	{
		if (sampling.get())
		{
			points = sampling->getPoints();
			numSamples = sampling->getCumNumberSamples();
		};
	};

	void generate(int pathIdx)
	{
		if (pathIdx == 0)
			pointIndex = 0;

		MCRandomNoCache::generate(pathIdx);

		const DoubleMatrix & temp = MCRandomNoCache::getRandomNumbers();

		int i;
		for (i = 0; i < sampleSize; i++)
		{
			buffer[0][i] = temp[i][0];
		}

		if (!(sampling.get()))
			return;

		while (pathIdx >= (*numSamples)[pointIndex])
		{
			pointIndex++;
			if (pointIndex >= numSamples->size())
			{
				// something's gone wrong
				throw ModelException("Drawing too many paths");
			};
		};

		// doctor the first one based on pathIdx

		buffer[0][0] = (*points)[pointIndex];
	};

	virtual const DoubleMatrix & getRandomNumbers() const {
		return buffer;
	};
};

///////////////////////////
IModel::WantsRiskMapping MCPathConfigCCM::wantsRiskMapping() const
{
	return IModel::WantsRiskMapping(true);
};

int
MCPathConfigCCM::storagePerPath(IMCProduct* product) const
{
	static const string routine = "MCPathConfigCCM::storagePerPath";

	throw ModelException("Not to be used", routine);

	return 0;
}


MCPathGeneratorSP
MCPathConfigCCM::pastPathGenerator(const IMCProduct* product)
{
	const MCProductClient* prodClient = dynamic_cast<const MCProductClient*>(product);
	MCPathGeneratorSP temp = MCPathGeneratorSP(new PathGen(prodClient)); // default has no past
	return temp;
}

IRandomSP
MCPathConfigCCM::getRandomGenerator() const
{
    if (!initRandState){
        randomGen->init();
    } else {
        randomGen->setState(initRandState.get());
    }
    return randomGen;
}


void MCPathConfigCCM::validatePop2Object()
{
	ccmCopula = CCMCopulaModelSP(new CCMCopulaModel());

	//if (!(sampling.get()))
	//	sampling = IMFSamplingSP(new GaussSampling());
};


MCPathGeneratorSP
MCPathConfigCCM::futurePathGenerator(
	int                                cachingMode,
	int                                numPaths,
	const MCPathGeneratorSP&           pastPathGenerator,
	const IMCProduct*                  prod,
	Control*                           control,
	Results*                           results,
	DateTimeArray&                     simDates ) //from the parent, IMCPathConfig
{
	try
	{
		bool cachingRequested = false; // JCP
		ccmCopula->Update();

		// JCP do i need to do something with the pastPathGenerator
		// JCP does the IMCProduct have number of markets (for RFL)? other sampling information?
		// JCP assume 1 market for now
		// JCP need to reset the randomCache before using it in the path generator

		const MCProductClient* prodClient = dynamic_cast<const MCProductClient*>(prod);
		MCPathGeneratorSP temp =
				MCPathGeneratorSP(
					new MCPathConfigCCM::PathGen(
								ccmCopula.get(),
								this,
								prodClient,
								numPaths,
								1));

		return temp;
		// JCP need to add number of factors and other sampling parameters to the interface properly
	}
	catch(exception& e)
	{
		throw ModelException(e, "MCPathConfigCCM::futurePathGenerator");
	}
}

MarketDataFetcherSP MCPathConfigCCM::marketDataFetcher() const
											//from the parent, IMCPathConfig
{
	return MarketDataFetcherSP(new MarketDataFetcher());
}

void
MCPathConfigCCM::getMarket(	const IModel*      model,
							const MarketData*  market,
							IInstrumentCollectionSP instrument)//from the parent, IMCPathConfig
{
	try
	{
		//Nothing happens - here might remove it later
	}
	catch(exception& e)
	{
		throw ModelException(e, "MCPathConfigCCM::getMarket");
	}
}


void MCPathConfigCCM::getComponentMarketData(
    const CModel*         model,
    const MarketData*     market,
    MarketObjectSP        mo,
    MarketDataFetcher* mdf) const
{
	static const string routine = "MCPathConfigCCM::getMarket";

	throw ModelException("Not to be used", routine);
};

void
MCPathConfigCCM::generateEngineTimeline(
	const DateTimeArray& productTimeline,
	DateTimeArray& engineTimeline,	//output
	IntArray& fromEngineIndexToProductIndex //output
	) const
{
	//start and end dates of product and engine timeline would match
	engineTimeline.clear();
	fromEngineIndexToProductIndex.clear();

	if (productTimeline.size() == 0)
		return;

	engineTimeline.push_back(productTimeline[0]);
	fromEngineIndexToProductIndex.push_back(0);

	bool checkForInsert = true;
	for (int i=1; i< productTimeline.size(); ++i)
	{
		if (timelineFreqObj.get())
		{
			while (checkForInsert)
			{
				DateTime date = timelineFreqObj->toDate(engineTimeline.back());
				if (date < productTimeline[i] )
				{
					engineTimeline.push_back(date);
					fromEngineIndexToProductIndex.push_back(i);
				}
				else
					checkForInsert = false;
			}
		}
		checkForInsert = true;
		engineTimeline.push_back(productTimeline[i]);
		fromEngineIndexToProductIndex.push_back(i);
	}
}

MarketObjectSP MCPathConfigCCM::getMarket(
		const CModel        *model,
        const MarketData    *market,
		const string        &name,
        const CClassConstSP &type) const
{

	static const string routine = "MCPathConfigCCM::getMarket";

	throw ModelException("Not to be used", routine);
}


void MCPathConfigCCM::PathGen::generatePath(int pathIdx)
{
	static const string routine = "MCPathConfigCCM::generatePath";
	randomGen->generate(pathIdx);
	simulation.computeSample(pathIdx);

	//we update the portfolioName state variables
	const DateTimeArray&	defTimeline = simulation.timeline();

	const IBaseSimulation::ISimulationState& simState
								= simulation.simulationState();

	const DoubleArray&		recoveries = simState.getRecoveries();
	const IntArray&			defIndices = simState.getDefaultIndices();

	int num = simState.getNumAssets();

	if (num!=(int)(portfolioNameDBase.size()))
	{
		stringstream message;
		message << "Number of Names in the model donot match the number of "
		        << "results from the MC engine";
		throw ModelException(routine, message.str());
	};

	for (int i=0; i<num; ++i)
	{
		int defaultTimeIndex = defIndices[i];

		//asset has defaulted
		if ((defaultTimeIndex < defTimeline.size()-1) && (defaultTimeIndex > -1))
		{
			//ASarma: must revisit this +1 business;
			(*portfolioNameDBase[i])[0].defaultTimeIndex =
				defaultTimeIndex+1;

			(*portfolioNameDBase[i])[0].productTimelineIndex =
				fromEngineIndexToProductIndex[defaultTimeIndex+1];
		}
		else
			(*portfolioNameDBase[i])[0].defaultTimeIndex = -1;

		(*portfolioNameDBase[i])[0].recovery  = recoveries[i];
	}

	//we do not update the SVDiscFactor as we use deterministic discount factors
}

IStateVariableSP
MCPathConfigCCM::PathGen::create(const IStateVariableGen *svGen)
{
	static const string routine = "MCPathConfigCCM::PathGen::create";

	try
	{
		return svDBase.find(svGen);
	}
	catch (exception& e)
	{
		throw ModelException(e, routine);
	}
}

//used by the past path generator
MCPathConfigCCM::PathGen::PathGen(const MCProductClient*   prodClient)
{
    StateVariableCollectorSP svCollector(new StateVariableCollector());

	//collect SVGens for credit loss configs
	prodClient->collectStateVars(svCollector);

    IElemStateVariableGenArray stateVarGenArray
		= svCollector->getElemStateVarGens();

	std::vector<const PortfolioName::IndexedSVGen*>
		jsv(filterStateVars<PortfolioName::IndexedSVGen>(stateVarGenArray));

	for (unsigned int iVar = 0; iVar < jsv.size(); ++iVar)
		svDBase.append(
			jsv[iVar],
			jsv[iVar]->createNewSV(this));

	//collect SVGens for discount factors
	std::vector<const SVGenDiscFactor*>
		jsv1(filterStateVars<SVGenDiscFactor>(stateVarGenArray));

	for (unsigned int j=0; j< jsv1.size(); ++j)
	{
		SVDiscFactorSP stateVariable(jsv1[j]->determinsticSV(false));
		svDBase.append(
				jsv1[j],
				stateVariable);
	}
}

class MCPathBaseCCM : public DependenceMakerGauss::Support,
					  virtual public IMCRandom::Callbacks
{
	int numRands;
	CDoubleMatrixSP matrix;

public:

	MCPathBaseCCM(int nfactors) : numRands(nfactors)
	{
		int i;
		matrix = CDoubleMatrixSP(new DoubleMatrix(nfactors, nfactors));

		for (i = 0; i < numRands; i++)
			(*matrix.get())[i][i] = 1.0;

	};

	CDoubleMatrixConstSP getGaussData() const
	{
	    return matrix;
	};

	virtual void configureAntithetics() {
		// All done by MCRandom
		return;
	}

	/** Configures pathGen for nonAntithetics. We need a better name */
	virtual void configureNonAntithetics() {
		// All done by MCRandom
		return;
	};
};

MCPathConfigCCM::PathGen::PathGen(
	CCMCopulaModel *cop,
	MCPathConfigCCM *mcPathConfig,
	const MCProductClient*   prodClient,
	int nSamples,
	int nFactors)
{
	try
	{
		//generate engine timeline
		DateTimeArray productTimeline =
			prodClient->getSimSeries()->getAllDates(); //ASarma - revisit
		DateTimeArray engineTimeline;
		mcPathConfig->generateEngineTimeline(  productTimeline,
											engineTimeline,
											fromEngineIndexToProductIndex);
		int engineTimelineSize = engineTimeline.size();

		//update the copulaModel and the engine
		cop->valueDate = prodClient->getToday();

		// collect SVGens
		StateVariableCollectorSP svCollector(new StateVariableCollector());
		prodClient->collectStateVars(svCollector);
		IElemStateVariableGenArray stateVarGenArray =
			svCollector->getElemStateVarGens();

		// collect SVGens for discount factors
		std::vector<const SVGenDiscFactor*>
			jsv1(filterStateVars<SVGenDiscFactor>(stateVarGenArray));

		for (unsigned int j=0; j< jsv1.size(); ++j)
		{
			SVDiscFactorSP stateVariable(jsv1[j]->determinsticSV(false));
			svDBase.append(
					jsv1[j],
					stateVariable);
		}

		// collect SVGens for CreditLossConfigs
		std::vector<const ICreditLossConfig::IIndexedSVGen*>
			jsv(filterStateVars<ICreditLossConfig::IIndexedSVGen>(stateVarGenArray));

		// Each PortfolioNameSV has a part SubResult that is registered with the
		// engine.
		// While for a credit (defined uniquely by it credit) we may have many
		// PortfolioName (in difference CDO2 tranches for instance), they
		// will be triggered to default by a single credit event

		// The following map would store the SubResult that may be shared by all
		// PortfolioNames corresponding to a credit asset.
		map<string, PortfolioNameIndexedSVSubResultArraySP>
			portfolioNameIndexedSVSubResultsMap;

		// We also create a map of CreditAsset that need to be passed to the
		// Engine
		map<string, CreditAssetConstSP> creditAssetMap;

		// We also create a map of doubles that need to be passed to the Engine
		// as well
		map<string, double> portfolioNameBetaMap;

		// From the SVGens that we collected before, we now create the
		// corresponding SVs
		for (unsigned int iVar = 0; iVar < jsv.size(); ++iVar)
		{
			ICreditLossConfigIndexedSVSP stateVariable =
				jsv[iVar]->createNewSV(this);

			PortfolioNameIndexedSVSP portfolioNameIndexedSV =
				PortfolioNameIndexedSVSP(
					dynamic_cast<PortfolioName::IndexedSVGen::SV*>
						(stateVariable.get()));

			CreditAssetConstSP creditAsset =
				portfolioNameIndexedSV->getAsset();
			double beta = portfolioNameIndexedSV->getBeta();

			// If the asset has a past default date, we convert it to a product
			// timeline and enginetimeline index and update the SV

			// if past default has occurred
			DateTimeConstSP defaultSettDate =
				portfolioNameIndexedSV->getDefaultSettDate();

			if (defaultSettDate.get())
			{
				//convert it to the engine timeline index
				int engineIndex = max(
					0,
					defaultSettDate->findNearest(engineTimeline));

				int productIndex = fromEngineIndexToProductIndex[engineIndex];

				portfolioNameIndexedSV->initializePath(
					CIntConstSP(
						CInt::create(engineIndex)),
						productIndex,
						engineTimelineSize);
			}
			else
				portfolioNameIndexedSV->initializePath(
					CIntConstSP(),
					-1,
					engineTimelineSize);

			//populate the map with the SubResults
			PortfolioNameIndexedSVSubResultArraySP
				portfolioNameIndexedSVSubResults;

			string name =
				portfolioNameIndexedSV->getName();

			//fetch the sub result from the map if it exists else create a new
			// one
			//the int is the index of the first occurance of the
			// PortfolioNameSV
			map<string, PortfolioNameIndexedSVSubResultArraySP>::iterator iter=
				portfolioNameIndexedSVSubResultsMap.find(name);

			if (iter == portfolioNameIndexedSVSubResultsMap.end())
			{
				portfolioNameIndexedSVSubResults =
					PortfolioNameIndexedSVSubResultArraySP(
						// JCP was 1 should it be engineTimelineSize-1
						new PortfolioNameIndexedSVSubResultArray(1));


				//create a new entry in the map
				portfolioNameIndexedSVSubResultsMap[name] =
					portfolioNameIndexedSVSubResults;
			}
			else
				//fetch the existing one
				portfolioNameIndexedSVSubResults = (iter->second);

			portfolioNameIndexedSV->setSubResults(
						portfolioNameIndexedSVSubResults);

			svDBase.append(
				jsv[iVar],
				stateVariable);//here state variable is portfolioNameIndexedSV

			//update the creditAsset in the map
			map<string, CreditAssetConstSP>::iterator iter2 =
				creditAssetMap.find(name);
			if (iter2 == creditAssetMap.end())
				creditAssetMap[name] = creditAsset;

			//update the betas in the map
			map<string, double>::iterator iter1 =
				portfolioNameBetaMap.find(name);
			if (iter1 == portfolioNameBetaMap.end())
				portfolioNameBetaMap[name] = beta;
		}

		// The map created above gave us the unique list of names.
		// From the map, we obtain
		// (i) a vector of PortfolioNameIndexedSVSubResults that resides in
		// this class
		// (ii)a vector of CreditAssetArrayConst that is passed to the engine
		// (iii) a vector of betas that is passed to the engine


        // convert the portfolioNameIndexedSVSubResults Map into a vector
		int num = portfolioNameIndexedSVSubResultsMap.size();
		portfolioNameDBase.clear();
		portfolioNameDBase.resize(num);

		map<string, PortfolioNameIndexedSVSubResultArraySP>::iterator iter =
			portfolioNameIndexedSVSubResultsMap.begin();
		int i=0;
		for (; iter!= portfolioNameIndexedSVSubResultsMap.end(); ++iter)
		{
			portfolioNameDBase[i] = (iter->second);
			++i;
		}

		//convert the credit asset map into a vector
		CreditAssetWrapperArray creditAssetWrapperArray;
		map<string, CreditAssetConstSP>::iterator iter2 =
			creditAssetMap.begin();
		for (; iter2 != creditAssetMap.end(); ++iter2)
		{
			CreditAsset* creditAsset = const_cast<CreditAsset*>
				((iter2->second).get());

			creditAssetWrapperArray.push_back(CreditAssetWrapperSP(
				new CreditAssetWrapper(creditAsset)));
		}

		//convert the beta map into a vector
		vector<double> betaVector;
		map<string, double>::iterator iter1 =
			portfolioNameBetaMap.begin();

		for (; iter1 != portfolioNameBetaMap.end(); ++iter1)
			betaVector.push_back(iter1->second);

		//set the CreditAssetArray to the engine(copula)
		cop->addNames(creditAssetWrapperArray);

		//set the betas  to the engine(copula)
		cop->updateBetas(betaVector);

		//update the engine

		simulation = cop->simulation(
			engineTimeline,
			nSamples);

		simulation.UpdateSim();

		int size = simulation.getNumberRands();

        CDoubleMatrix corrMatrix(size, size);

        for (i = 0; i < size; i++)
			corrMatrix[i][i] = 1.0;

		// ugly hack to config the default sampling (Gaussian) JCP
		GaussSampling* gs = dynamic_cast<GaussSampling*>(mcPathConfig->sampling.get());

		if (gs)
			gs->config(nSamples);

		MCPathBaseCCM support(simulation.getNumberRands());

		randomGen = IMCRandomSP( new StratMCRandomNoCache(
						MCRandomNoCache(
							0,
							DependenceSP(new Gauss(corrMatrix)),
							mcPathConfig->getRandomGenerator(),
							false, // JCP
							1,
							size,
							0),
						mcPathConfig->sampling,
						nSamples,
						size)
					);

		//setting the random generator to the simulation object
		simulation.setRandoms(randomGen);
	}
	catch(exception& e)
	{
		throw ModelException(e, "MCPathConfigCCM::constructor");
	}
}

//from the parent, IMCPathConfig
bool MCPathConfigCCM::vegaMatrixSupported() const
{
	return false;
}

//from the parent, IMCPathConfig
void
MCPathConfigCCM::saveRandomGeneratorState()
{
	initRandState = IRandom::StateSP(randomGen->getState());
}

//from the parent, IMCPathConfig
LRGenerator*
MCPathConfigCCM::createLRGenerator(
	const MCPathGeneratorSP& pathGen,
	int                      nbIter,
	int                      nbSubSamples)
{
	throw ModelException(
			"Method not defined",
			"MCPathConfigCCM::createLRGenerator");

	return 0;
}

//from the parent, IMCPathConfig
ObjectArraySP
MCPathConfigCCM::getParameters()
{
	throw ModelException(
			"Method not defined",
			"MCPathConfigCCM::getParameters");
}

CClassConstSP const MCPathConfigCCM::TYPE =
	CClass::registerClassLoadMethod(
		"MCPathConfigCCM",
		typeid(MCPathConfigCCM),
		load);

void
MCPathConfigCCM::load(CClassSP& clazz)
{
    clazz->setPublic();

	REGISTER(MCPathConfigCCM, clazz);

	SUPERCLASS(CObject);

	IMPLEMENTS(IMCPathConfig);

	FIELD(randomGen, "Random Number Generator");

	FIELD(timelineFreq, "Frequency for the timeline = 1D 1W 1M 3M etc");
    FIELD_MAKE_OPTIONAL (timelineFreq);

	FIELD(sampling, "Sampling scheme");
	FIELD_MAKE_OPTIONAL(sampling);

	FIELD_NO_DESC       (timelineFreqObj);
    FIELD_MAKE_TRANSIENT(timelineFreqObj);

	EMPTY_SHELL_METHOD(defaultMCPathConfigCCM);
}

bool MCPathConfigCCMLoad() {
    return MCPathConfigCCM::TYPE != NULL;
}

// not sure if this canbe avoided by deriving from another class,
// not MCPathConfig
bool MCPathConfigCCM::carefulRandoms() const
{
	return true;
};


MCPathGeneratorSP MCPathConfigCCM::makePathGenerator(
	bool                               cachingRequested,
    int                                numPaths,
    const MCPathGeneratorSP&           pastPathGenerator,
    const IMCProduct*                   prod,
    Control*                           control,
    Results*                           results,
    DateTimeArray&                     simDates)
{
	const MCProductClient* prodClient = dynamic_cast<const MCProductClient*>(prod);

	// this will have the relevant simulation parameters

	MCPathGeneratorSP temp =
		MCPathGeneratorSP (new MCPathConfigCCM::PathGen(
				ccmCopula.get(),
				this,
				prodClient,
				numPaths,
				1));

	return temp;
};


IObject*
MCPathConfigCCM::defaultMCPathConfigCCM()
{
	return new MCPathConfigCCM(TYPE);
}

bool MCPathConfigCCMLinkIn(){
    return (MCPathConfigCCM::TYPE != NULL);
}

DRLIB_END_NAMESPACE
