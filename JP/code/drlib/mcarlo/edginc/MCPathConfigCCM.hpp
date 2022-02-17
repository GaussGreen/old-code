//----------------------------------------------------------------------------
//
// Group       : CH Quantitative Research
//
// Filename    : MCPathConfig.hpp
//
// Description : Represents Path Config that would generate PathGenerators that
//               price an instrument using CDO
//
// Date        : Aug 2006
//
//----------------------------------------------------------------------------

#ifndef DR_MCPATHCONFIGCCM_HPP
#define DR_MCPATHCONFIGCCM_HPP

#include "edginc/MCRandom.hpp"
#include "edginc/MCPathConfig.hpp"
#include "edginc/CCMCopulaModel.hpp"
#include "edginc/SingleCreditAsset.hpp"
#include "edginc/PortfolioName.hpp"
#include "edginc/GaussSampling.hpp"
#include "edginc/AdaptiveSampling.hpp"



DRLIB_BEGIN_NAMESPACE


// after MCPathConfigLN
class MCARLO_DLL MCPathConfigCCM :	public CObject,
									public IMCPathConfig
{

private:

	CCMCopulaModelSP	ccmCopula;
	IRandomSP			randomGen;
	IRandom::StateSP    initRandState; // transient, not registered $unregistered
	string				timelineFreq;

//transient
	MaturityPeriodConstSP timelineFreqObj;

public:

	static CClassConstSP const TYPE;

	// a path generator, defined in the .cpp. After MCPathConfigLN::Generator
	class PathGen;
	friend class PathGen;

	DECLARE(PathGen);

	MCPathConfigCCM(CCMCopulaModelSP copula, CClassConstSP clazz = TYPE) :
		CObject(clazz),
		ccmCopula(copula),
		sampling(IMFSamplingSP(NULL))
		{};

	MCPathConfigCCM(CClassConstSP clazz = TYPE) :
		CObject(clazz),
		ccmCopula(CCMCopulaModelSP( new CCMCopulaModel())),
		timelineFreq(""),
		sampling(IMFSamplingSP(NULL))
		{};

	MCPathGeneratorSP makePathGenerator(
        bool                               cachingRequested,
        int                                numPaths,
        const MCPathGeneratorSP&           pastPathGenerator,
        const IMCProduct*                  prod,
        Control*                           control,
        Results*                           results,
        DateTimeArray&                     simDates
	);

	virtual ~MCPathConfigCCM() {}

	virtual int storagePerPath(IMCProduct* product) const; //from the parent, IMCPathConfig

	virtual MCPathGeneratorSP pastPathGenerator(const IMCProduct*);//from the parent, IMCPathConfig

	IRandomSP	getRandomGenerator() const;

	// would also take in sampling information, including the list of dates
	// in as much as the product has model information, that would be in the product as well
	// would this also discern whether the path generator we create does indexDates or plain dates?
    virtual MCPathGeneratorSP futurePathGenerator(	int                                cachingMode,
													int                                numPaths,
													const MCPathGeneratorSP&           pastPathGenerator,
													const IMCProduct*                  prod,
													Control*                           control,
													Results*                           results,
													DateTimeArray&                     simDates ); //from the parent, IMCPathConfig



	//from the parent, IMCPathConfig
	virtual MarketDataFetcherSP marketDataFetcher() const;

	//from the parent, IMCPathConfig
	virtual void validatePop2Object();

	//from the parent, IMCPathConfig
	virtual void getMarket(const IModel*      model,
                        const MarketData*  market,
                        IInstrumentCollectionSP instrument);

	virtual MarketObjectSP getMarket(
		const CModel        *model,
        const MarketData    *market,
		const string        &name,
        const CClassConstSP &type) const ;//from the parent, IMCPathConfig

	virtual bool vegaMatrixSupported() const;//from the parent, IMCPathConfig

	virtual void saveRandomGeneratorState();//from the parent, IMCPathConfig

	//from the parent, IMCPathConfig
	virtual LRGenerator* createLRGenerator(
		const MCPathGeneratorSP& pathGen,
		int                      nbIter,
		int                      nbSubSamples);

	virtual ObjectArraySP getParameters();//from the parent, IMCPathConfig

	bool carefulRandoms() const;

	CCMCopulaModelSP getCopula() const {return ccmCopula;};

	// here we catch all calls for market data, and build the models'
	// (ccmCopula) list of names
	void getComponentMarketData(
		const CModel*         model,
		const MarketData*     market,
		MarketObjectSP        mo,
		MarketDataFetcher*    mdf) const;

	// Generates an engineTimeline from the product timeline by making sure
	// that no two adjacent product timeline dates differ by the timeline
	// frequency supplied.
	virtual void generateEngineTimeline(
		const DateTimeArray& productTimeline, // input
		DateTimeArray& engineTimeline,	//output
		IntArray& fromEngineIndexToProductIndex //output
		) const;

protected:

	static void load(CClassSP& clazz);

	static IObject* defaultMCPathConfigCCM();

	IMFSamplingSP sampling;

public:

	virtual IModel::WantsRiskMapping wantsRiskMapping() const;

	IMFSamplingSP getSampling() {return sampling;};

	void setSampling(IMFSamplingSP newSampling)
	{
		sampling = newSampling;
	};

};

typedef std::map<string, ICreditLossConfigIndexedSVSP> CreditLossConfigIndexedSVMap;

// Path generator, modelled after MCPathConfigLN::PathGenSpot
class MCPathConfigCCM::PathGen :  public MCPathGenerator,
						          public virtual VirtualDestructorBase
{
public:

	friend class MCPathConfigCCM;

	// takes in a CCMCopulaModel and some simulation parameters,
	// initializes a simulator

	const CCMCopulaModel::Simulation &Simulation() const {return simulation;};

	PathGen(
			CCMCopulaModel *cop,
			MCPathConfigCCM *mcPathConfig,
			const MCProductClient*   prodClient,
			int nSamples,
			int nFactors);

	PathGen(const MCProductClient*   prodClient);

	virtual void generatePath(int nPath);

	bool hasPast() const {return false;};

	bool doingPast() const {return false;};

	IStateVariableSP create(const IStateVariableGen *svGen);

	// for compatibility with MC

    virtual int NbSimAssets() const { return 0; }

    virtual const double* Path(int iAsset, int iPath) const { return 0; }

    virtual double refLevel(int iAsset, int iPath) const { return 0; }

    virtual double maxDriftProduct(int iAsset) const { return 0; }

    virtual int begin(int iAsset) const { return 0; }

    virtual int end(int iAsset) const { return 0; }

	int getPathIndex() const {return nowPathIdx;};

private:
	CCMCopulaModel::Simulation		simulation;
	MCPathConfigCCM*				mcPathConfig;
	
	/* This vector will be the same size as the engine timeline.
	   It will be populated by an int which would be the corresponding product timeline index;
	*/
	IntArray						fromEngineIndexToProductIndex; 


	/** this vector has been created so that there is a 1:1 correspondence of names with the PathGenerator
		and we avoid a lookup inside a data structure like a map
	*/
	vector<PortfolioNameIndexedSVSubResultArraySP>	portfolioNameDBase;
    StateVarDBase                    svDBase;        //!< Collection of Generators + statevars
    int                              nowPathIdx;     //!< Current path idx

protected:

	IMCRandomSP               randomGen;

};

DRLIB_END_NAMESPACE

#endif
