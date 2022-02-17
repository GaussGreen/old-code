//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CollarGamma.cpp
//
//   Description : Average Collar with Gamma adjustement
//
//   Author      : Bruno O Melka
//
//   Date        : 10 Oct 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Generic1Factor.hpp"
#include "edginc/ModelException.hpp"
#include "edginc/Maths.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/Format.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/OutputRequestUtil.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/Asset.hpp"
#include "edginc/SVGenSpot.hpp"
#include "edginc/CashSettlePeriod.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/MCPathBase.hpp"
#include "edginc/SampleList.hpp"
#include "edginc/Black.hpp"
#include "edginc/EquityCache.hpp"


DRLIB_BEGIN_NAMESPACE

//** set vol  */
static void setVol(CAsset* asset, double vol)
{
    EquityCache* assetCache = dynamic_cast<EquityCache*>(asset);
    if (assetCache) // this is american case which has cache for now
        assetCache->setVol(vol);
    else
    {
        // temp solution for euro env
        EquityBase* eq = dynamic_cast<EquityBase*>(asset);
        VolLevel shift(vol);
        VolLevel::Shift* tweak = dynamic_cast<VolLevel::Shift*>(const_cast<CVolBase*>(eq->getVol().get()));
        tweak->sensShift(&shift);
    }
}

/** Compute the Black Scholes Gamma foer an average collar*/
static double blackCollarGamma (double spot,
								double fwd,
								double strikeCall,
								double strikePut,
								double variance) {

	if (!Maths::isPositive(variance) || !Maths::isPositive(fwd)) {
		return 0.0;
	}
	else {
		double sqrtVar = sqrt(variance);            
		double nCall = Maths::isPositive(strikeCall) ? N1Density((log(fwd / strikeCall) + (0.5 * variance)) / sqrtVar) : 0.0;
		double nPut = Maths::isPositive(strikePut) ? N1Density((log(fwd / strikePut) + (0.5 * variance)) / sqrtVar) : 0.0;
		return (nPut - nCall) / (spot * sqrtVar);
	}
}

/** Compute inputs for gamma calculation */
static double avgCollarGamma(double strikeCall,
							   double strikePut,
							   const  SampleListSP avgOut,
							   const  DateTime& valueDate,
							   CAsset* asset,
							   double vol,
							   double weight) {
						 
	setVol(asset, vol);
    ATMVolRequestSP atm(new ATMVolRequest());
    CVolProcessedBSSP volBS(asset->getProcessedVol(atm.get()));
	double avgFwd = avgOut->futureSampleSum(asset, valueDate);
    double sumSoFar = avgOut->sumToDate(valueDate);
    double avgVariance = avgOut->averageVariance(volBS.get(), valueDate, false);
	return blackCollarGamma(sumSoFar * weight, avgFwd, strikeCall - sumSoFar, strikePut - sumSoFar, avgVariance);
}


/** CollarGamma class */ 
class CollarGamma: public Generic1Factor,
                public virtual LastSensDate,
                public virtual MonteCarlo::IIntoProduct {
public:
    static CClassConstSP const TYPE; 

    virtual void validatePop2Object(){
        static const string method = "CollarGamma::validatePop2Object";
    }

    virtual void Validate() {
        static const string method = "CollarGamma::Validate";
    }

	DateTime endDate(const Sensitivity* sensControl) const {
		return monitorDates[monitorDates.size() - 1];
	}

	DateTime getValueDate() const {
		return valueDate;
	}

    /** Implementation of MonteCarlo::IntoProduct interface */
    IMCProduct* createProduct(const MonteCarlo* model) const; // see below

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CollarGamma, clazz);
        SUPERCLASS(Generic1Factor);
        IMPLEMENTS(LastSensDate);
        IMPLEMENTS(MonteCarlo::IIntoProduct);
        EMPTY_SHELL_METHOD(defaultCollarGamma);
		FIELD(strikePut, "strikePut");
		FIELD(strikeCall, "strikeCall");
		FIELD(monitorDates, "averaging schedule");
		FIELD(pastValues, "past values");
		FIELD(gammaScaling, "use Gamma scaling or not");
		FIELD(gammaWidth, "gammaWidth");
		FIELD(gammaRange, "gammaRange");
		FIELD(gammaStart, "gammaStart");
    }
	 
    static IObject* defaultCollarGamma(){
        return new CollarGamma();
    }

private:

    CollarGamma():Generic1Factor(TYPE) {}; 
    CollarGamma(const CollarGamma& rhs);
    CollarGamma& operator=(const CollarGamma& rhs);

	class			MCSV;
	friend class	MCSV;
	class			MC;
	friend class	MC;

	double			strikePut;		// strike Put
	double			strikeCall;		// strike Call
	DateTimeArray	monitorDates;	// averaging schedule
	IPastValuesSP   pastValues;		// all historic spots

	bool			gammaScaling;	// use Gamma scaling or not
	double			gammaWidth;
	double			gammaRange;
	double			gammaStart;

};

CClassConstSP const CollarGamma::TYPE = CClass::registerClassLoadMethod(
    "CollarGamma", typeid(CollarGamma), CollarGamma::load);

/* MC product class for CollarGamma */
class CollarGamma::MC: public IMCProduct, virtual public IMCProductLN {

private:

	DateTimeArray	simulationDates;
	DoubleArray		values;
	DoubleArray		weights;
	DoubleArray		newPath;
	
	double			averageSoFar;
	double			strikePut;				
	double			strikeCall;
	double			discount;
	double			gammaWidth;
	double			gammaRange;
	double			gammaStart;
	bool			gammaScaling;
	CAssetSP		asset;

public:
    // for the LogNormal path generator
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGenerator,
                                    int                     iAsset) const {
        CVolRequestLNArray   reqarr(1); // one interp level/path per asset here
        reqarr[0] = CVolRequestLNSP(new ATMVolRequest());     
        return reqarr;
    }

	/** Use this opportunity to do any model driven initialisation
        of the instrument. e.g closed from barrier adjustments */
    virtual void initialiseLN(const MCPathGenerator* pathGenerator) const {}

	/* override the initial function. All discounting is done in payoff */
	virtual double pvFromPaymentDate() const{
        return discount;
    }

	/** compute total vol to step */
	void computeTotalVols(const DoubleMatrix& vols, DoubleArray& totalVols) {
		const double* pathVols = vols[0];
		double totalVar = 0.0;
		int num = vols.numRows();
		for (int i = 0; i < num; i++) {
			totalVar += pathVols[num - i - 1] * pathVols[num - i - 1];
			totalVols[num - i - 1] = sqrt(totalVar);
		}
	}

	/** estimate gammas from initial path */
	void computeGammas(DoubleArray& gammas, const IPathGenerator* pathGen, const DoubleArray& totalVols) {

		DoubleArray newValues(values);
		newValues.resize(simulationDates.size());

		for (int i = pathGen->begin(0) + 1, j = 1; i < pathGen->end(0); i++, j++) {
			newValues[i] = pathGen->Path(0,0)[i];
			SampleListSP newAvgOut(new SampleList(simulationDates, newValues, weights));
			gammas[i] = avgCollarGamma(strikeCall,
										strikePut, 
										newAvgOut, 
										simulationDates[i], 
										asset.get(), 
										totalVols[j], 
										weights.size() / double(i+1));
		}
	}

	/** compute new vols with gamma adjustment */
	void computeAdjustedVols(DoubleArray& newVols, const DoubleMatrix& vols, const DoubleArray& gammas, int start) {
		const double* pathVols = vols[0];
		for (int i = start; i < gammas.size(); i++) {
			double factor = (gammas[i] > 0) ? Maths::max(0.0, gammas[i] - gammaStart) : Maths::min(0.0, gammas[i] + gammaStart);
			factor = fabs(factor / gammaWidth);
			double scaling = 1.0/gammaRange + (1.0 - 1.0/gammaRange) * exp(-factor);
			scaling = (gammas[i] > 0) ? scaling : 1.0/scaling;
			newVols[i] = pathVols[i - start] * scaling;
		}
	}

	/** resimulate a new path with new vols */
	void reSimulatePath(DoubleArray& newPath, const DoubleMatrix& drifts, const DoubleMatrix& randoms, const DoubleArray& newVols, int start) {

		const double* pathDrift = drifts[0];
		const double* pathRandoms = randoms[0];

		for (int i = 0, modStep = start; modStep < newVols.size(); i++, modStep++) {
			double sigmaDtDz = newVols[modStep] * pathRandoms[i];           
			newPath[modStep] = (i == 0) ? values[start] : newPath[modStep-1];
			newPath[modStep] *= pathDrift[i] * exp(sigmaDtDz);
		}
	}

    /** equivalent to InstIntoMCProduct. Need to call parent's constructor */
    MC(const CollarGamma*          inst,
       const SimSeriesSP&       simSeries,
       InstrumentSettlementSP   instSettle):
	IMCProduct(inst->asset.get(),
				inst->valueDate,
				inst->discount.get(),
				IRefLevelSP(IRefLevel::Util::makeZero(inst->valueDate)), // fix!
				simSeries, 
				inst->pastValues,
				instSettle.get(),
				inst->monitorDates[inst->monitorDates.size() - 1]),
	simulationDates(inst->monitorDates),
	weights(simulationDates.size(), 1.0 / simulationDates.size()),
	newPath(simulationDates.size()),
	averageSoFar(0.0),
	strikePut(inst->strikePut),
	strikeCall(inst->strikeCall),
	discount(inst->discount->pv(inst->valueDate, simulationDates[simulationDates.size() - 1])),
	gammaWidth(inst->gammaWidth),
	gammaRange(inst->gammaRange),
	gammaStart(inst->gammaStart),
	gammaScaling(inst->gammaScaling),
	asset(dynamic_cast<CAsset*>(inst->asset.get()->clone())) {

		int numPast = inst->pastValues->getNbPastValues(0);
		values = DoubleArray(inst->pastValues->getPastValues(simulationDates, 0, simulationDates[numPast - 1]));
	}

    /** Called within the simulation loop */
    virtual void payoff(const IPathGenerator*	pathGen,
                        IMCPrices&					prices ) {

		// if gamma scaling, recompute the whole path with adjusted vols
		if (gammaScaling && !pathGen->doingPast()) {
			MCPathBase::Generator* gen = dynamic_cast<MCPathBase::Generator*>(const_cast<IPathGenerator*>(pathGen) );
			if(gen) {
				const MCGammaAdjustedPathGenerator* pathLN = dynamic_cast<const MCGammaAdjustedPathGenerator*>(gen->getPathBase().get());
				if (pathLN) {

					const DoubleMatrix& vols = pathLN->getVols(0);
					const DoubleMatrix& drifts = pathLN->getDrifts(0);
					const DoubleMatrix& randoms = pathLN->getRandoms(0);
				
					DoubleArray totalVols(vols.numRows());
					computeTotalVols(vols, totalVols);
					DoubleArray gammas(pathGen->end(0), 0.0);
					computeGammas(gammas, pathGen, totalVols);
					DoubleArray newVols(pathGen->end(0));
					computeAdjustedVols(newVols, vols, gammas, pathGen->begin(0));
					reSimulatePath(newPath, drifts, randoms, newVols, pathGen->begin(0));
				}
			}
		}

		// compute averages
		double average = averageSoFar;
		for (int iStep = pathGen->begin(0); iStep < pathGen->end(0); iStep++) {
			double spot = (!pathGen->doingPast() && gammaScaling) ? newPath[iStep] : pathGen->Path(0,0)[iStep];
			average = (iStep == 0) ? spot : (average * (iStep - 1) + spot) / iStep;
		}

		// compute average for past values
		if (pathGen->doingPast()) {
			averageSoFar = average;
		}

		// compute value of option.
		if (!pathGen->doingPast() || !hasFuture()) {
			double myPayoff = Maths::max((strikePut - average), 0.0) - Maths::max((average - strikeCall), 0.0); 
			prices.add(myPayoff);
		}
	}

};

/* MC product class for CollarGamma with state variables */
class CollarGamma::MCSV: public MCProductClient, virtual public IMCProductLN{

private:

    SVGenSpot::IStateVarSP		assetSV;		// asset state variable
    SVGenSpotSP				assetGen;		// generator for asset

	double	averageSoFar;
	double	strikePut;				
	double	strikeCall;
	double  discount;

protected:
    /** Override default method on IMCProduct. This method is called every time
        the path generator is changed (which is, at the moment, when the
        past path generator is created, and then when the future path
        generator is created  */
    virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen){
        assetSV = assetGen->getSpotSV(newPathGen);
    }

public:
    // for the LogNormal path generator
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGenerator,
                                    int                     iAsset) const {
        CVolRequestLNArray   reqarr(1); // one interp level/path per asset here
        reqarr[0] = CVolRequestLNSP(new ATMVolRequest()); 
        return reqarr;
    }

	/** Use this opportunity to do any model driven initialisation
        of the instrument. e.g closed from barrier adjustments */
    virtual void initialiseLN(const MCPathGenerator* pathGenerator)const {}

	/* override the initial function. All discounting is done in payoff */
	virtual double pvFromPaymentDate() const{
        return discount;
    }

	/** Appends 'true' (ie non derived) state variable generators
        required to the supplied collector.*/
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const{
        svCollector->append(assetGen.get()); 
    }

    /** equivalent to InstIntoMCProduct. Need to call parent's constructor */
    MCSV(const CollarGamma*          inst,
       const SimSeriesSP&       simSeries,
       InstrumentSettlementSP   instSettle):
	MCProductClient(inst->asset.get(),
					inst->valueDate,
					inst->discount.get(),
					IRefLevelSP(IRefLevel::Util::makeZero(inst->valueDate)), // fix!
					simSeries, // fix!
					inst->pastValues,
					instSettle.get(),
					inst->monitorDates[inst->monitorDates.size() - 1]),
	averageSoFar(0.0),
	strikePut(inst->strikePut),
	strikeCall(inst->strikeCall),
	discount(inst->discount->pv(inst->valueDate, inst->monitorDates[inst->monitorDates.size() - 1])),
	assetGen(new SVGenSpot(1, inst->monitorDates))	 { }

    /** Called within the simulation loop */
    virtual void payoff(const IPathGenerator*  	pathGen,
                        IMCPrices&					prices ) {

        const SVPath& pathAsset = assetSV->path(0);
		
		// compute averages
		double average = averageSoFar;
		for (int iStep = pathAsset.begin(); iStep < pathAsset.end(); iStep++) {
			double spot = pathAsset[iStep];
			average = (iStep == 0) ? spot : (average * (iStep - 1) + spot) / iStep;
		}

		// compute average for past values
		if(pathGen->doingPast()) {
			averageSoFar = average;
		}

		// compute value of option.
		if (!pathGen->doingPast() || !hasFuture()) {
			double myPayoff = Maths::max((strikePut - average), 0.0) - Maths::max((average - strikeCall), 0.0); 
			prices.add(myPayoff);
		}
    }

};

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* CollarGamma::createProduct(const MonteCarlo* model) const {
    // the simSeries passed to the IMCProduct is redundant
    // well not quite. The MC wants to know the last sim date
    SimSeriesSP simSeries(new SimSeries(1)); /* create empty one */
    simSeries->addDates(monitorDates);
    InstrumentSettlementSP instSettle(new CashSettlePeriod(0));
	if (model->stateVarUsed()){
		return new MCSV(this, simSeries, instSettle);
	}
	return new MC(this, simSeries, instSettle);
}

// for class loading 
bool CollarGammaLoad() {
    return (CollarGamma::TYPE != 0);
}

DRLIB_END_NAMESPACE
