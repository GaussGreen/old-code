//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PickAverage.cpp
//
//   Description : Option on Average in which the samples are picked according
//                 to a set of rules (Pick Types) specified:
//					PickType 1 -- if the sample is ABOVE_INI_FLOOR
//					PickType 2 -- if the sample is ABOVE_PREV_SAMPLE
//					PickType 3 -- if the sample is ABOVE_PREV_PICKED
//					PickType 4 -- if the sample is ABOVE_AVG_SO_FAR
//
//   Date        : May 2003
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/Maths.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/Reprice.hpp"
#include "edginc/MonteCarlo.hpp"

// For PickRefTypes
#define	  INI_FLOOR		1
#define	  PREV_SAMPLE	2
#define	  PREV_PICKED	3
#define	  AVG_SO_FAR	4


DRLIB_BEGIN_NAMESPACE

/** PickAverage product - option based upon the average samples picked */
class PickAverage: public GenericNFBase, 
                   virtual public IMCIntoProduct {
protected:
    /// fields ////////
    int							avgPickType;		// 
    double						initialFloor; 
    bool						isPickAbove;			// pick above or below

    DateTimeArray				averageOutDates;
	DoubleArray					weights;
	string						weightType;
	bool						isCall;
	double						strike;

public:
    static CClassConstSP const TYPE;
    friend class PickAverageMC;

    // validation
    void validatePop2Object(){
        static const string routine("PickAverage::validatePop2Object");
        GenericNFBase::validatePop2Object();

        // check that we've got as many weights as there are assets 
        // and if % weights that they sum to 100%
        if (weightType=="P") {
            AssetUtil::checkWeights(weights, assets->NbAssets());
        } else if (weightType=="U") {
            if (assets->NbAssets() != weights.size()){
                throw ModelException(routine,
                                     "Different number of assets ("+
                                     Format::toString(assets->NbAssets())+")"
                                     " to weights ("+
                                     Format::toString(weights.size())+")");
            }
        } else {
            throw ModelException(routine, "PickAverage: weightType must be U or P,"
                                 " but " + weightType + " given");
        }
        if (Maths::isNegative(strike)){
            throw ModelException(routine, "strike ("+
                                 Format::toString(strike)+") is negative");
        }
        // validate dates are not empty - order is handled by SimSeries
        if (averageOutDates.empty()) {
            throw ModelException(routine, "No averageOutDates supplied!");
        }

        if( !(avgPickType == INI_FLOOR || avgPickType == PREV_SAMPLE ||
              avgPickType == PREV_PICKED || avgPickType == AVG_SO_FAR) ){
            throw ModelException(routine, "avgPickType must be 1-4!");
        }
    }
   
    // returns the array of all sampling dates the instrument will ever need
    // excluding (possibly) the ref level dates (handled in GenericNFBase)
    const DateTimeArray samplingDates() const {
        return averageOutDates;
    }

    /** Implementation of MonteCarlo::IntoProduct interface - the
        implementation of this is below */
    virtual IMCProduct* createProduct(
        const MonteCarlo* model) const; // see below

private:
    PickAverage(): GenericNFBase(TYPE) {} // for reflection
    PickAverage(const PickAverage& rhs); // not implemented
    PickAverage& operator=(const PickAverage& rhs); // not implemented

    static IObject* defaultPickAverage(){
        return new PickAverage();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(PickAverage, clazz);
        SUPERCLASS(GenericNFBase);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultPickAverage);
        FIELD(avgPickType,       "Avg Pick Type:  1 - 4");
        FIELD(initialFloor,      "as % of the Basket Reference Level");
        FIELD(averageOutDates,   "Avg Dates");
        FIELD(weightType,		"Basket weight type:  P or U");
        FIELD(weights,			"Basket weights");
        FIELD(isCall,			"TRUE (Call), FALSE (Put)");
        FIELD(isPickAbove,		"TRUE (pick above), FALSE (pick below)");
        FIELD(strike,			"strike");
        
		clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

/* MC product class for PickAverage - with isRefAvgOut false */
class PickAverageMC : public IMCProduct,
					  virtual public IMCProductLN {
private:
    const PickAverage*      inst;         // reference to original instrument
    int                     nbAssets;     // convenient
    bool                    isUnitBasket;
	double					AboveOrBelow; // +1 if pick above; -1 if below

	// historical values from payoff
	double					baskRef;		 // initial basket ref
	double					baskSumOutSoFar;
    double                  comparisonLevelSoFar; // the reference level for pick decision
	long					nbPickedSoFar;

protected:

public:   
    /** equivalent to InstIntoMCProduct. Need to call parent's constructor */
    PickAverageMC(const PickAverage*            inst,
	              const SimSeriesSP&       simSeries):
        IMCProduct(inst->assets.get(),
                  inst->valueDate,
                  inst->discount.get(),
                  inst->refLevel,
                  simSeries,
                  inst->pastValues,
                  inst->instSettle.get(),
                  simSeries->getLastDate()),
	    inst(inst),
        nbAssets(getNumAssets()),
		isUnitBasket(inst->weightType == "U"),
		AboveOrBelow(inst->isPickAbove ? +1.0 : -1.0),
        baskRef(0.0),
		baskSumOutSoFar(0.0),
		comparisonLevelSoFar(inst->initialFloor),
		nbPickedSoFar(0) {
		}

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen)const{
        // empty
    }

    /** Called within the simulation loop */
    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&                prices) {
        static const string routine("PickAverageMC::payoff");
        try {
            int		beginIdx = pathGen->begin(0); // same for all assets
            int		endIdx   = pathGen->end(0);
			
			long	nbPicked		= nbPickedSoFar;
			double	baskSumOut		= baskSumOutSoFar;
			double	comparisonLevel = comparisonLevelSoFar;

			int		iAsset			= 0;

            // Set Basket Ref for "U" type.  Only do once
			if ((Maths::isZero(baskRef)) && (isUnitBasket)) {
				for (iAsset=0; iAsset<nbAssets; iAsset++) {
					baskRef += inst->weights[iAsset] * pathGen->refLevel(iAsset, 0/*iPath*/);
				}
			}
				
            // 
            for (int iStep=beginIdx; iStep<endIdx; iStep++) {
				// Calculate Basket level, in %
				double baskLevel = 0.0;
				if (isUnitBasket) {
					// "U type, baskLevel is in % of baskRef
					for (iAsset=0; iAsset<nbAssets; iAsset++) {
						baskLevel += inst->weights[iAsset] * pathGen->Path(iAsset, 0/*iPath*/)[iStep];
					}
					baskLevel /= baskRef;
				}
				else {
					// "P" type
					for (iAsset=0; iAsset<nbAssets; iAsset++) {
						baskLevel += inst->weights[iAsset] * pathGen->Path(iAsset, 0/*iPath*/)[iStep] 
									/ pathGen->refLevel(iAsset, 0/*iPath*/);
					}
				}

				// Check picks
				if (Maths::isPositive(AboveOrBelow * (baskLevel - comparisonLevel))) {
					// The sample is picked
					baskSumOut += baskLevel;
					nbPicked ++;
					
					// Reset comparisonLevel for ABOVE_PREV_PICKED or ABOVE_AVG_SO_FAR
					if (inst->avgPickType == PREV_PICKED) {
						comparisonLevel = baskLevel;
					}
					else if (inst->avgPickType == AVG_SO_FAR) {
						comparisonLevel = baskSumOut / nbPicked;
					}
				}

				// Reset comparisonLevel for PREV_SAMPLE, 
				// reguardless whether or not the sample is picked
				if (inst->avgPickType == PREV_SAMPLE) {
					comparisonLevel = baskLevel;
				}
				// For ABOVE_INI_FLOOR, do nothing
			
			}

			if (pathGen->doingPast()) {  // preserve values
				nbPickedSoFar = nbPicked;
				baskSumOutSoFar = baskSumOut;
				comparisonLevelSoFar = comparisonLevel;
			}
            
			// Price:  Perf set to strike when nbPicked = 0
			if (nbPicked > 0) {
				double perf = (nbPicked == 0) ? inst->strike : (baskSumOut / nbPicked);
				double myPayoff = inst->isCall ? (perf - inst->strike): (inst->strike - perf);
	            prices.add(inst->notional * prices.maxWithZero(myPayoff));
			} else {
	            prices.add(0.0);
			}
        }
        catch (exception& e) {
            throw ModelException(e, routine);
        }
    }

    // for the LogNormal path generator
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
                                    int                     iAsset) const {
		static const string routine = "PickAverage::getVolInterp";
        CVolRequestLNArray reqarr(1); // one interp level/path per asset here

        const IRefLevel* refLevel = getRefLevel();
        const DateTime& startDate = refLevel->getAllDates().front();
        const DateTime& today = getToday();
        bool  fwdStarting = startDate.isGreater(today);
        double interpLevel;

        if (fwdStarting){
            interpLevel = inst->strike;
        } else {
            /* not forward starting - some samples have fixed already
            (this includes averaging in) */
            int numDates = inst->averageOutDates.size();
            int numRemaining = 
                today.numFutureDates(inst->averageOutDates);
            // moneyness is from basket levels
            interpLevel = (numDates * inst->strike - baskSumOutSoFar)/ numRemaining;
            interpLevel *= pathGen->refLevel(iAsset, 0);
		}              

		const SimSeries* simSeries = getSimSeries();
        reqarr[0] = CVolRequestLNSP(new LinearStrikeTSVolRequest(interpLevel,
                                                                     startDate,
                                                                     simSeries->getLastDate(),
                                                                     fwdStarting));
        return reqarr;
    }
};

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* PickAverage::createProduct(const MonteCarlo* model) const {
    // we need to create a SimSeries object which says which assets need
    // which dates to be simulated
    SimSeriesSP simSeries(new SimSeries(assets->NbAssets())); /* create empty
                                                                 one */
    simSeries->addDates(averageOutDates);
    return new PickAverageMC(this, simSeries);
}

CClassConstSP const PickAverage::TYPE = CClass::registerClassLoadMethod(
    "PickAverage", typeid(PickAverage), PickAverage::load);

// * for class loading (avoid having header file) */
bool PickAverageLoad() {
    return (PickAverage::TYPE != 0);
}

DRLIB_END_NAMESPACE
