//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : LadderAverage.cpp
//
//   Description : Ladder Average option with partial lock-in
//
//   Date        : October 2003
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
DRLIB_BEGIN_NAMESPACE

/** LadderAverage product - option based upon the average samples picked */
class LadderAverage: public GenericNFBase, 
                   virtual public IMCIntoProduct {
protected:
    /// fields ////////
    DateTimeArray				monitoringDates;
	DateTimeArray				averageOutDates;
	DoubleArray					weights;
	string						weightType;
	
	double						strike;
	double						cap;
	DoubleArray					ladderLevels;
	DoubleArray					lockInAmount;

public:
    static CClassConstSP const TYPE;
    friend class LadderAverageMC;

    // validation
    void validatePop2Object(){
        static const string routine("LadderAverage::validatePop2Object");
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
            throw ModelException(routine, "LadderAverage: weightType must be U or P,"
                                 " but " + weightType + " given");
        }
        
		// validate dates are not empty - order is handled by SimSeries
        if (averageOutDates.empty()) {
            throw ModelException(routine, "No AverageOut Dates supplied!");
        }
        if (monitoringDates.empty()) {
            throw ModelException(routine, "No Monitoring Dates supplied!");
        }

		// strike/cap must be positive
        if (Maths::isNegative(strike)) {
            throw ModelException(routine, "strike ("+
                                 Format::toString(strike)+") is negative");
        }

		// check the ladder levels
        if ((ladderLevels.empty()) || (lockInAmount.empty())) {
            throw ModelException(routine, "No ladder level or lock-in amount supplied!");
        }
        if ((ladderLevels.size()) != (lockInAmount.size())) {
            throw ModelException(routine, "The number of lock-in amounts must be the same as that of ladder levels!");
        }
		// 1st ladder level must be higher than the strike
        if (Maths::isPositive(strike - ladderLevels[0])) {
            throw ModelException(routine, "Ladder levels should be higher than strike.");
        }
		// Ladder Levels in increasing order
		for (int iItem = 1; iItem < ladderLevels.size(); iItem++) {
	        if (Maths::isPositive(ladderLevels[iItem-1] - ladderLevels[iItem])) {
		        throw ModelException(routine, "Ladder levels must be in increasing order.");
			}
	        if (Maths::isPositive(lockInAmount[iItem-1] - lockInAmount[iItem])) {
		        throw ModelException(routine, "Lock-In Amount must be in increasing order.");
			}
		}
        
		// Cap is above the 
		if (Maths::isNegative(cap - strike - lockInAmount[ladderLevels.size() - 1])) {
            throw ModelException(routine, "cap level must be more than (strike + highest lockIn amount).");
        }
    }
   
    // returns the array of all sampling dates the instrument will ever need
    // excluding (possibly) the ref level dates (handled in GenericNFBase)
    const DateTimeArray samplingDates() const {
        return DateTime::merge(averageOutDates, monitoringDates);
    }

    /** Implementation of MonteCarlo::IntoProduct interface - the
        implementation of this is below */
    virtual IMCProduct* createProduct(
        const MonteCarlo* model) const; // see below

private:
    LadderAverage(): GenericNFBase(TYPE) {} // for reflection
    LadderAverage(const LadderAverage& rhs); // not implemented
    LadderAverage& operator=(const LadderAverage& rhs); // not implemented

    static IObject* defaultLadderAverage(){
        return new LadderAverage();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(LadderAverage, clazz);
        SUPERCLASS(GenericNFBase);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultLadderAverage);
        FIELD(monitoringDates,   "Monitoring Dates");
        FIELD(averageOutDates,   "AverageOut Dates");
        FIELD(weightType,		"Basket weight type:  P or U");
        FIELD(weights,			"Basket weights");
        FIELD(strike,			"strike");
        FIELD(cap,				"cap");
        FIELD(ladderLevels,		"Ladder Levels");
        FIELD(lockInAmount,		"Lock-In amount");

		clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

/* MC product class for LadderAverage */
class LadderAverageMC: public IMCProduct,
					   virtual public IMCProductLN {
private:
    const LadderAverage*    inst;			// reference to original instrument
    int                     nbAssets;		// convenient
    bool                    isUnitBasket;

	double					baskRef;		// initial basket ref
	int						nbAverages;		// no. of averages
	int						nbLadders;		// no. of ladders
	int						nbSteps;		// no. of sim steps (incl. average and monitoring)

    // for mapping
	IntArray                averageMap;   // to track averaging dates
    IntArray                monitorMap;   // to track monitoring dates

	// historical values from payoff
	double					baskSumOutSoFar;
    int		                nbLadderReachedSoFar;

	// working variables
	DoubleArray				baskLevels;

protected:

public:   
    /** equivalent to InstIntoMCProduct. Need to call parent's constructor */
    LadderAverageMC(const LadderAverage*     inst,
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
        baskRef(0.0),
		nbAverages(inst->averageOutDates.size()),
		nbLadders(inst->ladderLevels.size()),
		nbSteps(simSeries->numDates(0)),
		baskSumOutSoFar(0.0),
		nbLadderReachedSoFar(0),
		baskLevels(nbSteps,0.0) {
			// mapping
			bool isTrivial;
			averageMap = DateTime::createMapping(simSeries->getAllDates(),
				                             inst->averageOutDates,
					                         isTrivial);
	        monitorMap = DateTime::createMapping(simSeries->getAllDates(),
		                                     inst->monitoringDates,
			                                 isTrivial);
		}

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen)const{
        // empty
    }

    /** Called within the simulation loop */
    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&                prices) {
        static const string routine("LadderAverageMC::payoff");
        try {
            int		beginIdx = pathGen->begin(0); // same for all assets
            int		endIdx   = pathGen->end(0);
			
			long	nbLadderReached = nbLadderReachedSoFar;
			double	baskSumOut		= baskSumOutSoFar;

			int		iAsset			= 0;
			int		iStep			= 0;
			double	finalStrike		= inst->strike; // to be updated based on final ladder reached
			double	finalLockIn		= 0.0;			// to be updated based on final ladder reached
			double	finalAvg		= 0.0;

			// Set Basket Ref for "U" type.  Only do once
			if ((Maths::isZero(baskRef)) && (isUnitBasket)) {
				for (iAsset=0; iAsset<nbAssets; iAsset++) {
					baskRef += inst->weights[iAsset] * pathGen->refLevel(iAsset, 0/*iPath*/);
				}
			}
				
            // basket levels
            for (iStep = beginIdx; iStep < endIdx; iStep++) {
				// Calculate Basket level, in %
				baskLevels[iStep] = 0.0;
				if (isUnitBasket) {
					// "U type, baskLevel is in % of baskRef
					for (iAsset=0; iAsset<nbAssets; iAsset++) {
						baskLevels[iStep] += inst->weights[iAsset] * pathGen->Path(iAsset, 0/*iPath*/)[iStep];
					}
					baskLevels[iStep] /= baskRef;
				}
				else {
					// "P" type
					for (iAsset=0; iAsset<nbAssets; iAsset++) {
						baskLevels[iStep] += inst->weights[iAsset] * pathGen->Path(iAsset, 0/*iPath*/)[iStep] 
											 / pathGen->refLevel(iAsset, 0/*iPath*/);
					}
				}
			}

	        // baskSumOut
		    for (iStep = beginIdx + averageMap[beginIdx]; iStep < endIdx; 
			     iStep++, iStep += averageMap[iStep]) {
				baskSumOut += baskLevels[iStep];
	        }
		    
	        // check ladders
		    for (iStep = beginIdx + monitorMap[beginIdx]; iStep < endIdx; 
			     iStep++, iStep += monitorMap[iStep]) {
				while ((nbLadderReached < nbLadders) &&
					   (Maths::isPositive(baskLevels[iStep] - inst->ladderLevels[nbLadderReached]))) {
					nbLadderReached++;
				} 
	        }

			if (pathGen->doingPast()) {  
				// information updated by past values
				nbLadderReachedSoFar = nbLadderReached;
				baskSumOutSoFar = baskSumOut;
			}
			else
			{
				// Price: 
				if (nbLadderReached > 0) {
					finalStrike += inst->lockInAmount[nbLadderReached-1];
					finalLockIn  = inst->lockInAmount[nbLadderReached-1];
				}
				finalAvg = baskSumOut / nbAverages;
				double myPayoff = finalLockIn + prices.maxWithZero(finalAvg - finalStrike)
											  - prices.maxWithZero(finalAvg - inst->cap);
			    prices.add(inst->notional * myPayoff);
			}
        }
        catch (exception& e) {
            throw ModelException(e, routine);
        }
    }

    // for the LogNormal path generator
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
                                    int                     iAsset) const {
		static const string routine = "LadderAverage::getVolInterp";
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
IMCProduct* LadderAverage::createProduct(const MonteCarlo* model) const {
    // we need to create a SimSeries object which says which assets need
    // which dates to be simulated
    SimSeriesSP simSeries(new SimSeries(assets->NbAssets())); /* create empty
                                                                 one */
    simSeries->addDates(averageOutDates);
    simSeries->addDates(monitoringDates);
    return new LadderAverageMC(this, simSeries);
}

CClassConstSP const LadderAverage::TYPE = CClass::registerClassLoadMethod(
    "LadderAverage", typeid(LadderAverage), LadderAverage::load);

// * for class loading (avoid having header file) */
bool LadderAverageLoad() {
    return (LadderAverage::TYPE != 0);
}

DRLIB_END_NAMESPACE
