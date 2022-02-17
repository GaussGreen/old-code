//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : BasketOfLookBacks.cpp
//
//   Description : Option on Basket of Lookbacks
//
//   Date        : Novemebr 2003
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/Maths.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/CliquetVolRequest.hpp"
#include "edginc/IAggregate.hpp"

DRLIB_BEGIN_NAMESPACE

/*****************************************************************************/

// BasketOfLookBacks product
class BasketOfLookBacks: public GenericNFBase, 
                       virtual public IMCIntoProduct{
private: 
    // Simulation
	DateTimeArray					AvgOutDates;
	DateTimeArray					couponDates;			// couponDates must be subset of monitoringDate
	// Baskets
    IAggregateMakerSP				assetBasket;
    IAggregateMakerSP				timeBasket;
    IDoubleArrayModifierMakerSP		timeBasketComponents;
	// Lookback
	bool							isLookbackOnBasket;		// TRUE = lookback on assetBasket on coupon date
															// FALSE = loolback on asset levels
	bool							isLookbackForMax;		// TRUE = lookback for maximum
															// FALSE = lookback for minimum
	bool							isLookbackFromStart;	// TRUE = lookback from start to current coupon date (incl)
															// FALSE = lookback from previous coupon date (excl) to current coupon date (incl)
	// Overall Call Option
    IDoubleArrayModifierMakerSP		overallOption;				
    bool							isCliquetStyle;
    bool							isRefPrevExtreme;		// only applicable, if isCliquetStyle = TRUE and isLookbackOnBasket = FALSE

public:
    static CClassConstSP const TYPE;
    friend class BasketOfLookBacksMC;

    // Validation
    void validatePop2Object(){
        static const string routine = "BasketOfLookBacks::validatePop2Object";
        GenericNFBase::validatePop2Object();
        if (AvgOutDates.empty()) {
            throw ModelException(routine, "No AvgOutDates given!");
        }
        if (couponDates.empty()) {
            throw ModelException(routine, "No couponDates given!");
        }
        // couponDates must be subset of AvgOutDates
        if (!DateTime::isSubset(AvgOutDates, couponDates)) {
            throw ModelException(routine, "couponDates must be a subset of AvgOutDates");
        }
 		// require that couponDates do not start before AvgOutDates
        const DateTime& firstMon = AvgOutDates[0];
        const DateTime& firstCpn = couponDates[0];
        if (firstCpn < firstMon) {
            throw ModelException(routine, "Cannot have first coupon date " + firstCpn.toString() + 
                                 " before first monitoring date " + firstMon.toString());
        }
        // require that AvgOutDates do not end after couponDates
        const DateTime& lastMon = AvgOutDates[AvgOutDates.size()-1];
        const DateTime& lastCpn = couponDates[couponDates.size()-1];
        if (lastMon > lastCpn) {
            throw ModelException(routine, "Cannot have final monitoring date " + lastMon.toString() + 
                                 " after final coupon date " + lastCpn.toString());
        }

        // if isRefPrevExtreme=TRUE, then we must have isLookbackOnBasket=FALSE
		if (isCliquetStyle && isRefPrevExtreme && isLookbackOnBasket) {
			throw ModelException(routine, "Cannot have isRefPrevExtreme=TRUE, if isLookbackOnBasket=TRUE");
		}
		// if isRefPrevExtreme=TRUE, then we must have isCliquetStyle=TRUE
		if (!isCliquetStyle && isRefPrevExtreme) {
			throw ModelException(routine, "Cannot have isRefPrevExtreme=TRUE, if isCliquetStyle=FALSE");
		}
    } // End of Validation
   
    // returns the array of all sampling dates the instrument will ever need
    // excluding (possibly) the ref level dates (handled in GenericNFBase)
    const DateTimeArray samplingDates() const {
        return AvgOutDates;
    }

    // Implementation of MonteCarlo::IntoProduct interface (details are below)
    virtual IMCProduct* createProduct(const MonteCarlo* model) const; 

private:
    BasketOfLookBacks(): GenericNFBase(TYPE) {} // for reflection
    BasketOfLookBacks(const BasketOfLookBacks& rhs); // not implemented
    BasketOfLookBacks& operator=(const BasketOfLookBacks& rhs); // not implemented

    static IObject* defaultBasketOfLookBacks(){
        return new BasketOfLookBacks();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(BasketOfLookBacks, clazz);
        SUPERCLASS(GenericNFBase);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultBasketOfLookBacks);
        FIELD(timeBasket,             "comment timeBasket");
        FIELD(timeBasketComponents,   "comment timeBasketComponents");
        FIELD(assetBasket,            "comment assetBasket");
        FIELD(overallOption,          "comment overallOption");
        FIELD(AvgOutDates,		"comment AvgOutDates");
        FIELD(couponDates,			"comment couponDates");
        FIELD(isLookbackOnBasket,	"isLookbackOnBasket");
        FIELD(isLookbackForMax,		"isLookbackForMax");
		FIELD(isLookbackFromStart,	"isLookbackFromStart");
        FIELD(isCliquetStyle,        "isCliqueStyle");
		FIELD(isRefPrevExtreme,		"isRefPrevExtreme");
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

// MC product class for BasketOfLookBacks
class BasketOfLookBacksMC: public IMCProduct,
						   virtual public IMCProductLN {
private:
    const BasketOfLookBacks*	inst;		// reference to original instrument
    int							nbAssets;	// number of assets, assignment via constructor
	
	int							iCoupon;	    // number of coupon dates
	DoubleArray					extreme;        // min/max for each asset [nbAssets]
    DoubleArray					refLevel;	    // refLevel for each asset [nbAssets]
    DoubleArray					refLevelTemp;	// refLevel for each asset [nbAssets]
    SimpleDoubleArray			coupons;	    // info for every coupon date [nbCoupons]

	SimpleDoubleArray			assetComps;	// info for assets before performance modifier (length = nbAssets)
	
    // Historical Values -- assignment either zero or from run through past
    int                      iCouponSoFar;
	DoubleArray              extremeSoFar;	    // [nbAssets]
    DoubleArray              refLevelSoFar;     // [nbAssets]
    DoubleArray              refLevelSoFarTemp; // [nbAssets]
    SimpleDoubleArray        couponsSoFar;      // [nbCoupons]
    
    // Operational aggregation and performance calcs
    IAggregateSP             timeBasket;
    IDoubleArrayModifierSP   timeBasketComponents;
    IAggregateSP             assetBasket;
    IDoubleArrayModifierSP   overallOption;
    
	TrivialDoubleArray       BaskLookBk;	// just a double ... may be more natural way to phrase this?

    IntArray                 couponMap;		// [nbmonitoringDates], which monitoring dates match coupon dates

public:
    
    // equivalent to InstIntoMCProduct
    BasketOfLookBacksMC(const BasketOfLookBacks* inst,
						const SimSeriesSP& simSeries):
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
        refLevel(nbAssets, 0.0),
        refLevelTemp(nbAssets, 0.0),
        coupons(inst->couponDates.size(), 0.0),
        assetComps(nbAssets, 0.0),
        iCouponSoFar(0),
        extremeSoFar(nbAssets, 1000.0 - 2000.0 * inst->isLookbackForMax),
        refLevelSoFar(nbAssets, 0.0),
        refLevelSoFarTemp(nbAssets, 100000.0 - 200000.0 * inst->isLookbackForMax),
        couponsSoFar(inst->couponDates.size(), 0.0),
        BaskLookBk(0.0) {
        
		// which monitoring dates match coupon dates
        bool isTrivial;
        couponMap = DateTime::createMapping(inst->AvgOutDates,
                                            inst->couponDates,
                                            isTrivial);

        timeBasketComponents = IDoubleArrayModifierSP(inst->timeBasketComponents->getModifier(&coupons));
        timeBasket = IAggregateSP(inst->timeBasket->getAggregate(&coupons));
        
        assetBasket = IAggregateSP(inst->assetBasket->getAggregate(&assetComps));
        
		overallOption = IDoubleArrayModifierSP(inst->overallOption->getModifier(&BaskLookBk));
    }

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen)const{
        // empty
    }

    
	// called twice, if there are monitoring dates in the past; called once, otherwise
	void payoff(const IPathGenerator* pathGen, IMCPrices& prices) {
		
		// same simulation dates for all assets
		// note: if payoff fct is called twice, then input concerning pathGen->begin and pathGen->end is different!
		int    beginIdx = pathGen->begin(0);	
        int    endIdx   = pathGen->end(0);
        
		int    iAsset;
		
		// load historical values (only non-zero if there is some past)
        iCoupon = iCouponSoFar;		
        extreme = extremeSoFar;
        refLevel = refLevelSoFar;
        refLevelTemp = refLevelSoFarTemp;
        coupons = couponsSoFar;

		// determination of refLevel, depending on length of past and whether or not isCliquetStyle
        for (iAsset = 0; iAsset < nbAssets; iAsset ++) {
			if (!inst->isCliquetStyle || iCoupon == 0) {
				// if [non-cliquet] or if [cliquet and first coupon], then refLevel = pathGen->refLevel 
                refLevel[iAsset] = pathGen->refLevel(iAsset, 0);
            } else {
				// if [cliquet and second+ coupon], then refLevel = refLevelSoFar
                refLevel[iAsset] = refLevelSoFar[iAsset];
                refLevelTemp[iAsset] = refLevelSoFarTemp[iAsset];
            }
		}

		// run through all simulation dates
        for(int iStep=beginIdx; iStep<endIdx; iStep++) {
			bool isCouponDate = (couponMap[iStep]==0);

			// alternative 1: aggregate first, then take MAX/MIN
			if (inst->isLookbackOnBasket) {
				for(iAsset=0; iAsset<nbAssets; iAsset++) {
					assetComps[iAsset] = pathGen->Path(iAsset,0)[iStep] / refLevel[iAsset];
				}

				// aggregation
				double SumAcrossAssets = assetBasket->aggregate();

				// update value of aggregated MIN/MAX at each monitoring date
				if(inst->isLookbackForMax) {
					extreme[0] = Maths::max(extreme[0], SumAcrossAssets);
					} else {
					extreme[0] = Maths::min(extreme[0], SumAcrossAssets);
				}
				// save aggregated MIN/MAX at each coupon date
				if (isCouponDate) {
					coupons[iCoupon] = extreme[0]; // save value for coupon date
					iCoupon++;
					// reset extreme, if isLookbackFromStart = FALSE
					if(!inst->isLookbackFromStart) {extreme[0] = 1000.0 - 2000.0 * inst->isLookbackForMax;} 
					// reset refLevel, if isCliquetStyle = TRUE
					if(inst->isCliquetStyle) {
						for(iAsset=0; iAsset<nbAssets; iAsset++) {
							refLevel[iAsset] = pathGen->Path(iAsset,0)[iStep]; 
						}
					}
				} 
			} // end of alternative 1
			
			// alternative 2: take MIN/MAX first, then aggregate
			if (!inst->isLookbackOnBasket) {
				for(iAsset=0; iAsset<nbAssets; iAsset++) {
                    assetComps[iAsset] = pathGen->Path(iAsset,0)[iStep] / refLevel[iAsset];
                    
                    // update value of MIN/MAX for each asset, at each monitoring date
					if(inst->isLookbackForMax) {
						extreme[iAsset] = Maths::max(extreme[iAsset], assetComps[iAsset]);
                        refLevelTemp[iAsset] = Maths::max(refLevelTemp[iAsset], pathGen->Path(iAsset,0)[iStep]);
						} else {
						extreme[iAsset] = Maths::min(extreme[iAsset], assetComps[iAsset]);
                        refLevelTemp[iAsset] = Maths::min(refLevelTemp[iAsset], pathGen->Path(iAsset,0)[iStep]);
					}
                    // save value of MIN/MAX for each asset, at each coupon date
					if (isCouponDate) {
						// rename current MIN/MAX in order to apply assetBasket->aggregate();
						assetComps[iAsset] = extreme[iAsset];
						// reset extreme if isLookbackFromStart = FALSE
						if(!inst->isLookbackFromStart) {extreme[iAsset] = 1000.0 - 2000.0 * inst->isLookbackForMax;} 
						// reset refLevel, if isCliquetStyle = TRUE
						if(inst->isCliquetStyle) {
							if(!inst->isRefPrevExtreme) {
								refLevel[iAsset] = pathGen->Path(iAsset,0)[iStep];
							} else {
                                refLevel[iAsset] = refLevelTemp[iAsset];
							}
                        refLevelTemp[iAsset] = 100000.0 - 200000.0 * inst->isLookbackForMax;
                        }
					}
				}
				// Aggregation
				if (isCouponDate) {
					coupons[iCoupon] = assetBasket->aggregate();
					iCoupon++;
				}
			} // end of alternative 2
		} // end of run through all simulation dates

        // preserve values for past - before we modify/sort the coupons.
        if (pathGen->doingPast()){ 
            extremeSoFar = extreme;
            refLevelSoFar = refLevel;
            couponsSoFar = coupons;
            iCouponSoFar = iCoupon;
        }
        if (!pathGen->doingPast() || !hasFuture()) {
            // Compute a payoff, but only when we have a "complete" situation : either 
            // doingPast() and all is past, or !doingPast().
            // Now perform the time basket - apply perf modifiers
            timeBasketComponents->apply();
            // then make a basket of them
            BaskLookBk() = timeBasket->aggregate();
            overallOption->apply();
            prices.add(inst->notional * BaskLookBk()); 
        }
    }

    // for the LogNormal path generator
	// this is run between past and simulation runs for the future
	// this is run for each asset separately
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen, int iAsset) const {
		static const string routine = "BasketOfLookBacksMC::getVolInterp";
        CVolRequestLNArray reqarr(1);
        const DateTime&    startDate = getRefLevel()->getAllDates().front();
        const DateTime&    today = getToday();
        const DateTime&    lastSimDate = getSimSeries()->getLastDate();
        bool               fwdStarting = startDate.isGreater(today);
        double             interpLevel = inst->timeBasketComponents->getInterpLevel(iAsset);

        if (inst->isCliquetStyle) {
			// get hold of the future strike dates
            int numLiveCliqs = inst->couponDates.size() - iCouponSoFar; // number of remaining coupons
            if (numLiveCliqs<=0) {throw ModelException(routine, "No future coupons!?");}

			// determine the dates of the remaining coupons
            DateTimeArray liveCliqStartDates(numLiveCliqs);
            for (int iCliquet = 0; iCliquet < numLiveCliqs; iCliquet++){
                int iCoupon = iCouponSoFar+iCliquet-1;
                liveCliqStartDates[iCliquet] = iCoupon<0?startDate:inst->couponDates[iCoupon];
            }
            // same strike levels (percentage) for future coupon dates 
			// however, strike level of current coupon period has to be adjusted (absolute, instead of percentage)
            DoubleArray  strikes(numLiveCliqs, interpLevel);
            
			if (!fwdStarting){
                double refLevel = iCouponSoFar == 0 ? pathGen->refLevel(iAsset,0) : refLevelSoFar[iAsset];
                strikes[0] = refLevel * interpLevel;
            }
            
			reqarr[0] =  CVolRequestLNSP(new CliquetVolRequest(fwdStarting, 
                                                               liveCliqStartDates, 
                                                               lastSimDate,
                                                               strikes));
		} else { 
            if (!fwdStarting){
                interpLevel = interpLevel * pathGen->refLevel(iAsset, 0);
			}
			reqarr[0] = CVolRequestLNSP(new LinearStrikeTSVolRequest(interpLevel,
                                                                     startDate,
                                                                     lastSimDate,
                                                                     fwdStarting));
	}
    return reqarr;
    }
};

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* BasketOfLookBacks::createProduct(const MonteCarlo* model) const {

    // XXX Cliquet style resetting is not supported with implied - how enforce that?

    // we create a SimSeries object which says which assets need
    // which dates to be simulated
    SimSeriesSP simSeries(new SimSeries(assets->NbAssets())); /* create empty
                                                                 one */
    simSeries->addDates(AvgOutDates);
    return new BasketOfLookBacksMC(this, simSeries);
}

CClassConstSP const BasketOfLookBacks::TYPE = CClass::registerClassLoadMethod(
    "BasketOfLookBacks", typeid(BasketOfLookBacks), BasketOfLookBacks::load);

// * for class loading (avoid having header file) */
bool BasketOfLookBacksLoad() {
    return (BasketOfLookBacks::TYPE != 0);
}

DRLIB_END_NAMESPACE





