//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : AverageLocker.cpp
//
//   Description : A running lock-in average of the basket performance. At maturity, it is a call option
//                  with participation rate.
//
//   Date        : October 2005
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SimSeries.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/MCPathConfig.hpp"
#include "edginc/LegalTerms.hpp"
#include "edginc/IDoubleArray.hpp"
#include "edginc/IAggregate.hpp" // for basket performance

DRLIB_BEGIN_NAMESPACE

class AverageLocker: public GenericNFBase, 
                     virtual public IMCIntoProduct {
protected:

    DateTimeArray				monitoringDates;    // all dates after ref date to the end
    double						strike;        // the final strike
    double						participation;      // the participation in the call option
	double						lockThreshold;				// level at which basket goes above to get locked
	string						lockBasedMech;	// lock-in based on Average or Max So Far
    IAggregateMakerSP           assetBasket;   // pct/rainbow/product basket etc of ...

public:
    static CClassConstSP const TYPE;
    friend class AverageLockerMC;
    friend class AverageLockerSVMC;


	void Validate() {

		static const string method = "AverageLocker::Validate";

        // Call super class to get observation map populated
        GenericNFBase::Validate();

		try {
			//Have to do any refLevel activity after GetMarket
			// need to make sure the ref level averaging in has finished by the time we hit monitoring dates
			DateTimeArray overlap = refLevel->getFutureDates(monitoringDates[0]);
			if (!overlap.empty()) {
				throw ModelException(method, "Monitoring dates cannot start ("+ monitoringDates[0].toString()
					+") until the reference level has been set");
			}
		} catch (exception& e) {
			throw ModelException(e,method);
		}
	}

    // validation
    void validatePop2Object(){
        static const string method = "AverageLocker::validatePop2Object";
        GenericNFBase::validatePop2Object();
    
        // validate dates are not empty - order is handled by SimSeries
        if (monitoringDates.empty()) {
            throw ModelException(method, "No monitoring dates supplied!");
        }

		if (!Maths::isPositive(strike)){
            throw ModelException(method, "strike ("+Format::toString(strike)+") should be positive");
        }
		if (!Maths::isPositive(lockThreshold)){
            throw ModelException(method, "lockThreshold ("+Format::toString(lockThreshold)+") should be positive");
        }
		if (!(lockBasedMech=="Average"||lockBasedMech=="MaxSoFar")){
            throw ModelException(method, "lockBasedMech should be Average or MaxSoFar");
        }
	}

    // returns the array of all sampling dates the instrument will ever need
    // excluding (possibly) the ref level dates (handled in GenericNFBase)
    const DateTimeArray samplingDates() const {
        return monitoringDates;
    }

    /** Implementation of MonteCarlo::IntoProduct interface - the
        implementation of this is below */
    virtual IMCProduct* createProduct( const MonteCarlo* model) const; // see below

private:
    AverageLocker(): GenericNFBase(TYPE) {} // for reflection
    AverageLocker(const AverageLocker& rhs);     // not implemented
    AverageLocker& operator=(const AverageLocker& rhs); // not implemented

    static IObject* defaultAverageLocker(){
        return new AverageLocker();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(AverageLocker, clazz);
        SUPERCLASS(GenericNFBase);
        IMPLEMENTS(MonteCarlo::IIntoProduct);
        EMPTY_SHELL_METHOD(defaultAverageLocker);
        FIELD(monitoringDates, "all dates on which levels are monitored");
		FIELD(strike, "the final strike");
        FIELD(participation, "the participation in the call option");
		FIELD(lockThreshold, "level at which basket goes above to get locked");
		FIELD(lockBasedMech, "lock-in based on Average or MaxSoFar");
        FIELD(assetBasket, "How to aggregate asset perfs");
    }
};

// string constants for lock mechanism
#define AVERAGE_MECHANISM       "Average"
#define MAX_SO_FAR_MECHANISM    "MaxSoFar"

// interface for lock mechanism
class ILockMechanism : public virtual IObject {
public:
    static CClassConstSP const TYPE;

    // keep updating n, lockBsk in the loop
   	virtual	void updateLockLevel(double basket) = 0;

    //only used once for final payoff
	virtual	double getLockLevel() = 0;
};

CClassConstSP const ILockMechanism::TYPE = CClass::registerInterfaceLoadMethod(
    "ILockMechanism", typeid(ILockMechanism), 0);

typedef smartPtr<ILockMechanism> LockMechanismSP;

// Lock-in mechanism based on running average
class AverageLockMechanism: public CObject,
                            virtual public ILockMechanism {
public:
    static CClassConstSP const TYPE;

    AverageLockMechanism(double lockThreshold) : CObject(TYPE),
            n(0), lockBsk(0.0), threshold(lockThreshold) {}	
	
	void updateLockLevel(double basket){
		if(n == 0) {
			if(basket > threshold) {
				lockBsk = basket;
				n++;
			}
		} else {
			if(basket > lockBsk) {
				n++;
				lockBsk = Maths::max(lockBsk, ((n-1)*lockBsk + basket)/n);
			}
		}
	}

	double getLockLevel() {
        return lockBsk;
    }
		
private:
    AverageLockMechanism(const AverageLockMechanism &rhs);
    AverageLockMechanism& operator=(const AverageLockMechanism& rhs);

    static IObject* defaultAverageLockMechanism() {
        return new AverageLockMechanism(0.0);
    }

    static void load(CClassSP& clazz) {
        REGISTER(AverageLockMechanism, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ILockMechanism);
        EMPTY_SHELL_METHOD(defaultAverageLockMechanism);
        FIELD(threshold,   "threshold");
        FIELD(lockBsk,           "lockBsk");
        FIELD(n,      "n");
    }

	int n;
	double lockBsk; 
    double threshold;
};

CClassConstSP const AverageLockMechanism::TYPE = CClass::registerClassLoadMethod(
    "AverageLockMechanism", typeid(AverageLockMechanism), load);

//Lock-in mechanism based on Max so far
class MaxLockMechanism: public CObject,
                        virtual public ILockMechanism {
public:
    static CClassConstSP const TYPE;

	MaxLockMechanism(double lockThreshold) : CObject(TYPE),
            n(0), lockBsk(0.0), threshold(lockThreshold), sum(0.0) {}	
	
	void updateLockLevel(double basket){
		if(n == 0) {
			if(basket > threshold) {
				lockBsk = basket;
				n++;
				sum = lockBsk;
			}
		} else {
			if(basket > lockBsk) {
				n++;
				lockBsk = basket;
				sum = (sum*(n-1) + lockBsk) / n;
			}
		}			
	}

	double getLockLevel() {
        return sum;
    }
		
private:
    MaxLockMechanism(const MaxLockMechanism &rhs);
    MaxLockMechanism& operator=(const MaxLockMechanism& rhs);

    static IObject* defaultMaxLockMechanism() {
        return new MaxLockMechanism(0.0);
    }

    static void load(CClassSP& clazz) {
        REGISTER(MaxLockMechanism, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ILockMechanism);
        EMPTY_SHELL_METHOD(defaultMaxLockMechanism);
        FIELD(threshold,   "threshold");
        FIELD(lockBsk,           "lockBsk");
        FIELD(n,      "n");
        FIELD(sum, "sum");
    }

	int n;
	double lockBsk; 
	double threshold;
	double sum;
};

CClassConstSP const MaxLockMechanism::TYPE = CClass::registerClassLoadMethod(
    "MaxLockMechanism", typeid(MaxLockMechanism), load);

/* MC product class for AverageLocker */
class AverageLockerMC: public IMCProduct,
                    virtual public IMCProductLN,
                    virtual public IMCProductImplied {

public:
	    
	AverageLockerMC(const AverageLocker*      inst,
                 const DateTimeArray& monitorDates,
                 const SimSeriesSP&   simSeries):
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
        perf(nbAssets, 0.0)
		{
			assetBasket = IAggregateSP(inst->assetBasket->getAggregate(&perf));
			
			if (inst->lockBasedMech == AVERAGE_MECHANISM) {
				lockerSoFar = LockMechanismSP(new AverageLockMechanism(inst->lockThreshold));
			} else {
				lockerSoFar = LockMechanismSP(new MaxLockMechanism(inst->lockThreshold));
			}
		}

    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&                prices) {
        int beginIdx = pathGen->begin(0); // same for all assets
        int endIdx   = pathGen->end(0);
		
		LockMechanismSP locker = LockMechanismSP(copy(lockerSoFar.get()));
		
		for (int iStep = beginIdx; iStep < endIdx; iStep++) {    
            for(int iAsset=0; iAsset<nbAssets; iAsset++) {
                perf[iAsset] = pathGen->Path(iAsset, 0/*iPath*/)[iStep] 
                                / pathGen->refLevel(iAsset, 0/*iPath*/);
			}
        
			double basket = assetBasket->aggregate();

       		
			locker->updateLockLevel(basket);
		
		}

		if (doingPast()) {
			lockerSoFar = LockMechanismSP(copy(locker.get()));
		}

        if (!doingPast() || !hasFuture()) {
			double lockBsk = locker->getLockLevel();
			double value = inst->participation * Maths::max(0.0, lockBsk - inst->strike);
	        prices.add(inst->notional * value); 
		}
	}

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseImplied(const  IMCPathGenerator*  pathGen)const{}

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen)const{
        static const string routine = "AverageLockerMC::initialiseLN";
        throw ModelException(routine, "Methodology not supported");
    }

    // any old level so that MC implied works
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
                                    int                     iAsset) const {
        static const string method = "AverageLockerMC::getVolInterp";

        try {
            // one interp level per asset
            CVolRequestLNArray reqarr(1);
            const DateTime& today = getToday();
            const DateTime& startDate = refLevel->getAllDates().front();
            const DateTime& lastSimDate = getSimSeries()->getLastDate();
            bool  fwdStarting = startDate.isGreater(today);

            double interpLevel  = 1.0;

            reqarr[0] = CVolRequestLNSP(new LinearStrikeTSVolRequest(interpLevel,
                                                                     startDate,
                                                                     lastSimDate,
                                                                     fwdStarting));

            return reqarr;
        } catch (exception& e) {
            throw ModelException(e, method);
        }
    }

private:
    const AverageLocker*		inst;         // Instrument
    int                         nbAssets;     // convenient
    SimpleDoubleArray           perf;         // current perf per asset
    IAggregateSP				assetBasket;
	LockMechanismSP				lockerSoFar;
};

////////////////////////////////////////////////////////////////////////////////////////////////////// 

/* MC product class for AverageLockerSV */
class AverageLockerSVMC: public MCProductClient,
                    virtual public IMCProductLN,
                    virtual public IMCProductImplied {

public:
    /** Appends 'true' (ie non derived) state variable generators
        required to the supplied collector.*/
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const{
        static const string routine = "AverageLockerSVMC::collectStateVars";
        try{
            svCollector->append(spotGen.get());             // spot level
            svCollector->append(refLevelGen.get());         // reference level
            svCollector->append(dfGen.get());               // and a DiscFactor one
        } catch (exception& e) {
            throw ModelException(e, routine);
        }
    }

    /** Override default method on IMCProduct. This method is called every time
        the path generator is changed (which is, at the moment, when the
        past path generator is created, and then when the future path
        generator is created  */
    virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen) {
        static const string routine = "AverageLockerSVMC::pathGenUpdated";
        try{
            spotSV = spotGen->getSpotSV(newPathGen);
            refLevelSV = refLevelGen->getRefLevelSV(refLevelSV, newPathGen);
            dfSV = dfGen->getSVDiscFactor(newPathGen);
        } catch (exception& e) {
            throw ModelException(e, routine);
        }
    }

    AverageLockerSVMC(const AverageLocker*      inst,
                 const DateTimeArray& monitorDates,
                 const SimSeriesSP&   simSeries):
        MCProductClient(inst->assets.get(),
                        inst->valueDate,
                        inst->discount.get(),
                        inst->refLevel,
                        simSeries,
                        inst->pastValues,
                        inst->instSettle.get(),
                        simSeries->getLastDate()),
        inst(inst),
        nbAssets(getNumAssets()),
        perf(nbAssets, 0.0),
        assetBasket(inst->assetBasket->getAggregate(&perf)),
		spotGen(new SVGenSpot(simSeries)),
        refLevelGen(inst->refLevel->createStateVarGen(getMultiFactors(), getToday())),
        dfGen(new SVGenDiscFactor(inst->valueDate,inst->discount.getSP(),
                               inst->instSettle, simSeries->getLastDate())) {
			if (inst->lockBasedMech == AVERAGE_MECHANISM) {
				lockerSoFar = LockMechanismSP(new AverageLockMechanism(inst->lockThreshold));
			} else {
				lockerSoFar = LockMechanismSP(new MaxLockMechanism(inst->lockThreshold));
			}
        }

    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&                prices) {
        

        // Same begin & end for all assets, so read from the first
        const SVPath& path = spotSV->path(0);
        int    beginIdx = path.begin();
        int    endIdx   = path.end();

		LockMechanismSP locker = LockMechanismSP(copy(lockerSoFar.get()));
		
		for (int iStep = beginIdx; iStep < endIdx; iStep++) {    
            for(int iAsset=0; iAsset<nbAssets; iAsset++) {
                const SVPath& path = spotSV->path(iAsset);
                perf[iAsset] = path[iStep] / refLevelSV->refLevel(iAsset);
			}
        
			double basket = assetBasket->aggregate();	

			locker->updateLockLevel(basket);
		}


		if (doingPast()) {
			lockerSoFar = LockMechanismSP(copy(locker.get()));
		}

        // needs the discount factor since this is the sv case
        if (!doingPast() || !hasFuture()) {
			double lockBsk = locker->getLockLevel();
			double value = inst->participation * Maths::max(0.0, lockBsk - inst->strike);
	        prices.add(inst->notional * value * dfSV->firstDF()); 
		}
	}
	

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseImplied(const  IMCPathGenerator*  pathGen)const{}

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen)const{
        static const string routine = "AverageLockerSVMC::initialiseLN";
        throw ModelException(routine, "Methodology not supported");
    }


    // any old level so that MC implied works
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
                                    int                     iAsset) const {
        static const string method = "AverageLockerSVMC::getVolInterp";

        try {
            // one interp level per asset
            CVolRequestLNArray reqarr(1);
            const DateTime& today = getToday();
            const DateTime& startDate = refLevelGen->getAllDates().front();
            const DateTime& lastSimDate = getSimSeries()->getLastDate();
            bool  fwdStarting = startDate.isGreater(today);

            double interpLevel  = 1.0;

            reqarr[0] = CVolRequestLNSP(new LinearStrikeTSVolRequest(interpLevel,
                                                                     startDate,
                                                                     lastSimDate,
                                                                     fwdStarting));

            return reqarr;
        } catch (exception& e) {
            throw ModelException(e, method);
        }
    }

private:
    const AverageLocker*		inst;         // Instrument
    int                         nbAssets;     // convenient
    SimpleDoubleArray                 perf;         // current perf per asset
    IAggregateSP              assetBasket;
	LockMechanismSP				lockerSoFar;

   
    // State variables and generators
    SVGenSpotSP                  spotGen;      // Generator for spot
    IRefLevel::IStateVarGenSP refLevelGen;  // Generator for ref level
    SVGenDiscFactorSP            dfGen;        // Generator for discount factors
    SVGenSpot::IStateVarSP       spotSV;       // Spot state variable
    IRefLevel::IStateVarSP    refLevelSV;   // Ref level state variable
    SVDiscFactorSP dfSV;         // Df state variable
};

//////////////////////////////////////////////////////////////////////////


/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* AverageLocker::createProduct(const MonteCarlo* model) const {
    static const string method = "AverageLocker::createProduct";

    try {
        int nbAssets = assets->NbAssets();
        
        // Create a SimSeries object which says which assets need
        // which dates to be simulated
        SimSeriesSP simSeries(new SimSeries(nbAssets));

        // do we crop these dates to only be future ones??
        simSeries->addDates(monitoringDates);
    
        if(model->stateVarUsed()) {
            return new AverageLockerSVMC(this, monitoringDates, simSeries);
        } else {
            // Otherwise, use old methodology
            return new AverageLockerMC(this, monitoringDates, simSeries);
        }
    } catch(exception& e) {
        throw ModelException(e, method);
    }
}

CClassConstSP const AverageLocker::TYPE = CClass::registerClassLoadMethod(
    "AverageLocker", typeid(AverageLocker), AverageLocker::load);

// * for class loading (avoid having header file) */
bool AverageLockerLoad() {
    return (AverageLocker::TYPE != 0);
}

DRLIB_END_NAMESPACE
