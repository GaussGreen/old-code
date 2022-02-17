//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FinalPerf.cpp
//
//   Description : Class to express the final performance:
//
//   Date        : January 2006
//
//
//----------------------------------------------------------------------------

#include <math.h>
#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Class.hpp"
#include "edginc/Format.hpp"

#include "edginc/IFinalPerf.hpp"
#include "edginc/Model1F.hpp"
#include "edginc/Barrier.hpp"

//to remove, need it here only for have the FP_MIN
#include "edginc/Tree1f.hpp"

DRLIB_BEGIN_NAMESPACE

void IFinalPerfMaker::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER_INTERFACE(IFinalPerfMaker, clazz);
    EXTENDS(IObject);
}    

CClassConstSP const IFinalPerfMaker::TYPE = CClass::registerInterfaceLoadMethod(
    "IFinalPerfMaker", typeid(IFinalPerfMaker), load);

//---------------------------------------------------------------------//
// Standard Vanilla Performance.  
class VanillaPerf : virtual public IFinalPerf {
public:
    VanillaPerf(VanillaPerfMaker*  maker,
                DateTimeArray      timeLine,
                DateTime           valueDate,
                double             refLevel,
                const CAsset*        asset,
                const YieldCurveConstSP   discount,
                const InstrumentSettlementSP instSettle):
                 maker(maker), timeLine(timeLine), refLevel(refLevel), 
                 perfs(1, 0.0){// only need a double, not necessary to be array.
        
        int numSteps = timeLine.size();

        // set up generic performance at maturity.
        performance = IDoubleArrayModifierSP(maker->matPerformance->getModifier(&perfs));        

        // set up settle pv
        matDateTime = maker->matDateTime;
        settlePV = discount->pv(matDateTime, instSettle->settles(matDateTime, asset));
    };

    virtual double getIntrinsic(double s, double k){
        perfs[0] = s/refLevel;
        performance->apply(0);
        double value = perfs[0];
        return value;
    };

    // avgValue return already pv adjusted.  
    virtual double getValue(int iPrice, int simDateIdx, double s, double k){
        perfs[0] = s/refLevel;
        performance->apply(0);
        double value = perfs[0];
        return value * settlePV;
    };

private:
    VanillaPerfMakerSP maker;
    DateTime           matDateTime;
    DateTimeArray      timeLine;
    double             refLevel;

    IDoubleArrayModifierSP performance;
    SimpleDoubleArray      perfs;    

    double              settlePV;
};

typedef refCountPtr<VanillaPerf> VanillaPerfSP;


VanillaPerfMaker::VanillaPerfMaker(IDoubleArrayModifierMakerSP  matPerformance,
                                    DateTime               matDateTime) : CObject(TYPE),
                matPerformance(matPerformance), matDateTime(matDateTime){};

// 
IFinalPerfSP VanillaPerfMaker::getFinalPerf(const CAsset*        asset,
                                            const DateTime       valueDate,
                                            const double         refLevel,
                                            const YieldCurveConstSP   discount,
                                            const InstrumentSettlementSP instSettle,
                                            const InstrumentSettlementSP premiumSettle,
                                            string                       ccyTreatment,
                                            const DateTimeArray& simDates)
{
    return VanillaPerfSP(new VanillaPerf(this, simDates, valueDate, refLevel, asset, discount, instSettle));	
}

// validation & and make aop for the calling from the public interface.    
void VanillaPerfMaker::validatePop2Object(){
    static const string routine = "VanillaPerfMaker::validatePop2Object";    
}

DateTimeArray VanillaPerfMaker::getMatDates(){
    DateTimeArray mats(1,matDateTime);
    return mats;
}

DateTime VanillaPerfMaker::getMatStartDate(){
    return matDateTime;
};

bool VanillaPerfMaker::sensShift(Theta* shift, CAsset* asset, DateTime valueDate){
    return true;
}

class VanillaPerfMakerHelper{
public:
    static IObject* defaultVanillaPerfMaker(){
        return new VanillaPerfMaker();
    }
    // necessary to be all optional so as to allow empty in IMS.
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(VanillaPerfMaker, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IFinalPerfMaker);
        EMPTY_SHELL_METHOD(defaultVanillaPerfMaker);
        FIELD(matDateTime,   "maturity Date");
        FIELD_MAKE_OPTIONAL(matDateTime);
        FIELD(matPerformance, "Overall Perf Modifier");
        FIELD_MAKE_OPTIONAL(matPerformance);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

CClassConstSP const VanillaPerfMaker::TYPE =
CClass::registerClassLoadMethod("VanillaPerfMaker", 
                                typeid(VanillaPerfMaker), VanillaPerfMakerHelper::load);

//---------------------------------------------------------------------//
// Standard Average Out Performance.  AvgOutPerf itself is already exposed
// by EGBond or other instruemnts.  It can be one of IFinalPerf, here.
class AveragePerf : virtual public IFinalPerf {
public:
    AveragePerf(AveragePerfMaker*  maker,
                DateTimeArray      timeLine,
                const DateTime      valueDate,
                const CAsset*        asset,
                const InstrumentSettlementSP instSettle,
                const InstrumentSettlementSP premiumSettle,
                string                       ccyTreatment,
                const YieldCurveConstSP   discount):
         maker(maker), timeLine(timeLine),isExpired(false),valueDate(valueDate),
         discount(discount), instSettle(instSettle),asset(asset) {
             // make own AvgOutPerfSP.
             myAOP = maker->getAvgOutPerf();
             myAOP->setInstInfo(instSettle,
                                premiumSettle,
                                asset,
                                ccyTreatment,
                                discount,
                                valueDate);

             // keep the payoff timing index.
             matDateTime = maker->getMatStartDate();             
             if (valueDate>matDateTime)
                isExpired = true;
             else
                payTiming = matDateTime.find(timeLine);
    };

    virtual double getIntrinsic(double s, double k){
        double value = avgValue(k);
        value /= discount->pv(instSettle->settles(myAOP->avgOut->getLastDate(), asset));
        return value;
    };

    // avgValue return already pv adjusted.  
    virtual double getValue(int iPrice, int simDateIdx, double s, double k){
        // ignore iPrice.
        if(isExpired)
            return avgValue(s);
        else if (simDateIdx == payTiming)            
            return avgValue(s, k, matDateTime);
        else
            return 0.0;
    };

    // return average out closed form price.
    double avgValue(double const spot, double SpotRef, const DateTime& when) const {
        double s = Maths::max  (spot,0.00001);
        double strikeScale = SpotRef/s;
        double notl = s/SpotRef;
        double refLevel = s;
        double value;
        value = myAOP->AvgValue(strikeScale, refLevel, notl, when);
        return value;
    }

    // return average out price, only after all simulation date.
    double avgValue(double SpotRef) const {
        double value;
        value = myAOP->AvgValue(SpotRef, SpotRef, 1.0, valueDate);   
        return value;
    }

private:
    AveragePerfMakerSP maker;
    AvgOutPerfSP        myAOP;
    DateTime           matDateTime;
    DateTimeArray      timeLine;
    DateTime           valueDate;
    int payTiming;
    bool isExpired;
    YieldCurveConstSP   discount;
    InstrumentSettlementSP instSettle;
    const CAsset*        asset;
};

typedef refCountPtr<AveragePerf> AveragePerfSP;


AveragePerfMaker::AveragePerfMaker(AvgOutPerfMakerSP aopMaker)
                                    :CObject(TYPE),aopMaker(aopMaker)
{
	// Construct
    aop = AvgOutPerfSP(new AvgOutPerf(aopMaker));
}

// 
IFinalPerfSP AveragePerfMaker::getFinalPerf(const CAsset*        asset,
                                            const DateTime       valueDate,
                                            const double         refLevel,
                                            const YieldCurveConstSP   discount,
                                            const InstrumentSettlementSP instSettle,
                                            const InstrumentSettlementSP premiumSettle,
                                            string                       ccyTreatment,
                                            const DateTimeArray& simDates)
{
    return AveragePerfSP(new AveragePerf(this, simDates, valueDate, asset, 
                                         instSettle, premiumSettle, ccyTreatment, 
                                         discount));	
}

// validation & and make aop for the calling from the public interface.    
void AveragePerfMaker::validatePop2Object(){
    static const string routine = "AveragePerfMaker::validatePop2Object";
    aop = AvgOutPerfSP(new AvgOutPerf(aopMaker));
}


AvgOutPerfSP AveragePerfMaker::getAvgOutPerf()
{
    if (!aop.get())
        aop = AvgOutPerfSP(new AvgOutPerf(aopMaker));   // need to remake for tweeks.
	return aop;
}

DateTimeArray AveragePerfMaker::getMatDates(){
    if (!aop.get())
        aop = AvgOutPerfSP(new AvgOutPerf(aopMaker));   // need to remake for tweeks.
    return aop->avgOut->getDates();
}

// return average out start date for tree.
// Finding out the end of tree.  
// if there's just 1 avg out date (which must be maturity).
// Otherwise, set the end of the day before avg out starts
DateTime AveragePerfMaker::getMatStartDate()
{
    DateTimeArray avgDates = getMatDates();
    if (avgDates.size() == 1)
        return avgDates[0];
    else                 
        return DateTime(avgDates[0].getDate()-1, DateTime::END_OF_DAY_TIME);	    
};


// 
bool AveragePerfMaker::sensShift(Theta* shift, CAsset* asset, DateTime valueDate)
{
    if (!aop.get())
        aop = AvgOutPerfSP(new AvgOutPerf(aopMaker));   // need to remake for tweeks.
	aop->avgOut->roll(shift->getUtil(valueDate), 0, asset);
    return true;
}

class AveragePerfMakerHelper{
public:
    static IObject* defaultAveragePerfMaker(){
        return new AveragePerfMaker();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(AveragePerfMaker, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IFinalPerfMaker);
        EMPTY_SHELL_METHOD(defaultAveragePerfMaker);
        FIELD(aopMaker,          "AvgOutPerfMaker");
        FIELD_MAKE_OPTIONAL(aopMaker);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

CClassConstSP const AveragePerfMaker::TYPE =
CClass::registerClassLoadMethod("AveragePerfMaker", 
                                typeid(AveragePerfMaker), AveragePerfMakerHelper::load);


//---------------------------------------------------------------------//
// Knock-In performance.
class KnockInPerf : virtual public IFinalPerf {
public:
    KnockInPerf(KnockInPerfMaker*  maker,
                DateTimeArray      timeLine,
                DateTime           valueDate,
                double             refLevel,
                const CAsset*        asset,
                const YieldCurveConstSP   discount,
                const InstrumentSettlementSP instSettle):
                 maker(maker), timeLine(timeLine), refLevel(refLevel), 
                 perfs(1, 0.0){// only need a double, not necessary to be array.
        
        int numSteps = timeLine.size();
        // initialization of array per steps.
        isBarStep.resize(numSteps,false);

        // setting up barrier information
        DateTimeArray barDates = maker->barrier->getDates();
        barLevels = DoubleArraySP(new DoubleArray(numSteps+1,-1.0));
        barMap = Barrier::setStepLevel(maker->barrier, refLevel, timeLine, *barLevels);
        
        int i=0;
        for (i = barMap[i]; i < numSteps; i++, i += barMap[i])
            isBarStep[i] = true;

        // copy the original to adjust barrier
        adjBarLevels = barLevels;

        // set up generic performance at maturity.
        performance = IDoubleArrayModifierSP(maker->matPerformance->getModifier(&perfs));        

        // set up settle pv
        matDateTime = maker->matDateTime;
        settlePV = discount->pv(matDateTime, instSettle->settles(matDateTime, asset));
    };

    virtual double getIntrinsic(double s, double k){
        double value = 0.0;
        if (maker->isBreached){
            perfs[0] = s/refLevel;
            performance->apply(0);
            value = perfs[0];
        }
        return value;
    };

    virtual double getValue(int iPrice, int simDateIdx, double s, double k){
        double value;
        if (iPrice == 2){
            perfs[0] = s/refLevel;
            performance->apply(0);
            value = perfs[0];
        }
        else 
            value = 0.0;     // not knocked-in
        return value * settlePV;
    };

    virtual void upDatePerf(CTree1f* tree1f, int step){
        if (isBarStep[step] && !maker->intraDayMonitor && !maker->isBreached){
            vector<double> vol;            
            tree1f->GetStepVol(step, vol, &(*barLevels)[step], 0, 0);                                
            Barrier::BarrierAdjustment(vol[0], false, (*adjBarLevels)[step]);          // always down in!!
        }
    }

    // return true if it should have insert node level.
    virtual bool insNodeLevel(int simDateIdx, double& insLvl){
        insLvl = (*adjBarLevels)[simDateIdx];
        return true;
    }

    virtual double getValue0(int simDateIdx, int jNode, int iStart, int iEnd, const vector< double * > & p){
        if (simDateIdx == 0)
            return  p[iStart][jNode];
        else
            return 0.0;
    }

    // update the tree prices by judging KI monitor.
    virtual void upDateValues(int simDateIdx, const double s, double & p1, const double p2){
        if (   s < (*adjBarLevels)[simDateIdx] * (1.0+FP_MIN)
            || maker->isBreached  ){
            p1 = p2;
        }
    }

private:
    KnockInPerfMakerSP maker;
    AvgOutPerfSP        myKIP;
    DateTime           matDateTime;
    DateTimeArray      timeLine;
    double             refLevel;

    DoubleArraySP      barLevels;
    DoubleArraySP      adjBarLevels;
    
    IntArray   barMap;
    BoolArray  isBarStep;

    IDoubleArrayModifierSP performance;
    SimpleDoubleArray      perfs;    

    double              settlePV;
};

typedef refCountPtr<KnockInPerf> KnockInPerfSP;


IFinalPerfSP KnockInPerfMaker::getFinalPerf(const CAsset*        asset,
                                            const DateTime       valueDate,
                                            const double         refLevel,
                                            const YieldCurveConstSP   discount,
                                            const InstrumentSettlementSP instSettle,
                                            const InstrumentSettlementSP premiumSettle,
                                            string                       ccyTreatment,
                                            const DateTimeArray& simDates)
{
    return KnockInPerfSP(new KnockInPerf(this, simDates, valueDate, refLevel, asset, discount, instSettle));	
}


DateTimeArray KnockInPerfMaker::getMatDates(){
    DateTimeArray mats(1,matDateTime);
    return mats;
}

DateTime KnockInPerfMaker::getMatStartDate(){
    return matDateTime;
};

bool KnockInPerfMaker::sensShift(Theta* shift, CAsset* asset, DateTime valueDate){
    return true;
}

int KnockInPerfMaker::addInsNode(){
	return 1;       // need one insert node. 
}

// need one more price layer.
int KnockInPerfMaker::addNumPrices(){
	return 2;       // One is plain(or KO), the other is Knocked-In price.
}

bool KnockInPerfMaker::hasCritDates(DateTimeArray& critDates){
    if (barrier->getInterp() == Schedule::INTERP_NONE){
        const DateTimeArray& barDates = barrier->getDates();
        for (int i=0; i<barrier->length(); i++){
            // 2 days gap seem to be best for convergence
            critDates.push_back(barDates[i].rollDate(-2)); 
            critDates.push_back(barDates[i]);
        }
    } else {
        critDates.push_back(barrier->firstDate());
        critDates.push_back(barrier->lastDate());
    }

    return true;
}

bool KnockInPerfMaker::useEcoBar(){
    if (ecoBarrier.get()){
        barrier = ecoBarrier;
        return true;
    }else 
        return false;
}

ScheduleSP KnockInPerfMaker::getBarrier(bool isEconomic, bool &isContinuous) const{
    ScheduleSP bar = (isEconomic && ecoBarrier.get()) ? ecoBarrier: barrier;
    isContinuous = intraDayMonitor;
    return bar;
}

bool KnockInPerfMaker::isEnd(DateTime valDate,double spot,DateTime* endDate, double& barLvl){
    bool isHit = false; 
    if (isBreached){
        *endDate = breachDate;
        barLvl = barrier->interpolate(breachDate);        
        isHit = true;
    }else if (intraDayMonitor || valDate.getTime() == barrier->firstDate().getTime()){
        barLvl = barrier->interpolate(valDate);
        if (spot < barLvl){
            *endDate = valDate;
            isHit = true;
        }
    }
    return isHit;
}

class KnockInPerfMakerHelper{
public:
    static IObject* defaultKnockInPerfMaker(){
        return new KnockInPerfMaker();
    }

    // necessary to be all optional so as to allow empty in IMS.
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(KnockInPerfMaker, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IFinalPerfMaker);
        EMPTY_SHELL_METHOD(defaultKnockInPerfMaker);
        FIELD(barrier,          "Knock-In barrier");
        FIELD_MAKE_OPTIONAL(barrier);
        FIELD(ecoBarrier,          "Knock-In barrier");
        FIELD_MAKE_OPTIONAL(ecoBarrier);
        FIELD(intraDayMonitor, "intra-day barrier. false=>look at only specified timing of the day");
        FIELD_MAKE_OPTIONAL(intraDayMonitor);
        FIELD(matDateTime,   "maturity Date");
        FIELD_MAKE_OPTIONAL(matDateTime);
        FIELD(matPerformance, "Overall Perf Modifier");
        FIELD_MAKE_OPTIONAL(matPerformance);
        FIELD(isBreached,  "is barrier already breached?");
        FIELD_MAKE_OPTIONAL(isBreached);
        FIELD(breachDate,  "the breach date");
        FIELD_MAKE_OPTIONAL(breachDate);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

CClassConstSP const KnockInPerfMaker::TYPE =
CClass::registerClassLoadMethod("KnockInPerfMaker", 
                                typeid(KnockInPerfMaker), KnockInPerfMakerHelper::load);

/*********************************************************************/

/* The ultimate wrapping of IRebateMaker's mainly for use in Pyramid 
 */
#define FNL_PERF_TYPE_AVG     "Average"
#define FNL_PERF_TYPE_KNOCKIN "KnockIn"
#define FNL_PERF_TYPE_VAN     "Vanilla"

class FinalPerfMakerWrapper : public CObject,
                           virtual public IFinalPerfMaker{
public: // how can I have this protected or private?
    string                   finalPerfMakerType;
    AveragePerfMakerSP       avgMaker;
    KnockInPerfMakerSP       knockInMaker;
    
private:
    IFinalPerfMakerSP  realMaker; // $unregistered

public:
    static CClassConstSP const TYPE;

    virtual IFinalPerfSP getFinalPerf(const CAsset*        asset,
                                      const DateTime       valueDate,
                                      const double         refLevel,
                                      const YieldCurveConstSP   discount,
                                      const InstrumentSettlementSP instSettle,
                                      const InstrumentSettlementSP premiumSettle,
                                      string                       ccyTreatment,
                                      const DateTimeArray& simDates){
        return realMaker->getFinalPerf(asset, valueDate, refLevel, discount, 
                                       instSettle, premiumSettle, ccyTreatment, 
                                       simDates);
    }

    // Implementation of the Theta shift interface
    virtual bool sensShift(Theta* shift, CAsset* asset, DateTime valueDate) {
        return realMaker->sensShift(shift, asset, valueDate);
    }

    virtual DateTimeArray getMatDates(){
        return realMaker->getMatDates();        
    }

    virtual DateTime getMatStartDate(){
        return realMaker->getMatStartDate();
    }

    virtual bool useEcoBar(){
        return realMaker->useEcoBar();
    };

    virtual bool isEnd(DateTime valDate, double spot, DateTime* endDate, double &barLvl) {
        return realMaker->isEnd(valDate, spot, endDate, barLvl);        
    };


    
    // validation
    void validatePop2Object(){
        static const string routine = "FinalPerfMakerWrapper::validatePop2Object";

        if (finalPerfMakerType.empty()){
            throw ModelException(routine, "Blank Rebate Maker specified!");
        }
        if (finalPerfMakerType==FNL_PERF_TYPE_AVG) {
            if (avgMaker.get()) {
                realMaker = avgMaker;
            } else {
                throw ModelException(routine, "Expected Average Perf Maker but none supplied!");
            }
        } else if (finalPerfMakerType==FNL_PERF_TYPE_KNOCKIN) {
            if (knockInMaker.get()) {
                realMaker = knockInMaker;
            } else {
                throw ModelException(routine, "Expected knock in perf Maker but none supplied!");
            }
        } else if (finalPerfMakerType==FNL_PERF_TYPE_VAN) {
            if (knockInMaker.get()) {
                realMaker = knockInMaker;
            } else {
                throw ModelException(routine, "Expected knock in perf Maker but none supplied!");
            }
        } else {
            throw ModelException(routine, "Unrecognised Rebate Maker " + finalPerfMakerType);
        }
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FinalPerfMakerWrapper, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IFinalPerfMaker);
        EMPTY_SHELL_METHOD(defaultFinalPerfMakerWrapper);
        FIELD(finalPerfMakerType, "Average or KnockIn");
        FIELD(avgMaker,  "Average Perf Maker");
        FIELD_MAKE_OPTIONAL(avgMaker);
        FIELD(knockInMaker,  "Knock-In perf Maker");
        FIELD_MAKE_OPTIONAL(knockInMaker);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
     
    // for reflection
    FinalPerfMakerWrapper(): CObject(TYPE){}

    FinalPerfMakerWrapper(CClassConstSP clazz): CObject(clazz){}

    static IObject* defaultFinalPerfMakerWrapper(){
        return new FinalPerfMakerWrapper();
    }

    /** Override clone method to make sure realMaker references
        the maker inside the clone and not the original */
    IObject* clone() const{
        // first clone all the registered fields
        IObject*  copy = CObject::clone();
        FinalPerfMakerWrapper* rmw = dynamic_cast<FinalPerfMakerWrapper*>(copy);
        if (!rmw){
            throw ModelException("FinalPerfMakerWrapper::clone"); // shouldn't happen
        }
        rmw->validatePop2Object();
        return copy;
    }

};

typedef smartPtr<FinalPerfMakerWrapper> FinalPerfMakerWrapperSP;

CClassConstSP const FinalPerfMakerWrapper::TYPE =
CClass::registerClassLoadMethod("FinalPerfMakerWrapper", 
                                typeid(FinalPerfMakerWrapper), load);



DRLIB_END_NAMESPACE

