//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : EqLinkCashFlow.hpp
//
//   Description   Coupon stream linked to Equity generalized performance
//
//
//   $Log: EqLinkCashFlow.cpp,v $
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/EqLinkCashFlow.hpp"
#include "edginc/UtilFuncs.hpp"

DRLIB_BEGIN_NAMESPACE

EqLinkCashFlow::EqLinkCashFlow(const IDoubleArrayModifierMakerSP eqLinkPerf,
                                const DateTimeArray observDates): CObject(TYPE), 
                                eqLinkPerf(eqLinkPerf), observDates(observDates){
            validatePop2Object();
};

EqLinkCashFlow::EqLinkCashFlow(const IDoubleArrayModifierMakerSP eqLinkPerf,
                                const DateTimeArray observDates,
                                const DateTimeArray paymentDates): CObject(TYPE), 
                                eqLinkPerf(eqLinkPerf), observDates(observDates), paymentDates(paymentDates){
            validatePop2Object();
};

void EqLinkCashFlow::validatePop2Object() {
    static const string method("EqLinkCashFlow::validatePop2Object");
    try {

        if (paymentDates.size() > 0){
            if (observDates.size() != paymentDates.size()){
                throw ModelException(method, "paymentDates size[" + Format::toString(paymentDates.size()) + 
                                             "] should be empty or same to observDates size[" 
                                             + Format::toString(observDates.size()) + "].");
            }
            for (int i=0; i<observDates.size(); i++){
                if (observDates[i] > paymentDates[i]){
                    throw ModelException(method, "paymentDates[" + paymentDates[i].toString() 
                                                   + "] should be equal or later than observDates size[" 
                                                   + observDates[i].toString() + "].");
                }
            }            
        }

//        SimpleDoubleArray         dummy(1,0.0);    
//        IDoubleArrayModifierSP    dummy_perf = IDoubleArrayModifierSP(eqLinkPerf->getModifier(&dummy));
//        if (PerfTypeNothingMaker::TYPE->isInstance(eqLinkPerf)){
//        if (PerfTypeNothingMaker::TYPE->isInstance(eqLinkPerf)){
//            throw ModelException(method, "PerfType = NONE is not allowed.");
//        }

    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}
    
DateTimeArray EqLinkCashFlow::getObservDates() const
{
	return observDates;
}

DateTimeArray EqLinkCashFlow::getPaymentDates() const
{
	return paymentDates;
}


IDoubleArrayModifierMakerSP EqLinkCashFlow::getEqLinkPerf()
{
	return eqLinkPerf;
}

// return KNOWN_CASH_FLOW....
// when the observation is not finished, it retuns 0 cash flow to meet MO's request.
CashFlowArray EqLinkCashFlow::getKnownCashFlow(CashFlowArray samples, 
                                                 DateTime valDate, 
                                                 double refLevel, 
                                                 const InstrumentSettlementSP  instSettle,
                                                 const CAssetWrapper    asset,
                                                 bool  isExcludeKnownToday){
    static const string method("EqLinkCashFlow::setKnownCashFlow");
    try {        
        int nObserv = observDates.size();
        CashFlowArray cfl(0);
        DateTimeArray sampleDates = CashFlow::dates(samples);
        if (sampleDates.size()>=observDates.size()){
            vector<int> index = DateTime::getIndexes(sampleDates, observDates);

            SimpleDoubleArray         perfs(nObserv,0.0);    // equity performances to be processed and aggregated for eqLinkCpn
            IDoubleArrayModifierSP    performance = IDoubleArrayModifierSP(eqLinkPerf->getModifier(&perfs));
            
            DateTime payDate;
            for (int i=0; i<observDates.size(); i++){
                if (paymentDates.size()>0)
                    payDate = paymentDates[i];
                else
                    payDate = instSettle->settles(samples[index[i]].date, asset.get());
//                if (  ( isExcludeKnownToday && observDates[i] <  valDate && valDate < payDate)
//                    ||(!isExcludeKnownToday && observDates[i] <= valDate && valDate < payDate) ){
                if (  ( isExcludeKnownToday && observDates[i] <  valDate)
                    ||(!isExcludeKnownToday && observDates[i] <= valDate) ){
                    if (valDate == observDates[i]) // use Spot if valDate is observDate Timing.
                        perfs[i] = asset->getSpot()/refLevel;
                    else
                        perfs[i] = samples[index[i]].amount/refLevel;
                    performance->apply(i);
                    double payment = perfs[i];                
                    cfl.push_back(CashFlow(payDate, payment));
                }else{ // return 0 value, but dates.
                    cfl.push_back(CashFlow(payDate, 0.0));
                }
            }    
        }
        return cfl;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


EqLinkCashFlow::EqLinkCashFlow():CObject(TYPE){
    paymentDates = DateTimeArray(0);
}; 

/** Invoked when Class is 'loaded' */
void EqLinkCashFlow::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    clazz->setDescription("CallableEquityKOSwap Equity Linked CashFlows");
    REGISTER(EqLinkCashFlow, clazz);
    SUPERCLASS(CObject);
    FIELD(eqLinkPerf,       "performace of each observation date");
    FIELD_MAKE_OPTIONAL(eqLinkPerf);
    FIELD(observDates,        "maturities of each eq link coupons");
    FIELD_MAKE_OPTIONAL(observDates);
    FIELD(paymentDates,       "coupon payment dates");
    FIELD_MAKE_OPTIONAL(paymentDates);
    EMPTY_SHELL_METHOD(defaultEqLinkCashFlow);
}

IObject* EqLinkCashFlow::defaultEqLinkCashFlow(){
    return new EqLinkCashFlow();
}

CClassConstSP const EqLinkCashFlow::TYPE = CClass::registerClassLoadMethod(
    "EqLinkCashFlow", typeid(EqLinkCashFlow), EqLinkCashFlow::load);

//////////////////////////////////////////////////////////////
//      public maker class for with KO info.  
//      No plan to expose to public interface. 
//      Just for CallableKOSwap, not for EGKBond nor any others.
//////////////////////////////////////////////////////////////
EqCpnKOMaker::EqCpnKOMaker(const EqLinkCashFlowSP eqCpn,
                                const DateTime  initialAccrueDate,
                                const KOStubRuleSP koRule,
                                bool inAdvance): CObject(TYPE), eqCpn(eqCpn),inAdvance(inAdvance),
                                koRule(koRule),initialAccrueDate(initialAccrueDate){
            eqCpn->validatePop2Object();
};

void EqCpnKOMaker::validatePop2Object() {
    static const string method("EqCpnKOMaker::validatePop2Object");
    try {
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}
    
DateTimeArray EqCpnKOMaker::getObservDates(){
	return eqCpn->getObservDates();
}

DateTimeArray EqCpnKOMaker::getPaymentDates(){
	return eqCpn->getPaymentDates();
}


IDoubleArrayModifierMakerSP EqCpnKOMaker::getEqLinkPerf()
{
	return eqCpn->getEqLinkPerf();
}

DateTime EqCpnKOMaker::getInitialAccrueDate()
{
	return initialAccrueDate;
}

KOStubRuleSP EqCpnKOMaker::getKOStubRule()
{
	return koRule;
}

CashFlowArray EqCpnKOMaker::getKnownCashFlow(CashFlowArray samples, 
                                                 DateTime valDate, 
                                                 double refLevel, 
                                                 const InstrumentSettlementSP  instSettle,
                                                 const CAssetWrapper    asset,
                                                 bool isExcludeKnownToday){
    return eqCpn->getKnownCashFlow(samples, valDate, refLevel, instSettle, asset, isExcludeKnownToday);
}

// return next payment date corresponding to hit Date.
// currently, "S" is not supported.  So always payment date is next pay date.
bool EqCpnKOMaker::getPayDate(const DateTime hitDate, const InstrumentSettlementSP  instSettle,
                                  DateTime* payDate) const
{
    static const string method("EqCpnKOMaker::getPayDate");
    try {
        DateTimeArray payDts = eqCpn->getPaymentDates();
        int iPayIdx;
        bool found = false;
        if (payDts.size() == 0){     
            // use instSettle.
            *payDate = instSettle->settles(hitDate, 0);            
            found = true;            
        }else{
            iPayIdx = Neighbour(hitDate, payDts,0,payDts.size()-1,1);
            if (iPayIdx>=0){
                *payDate = payDts[iPayIdx];
                found = true;
            }
        }
        return found;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

EqCpnKOMaker::EqCpnKOMaker():CObject(TYPE), inAdvance(inAdvance), 
                                    initialAccrueDate(DateTime(0,0)){}; 

/** Invoked when Class is 'loaded' */
void EqCpnKOMaker::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    clazz->setDescription("Equity Linked CashFlows with KO Info");
    REGISTER(EqCpnKOMaker, clazz);
    SUPERCLASS(CObject);
    FIELD(eqCpn,"Equity Linked Cash Flow");
    FIELD_MAKE_OPTIONAL(eqCpn);
    FIELD(initialAccrueDate, "initial accrue start date.");
    FIELD_MAKE_OPTIONAL(initialAccrueDate);
    FIELD(koRule, "ko rule for early terminate event");
    FIELD_MAKE_OPTIONAL(koRule);
    FIELD(inAdvance,"inAdvance = true, then if KO & obs is same date, "
                    "coupon would be cancelled even for stubRule = 'B'.");
    FIELD_MAKE_TRANSIENT(inAdvance);
    EMPTY_SHELL_METHOD(defaultEqCpnKOMaker);
}

IObject* EqCpnKOMaker::defaultEqCpnKOMaker(){
    return new EqCpnKOMaker();
}

CClassConstSP const EqCpnKOMaker::TYPE = CClass::registerClassLoadMethod(
    "EqCpnKOMaker", typeid(EqCpnKOMaker), EqCpnKOMaker::load);


//////////////////////////////////////////////////////////////
//      Internal Class
//      EqKOLeg (with KO rule)  
//////////////////////////////////////////////////////////////

EqCpnKOLeg::EqCpnKOLeg():isNull(true),inAdvance(false),payIndex(0){};

EqCpnKOLeg::EqCpnKOLeg(const EqLinkCashFlowSP elc,
                         SimpleDoubleArray& perfs,
                         DateTime valueDate,
                         KOStubRuleSP koRule,
                         bool inAdvance,
                         DateTime initialAccrueDate): // first accrue start Date. 
                            elc(elc), isNull(false),valueDate(valueDate),
                            isAlreadyApplied(false),inAdvance(inAdvance),
                            koRule(koRule),initialAccrueDate(initialAccrueDate),payIndex(0){
             performance = IDoubleArrayModifierSP(elc->getEqLinkPerf()->getModifier(&perfs));        
             couponPayments = PerfContainer(&perfs);
             remainedValue = 0.0;
}

// constructor wirh EqCpnKOMaker
EqCpnKOLeg::EqCpnKOLeg(const EqCpnKOMakerSP elcKO,
                         SimpleDoubleArray& perfs,
                         DateTime valueDate):
                            elc(elcKO->eqCpn), isNull(false),valueDate(valueDate),
                            isAlreadyApplied(false),koRule(elcKO->koRule), inAdvance(elcKO->inAdvance),
                            initialAccrueDate(elcKO->initialAccrueDate),payIndex(0){
             performance = IDoubleArrayModifierSP(elc->getEqLinkPerf()->getModifier(&perfs));        
             couponPayments = PerfContainer(&perfs);
             remainedValue = 0.0;
}

double EqCpnKOLeg::getPerf(const int idx){
    double value = 0.0;
    if (idx>=0){
        performance->apply(idx);
        value = couponPayments.getPerf(idx);
    }
    return value;
}

// this is for Tree base
double EqCpnKOLeg::getValue(int step){
    static const string method = "EqCpnKOLeg::getValue";
    try{    
        double value = 0.0;
        // validation
        if (obsIndex.size() <= step)
            throw ModelException(method, "not correctly set up obsIndex in EqCpnKOLeg (internal object error).");
        
        int iObs = obsIndex[step];
        
        // validation
        if (discFacts.size() <= iObs)
            throw ModelException(method, "not correctly set up discFacts in EqCpnKOLeg (internal object error).");

        // return value only when step is observDate.  Otherwise, return value=0.0.
        if (monStepMap[step] == 0){
            value = getPerf(iObs) * discFacts[iObs];
        }
        return value;        
    }
    catch (exception& e)
    {            
        throw ModelException(e, method);
    }
}


//class getValue_oper : public SliceMarker< getValue_oper >
EqCpnKOLeg::getValue_oper::getValue_oper(const TreeSlice & spot,
                      int step,
                      double refLevel,
                      double &eqPerfAddress,
                      EqCpnKOLeg & eqCpn)
                :
                spot( spot ),
                step( step ),
                refLevel( refLevel ),
                eqCpn( eqCpn )
{
    pPerf = &eqPerfAddress;
}

double EqCpnKOLeg::getValue_oper::apply( double s ) const{
        *pPerf = s/refLevel;
        return eqCpn.getValue(step);
}


// This is for future, to set up "S" rule.  Now, it just return 1 (B) or 0 (N).
// When the time step is observation date, then
// "B" : Do not cancel the cpn
// "N" : Cancel the cpn.
// Not Only exact time point, the same date is also considered as "KO" for "factors"
// because of Tree usage.
void EqCpnKOLeg::setKOFactor(const DateTimeArray timeLine, 
                             const vector<bool>& stepIsKO, 
                             vector<double>& factors){
    static const string method = "EqCpnKOLeg::setKOFactor";
    try{
        int numSize = timeLine.size();
        int i;
        factors.clear();

        if (koRule->getKOStubRule() == "B")
            factors.resize(numSize+1, 1.0);         // always 1.0 for "B"
        else if (koRule->getKOStubRule() == "N")
            factors.resize(numSize+1, 0.0);         // always 0.0 for "N"
        else
            throw ModelException(method, "koStubRule should be 'N' or 'B'.");

        if (koRule->getpayTimingRule() != "AsSchedule")   
            throw ModelException(method, "payTimingRule should be 'AsSchedule'.");

        // for "B", all dates will have 1.0 (no cancel), 
        // if the time point is strictly later than initialAccrueDate
        // or could be 0 when inAdvance = true
        if (koRule->getKOStubRule() == "B"){        
            // modify the KO factor for initialAccrueDate
            i=0;
            while (timeLine[i] <= initialAccrueDate){    // same date of initDate is exclusive.
                factors[i] = 0.0;
                i++;
            }

            // modify the KO factor for 'inAdvance
            if (inAdvance){
                // validate the array size.
                if (timeLine.size() != monStepMap.size()-1){
                    throw ModelException(method, "Internal Error. timeLine size[="
                        + Format::toString(timeLine.size()) 
                        + "] should be equal to monStepMap size - 1 [="
                        + Format::toString(monStepMap.size()-1) + "].");
                }
                if (timeLine.size() != stepIsKO.size()){
                    throw ModelException(method, "Internal Error. timeLine size[="
                        + Format::toString(timeLine.size()) 
                        + "] should be equal to stepIsKO size[="
                        + Format::toString(stepIsKO.size()) + "].");
                }

                int iStepMon = 0;
                for(int iStep = 0; iStep<timeLine.size(); iStep++, iStep += monStepMap[iStep] ){
                    iStepMon = iStep;
                    // give 0 for all time points on the KO date && monitoring (observation) date.
                    while (iStep>0 && stepIsKO[iStep] && timeLine[iStep].equals(timeLine[iStepMon], false) ){
                        factors[iStep] = 0.0;
                        iStep--;
                    }                    
                    iStep = iStepMon;
                }
            }
        }
    }
    catch (exception& e)
    {            
        throw ModelException(e, method);
    }

}

// set flag whether the date step is monitoring date or not.
IntArray EqCpnKOLeg::setMonStep(const DateTimeArray timeLine){
    static const string method = "EqCpnKOLeg::setMonStep";
    try{    
        int numSize = timeLine.size();
        int i;
        DateTimeArray observDates = elc->getObservDates();
        // trim the observDates for tree products.  
        //  (MC always has past in timeLine, but not in Tree)
        DateTimeArray trimObsDts;
        int iObs = 0;
        for (i=0; i<observDates.size(); i++){
            if (observDates[i] >= timeLine[0])
                trimObsDts.push_back(observDates[i]);
            else
                iObs ++;        //skip the past observ Date.

            if (observDates[i] > timeLine[numSize-1])
                throw ModelException(method, "The last simulation date["
                                             + timeLine[numSize-1].toString() 
                                             + "] should be later than or equal to \n"
                                             + "the last of eqLinkCpn observation date["
                                             + observDates[i].toString()
                                             + "].");
        }


//        isMonStep.resize(numSize+1,false);
//        int j = timeLine[0].findUpper(observDates);  ;
//        for (i=0; i< numSize && j<observDates.size(); i++){
//            if (timeLine[i] == observDates[j]){
//                isMonStep[i] = true;
//                j++;
//            }
//        }
//        if ( j != observDates.size())
//            throw ModelException(method, "The last simulation date (timeLine) is before the last of eqLinkCpn observation date.");

        // create the observatoin map
        bool isTrivial;
        monStepMap = DateTime::createMapping(timeLine,trimObsDts,isTrivial);

        // make a IndexArray to know the position of iObs against iStep.
        // currently, obsIndex has previous observation index if the time step is not observation point.
        // this is strange.....
        obsIndex.clear();
        obsIndex.resize(monStepMap.size()-1);
        for (int iStep = 0; iStep < monStepMap.size()-1;  iStep++){
            obsIndex[iStep] = iObs;
            if (monStepMap[iStep] == 0)
                iObs++;
        }

        return monStepMap;
    }
    catch (exception& e)
    {            
        throw ModelException(e, method);
    }

}


// set up the discount factor.  from payment date to observ date (mainly for tree).
// check whether there is paydates or not, and judge use instSettle or not.
void EqCpnKOLeg::setDiscFacts(const DateTimeArray timeLine, 
                              const YieldCurveConstSP discount,
                              const InstrumentSettlementSP  instSettle,
                              const CAssetWrapper asset){
    static const string method = "EqCpnKOLeg::setDiscFacts";
    try{    
        discFacts.clear();
        DateTimeArray cpnDts = elc->paymentDates;
        DateTimeArray obsDts = elc->observDates;
        if (cpnDts.size() > 0)
            setDiscFacts(timeLine, discount);
        else{            
            // make disc factors array on observe dates.
            discFacts.resize(obsDts.size());
            for (int i=0;i<obsDts.size();i++){
                discFacts[i] = instSettle->pvAdjust(obsDts[i],
                                                    discount.get(), 
                                                    asset.get());
            }
        }
    }
    catch (exception& e)
    {            
        throw ModelException(e, method);
    }
}

// set up the discount factor.  from payment date to observ date (mainly for tree).
// only available when payment dates are given.
void EqCpnKOLeg::setDiscFacts(const DateTimeArray timeLine, const YieldCurveConstSP discount){
    static const string method = "EqCpnKOLeg::setDiscFacts";
    try{    
        discFacts.clear();
        DateTimeArray cpnDts = elc->paymentDates;
        DateTimeArray obsDts = elc->observDates;
        if (cpnDts.size() != obsDts.size())
            throw ModelException(method, "number of paymentDates and observDates are should be same!!");

        // make disc factors array on observe dates.
        discFacts.resize(obsDts.size());
        for (int i=0;i<obsDts.size();i++){
            discFacts[i] = discount->pv(obsDts[i], cpnDts[i]);
        }
    }
    catch (exception& e)
    {            
        throw ModelException(e, method);
    }
}


// set up the payment date by using instSettle.  Ignore the EqCashFlow's paymentDates
void EqCpnKOLeg::setPayDates(const DateTimeArray timeLine,const InstrumentSettlementSP  instSettle) 
{
    ///////////  NOT TESTED!!! ////////////////////
    throw ModelException("EqCpnKOLeg::setPayDates", "not Tested. ");
    payDates.resize(0); // clear
    for (int i=0; i<timeLine.size(); i++){
        DateTime payDate = instSettle->settles(timeLine[i], 0);
        payDates.push_back(payDate);
    }    	
}

// set up the payment by using EqCashFlow's paymentDates
void EqCpnKOLeg::setPayDatesAndDiscFacts(const DateTimeArray timeLine, const YieldCurve* discount, const DateTime toDate) 
{
    static const string method = "EqCpnKOLeg::setPayDatesAndDiscFacts";
    try{    
        // need validation about monStepMap
        //if (monStepMap.size()>0)

        payDates.resize(0); // clear
        DateTimeArray cpnDts = elc->paymentDates;
        int i;
        // first, make disc factors array on coupon payment dates.
        discFacts.resize(cpnDts.size());
        for (i=0;i<cpnDts.size();i++){
            discFacts[i] = discount->pv(toDate, cpnDts[i]);
        }
        // second, make payment dates array corresponding simulation date.
        // currently, all payment occurs at next coming payment dates.
        // need to work to cover more various payment case.
        // NB :  not related to observDate!!
        int iCpn = 0;
        for (i=0; i<timeLine.size(); i++){
            int j = timeLine[i].findUpper(cpnDts);
            if (j == cpnDts.size()){
                throw ModelException("EqCpnKOLeg::setPayDatesAndDiscFacts",
                                    "simulatoin dates [" + timeLine[i].toString() + 
                                    "] is later than last coupon payment dates [" +
                                    cpnDts.back().toString() +
                                    "]. check the date and timing.");
            }
            payDates.push_back(cpnDts[j]);
        }   

        // finally, make a unique payment dates array and mapping from observeDates to new unique array.
        DateTimeArray uniqPayDates = elc->getPaymentDates();
        DateTime::removeDuplicates(uniqPayDates,false);
        payIndex.clear();
        payIndex.resize(uniqPayDates.size());
        payIndex[0] = 0;
        for (i=1; i<cpnDts.size(); i++){
            if (cpnDts[i-1] == cpnDts[i])
                payIndex[i] = payIndex[i-1];
            else
                payIndex[i] = payIndex[i-1]+1;
        }
        if (payIndex.size() != uniqPayDates.size())
            throw ModelException(method, "model error.  iCpn shouldn't be over the cpnDts size.");
    }
    catch (exception& e)
    {            
        throw ModelException(e, method);
    }
}

// for MC product, calculate the range coupon value.
double EqCpnKOLeg::getValue(const int   startStep,            //(I) start step of path
                            const int   endStep,              //(I) calculate up to endStep                             
                            const int   nbAssets){            
    static const string method = "EqCpnKOLeg::getValue";
    try{
        double value = 0.0; //remainedValue;        // start from the not yet paid coupon.
        int iStep;
//        int nbAssets=pathGen->NbSimAssets();

        // add the range coupon before endStep        
        //iStep = pathGen->begin(0) + monStepMap[pathGen->begin(0)];
        iStep = startStep;
        int iObs = obsIndex[iStep];
        int iPaySize = elc->paymentDates.size();
        for(; iStep<endStep; iStep += monStepMap[iStep+1] + 1, iObs++){
            if (iStep+monStepMap[iStep] == monStepMap.size()-1){
                // it reach to the end of monStepMap (timeLine +1).  
                // It means the final timeLine is not observation date (e.g. inAdvance)
                // and no more equity observation remained.
                break;
            }
            if (!isAlreadyApplied)
                performance->apply(iObs);
            if (elc->paymentDates[iObs]>valueDate)
                value += couponPayments.getPerf(iObs) * discFacts[iObs];
        }      
        value += remainedValue;
        return value;
    }
    catch (exception& e)
    {            
        throw ModelException(e, method);
    }
}

// for MC product, calculate the range coupon value.
// <NB> : discFactsSV is not guranteed to be same time step to paymentDate!!
// This is used from only RYM, and RYM mades discFactsSV from this class's paymentDates.
// Maybe we'd like to SVGen and SV in this class, or pass the DiscFact array as doubleArray.
double EqCpnKOLeg::getValue(const int   startStep,            //(I) start step of path
                            const int   endStep,              //(I) calculate up to endStep                             
                            const int   nbAssets,
                            SVDiscFactorSP discFactsSV){            
    static const string method = "EqCpnKOLeg::getValue";
    try{
        double value = 0.0; //remainedValue;        // start from the not yet paid coupon.
        int iStep;
//        int nbAssets=pathGen->NbSimAssets();

        // add the range coupon before endStep        
        iStep = startStep;
        for(; iStep<endStep; iStep++, iStep += monStepMap[iStep]){
            int iObs = obsIndex[iStep];        
            if (!isAlreadyApplied)
                performance->apply(iObs);
            if (elc->paymentDates[iObs]>valueDate){
                value += couponPayments.getPerf(iObs) * discFactsSV->path()[payDateIndex(iStep)];
            }
        }      
        value += remainedValue;
        return value;
    }
    catch (exception& e)
    {            
        throw ModelException(e, method);
    }
}

// return the payment date corresponding to KO dates)
DateTimeArray EqCpnKOLeg::getPayDatesOnTimeLine() const
{
    return payDates;
}

// return the payment date of elc, after removing duplicate of input.
DateTimeArray EqCpnKOLeg::getPayDates() const
{
    DateTimeArray uniqPayDates = elc->getPaymentDates();
    DateTime::removeDuplicates(uniqPayDates,false);
    return uniqPayDates;
}
// return index of payment dates array (after removed duplicated one),
// which corresponding to the timeLine[step]
int EqCpnKOLeg::payDateIndex(int step)
{
    int iObs = obsIndex[step];
    return payIndex[iObs];
}

// return payment date,
// which corresponding to the timeLine[step]
DateTime EqCpnKOLeg::payDateOnKO(int step)
{
    int payIdx = obsIndex[step];
    return payDates[payIdx];
}

                
// for MC version
// for the coupon observed are going to use discount factor of today, not state variable!!
// for SRM, may need to review!!
void EqCpnKOLeg::makeKnownCashFlow(const int   startStep,
                                   const int   endStep,  
                                   const int   nbAssets,
                                   const bool  hasFuture){          //(I) if hasFuture = false, no need to calc Ramained Value
    static const string method = "EqLinkCashFlow::makeKnownCashFlow";
    try{
        double value =0.0;    
        int iStep, j=0;
        //int nbAssets=pathGen->NbSimAssets();
        
        // initialization
        knownCFL = CashFlowArray(0);
        
        // add the range coupon before endStep
        // not payDates, whose payDate are corresponding to KO event,
        // but use the original pay dates of ELC.  
        iStep = startStep; //pathGen->begin(0) + monStepMap[pathGen->begin(0)];
        int iObs = obsIndex[iStep];
        for(; iStep<endStep; iStep += monStepMap[iStep+1] + 1, iObs++){
            performance->apply(iObs);
            value = couponPayments.getPerf(iObs);
            knownCFL.push_back(CashFlow(elc->paymentDates[iObs], value));
            if (elc->paymentDates[iObs]>valueDate)
                remainedValue += value * discFacts[iObs];
        }            
        
        // if option terminated (hasFuture==false), pathGen->begin(0) returns 0.
        // rather than endIdx.  So, it will calculate all remained coupons at getValues.
        // Thus, we need to skip "performance->apply" , otherwise
        // it will apply twice (i.e. getValue used the underlying = performed coupon)  
        if (!hasFuture){
            remainedValue = 0.0;        
            isAlreadyApplied = true;
        }
        CashFlow::aggregate(knownCFL);
    }
    catch (exception& e)
    {            
        throw ModelException(e, method);
    }
}

    // return a copy of knownCFL
CashFlowArraySP EqCpnKOLeg::getKnownCashFlows()
{
    CashFlowArraySP cfl(new CashFlowArray(0));
    for (int i=0; i<knownCFL.size();i++)
        cfl->push_back(knownCFL[i]);
    return cfl;
}

// a class to get the genereric performance results.
EqCpnKOLeg::PerfContainer::PerfContainer():components(0){};
EqCpnKOLeg::PerfContainer::PerfContainer(IDoubleArray* components):components(components) {};

double EqCpnKOLeg::PerfContainer::getPerf() {
    double val = 0.0;
    for(int i=0; i<components->size(); i++) {
        val += (*components)[i];
    }
    return val;
}
double EqCpnKOLeg::PerfContainer::getPerf(int index) {
    return (*components)[index];
}

DRLIB_END_NAMESPACE
