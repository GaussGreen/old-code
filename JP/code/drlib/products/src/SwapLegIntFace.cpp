//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SwapLegIntFace.cpp
//
//   Description   Swap Leg Inter Face, to make libor and fixed cash flow.
//
//
//	 Revision 1.37 2006/02/24 klaw
//	 Revised definition of Stub rule "B"
//
//   $Log: SwapLegIntFace.cpp,v $
//   Revision 1.36  2005/05/05 09:35:59  snesbitt
//   Fix some unix compilation problems.
//   Small changes after code review
//
//   Revision 1.35  2005/05/03 16:44:22  snesbitt
//   Add state var support
//   I suspect a memory leak
//
//   Revision 1.34  2005/05/03 12:53:29  snesbitt
//   Fix loop upper bound in getKOPV to avoid ABRs
//
//   Revision 1.33  2005/05/03 09:15:46  qhou
//   fix mem access in BarrierPay when assign payPerStep
//
//   Revision 1.32  2005/04/29 10:26:08  qhou
//   fix bug -- need to use DCC to calc fix leg coupon in barrierPay
//
//   Revision 1.31  2005/04/27 10:12:50  qhou
//   revise/fix Libor/Fix leg barrierPay preprocess to link mon date with coupon periods
//   and handle properly known cf
//
//   Revision 1.30  2005/04/26 08:57:01  Kkitazaw
//   Fixing getCashFlowOnKO and related to koStub = "S".
//
//   Revision 1.29  2005/04/25 04:44:17  qhou
//   add relevantMonDts function. cpn calc for stub N/B in preprocess
//
//   Revision 1.28  2005/04/22 09:42:00  Kkitazaw
//   Added dcc factor for FixedLeg getCashFlowArray.
//
//   Revision 1.27  2005/04/20 12:11:36  mrobson
//   Code was using member function as a bool instead of actually calling the
//   function
//
//   Revision 1.26  2005/04/20 08:40:19  Kkitazaw
//   renewal KOSettle and getKOPV to do koStubRule = "S".
//
//   Revision 1.25  2005/04/18 10:57:59  qhou
//   add createBarrierPay for libor/fixed legs
//
//   Revision 1.24  2005/04/04 13:01:03  jmusset
//   Removed myKibor and myDiscCurve
//
//   Revision 1.23  2005/03/30 08:49:51  aswain
//   redo KOSettle again
//
//   Revision 1.22  2005/03/29 11:40:06  aswain
//   move KOSettle into source file
//
//   Revision 1.21  2004/11/02 05:16:23  Kkitazaw
//   Add myLibor, myDiscCurve, koSettle, knownCashFlow.
//   Add makeKnownCashFlow, setKOSettle, getKnownCashFlow(), checkMarket, setCouponCurve,
//   getPVsAlongKOTimeStep, getNextKOCoupnPV.
//   Add KOSettle class.
//
//   Revision 1.20  2004/10/13 09:51:28  Kkitazaw
//   Remove valDate from makePayStream.
//    Added LiborLeg::getKOPV() and getCashFlowArray which used internal valueDate.
//   Added validation. getKOPV() "S" is not allowed.  Fixing rate cannot be zero in the past.
//
//   Revision 1.19  2004/09/10 18:58:17  snesbitt
//   Adding StreamRebate to IRebate
//
//   Revision 1.18  2004/09/06 17:00:36  jmusset
//   Fixing Memory Leakn in LiborLeg:getKOPV (with Keiji).
//
//   Revision 1.17  2004/08/22 23:49:15  Kkitazaw
//   getKOPV, fixed the memory bug.
//
//   Revision 1.14  2004/08/09 10:37:57  Kkitazaw
//   change FixedLeg::getNextCoupon, to remove accrued so far case. Fixed getKOPV (Not Tested, Yet!)
//
//   Revision 1.13  2004/07/15 10:03:59  jmusset
//   Used Curve rate instead of discount curve to create the FloatRates object
//
//   Revision 1.12  2004/07/13 10:02:04  jmusset
//   Added coupon curve and value date fields to LiborLeg class
//   Moved some code from the hpp to the cpp
//
//   Revision 1.11  2004/06/14 13:33:20  kbeguir
//   added condition to handle non existent LiborLeg in validatePop2Object()
//
//   Revision 1.10  2003/12/03 06:21:25  Kkitazaw
//   Added getKOPV function for both , and getAccrued, getNextCoupon for Fixed Leg.
//
//   Revision 1.9  2003/07/23 04:15:52  Kkitazaw
//   Tiny Change in error message in Validate.
//
//   Revision 1.8  2003/06/04 10:47:03  Kkitazaw
//   Added new flavour getPV for FixedLeg, too.
//
//   Revision 1.7  2003/05/13 14:07:01  aswain
//   LiborLeg::getPV now includes payments on fromDate provided they're after valueDate
//
//   Revision 1.6  2003/04/04 10:07:55  aswain
//   added new flavour of getPV + new getAccrued
//
//   Revision 1.5  2002/12/19 11:28:58  kkitazaw
//   Make Compouding Optional.
//
//   Revision 1.4  2002/12/19 09:20:38  kkitazaw
//   Fixed warning bug.
//
//   Revision 1.3  2002/10/25 12:09:46  mrobson
//   Now compiles under gcc
//
//   Revision 1.2  2002/10/25 10:32:49  kkitazaw
//   Remove prefix of LiborLeg and FloatLeg .
//
//   Revision 1.1  2002/10/24 09:36:36  kkitazaw
//   Contains FixedLeg and LiborLeg as cash flow interface for ExtendableNote, StepDownBond and so on.
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#define QLIB_SWAPLEGINTFACE_CPP
#include "edginc/SwapLegIntFace.hpp"
#include "edginc/Maths.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/BadDayConventionFactory.hpp"
#include "edginc/BarrierUtil.hpp"
#include "edginc/PhysicalDelivery.hpp"
#include "edginc/RateConversion.hpp"
#include "edginc/CompoundBasis.hpp"

DRLIB_BEGIN_NAMESPACE


//------------------------------------------------------------------//
//  LiborLeg : The interface for Libor Leg  (deterministic rate)
//------------------------------------------------------------------//

/** overrides default */
void LiborLeg::validatePop2Object(){

    static const string method = "LiborLeg::validatePop2Object";
    
    try
    { 
        int numItems = PayDates.size();
        bool allZero = ((numItems==0) && (RefixDates.size()==0) && (Spreads.size()==0) && (Weights.size()==0) && (AccrualDates.size()==0));
        
        if ((numItems != RefixDates.size() ||
             numItems != Fixings.size() ||
             numItems != Spreads.size() ||
             numItems != Weights.size() ||
             numItems != (AccrualDates.size()-1)) && !allZero)
        {
            throw ModelException(method,
                                 "The number of pay dates (size = " + Format::toString(PayDates.size()) +
                                 "), \nrefix dates (size =" + Format::toString(RefixDates.size()) +
                                 "), \nfixings (size =" + Format::toString(Fixings.size()) + 
                                 "), \nspreads (size =" + Format::toString(Spreads.size()) + 
                                 "), \nWeights (size =" + Format::toString(Weights.size()) + 
                                 ") must be equal.  \nThe number of accrualDates (size =" + Format::toString(AccrualDates.size()) +
                                 ") must be one greater than this. \n");
        }
        
        for (int i=1; i<numItems; i++)
        {// checke cronical date array or not
            if (RefixDates[i-1]>RefixDates[i])
                throw ModelException(method, "RefixDates (" + Format::toString(i) + ") [ " + RefixDates[i-1].toString() + 
                                     " ] should be ealier than next date [ " + RefixDates[i].toString() + " ] \n");
            if (PayDates[i-1]>PayDates[i])
                throw ModelException(method, "PayDates (" + Format::toString(i) + ") [ " + PayDates[i-1].toString() + 
                                     " ] should be ealier than next date [ " + PayDates[i].toString() + " ] \n");
            if (AccrualDates[i-1]>AccrualDates[i])
                throw ModelException(method, "AccrualDates (" + Format::toString(i) + ") [ " + AccrualDates[i-1].toString() + 
                                     " ] should be ealier than next date [ " + AccrualDates[i].toString() + " ] \n");
        }
        
        if (Compounding.size()>0)
        {
            if (numItems != Compounding.size())
            {
                throw ModelException(method,
                                     "The number of Compounding (size =" + Format::toString(Compounding.size()) + 
                                     ") must be equal to " + Format::toString(numItems) + ".  ");
            }
            if (Compounding[numItems-1])
            {
                throw ModelException(method,
                                     "The last period cannot be compounded.");
            }
        }

        // store the original spreads
        orgSpreads = CDoubleArray(Spreads);
    }
    catch (exception& e) 
    {
        throw ModelException(e, method);
    }
    
}

// helpers
class LiborLegHelper {
public:
    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(LiborLeg, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(Theta::Shift);
        EMPTY_SHELL_METHOD(defaultLiborLeg);
        
        FIELD(RefixDates,   "refix dates for IR leg");
        FIELD(AccrualDates, "accrual dates for Libor");
        FIELD(PayDates,     "payment dates for Libor");
        FIELD(Fixings,      "fixings for Libor");
        FIELD(Spreads,      "spreads for Libor");
        FIELD(Weights,      "weighes for Libor");
        FIELD(Notionals,    "notionals for Libor");
        FIELD(PayDCC,       "The pay day count convention for Libor");
        FIELD(Compounding,  "compoundings for Libor");
        FIELD_MAKE_OPTIONAL(Compounding);
        FIELD(CompoundFlat, "compound flat");
        FIELD_MAKE_OPTIONAL(CompoundFlat);
        FIELD(RateDCC,      "The rate day count convention for Libor");
        FIELD(RateType,     "rate Type (e.g. 3M)");
        FIELD(BadDayConvention,  "bad day convention for libor leg");
        FIELD_MAKE_OPTIONAL(BadDayConvention);
        FIELD(couponCurve,  "Coupon curve");
        FIELD_MAKE_OPTIONAL(couponCurve);
        FIELD(valueDate,    "Value date");
        FIELD_MAKE_OPTIONAL(valueDate);
    }
    
    static IObject* defaultLiborLeg(){
        return new LiborLeg();
    }
};

CClassConstSP const LiborLeg::TYPE = CClass::registerClassLoadMethod(
    "LiborLeg", typeid(LiborLeg), LiborLegHelper::load);

// constructor
LiborLeg::LiborLeg(): CObject(TYPE)
{
    BadDayConvention = "NONE";
    CompoundFlat = false;
    hasKOSettle = false;
}

int LiborLeg::getSize()
{
    return PayDates.size();
}

// get the last coupon payment date
DateTime LiborLeg::getLastPayDate()
{
    return PayDates.back();
}

FloatRateSP LiborLeg::getFloatRate(const YieldCurve* discount){
    MaturityPeriodSP interval = MaturityPeriodSP(new MaturityPeriod(RateType));
    DayCountConventionSP liborRateDCC(DayCountConventionFactory::make(RateDCC));
    BadDayConventionConstSP badDayConvention(BadDayConventionFactory::make(BadDayConvention));
    HolidaySP hol(Holiday::weekendsOnly());
    /** if the coupon curve has been populated then we use it, 
        otherwise we use the discount curve */
    const YieldCurve* rateCurve = couponCurve.isEmpty()? 
        discount: couponCurve.get();

    FloatRateSP floater = FloatRateSP(new FloatRate(
                                          rateCurve,
                                          hol.get(),
                                          liborRateDCC.get(),
                                          badDayConvention.get(),
                                          interval.get(),
                                          interval.get(),
                                          false,
                                          0,
                                          &RefixDates));
    
    return floater;
}

// create a PayStream
PayStream* LiborLeg::makePayStream(const YieldCurve* discount) {
    static const string method("LiborLeg::makePayStream");
    try {
        DayCountConventionSP liborPayDCC(DayCountConventionFactory::make(PayDCC));
        BoolArray compound;
        for (int i = 0; i < Spreads.size(); i++) {
            if (!Compounding.empty()) {
                compound.push_back(Compounding[i]);
            }
            else {
                compound.push_back(false);
            }
        }
           
        FloatRateSP floater = this->getFloatRate(discount);

        PayStreamSP FloatPayStream(new PayStream(liborPayDCC.get(),
                                                 false,
                                                 floater.get(),
                                                 Notionals,              
                                                 RefixDates,
                                                 AccrualDates,
                                                 PayDates,
                                                 Fixings,
                                                 Spreads,
                                                 Weights,            
                                                 compound));  

        return FloatPayStream.release();
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// return the libor cash flow with internal value date and yield curve
// Need to call getMarket & setCouponCurve before.
CashFlowArrayConstSP LiborLeg::getCashFlowArray(){
    static const string method("LiborLeg::getCashFlwoArray");
    try {
        checkMarket();
        return getCashFlowArray(valueDate, couponCurve.get());
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }    
}

// return the libor cash flow.  The size is same to the paymentDate (input)
CashFlowArrayConstSP LiborLeg::getCashFlowArray(const DateTime& valueDate, 
                                                const YieldCurve* yield){    
    PayStreamSP FloatPayStream(makePayStream(yield));                                            
    CashFlowArrayConstSP libor = FloatPayStream->cashFlowList(valueDate, 0);
    return libor;
}

// return Cash Flow, but only cash flows later fromDate.
// Also, it's return the copy of cash flow
CashFlowArraySP LiborLeg::getCashFlowArray(const DateTime& valueDate, 
                                           const DateTime& fromDate, 
                                           const YieldCurve* yield){
    CashFlowArrayConstSP temp = getCashFlowArray(valueDate, yield);
    CashFlowArraySP libor (new CashFlowArray());
    for (int i=0; i<(*temp).size();i++){
        if((*temp)[i].date>fromDate)
            libor->push_back((*temp)[i]);
    }
    return libor;
}

// get PV from ValueDate to End, as of valueDate.
double LiborLeg::getPV(const DateTime& valueDate, const YieldCurve* discount){
    double result = 0.0;
    // calculate Libor Leg.	
    CashFlowArrayConstSP libor = getCashFlowArray(valueDate, 
                                    couponCurve.isEmpty()? discount: couponCurve.get());
    for (int i = 0; i < libor->size(); i++)
    {// pv the cash flow, include strickly future flows
        if ((*libor)[i].date > valueDate)
            result += (*libor)[i].amount * discount->pv(valueDate, (*libor)[i].date);
    }
    return result; 
}

// get pv of cashflows paying between fromDate & toDate as of fromDate
double LiborLeg::getPV(const DateTime&   valueDate,
                       const DateTime&   fromDate, 
                       const DateTime&   toDate, 
                       const YieldCurve* discount) {
    static const string method("LiborLeg::getPV");
    try {
        double pv = 0.0;
        CashFlowArrayConstSP libor = getCashFlowArray(valueDate, 
                                        couponCurve.isEmpty()? discount: couponCurve.get());
        for (int i = 0; i < libor->size(); i++) {
            DateTime paydate = (*libor)[i].date;
            if ((paydate>=fromDate && paydate<=toDate) && (paydate>valueDate)) {
                pv += (*libor)[i].amount * discount->pv(fromDate, paydate);
            }
        }
        return pv;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// compute accrued interest as of given date
double LiborLeg::getAccrued(const DateTime&   valueDate,
                            const DateTime&   when, 
                            const YieldCurve* yield) {
    static const string method("LiborLeg::getAccrued");
    try {
        PayStreamSP FloatPayStream(makePayStream(yield));                                          
        return FloatPayStream->getAccruedInterest(valueDate, when, 0);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// next coupon payment date
bool LiborLeg::getNextPayDate(const DateTime&   when, 
                              DateTime& nextPayDate){
    static const string method("LiborLeg::getNextPayDate");
    try {
        int k = Neighbour(when, PayDates, 0, PayDates.size()-1, 1);
        if (k<0){
            return false;               
        }      
        nextPayDate = PayDates[k];        
        return true;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// prepare KNOWN_CASH_FLOW
void LiborLeg::makeKnownCashFlow(){             
    static const string method = "LiborLeg::makeKnownCashFlow";
    try{
        // checke coupon curve is set from getMarket
        checkMarket();
        
        //initialization.
        knownCFL = CashFlowArraySP(new CashFlowArray(0));   
        
        // all the known cash flow
        if (RefixDates.size()>0) {
            CashFlowArrayConstSP  temp = getCashFlowArray(valueDate, couponCurve.get());
        
            for (int i=0; i<(*temp).size();i++){
                knownCFL->push_back((*temp)[i]);
            }
         }   
    }
    catch (exception& e)
    {            
        throw ModelException(e, method);
    }
};

CashFlowArraySP LiborLeg::getKnownCashFlows(){
    return knownCFL;
};

// tweak spreads for client valuation purpose
// isTweak = true, shift the spreads.  false recover the original values.
void LiborLeg::tweakSpreads(bool isTweak)
{
    int i;
    int nsize = Spreads.size();
    if (isTweak){
        for (i=0;i<nsize;i++){
            Spreads[i] = orgSpreads[i] + tweak_spread_size;  
        }
    }else{
        // store the original spreads
        for (i=0;i<nsize;i++){
            Spreads[i] = orgSpreads[i];
        }
    }
}

void LiborLeg::setFixingforThetaShift(const DateTime& valueDate, 
                                      const YieldCurve* yield,
                                      const DateTime& rollDate){
    DayCountConventionSP liborPayDCC(DayCountConventionFactory::make(PayDCC));
    
    FloatRateSP floater = this->getFloatRate(yield);

    DoubleArray fixings;

    for (int i=0; i<RefixDates.size(); i++)
    {
        if (rollDate.getDate() >= RefixDates[i].getDate() && valueDate.getDate() <= RefixDates[i].getDate())
            if (!(valueDate == RefixDates[i]) || Maths::isZero(Fixings[i]))
                Fixings[i] = floater->calculateRate(RefixDates[i]);
    }
}

/** populate couponCurve from product, if it doesn't have own Coupon Curve*/
void LiborLeg::setCouponCurve(const YieldCurveWrapper instCouponCurve)
{
    if (couponCurve.isEmpty()) {
        couponCurve.setObject(MarketObjectSP(copy(instCouponCurve.get())));
//        couponCurve->setProjectionCurve(true);  // use projection!
    }
}
                                                                    
/** populate from market cache */
void LiborLeg::getMarket(const IModel* model, const MarketData* market)
{
    market->GetReferenceDate(valueDate);
    if (!couponCurve.isEmpty()) {
        couponCurve.getData(model, market);
//        couponCurve->setProjectionCurve(true);  // use projection!
    }
}

// check market data is already set.
void LiborLeg::checkMarket(){
    static const string method = "LiborLeg::checkMarket";
    try  {
        if (couponCurve.isEmpty()) 
            throw ModelException(method,"couponCurve is not set.  Need to call setCouponCurve before.");
        if (valueDate.empty())
            throw ModelException(method,"valueDate is not initilized.  Need to call getMarket before.");
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Implementation of the Theta Shift interface */
bool LiborLeg::sensShift(Theta* shift) 
{
    static const string method = "LiborLeg::sensShift";
    try  {
        if (couponCurve.get()) {
            const DateTime& newDate = shift->rollDate(valueDate);
            setFixingforThetaShift(valueDate,
                                   couponCurve.get(),
                                   newDate);
            return true; // our components have theta type sensitivity
        }
        else {
            return false;
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }    
}

// -------------------------------------------------------------------
// LiborLeg State Var support
class LiborLeg::LiborLegSV::Imp {
public:
    Imp(LiborLegSP                    libor,
        PayStreamSP                   paystream,
        YieldCurveConstSP             discount,
        const DateTimeArray&          matDates,  // convenient
        vector<SVGenExpectedDiscFactorSP> expDfGens,
        SVGenDiscFactorSP                paymentDfGen,
        IStateVariableGen::IStateGen* pathGen):
        libor(libor),
        paystream(paystream),
        cfl(libor->RefixDates.size()),
        liborRateDCC(DayCountConventionFactory::make(libor->RateDCC)),
        fixings(libor->Fixings), // this fills in past fixings
        notionals(0),
        matDates(matDates),
        expDfGens(expDfGens),
        paymentDfGen(paymentDfGen) {

        if (!libor->Notionals.empty()) {
            notionals = &libor->Notionals;
        }

        // Prepare the SVs
        update(pathGen);
    }

    virtual ~Imp() {};

    virtual bool doingPast() const {
        return isDoingPast;
    }

    //// update state variables (e.g. going from past to future)
    void update(IStateVariableGen::IStateGen* pathGen){
        expDfSVs.clear();
        for(size_t i=0; i<expDfGens.size(); i++) {
            expDfSVs.push_back(expDfGens[i]->getSVExpectedDiscFactor(pathGen));
        }
        paymentDfSV  = paymentDfGen->getSVDiscFactor(pathGen);
        isDoingPast  = expDfSVs.size()>0?expDfSVs[0]->doingPast():true; // all may be past
    }

    // I'm not terribly convinced by the "valueDate" parameter here.
    const CashFlowArray* getCashFlowArray(const DateTime& valueDate) const {
        // Calculate fixings. 
        int offset = fixings.size() - expDfSVs.size();
        for(int i=0; i<fixings.size(); i++) {
            const DateTime& indexStartDate = libor->RefixDates[i];
            if (indexStartDate.isGreater(valueDate)) {
                double disc = expDfSVs[i-offset]->path()[0];
                fixings[i] = RateConversion::discountToRate(disc,
                                                            indexStartDate,
                                                            matDates[i],
                                                            liborRateDCC.get(),
                                                            CompoundBasis::SIMPLE);
            } 
        }
        // Turn them into cash flows
        paystream->cashFlowList(valueDate,
                                notionals,
                                fixings,
                                cfl);
        return &cfl;
    }

    double getPV(const DateTime& valueDate) const {
        const CashFlowArray* cfa = getCashFlowArray(valueDate);
        double result = 0.0;
        for (int j=0; j<cfa->size(); j++)
        {
            if ((*cfa)[j].date > valueDate)
    		result += (*cfa)[j].amount * paymentDfSV->path()[j]; // ? pv to "valueDate" but that can change???
        }
        return result;
    }

    // PV factor from refix date to today (valueDate?)
    double payDatePV(int idx) const {
        return paymentDfSV->path()[idx];
    }
    
    // fields
    LiborLegSP           libor;
    PayStreamSP          paystream;
    mutable CashFlowArray cfl;
    DayCountConventionSP liborRateDCC;
    mutable DoubleArray  fixings;
    const DoubleArray*   notionals;
    const DateTimeArray& matDates;
    bool                 isDoingPast;

    // SV stuff
    vector<SVGenExpectedDiscFactorSP>    expDfGens;  
    SVGenDiscFactorSP                    paymentDfGen;
    vector<SVExpectedDiscFactorSP> expDfSVs;
    SVDiscFactorSP                 paymentDfSV;
};

LiborLeg::LiborLegSV::LiborLegSV(LiborLegSP                    libor,
                                 PayStreamSP                   paystream,
                                 YieldCurveConstSP             discount,
                                 const DateTimeArray&          matDates,  // convenient
                                 vector<SVGenExpectedDiscFactorSP> expDfGens,
                                 SVGenDiscFactorSP                paymentDfGen,
                                 IStateVariableGen::IStateGen* pathGen) :
    me(new Imp(libor, paystream, discount, matDates, 
               expDfGens, paymentDfGen, pathGen)) {};

LiborLeg::LiborLegSV::~LiborLegSV() {}

bool LiborLeg::LiborLegSV::doingPast() const {
    return me->isDoingPast;
}

void LiborLeg::LiborLegSV::update(IStateVariableGen::IStateGen* pathGen) {
    me->update(pathGen);
}

PayStreamSP LiborLeg::LiborLegSV::getPayStream() const {
    return me->paystream;
}

const CashFlowArray* LiborLeg::LiborLegSV::getCashFlowArray(const DateTime& valueDate) const {
    return me->getCashFlowArray(valueDate);
}

double LiborLeg::LiborLegSV::getPV(const DateTime& valueDate) const {
    return me->getPV(valueDate);
}

double LiborLeg::LiborLegSV::payDatePV(int idx) const {
    return me->payDatePV(idx);
}

        // return libor leg 
smartPtr<LiborLeg> LiborLeg::LiborLegSV::getLibor() const{
    return me->libor;
}

class LiborLeg::LiborLegSVGen::Imp {
public:
    Imp(LiborLegSP              libor,
        YieldCurveConstSP       discount):
        libor(libor),
        discount(discount),
        floatRate(libor->getFloatRate(discount.get())),
        payStream(libor->makePayStream(discount.get())),
        refixDates(libor->RefixDates),
        matDates(libor->RefixDates.size()) {

        YieldCurveConstSP rateCurve = libor->couponCurve.isEmpty()? 
            discount: libor->couponCurve.getSP();
        expDfGens.clear();
        for(int i=0; i<refixDates.size(); i++) {
            matDates[i] = floatRate->findMatDate(refixDates[i]);
            if (refixDates[i] > libor->valueDate) {
                // only future ones
                expDfGens.push_back(SVGenExpectedDiscFactorSP(
                                        new SVGenExpectedDiscFactor(refixDates[i], // when to compute
                                                                 refixDates[i], // pv from here
                                                                 rateCurve,
                                                                 DateTimeArray(1, matDates[i]),
                                                                 false))); // compute log
            }
        }

        // Use this to get the cf dates
        CashFlowArrayConstSP cfa = payStream->cashFlowList(libor->valueDate, 0);
        DateTimeArray payDates = CashFlow::dates(*(cfa.get()));
        paymentDfGen = SVGenDiscFactorSP(new SVGenDiscFactor(libor->valueDate,
                                                       rateCurve,
                                                       payDates));
    }

    virtual ~Imp() {};

    LiborLeg::LiborLegSVSP getLiborLegSV(LiborLeg::LiborLegSVSP        oldStateVar,
                                         IStateVariableGen::IStateGen* pathGen) const {
        if (oldStateVar.get()) {
            // just update
            dynamic_cast<LiborLegSV&>(*oldStateVar).update(pathGen);
            return oldStateVar;
        }
        return LiborLegSVSP(new LiborLegSV(libor, payStream, discount, matDates, 
                                           expDfGens, paymentDfGen, 
                                           pathGen));
    }

    void collectStateVars(IStateVariableCollectorSP svCollector) const {
        for(size_t i=0; i<expDfGens.size(); i++) {
            svCollector->append(expDfGens[i].get());
        }
        svCollector->append(paymentDfGen.get());
    }

    // fields
    LiborLegSP        libor;
    YieldCurveConstSP discount;
    FloatRateSP       floatRate;
    PayStreamSP       payStream;
    DateTimeArray     refixDates; 
    DateTimeArray     matDates;
    vector<SVGenExpectedDiscFactorSP>    expDfGens; // StateVar generator for expected disc factors from refix dates to rate mat dates
    SVGenDiscFactorSP    paymentDfGen; // StateVar generator for disc factors from pay dates
};

LiborLeg::LiborLegSVGen::LiborLegSVGen(LiborLegSP        libor,
                                       YieldCurveConstSP discount):
    me(new Imp(libor, discount)){};

LiborLeg::LiborLegSVGen::~LiborLegSVGen() {}

/** Create the corresponding State Variable for this State
    Variable Generator (NB Implies one state variable per
    generator). The previous IStateVariableSP (may be null)
    should be passed in.  The return object may or may not be
    the same as oldStateVar.  */
IStateVariableSP LiborLeg::LiborLegSVGen::create(IStateVariableSP             oldStateVar,
                                                 IStateVariableGen::IStateGen* pathGen) const{
    // need to handle the case safely when oldStateVar.get() is 0
    if (oldStateVar.get()) {
        // just update
        dynamic_cast<LiborLegSV&>(*oldStateVar).update(pathGen);
        return oldStateVar;
    }
    return me->getLiborLegSV(LiborLegSVSP(   ), pathGen);
}

LiborLeg::LiborLegSVSP LiborLeg::LiborLegSVGen::getLiborLegSV(LiborLeg::LiborLegSVSP        oldStateVar,
                                                              IStateVariableGen::IStateGen* pathGen) const {
    return me->getLiborLegSV(oldStateVar, pathGen);
}

/** Appends 'true' (ie non derived) state variable generators
    required to the supplied collector. Implementations typically call
    IStateVariableCollector::append */
void LiborLeg::LiborLegSVGen::collectStateVars(
    IStateVariableCollectorSP svCollector) const {
    me->collectStateVars(svCollector);
}

LiborLeg::LiborLegSVGen* LiborLeg::createLiborLegSVGen(YieldCurveConstSP discount) const {
    return new LiborLegSVGen(LiborLegSP(copy(this)), discount);
}

// -------------------------------------------------------------------


//------------------------------------------------------------------//
//  FixedLeg : The interface for Fixed Leg  (deterministic rate)
//------------------------------------------------------------------//
void FixedLeg::validatePop2Object(){

    static const string method = "Fixed::validatePop2Object";

    try
    { 
        int numItems = PaymentDatesArray.size();
        if (numItems != AccrueStartDates.size() ||
            numItems != AccrueEndDates.size() ||
            numItems != CouponAmounts.size())
        {
            throw ModelException(method,
                                 "The number of pay dates (size = " + Format::toString(PaymentDatesArray.size()) +
                                 "), coupon amounts (size =" + Format::toString(CouponAmounts.size()) + 
                                 "), accrue start dates (size =" + Format::toString(AccrueStartDates.size()) +
                                 "), accrue end dates (size =" + Format::toString(AccrueEndDates.size()) + 
                                 ") must be equal.  \n");
        }
        int i; // MSVC broken        
        for (i=1; i<numItems; i++)
        {// checke cronical date array or not
            if (PaymentDatesArray[i-1]>PaymentDatesArray[i])
                throw ModelException(method, "PaymentDatesArray (" + Format::toString(i) + ") [ " + PaymentDatesArray[i-1].toString() + 
                                    " ] should be ealier than next date [ " + PaymentDatesArray[i].toString() + " ] \n");
            if (AccrueStartDates[i-1]>AccrueStartDates[i])
                throw ModelException(method, "AccrueStartDates (" + Format::toString(i) + ") [ " + AccrueStartDates[i-1].toString() + 
                                    " ] should be ealier than next date [ " + AccrueStartDates[i].toString() + " ] \n");
            if (AccrueEndDates[i-1]>AccrueEndDates[i])
                throw ModelException(method, "AccrueEndDates (" + Format::toString(i) + ") [ " + AccrueEndDates[i-1].toString() + 
                                    " ] should be ealier than next date [ " + AccrueEndDates[i].toString() + " ] \n");
        }
        for (i=0; i<numItems; i++)
        {
            if (AccrueStartDates[i]>AccrueEndDates[i])
                throw ModelException(method, "AccrueStartDates (" + Format::toString(i) + ") [ " + AccrueStartDates[i].toString() + 
                                    " ] should be ealier than AccrueEndDates [ " + AccrueEndDates[i].toString() + " ] \n");
        }
        for (i=0; i<numItems-1; i++)
        {
            if (AccrueEndDates[i]>AccrueStartDates[i+1])
                throw ModelException(method, "next AccrueStartDates (" + Format::toString(i+1) + ") [ " + AccrueStartDates[i+1].toString() + 
                                    " ] should be later than AccrueEndDates [ " + AccrueEndDates[i].toString() + " ] \n");
        }

    }
    catch (exception& e) 
    {
        throw ModelException(e, method);
    }

}

// helpers
class FixedLegHelper {
public:
    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(FixedLeg, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultFixedLeg);
                
        FIELD(AccrueStartDates, "FixedCoupon start date array");
        FIELD(AccrueEndDates, "FixedCoupon end date array");
        FIELD(PaymentDatesArray, "FixedCoupon payment dates array");
        FIELD(CouponAmounts, "fixed coupon amounts");
        FIELD(dcc, "fixed couopn day count");
        FIELD_MAKE_OPTIONAL(dcc);

    }

static IObject* defaultFixedLeg(){
        return new FixedLeg();
    }
};

CClassConstSP const FixedLeg::TYPE = CClass::registerClassLoadMethod(
    "FixedLeg", typeid(FixedLeg), FixedLegHelper::load);

// constructor
FixedLeg::FixedLeg(): CObject(TYPE){
    dcc = "NONE";
    hasKOSettle = false;
}

// get the last coupon payment date
DateTime FixedLeg::getLastPayDate()
{
    return PaymentDatesArray.back();
}

int  FixedLeg::getSize()
{
    return PaymentDatesArray.size();
}

bool FixedLeg::isDCC()
{
	return dcc!="NONE";
}

CashFlowArrayConstSP FixedLeg::getCashFlowArray(){
    //--Creating Fixed Cash Flow--//
    CashFlowArraySP fixedflow (new CashFlowArray());
    int num = CouponAmounts.size();
    double dccfact = 1.0;
    DayCountConventionSP daycount; 
    if (isDCC()) daycount = DayCountConventionSP(DayCountConventionFactory::make(dcc));    
    for (int i=0; i<num; i++){
        if (isDCC ()) dccfact = daycount->years(AccrueStartDates[i], AccrueEndDates[i]);
        fixedflow->push_back(CashFlow(PaymentDatesArray[i], CouponAmounts[i] * dccfact) );
    }
    return fixedflow;
}

//retern Fixed Cash Flow payment after baseDate //
CashFlowArraySP FixedLeg::getCashFlowArray(const DateTime& baseDate){
    CashFlowArraySP fixedflow(new CashFlowArray()) ;
    int num = CouponAmounts.size();
    double dccfact = 1.0;
    DayCountConventionSP daycount; 
    if (isDCC()) daycount = DayCountConventionSP(DayCountConventionFactory::make(dcc));
    for (int i=0; i<num; i++){
        if (baseDate<PaymentDatesArray[i]){
            if (isDCC()) dccfact = daycount->years(AccrueStartDates[i], AccrueEndDates[i]);
            fixedflow->push_back(CashFlow(PaymentDatesArray[i], CouponAmounts[i]*dccfact));
        }

    }
    return fixedflow;
}

// get PV from ValueDate to End, as of valueDate.
double FixedLeg::getPV(const DateTime& valueDate, const YieldCurve* discount){
	// add uncounted unconditional fix coupon
	double result = 0.0;
    CashFlowArrayConstSP cfl = getCashFlowArray();
    for (int j=0; j<cfl->size(); j++)
    {
        if ((*cfl)[j].date > valueDate){
            double amount = (*cfl)[j].amount;
            DateTime dddd = (*cfl)[j].date;
            double pvf = discount->pv(valueDate, (*cfl)[j].date);
    		result += (*cfl)[j].amount * discount->pv(valueDate, (*cfl)[j].date);
        }
    }
    return result;
}

// get pv of cashflows paying between fromDate & toDate as of fromDate
double FixedLeg::getPV(const DateTime&   valueDate,
                       const DateTime&   fromDate, 
                       const DateTime&   toDate, 
                       const YieldCurve* discount) {
    static const string method("FixedLeg::getPV");
    try {
        double pv = 0.0;
        CashFlowArrayConstSP cfl = getCashFlowArray();
        for (int i = 0; i < cfl->size(); i++) {
            DateTime paydate = (*cfl)[i].date;
            if ((paydate>=fromDate && paydate<=toDate) && (paydate>valueDate)) {
                pv += (*cfl)[i].amount * discount->pv(fromDate, paydate);
            }
        }
        return pv;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// compute accrued interest as of given date
double FixedLeg::getAccrued(const DateTime&   valueDate,
                            const DateTime&   when, 
                            const YieldCurve* discount) {
    static const string method("FixedLeg::getAccrued");
    try {                    
        double accrued = getAccrued(when)
                        * discount->pv(valueDate,when);         

        return accrued;    

    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// compute accrued coupon amount (no Discount)
double FixedLeg::getAccrued(const DateTime&   when) {
    static const string method("FixedLeg::getAccrued");
    try {
        // find the fee date immediately on or after this date
        int k = Neighbour(when, AccrueStartDates, 0, AccrueStartDates.size()-1, -1);
        if (k<0)
            return 0;
        else if (when >= PaymentDatesArray[k])
            return 0;

        if (isDCC()){
            DateTime lo = AccrueStartDates[k];
            DayCountConventionSP daycount = DayCountConventionSP(DayCountConventionFactory::make(dcc));
            return CouponAmounts[k] * daycount->years(lo, when);        
        }
        else
            return CouponAmounts[k];
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


// next coupon payment date
bool FixedLeg::getNextPayDate(const DateTime&   when, 
                              DateTime& nextPayDate){
    static const string method("FixedLeg::getNextPayDate");
    try {
        // find the payment date immediately on or after this date
        int k = Neighbour(when, PaymentDatesArray, 0, PaymentDatesArray.size()-1, 1);
        if (k<0){
            return false;
        }        
        nextPayDate = PaymentDatesArray[k];        
        return true;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


/***************************************************************/
// 
// support barrier pay
//
/***************************************************************/

class LiborLegBarrierPay : virtual public IBarrierPay
{
public:
    LiborLegBarrierPay(LiborLeg* leg, bool isOut, const string& koStubRule) 
        : leg(LiborLegSP(copy(leg))), isOut(isOut), koStubRule(koStubRule)
    {
        if( leg->Compounding.size() > 0 )
            throw ModelException("LiborLegBarrierPay", "compounding not supported yet");
    };
    ~LiborLegBarrierPay(){};

    void preprocess(DateTime valueDt, 
                    DateTimeArrayConstSP simDts, 
                    const DateTimeArray& monitorDates, 
                    const YieldCurve* discount);

    bool payOnce() const 
    { return false; }

    bool payIfHitOne() const 
    { return !isOut; }

    void payObsDates(DateTimeArray &obsDts, bool &hasContinuous) const
    {
        obsDts.clear();
        hasContinuous = false;
    }

    const DateTimeArray& relevantMonDates() const
    {
        return relevantMonDts;
    }

    double value(int step, DateTime valueDate, const MCPathGenerator* pathGen)
    {         
        return payPerStep[step]; 
    };
    
    void collect(int step, const MCPathGenerator* pathGenIn, double scalingFactor,
                 CashFlowArray& knownCFs, PhysicalDeliveryArray& phyDs)
    {
        int i=0;
        while( i<relevantMonDts.size() && (*simDates)[step]!=relevantMonDts[i])
            i++;
        if( i==relevantMonDts.size() )
            throw ModelException("LiborLegBarrierPay:collect", "Sim date " + (*simDates)[step].toString() 
                                + " not a monitor date");

        CashFlowArray &cf = cfs[i];
        for(i=0; i<cf.size(); i++)
            knownCFs.push_back(CashFlow(cf[i].date, scalingFactor * cf[i].amount));
    };
 
private:
    LiborLegSP      leg;
    bool            isOut;
    string          koStubRule;

    // transient member
    CashFlowCluster cfs;        // fv factor to maturity for each full/partial coupon pay
    DoubleArray     payPerStep;
    DateTimeArray   relevantMonDts;
    DateTimeArrayConstSP simDates;
};


#if defined(_MSC_VER)
#pragma optimize( "g", off)
#endif
void LiborLegBarrierPay::preprocess(DateTime valueDt, 
                                    DateTimeArrayConstSP simDts, 
                                    const DateTimeArray& monitorDates, 
                                    const YieldCurve* discount)
{
    static const string method("LiborLegBarrierPay::preprocess");

    simDates = simDts;

    PayStreamSP FloatPayStream(leg->makePayStream(discount));                                            
    CashFlowArrayConstSP libor = FloatPayStream->cashFlowList(valueDt, 0);        
    DoubleArray pvs;
    relevantMonDts.clear();
    cfs.clear();        

    // find first monitor date strictly after accrual start date
    int i=0, j=0;
    while( j<monitorDates.size() && monitorDates[j] <= leg->AccrualDates.front() )
        j++;
        
    // if all monitor before 1st accural start, insert last mon date to link to all cpns,
    // otherwise if swap/bond stub, insert last mon date before 1st accrue start to link to portion of 1st cpn 
    //      from 1st accrue start to the earlier of next monitor or 1st accrue end,
    // otherwise if none stub, insert last mon date before 1st accrue start to link to 1st cpn 
    //      if next mon date is after end of 1st accrue end

	//CHANGE
	//	if( j==monitorDates.size() || koStubRule == "S" || monitorDates[j] > leg->AccrualDates[1] ) //OLD
    if( j==monitorDates.size() || koStubRule == "S" || koStubRule == "B" || monitorDates[j] > leg->AccrualDates[1] )//NEW
    {
        // all(swap) or at least portion(none/bond) of 1st coupon period must have a past mon date linked to it
        if( j==0 )
            throw ModelException(method, "All or portion of libor payment on " + 
                                 libor->front().date.toString() + " has no monitor date associated with it.\n"
                                 "You cannot add no-contingent coupon, which has no barrier event.");
            
        relevantMonDts.push_back(monitorDates[j-1]);
        pvs.push_back(0.0);
        cfs.push_back(CashFlowArray(0));
    }
        
    DayCountConventionSP daycount(DayCountConventionFactory::make(leg->PayDCC));
    for(i=0; i<libor->size(); i++)
    {
        DateTime payDate = (*libor)[i].date;
        DateTime prevDate = leg->AccrualDates[i];
        double payLet, accrual0 = daycount->years(prevDate, leg->AccrualDates[i+1]);
        if( Maths::isZero(accrual0) ) continue;

          // if case "B", give accrual0 to pvs and cfs.back
		if( koStubRule == "B")	//NEW
		{						//NEW
            payLet = (*libor)[i].amount;	//NEW

                if( !Maths::isZero(payLet) ) {		//NEW	
                    if( payDate > valueDt)			//NEW	
                        pvs.back() += payLet  * discount->pv(valueDt, payDate);	//NEW
                    cfs.back().push_back(CashFlow(payDate, payLet));	//NEW
                }									//NEW
                prevDate = leg->AccrualDates[i+1];			//NEW
            }					//NEW
		
        while( j<monitorDates.size() && monitorDates[j] <= leg->AccrualDates[i+1] )
        {
      		if( koStubRule == "S" )
            {
                payLet = (*libor)[i].amount * daycount->years(prevDate, monitorDates[j])/accrual0;
                if( !Maths::isZero(payLet) ) {
                    if( payDate > valueDt)
                        pvs.back() += payLet  * discount->pv(valueDt, payDate);
                    cfs.back().push_back(CashFlow(payDate, payLet));
                }
                prevDate = monitorDates[j];
            }

			// if stub is N, the relevant monitor date is last one before/inclusive accrual (end) date
			// if case B, same as case N. ie. insert last mon date of this accural period into the
			// pv/cfs array to prepare for next accural period
            // if a coupon period contain no monitor date, it is conditional on the last past mon date

            if( koStubRule == "S" ||
                ( koStubRule == "N" && (j==(monitorDates.size()-1) || monitorDates[j+1] > leg->AccrualDates[i+1]) ) ||
				( koStubRule == "B" && (j==(monitorDates.size()-1) || monitorDates[j+1] > leg->AccrualDates[i+1]) ) ) //NEW
//                ( koStubRule == "B" && (j==0 || monitorDates[j-1] <= leg->AccrualDates[i]) ) ) //OLD
            {
                relevantMonDts.push_back(monitorDates[j]);
                pvs.push_back(0.0);
                cfs.push_back(CashFlowArray(0));
            }

			j++;
        }
            
        payLet = (*libor)[i].amount * daycount->years(prevDate, leg->AccrualDates[i+1])/accrual0;
        if( !Maths::isZero(payLet) ) {
            if( payDate > valueDt)
                pvs.back() += payLet * discount->pv(valueDt, payDate);
            cfs.back().push_back(CashFlow(payDate, payLet));
        }
    }
        
    payPerStep.resize(simDts->size());
    for(i=0, j=0; i<simDts->size(); i++)
    {
        payPerStep[i] = (j<relevantMonDts.size() && (*simDts)[i]==relevantMonDts[j])?pvs[j++]:0.0;
    }

};
#if defined(_MSC_VER)
#pragma optimize( "g", on)
#endif

IBarrierPaySP LiborLeg::createBarrierPay(bool isOut, const string& stubRule)
{
    return IBarrierPaySP(new LiborLegBarrierPay(this, isOut, stubRule));
}

class FixedLegBarrierPay : virtual public IBarrierPay
{
public:
    FixedLegBarrierPay(FixedLeg* leg, bool isOut, const string& koStubRule) 
        : leg(FixedLegSP(copy(leg))), isOut(isOut), koStubRule(koStubRule) 
    {
        if( koStubRule == "S" && leg->dcc == "NONE" )
            throw ModelException("FixedLegBarrierPay", "koStubRule S not allowed if DCC is NONE");
    };
    ~FixedLegBarrierPay(){};

    void preprocess(DateTime valueDt, DateTimeArrayConstSP simDts, const DateTimeArray& monitorDates, const YieldCurve* discount)
    {
        static const string method("FixedLegBarrierPay::preprocess");

        simDates = simDts;
        
        DoubleArray pvs;
        relevantMonDts.clear();
        cfs.clear();        
        
        // find first monitor date strictly after accrual start date
        int i=0, j=0;
        while( j<monitorDates.size() && monitorDates[j] <= leg->AccrueStartDates.front() )
            j++;
        
        // if all monitor before 1st accural start, insert last mon date to link to all cpns,
        // otherwise if swap/bond stub, insert last mon date before 1st accrue start to link to portion of 1st cpn 
        //      from 1st accrue start to the earlier of next monitor or 1st accrue end,
        // otherwise if none stub, insert last mon date before 1st accrue start to link to 1st cpn 
        //      if next mon date is after end of 1st accrue end
//        if( j==monitorDates.size() || koStubRule == "S" || monitorDates[j] > leg->AccrueEndDates.front() )

        if( j==monitorDates.size() || koStubRule == "S" ||koStubRule == "B" || monitorDates[j] > leg->AccrueEndDates.front() )
        {
            // all(swap) or at least portion(none/bond) of 1st coupon period must have a past mon date linked to it
            if( j==0 )
                throw ModelException(method, "All or portion of fix payment on " + 
                      leg->PaymentDatesArray.front().toString() + " has no monitor date associated with it .\n"
                      "You cannot add no-contingent coupon, which has no barrier event.");
            
            relevantMonDts.push_back(monitorDates[j-1]);
            pvs.push_back(0.0);
            cfs.push_back(CashFlowArray(0));
        }
        
        DayCountConventionSP daycount;
        if(leg->isDCC()) daycount = DayCountConventionSP(DayCountConventionFactory::make(leg->dcc));
        for(i=0; i<leg->CouponAmounts.size(); i++)
        {
            DateTime payDate = leg->PaymentDatesArray[i];
            DateTime prevDate = leg->AccrueStartDates[i];
            double payLet; 
			
            if( koStubRule == "B" )										//NEW
            {														//NEW
                payLet = leg->CouponAmounts[i];
                if (!!daycount)
                    payLet *= daycount->years(prevDate, leg->AccrueEndDates[i]); 	//NEW
                if( !Maths::isZero(payLet) ) {						//NEW
                    if( payDate > valueDt)							//NEW
                        pvs.back() += payLet  * discount->pv(valueDt, payDate);	//NEW
                    cfs.back().push_back(CashFlow(payDate, payLet));			//NEW
                }																//NEW
                prevDate = leg->AccrueEndDates[i];										//NEW
            }																	//NEW

            while( j<monitorDates.size() && monitorDates[j] <= leg->AccrueEndDates[i] )
            {
                if( koStubRule == "S" )
                {
                    payLet = leg->CouponAmounts[i] * daycount->years(prevDate, monitorDates[j]); // dayCount != 0
                    if( !Maths::isZero(payLet) ) {
                        if( payDate > valueDt)
                            pvs.back() += payLet  * discount->pv(valueDt, payDate);
                        cfs.back().push_back(CashFlow(payDate, payLet));
                    }
                    prevDate = monitorDates[j];
                }
                
                // if stub is N, the relevant monitor date is last one before/inclusive accrual (end) date
				// if case B, same as case N. ie. insert last mon date of this accural period into the
				// pv/cfs array to prepare for next accural period
				// if a coupon period contain no monitor date, it is conditional on the last past mon date
                
				if( koStubRule == "S" ||
                  ( koStubRule == "N" && (j==(monitorDates.size()-1) || monitorDates[j+1] > leg->AccrueEndDates[i]) ) ||
				  ( koStubRule == "B" && (j==(monitorDates.size()-1) || monitorDates[j+1] > leg->AccrueEndDates[i]) ))//NEW
                //  ( koStubRule == "B" && (j==0 || monitorDates[j-1] <= leg->AccrueStartDates[i]) ) )//OLD
                {
                    relevantMonDts.push_back(monitorDates[j]);
                    pvs.push_back(0.0);
                    cfs.push_back(CashFlowArray(0));
                }
                j++;
            }
            if (koStubRule == "B")
                payLet = 0.0;   // already added
            else
                payLet = leg->CouponAmounts[i] * (!daycount?1.0:daycount->years(prevDate, leg->AccrueEndDates[i]));
            if( !Maths::isZero(payLet) ) {
                if( payDate > valueDt)
                    pvs.back() += payLet * discount->pv(valueDt, payDate);
                cfs.back().push_back(CashFlow(payDate, payLet));
            }
        }
        
        payPerStep.resize(simDts->size());
        for(i=0, j=0; i<simDts->size(); i++)
        {
            payPerStep[i] = (j<relevantMonDts.size() && (*simDts)[i]==relevantMonDts[j])?pvs[j++]:0.0;
        }

    };

    bool payOnce() const 
    { return false; }

    bool payIfHitOne() const 
    { return !isOut; }

    void payObsDates(DateTimeArray &obsDts, bool &hasContinuous) const
    {
        obsDts.clear();
        hasContinuous = false;
    }
    
    const DateTimeArray& relevantMonDates() const
    {
        return relevantMonDts;
    }

    double value(int step, DateTime valueDate, const MCPathGenerator* pathGen)
    {         
        return payPerStep[step]; 
    };
    
    void collect(int step, const MCPathGenerator* pathGenIn, double scalingFactor,
                 CashFlowArray& knownCFs, PhysicalDeliveryArray& phyDs)
    {
        int i=0;
        while( i<relevantMonDts.size() && (*simDates)[step]!=relevantMonDts[i])
            i++;
        if( i==relevantMonDts.size() )
            throw ModelException("FixedLegBarrierPay:collect", "Sim date " + (*simDates)[step].toString() 
                                + " not a monitor date");

        CashFlowArray &cf = cfs[i];
        for(i=0; i<cf.size(); i++)
            knownCFs.push_back(CashFlow(cf[i].date, scalingFactor * cf[i].amount));
    };
 
private:
    FixedLegSP      leg;
    bool            isOut;
    string          koStubRule;

    // transient member
    CashFlowCluster cfs;        // fv factor to maturity for each full/partial coupon pay
    DoubleArray     payPerStep;
    DateTimeArray   relevantMonDts;
    DateTimeArrayConstSP simDates;
};

IBarrierPaySP FixedLeg::createBarrierPay(bool isOut, const string& stubRule)
{
    return IBarrierPaySP(new FixedLegBarrierPay(this, isOut, stubRule));
}

DRLIB_END_NAMESPACE
