//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FRBondVariable.cpp
//
//   Description : Read-only variable which provides bond-like values
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/FRController.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/Maths.hpp"
#include "edginc/FRFactory.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/SVGenExpectedDiscFactor.hpp"
#include "edginc/SVGenIRSwap.hpp"
#include "edginc/GetMarket.hpp"
#include "edginc/CashFlow.hpp"

DRLIB_BEGIN_NAMESPACE

static void checkAgainstStateVars(const string& method, FRController* frCtrl){
    if (frCtrl->stateVarUsed()){
        throw ModelException(method, "State Vars not supported yet");
    }
}

class FRBondVariable: public FR::LValueBase {
public:
    static CClassConstSP const TYPE;

    /** Returns the type that this variable represents */
    FRIfaces::VarType getType() const{
        return FRIfaces::doubleType;
    }

    /** returns a string that can be used to identify this object */
    virtual const string& getID() const{
        return id;
    }

    /** returns the input to define the variable in IMS */
    virtual StringArray getIMSInput() const {
        StringArray result(3);
        
        result[0] = "Bond type is not implemented in IMS";
        result[1] = id;
        result[2] = "Bond type is not implemented in IMS";

        return result;
    }

    virtual void validatePop2Object(){
        static const string method = "FRBondVariable::validatePop2Object";

        // initialise our id
        const_cast<string&>(id) = bondName; // validatePop2Object is privileged

        // convert the day count convention string to an object
        accruedDCC = DayCountConventionSP(
            DayCountConventionFactory::make(dayCountConvString));

        if (cashFlows->size()<1) {
            throw ModelException(method, "Require at least one cash flow!");
        }
        // validate that cashFlow dates are increasing 
        int i;
        for (i = 1; i < cashFlows->size(); i++){
            if ((*cashFlows)[i-1].date.isGreater((*cashFlows)[i].date)){
                throw ModelException(method, 
                                     "CashFlow dates are not increasing : [" + 
                                     Format::toString(i-1) +
                                     "] " + (*cashFlows)[i-1].date.toString() + 
                                     " > [" + Format::toString(i) + "] " +
                                     (*cashFlows)[i].date.toString());
            }
        }
        // validate that datedDate is before first cashFlowDate
        if (datedDate.isGreater((*cashFlows)[0].date)) {
            throw ModelException(method, "datedDate " + datedDate.toString() + 
                                 " must not be after first cashFlowDate " + 
                                 (*cashFlows)[0].date.toString());
        }

        switch (stubTypeString[0])
        {
        case 'n':
        case 'N':
            stubType = stubNone;
            break;
        case 's':
        case 'S':
            stubType = stubSwap;
            break;
        case 'b':
        case 'B':
            stubType = stubBond;
            break;
        default:
            throw ModelException(method, "Unrecognised stub type " + stubTypeString + 
                                 "\nExpect one of None, Swap, Bond");
            break;
        }

    }

    virtual FRIfaces::ILValue* getLValue(FRController* frCtrl) const{
        throw ModelException("FRBondVariable::getLValue", "Variable "+
                             id+" is read-only");
    }

private:
    typedef enum StubType{
        stubNone = 0,
        stubSwap,
        stubBond
    } StubType;

protected:
    /** creates FlexRule::IRValue which represents the value
        of this expression at a specific timepoint */
    virtual FRIfaces::IRValue* createRValue(int           index,
                                            FRController* frCtrl) const{
        checkAgainstStateVars("FRBondVariable:createRValue", frCtrl);
        // identify our ccyName
        const FRIfaces::IProductView* productView = frCtrl->productView();
        YieldCurveConstSP discount = productView->findYieldCurve(&ccyName);
            
        /* At each date will provide the PV of future complete coupons
           Plus possibly stubbed first coupon.  By working back from
           maturity to today we can do this in a single pass.  This
           should work for any sim date : before datedDate and after
           last cashFlow date */
        int        iCashFlow = cashFlows->size() - 1;
        double     price;
        double     dirtyPrice = 0;
        DateTime   dirtyDate = (*cashFlows)[iCashFlow].date;//also accrueEndDate
        DateTime   maturityDate = dirtyDate; // final cash flow date
        bool       isRedemptionIncluded = false;
        int        iPastDate = pastBondValues.get()?
            pastBondValues->size()-1 : -1;

        bool doSetFiller = !fillerSet;
        int numDates = productView->getSimDates()->size();
        const DateTime& today = discount->getToday();
        DateTimeArrayConstSP simDates(productView->getSimDates());
        for (int iSimDate = numDates-1; iSimDate >=0 ; iSimDate--)
        {
            const DateTime& valueDate = (*simDates)[iSimDate];

            // Careful with past - and rolling both "to today" and for theta 
            if (valueDate < today) {
                // If past rely on user to supply the bond price - locate
                // based on date to be safe
                while (iPastDate>=0 &&
                       (*pastBondValues)[iPastDate].date > valueDate) {
                    iPastDate--;
                }
                if (iPastDate>=0 &&
                    (*pastBondValues)[iPastDate].date == valueDate &&
                    Maths::isPositive((*pastBondValues)[iPastDate].amount)) {
                    price = (*pastBondValues)[iPastDate].amount;
                } else if (fillerSet) {
                    // For any rolled-past dates use the earliest
                    // price after orig valuation date
                    price = fillerPrice; 
                } else {
                    throw ModelException("FRBondVariable:createRValue", 
                                         "Past Bond Price missing for '" + 
                                         id+"' on "+valueDate.toString());
                }
            } else {
                // Redemption amount included (once only) on maturity
                // (final cash flow date)
                if (!isRedemptionIncluded &&
                    valueDate <= maturityDate) {
                    dirtyPrice = redemptionAmt;
                    isRedemptionIncluded = true;
                }

                // Capture cash flows future of this valueDate
                for(; iCashFlow>=0 &&
                        (*cashFlows)[iCashFlow].date >= valueDate;
                    iCashFlow--)
                {
                    dirtyPrice *= discount->pv((*cashFlows)[iCashFlow].date,
                                               dirtyDate);
                    dirtyDate = (*cashFlows)[iCashFlow].date;
                    dirtyPrice += (*cashFlows)[iCashFlow].amount;
                }
                // Treat any stub - iCashFlow now references "next
                // earlier" cash flow
                DateTime accrueStartDate = iCashFlow<0?
                    datedDate : (*cashFlows)[iCashFlow].date;
                double   accrued = 0;
                if (stubType!=stubNone &&
                    accrueStartDate < valueDate && 
                    dirtyDate >= valueDate) /* accrueEndDate >= valueDate means
                                               it's all past so ignore */
                { 
                    // By validating at least one cash flow we are ok
                    // referencing [iCashFlow+1]
                    accrued = (*cashFlows)[iCashFlow+1].amount *
                        accruedDCC->years(accrueStartDate, valueDate) /
                        accruedDCC->years(accrueStartDate, dirtyDate);

                    if (stubType==stubSwap) {
                        // swap convention pays accrued at next coupon date
                        price = dirtyPrice - accrued; 
                    } else {
                        price = dirtyPrice; // accrued handled later 
                    }
                } else {
                    price = dirtyPrice;
                }

                // Price at valueDate ... 
                price *= discount->pv(valueDate, dirtyDate);
            
                // Bond convention pays accrued at valueDate
                if (stubType==stubBond) {
                    price = dirtyPrice - accrued;
                }

                if (doSetFiller) {
                    /* Sorry about this.  The bond does not actually
                       have a "spot" without a yield curve, so this
                       allows us to pick the value at the earliest
                       date after original-today (since we know
                       pricing happens first) and use that at any
                       dates that are skipped (we only read this value
                       if we fail to find a supplied past closing -
                       see above). */
                    fillerPrice = price;
                    fillerSet = true;
                }
            }
            // avoid putting values in more than once - this can happen if
            // this routine is called again after it has failed (eg missing
            // sample). A calling method might ignore the error in the [vain]
            // hope that the particular expression will not be evaluated at
            // run time (eg consider 'a = i>3? value[i-3]: 0.0')
            if (!frCtrl->getRValue(this, iSimDate)){
                FRIfaces::IRValueSP rValue(new FR::RConstDouble(price));
                frCtrl->setRValue(iSimDate, this, rValue);
            }
        }
        return 0; // indicates that we've saved the value ourselves
    }

private:
    const string               id;      // transient - set upon construction
    const string               bondName;
    const string               ccyName;
    const string               stubTypeString;
    const DateTime             datedDate; // accrue start date for first coupon
    CashFlowArrayConstSP cashFlows;  // coupon amounts
    const double               redemptionAmt;
    const string               dayCountConvString;
    CashFlowArrayConstSP pastBondValues;

    // derived
    DayCountConventionSP accruedDCC; // $unregistered
    StubType             stubType; // $unregistered
    mutable double       fillerPrice; // $unregistered
    mutable bool         fillerSet; // $unregistered

    FRBondVariable(): LValueBase(TYPE), redemptionAmt(0), fillerPrice(0.0), fillerSet(false) {}

    static IObject* defaultFRBondVariable(){
        return new FRBondVariable();
    }
        
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FRBondVariable, clazz);
        SUPERCLASS(FR::LValueBase);
        EMPTY_SHELL_METHOD(defaultFRBondVariable);
        FIELD(id, "");
        FIELD_MAKE_TRANSIENT(id);
        FIELD(bondName, "bond name");
        FIELD(ccyName, "currency Name");
        FIELD(stubTypeString, "N S or B");
        FIELD(datedDate, "accrue start date for first cash flow");
        FIELD(cashFlows, "contiguous accrual periods");
        FIELD(redemptionAmt, "included at mat if no accruing");
        FIELD(dayCountConvString, "dayCountConvString");
        FIELD(pastBondValues, "pastBondValues");
        FIELD_MAKE_OPTIONAL(pastBondValues);
        clazz->setPublic(); // make visible to EAS/spreadsheet
        Addin::registerConstructor("FR_VAR_BOND",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a FR Bond Variable",
                                   TYPE);
    }
    
};

CClassConstSP const FRBondVariable::TYPE = CClass::registerClassLoadMethod(
    "FRBondVariable", typeid(FRBondVariable), load);

    
///////////////////////////////////////////////////////
// A more tractable version of the above. Also supports SV.
///////////////////////////////////////////////////////
class FRZeroBondVariable: public FR::LValueBase,
                          virtual public IGetMarket {
public:
    static CClassConstSP const TYPE;

    /** See comment in FRMMRateVariable */
    IObject* clone() const{
        return CObject::clone();
    }

    /** Returns the type that this variable represents */
    FRIfaces::VarType getType() const{
        return FRIfaces::doubleType;
    }

    /** returns a string that can be used to identify this object */
    virtual const string& getID() const{
        return id;
    }

    /** returns the input to define the variable in IMS */
    virtual StringArray getIMSInput() const {
        StringArray result(3);

        result[0] = "ZeroBond";
        result[1] = id;
        string part;
        if (!ccyName.getName().empty()) {
            part = ccyName.getName() + ",";
        }
        result[2] = part + maturityDate.toString();

        return result;
    }

    virtual void validatePop2Object(){
        static const string method = "FRZeroBondVariable::validatePop2Object";

        // initialise our id
        const_cast<string&>(id) = bondName; // validatePop2Object is privileged
    }

    virtual FRIfaces::ILValue* getLValue(FRController* frCtrl) const{
        throw ModelException("FRZeroBondVariable::getLValue", "Variable "+
                             id+" is read-only");
    }

    // IGetMarket
    virtual void getMarket(const IModel* model, const MarketData* market) {
        if (!ccyName.getName().empty()) {
            ccyName.getData(model, market);
        }
    }

    /** Class for holding expected discount factors along each point*/
    class ZeroBondSV: public FR::RValueDouble, 
                      public virtual FRIfaces::ILValueDouble,
                      public virtual FRController::IStateVarClient{
        struct MyRT{
            TGetValue*                func;
            TSetValue*                setFunc;
            SVExpectedDiscFactorSP discFactorSV;

            explicit MyRT(): 
                func(&getValue), setFunc(&setValue){}
            
            static void setValue(void* structure, double value){
                throw ModelException("ZeroBondSV::setValue", "Cannot assign a"
                                     " value to a const variable");
            }
            static double getValue(void* structure){
                MyRT* rt = (MyRT*)structure;
                double PV = rt->discFactorSV->firstDF();
                return PV;
            }
            
            void* operator new(size_t size){
                return FR::MemMgr::alloc(size);
            }
            void operator delete(void *ptr){
                FR::MemMgr::dealloc(ptr);
            }
        };
        SVGenExpectedDiscFactorSP    discFactorGen;      //!< Generator for spot
        MyRT*                     rt;
    public:
        /** Populates the collector with the state variables required by the
            various assets */
        virtual void collectStateVars(IStateVariableCollectorSP svCollector) const{
            svCollector->append(discFactorGen.get()); // ask for spot SV
        }

        /** To be called when the path generator changes (currently before doing
            the past, and before doing the future) */
        virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen){
            rt->discFactorSV = discFactorGen->getSVExpectedDiscFactor(newPathGen);
        }
        
        // simple constructor
        ZeroBondSV(SVGenExpectedDiscFactorSP discFactorGen): 
            discFactorGen(discFactorGen), 
            rt(new MyRT()){}
        
        ~ZeroBondSV(){
            delete rt;
        }
        // get the variable expressed as a double
        virtual double getValue() {
            return MyRT::getValue(rt);
        }
        virtual FRIfaces::IRValueDouble::RT* getRT(){
            return (FRIfaces::IRValueDouble::RT*)rt;
        }
        // 'set' the variable using a double - this throws an exception
        virtual void setValue(double value){
            rt->setFunc(rt, value);
        }
        void set(const IObjectConstSP& object){
            rt->setFunc(rt, 0.0); // this throws an exception
        }
        virtual void setReset(char* reset){} // Path object handles this
    };

protected:
    /** creates FlexRule::IRValue which represents the value
        of this expression at a specific timepoint */
    virtual FRIfaces::IRValue* createRValue(int           index,
                                            FRController* frCtrl) const{
        static const string method = "FRZeroBondVariable::createRValue";

        const FRIfaces::IProductView* productView = frCtrl->productView();
        int                numDates = productView->getSimDates()->size();
        const DateTime&    today = productView->getValueDate();
        DateTimeArrayConstSP simDates(productView->getSimDates());
        YieldCurveConstSP  yc = !!ccyName ? ccyName.getSP() : productView->findYieldCurve();
        // sanity check
        if (!yc.get()) {
            throw ModelException(method, 
                                 "No market data available for YieldCurve " + 
                                 !!ccyName ? ccyName.getName() : "<discount>");
        }        
        if (frCtrl->stateVarUsed()){
            const DateTime& calcDate = (*simDates)[index];
            if (calcDate > maturityDate) {
                return new FR::RConstDouble(0.0);
            } else if (calcDate > today) {
                //// can just create variable for this specific time point
                // XXX! SRM currently requires same dates for rateStart & calc in SVGenDiscFactor
                // DateTime rateStartDate(yc->settles(calcDate));
                DateTime rateStartDate(calcDate);
                SVGenExpectedDiscFactorSP dfGen(
                    new SVGenExpectedDiscFactor(calcDate, // when to compute
                                             rateStartDate, // pv from here
                                             yc,
                                             DateTimeArray(1, maturityDate),
                                             false)); // compute PV not log

                ZeroBondSV* rValue = new ZeroBondSV(dfGen);
                /* register with controller */
                frCtrl->addToStateVars(rValue, false /* do not free - we are
                                                      returning this object */);
                return rValue;
            } else {
                // It's important that we actually set a value for each date
                // Unfortunately we don't know what it is (as historic)
                // Creating this FR::LValueDouble results in the value
                // throwing an exception if it is used (and not set)
                return new FR::LValueDouble(id.c_str());
            }
        }
        for (int iSimDate = 0; iSimDate < numDates; iSimDate++) {
            // avoid putting values in more than once - this can happen if
            // this routine is called again after it has failed (eg missing
            // sample). A calling method might ignore the error in the [vain]
            // hope that the particular expression will not be evaluated at
            // run time (eg consider 'a = i>3? value[i-3]: 0.0')
            if (!frCtrl->getRValue(this, iSimDate)){
                const DateTime& calcDate = (*simDates)[iSimDate];
                if (calcDate > maturityDate) {
                    FRIfaces::IRValueSP rValue(new FR::RConstDouble(0.0));
                    frCtrl->setRValue(iSimDate, this, rValue);
                } else if (calcDate > today) {
                    double pv = yc->pv(calcDate, maturityDate);
                    FRIfaces::IRValueSP rValue(new FR::RConstDouble(pv));
                    frCtrl->setRValue(iSimDate, this, rValue);
                } else {
                    // It's important that we actually set a value for each date
                    // Unfortunately we don't know what it is (as historic)
                    // Creating this FR::LValueDouble results in the value
                    // throwing an exception if it is used (and not set)
                    FRIfaces::IRValueSP rValue(
                        new FR::LValueDouble(id.c_str()));
                    frCtrl->setRValue(iSimDate, this, rValue);
                }
            }
        }
        return 0; // indicates that we've saved the value ourselves
    }

private:
    const string               id;      // transient - set upon construction
    const string               bondName;
    YieldCurveWrapper          ccyName; //e.g. "currency-GBP-ERB". Not const, note.
    const DateTime             maturityDate;

    FRZeroBondVariable(): LValueBase(TYPE) {}

    FRZeroBondVariable(const string    name,
                       const DateTime& maturityDate): 
        LValueBase(TYPE), bondName(name),
        maturityDate(maturityDate) {
        validatePop2Object();
    }

    FRZeroBondVariable(const string    name,
                       const string    ccyName,
                       const DateTime& maturityDate): 
        LValueBase(TYPE), bondName(name), ccyName(ccyName),
        maturityDate(maturityDate) {
        validatePop2Object();
    }

    static FR::LValueBase* create(const string& varName, const StringArray& value) {
        static const string routine = "FRZeroBondVariable::create";
        if (value.size()!=1) {
            throw ModelException(routine, 
                                 "Per sim date initialisation not appropriate here");
        }
        if (value[0].empty()) {
            throw ModelException(routine, "Value should contain ccyName, maturity date.");
        }

        const char* pos = value[0].c_str();
        string params[3];
        int paramId = 0;
        while (*pos) {
            char   c;
            // skip whitespace
            while ((c = *pos) == ' ' || c == '\t'){
                pos++;
            }
            if (c == ','){
                // skip over our chosen separator and continue
                pos++;
                // count params
                paramId++;
                if (paramId>1) {
                    throw ModelException(routine, "Too many commas with "
                                         "variable " + varName + 
                                         " and value of " + value[0]);
                }
            } else {
                // next contiguous block is the text we want
                const char* varBegin = pos;
                do{
                    pos++; /* Get another character. */
                    c = *pos;
                } while (c != '\0' && c != ',' && c != ' ' && c != '\t');
                const char* varEnd = pos;
                // done with this one
                if (params[paramId].empty()) {
                    // first time through get the date
                    params[paramId] = string(varBegin, varEnd - varBegin);
                } else {
                    // if there's a next time, it'll be the time
                    params[paramId+1] = string(varBegin, varEnd - varBegin);
                }                    
            } 
        }

        if (paramId > 1) {
            throw ModelException(routine, "For variable " + varName + 
                                 " require [ccyName,] maturityDate "
                                 "in that order. ccyName is optional. "
                                 " Given " + value[0]);
        }
        string timeStr;
        DateTime matDate;
        FR::LValueBase* zb;
        if (paramId == 1) {
            // convert params[1 & 2] into a DateTime
            if (!params[2].empty()) {
                timeStr = params[2];
            } else {
                timeStr = DateTime::END_OF_DAY;
            }
            matDate = DateTime(params[1], timeStr);
            
            zb = new FRZeroBondVariable(varName,
                                        params[0],
                                        matDate);
        } else {
            // no ccyName given so interpret Date/Time only
            if (!params[1].empty()) {
                timeStr = params[1];
            } else {
                timeStr = DateTime::END_OF_DAY;
            }
            try {
                matDate = DateTime(params[0], timeStr);
            } catch (exception& e) {
                throw ModelException(e, routine,
                                     "If a single parameter is supplied it must "
                                     "represent a DateTime");
            }
            zb = new FRZeroBondVariable(varName,
                                        matDate);
        }
        return zb;
    }

    static IObject* defaultFRZeroBondVariable(){
        return new FRZeroBondVariable();
    }
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FRZeroBondVariable, clazz);
        SUPERCLASS(FR::LValueBase);
        IMPLEMENTS(IGetMarket);
        EMPTY_SHELL_METHOD(defaultFRZeroBondVariable);
        FIELD(id, "");
        FIELD_MAKE_TRANSIENT(id);
        FIELD(bondName, "bond name");
        FIELD(ccyName, "currency Name - optional, default = inst ccy");
        FIELD_MAKE_OPTIONAL(ccyName);
        FIELD(maturityDate, "bond maturity");
        clazz->setPublic(); // make visible to EAS/spreadsheet
        Addin::registerConstructor("FR_VAR_ZERO_BOND",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a FR Zero Bond Variable",
                                   TYPE);
        FRVarFactory::registerCreateMethod(FRZeroBondVariable::create, "ZeroBond");
    }
    
};

CClassConstSP const FRZeroBondVariable::TYPE = CClass::registerClassLoadMethod(
    "FRZeroBondVariable", typeid(FRZeroBondVariable), load);
    
    
///////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////

class FRMMYieldVariable: public FR::LValueBase{
public:
    static CClassConstSP const TYPE;

    /** Returns the type that this variable represents */
    FRIfaces::VarType getType() const{
        return FRIfaces::doubleType;
    }

    /** returns a string that can be used to identify this object */
    virtual const string& getID() const{
        return id;
    }

    /** returns the input to define the variable in IMS */
    virtual StringArray getIMSInput() const {
       StringArray result(3);
        
        result[0] = "MMYield type is not implemented in IMS";
        result[1] = id;
        result[2] = "MMYield type is not implemented in IMS";

        return result;
    }

    virtual void validatePop2Object(){
        static const string method = "FRMMYieldVariable::validatePop2Object";

        // initialise our id
        const_cast<string&>(id) = mmYieldName; // validatePop2Object is privileged

        if (startDateOffset >= endDateOffset) {
            throw ModelException(method, 
                                 "startDateOffset (" + 
                                 Format::toString(startDateOffset) +
                                 ") must be < endDateOffset (" +  
                                 Format::toString(endDateOffset) + ")");
        }
        
        if (mmBasis != 360 &&
            mmBasis != 365) {
            throw ModelException(method, 
                                 "mmBasis (" + 
                                 Format::toString(mmBasis) +
                                 ") should be 360 or 365.");
        }

        // validate that any past values dates are increasing 
        if (pastValues.get()) {
            for (int i = 1; i < pastValues->size(); i++){
                if ((*pastValues)[i-1].date.isGreater((*pastValues)[i].date)){
                    throw ModelException(method, 
                                         "Past Value dates are not increasing : [" + 
                                         Format::toString(i-1) +
                                         "] " + (*pastValues)[i-1].date.toString() + 
                                         " > [" + Format::toString(i) + "] " +
                                         (*pastValues)[i].date.toString());
                }
            }
        }
    }

    virtual FRIfaces::ILValue* getLValue(FRController* frCtrl) const{
        throw ModelException("FRMMYieldVariable::getLValue", "Variable "+
                             id+" is read-only");
    }

protected:
    /** creates FlexRule::IRValue which represents the value
        of this expression at a specific timepoint */
    virtual FRIfaces::IRValue* createRValue(int           index,
                                            FRController* frCtrl) const{
        static const string method = "FRMMYieldVariable::createRValue";
        checkAgainstStateVars("FRMMYieldVariable:createRValue", frCtrl);

        // identify our ccyName
        const FRIfaces::IProductView* productView = frCtrl->productView();
        YieldCurveConstSP  discount = productView->findYieldCurve(&ccyName);
        int                numDates = productView->getSimDates()->size();
        int                numPast = pastValues.get()? pastValues->size() : 0;
        int                iPast = 0;
        bool               doSetFiller = !fillerSet;
        
        // Final endDateOffset # dates cannot define this
        if (index >= numDates - endDateOffset){
            throw ModelException(method, "Request for MM_YIELD on last "+
                                 Format::toString(endDateOffset) +" day(s) of "
                                 "simulation. It is not defined in this "
                                 "interval");
        }
        const DateTime& today = discount->getToday();
        DateTimeArrayConstSP simDates(productView->getSimDates());
        for (int iSimDate = 0; iSimDate < numDates - endDateOffset; iSimDate++)
        {
            // avoid putting values in more than once - this can happen if
            // this routine is called again after it has failed (eg missing
            // sample). A calling method might ignore the error in the [vain]
            // hope that the particular expression will not be evaluated at
            // run time (eg consider 'a = i>3? value[i-3]: 0.0')
            if (!frCtrl->getRValue(this, iSimDate)){
                const DateTime& startDate =
                    (*simDates)[iSimDate+startDateOffset];
                const DateTime& endDate = (*simDates)[iSimDate+endDateOffset];
                double   yield;

                // Careful with past - and rolling both "to today" and for theta
                if (startDate < today) {
                    // If past rely on user to supply the yield - locate
                    // based on date to be safe. 
                    DateTime simDate = (*simDates)[iSimDate];
                    for( ; iPast < numPast && 
                             (*pastValues)[iPast].date < simDate; iPast++) {
                        // empty - just finding index
                    }
                    // now iPast refers to date on or after simDate
                    if (iPast < numPast &&
                        (*pastValues)[iPast].date == simDate &&
                        Maths::isPositive((*pastValues)[iPast].amount)) {
                        // given past value to use
                        yield = (*pastValues)[iPast].amount;
                    } else if (fillerSet) {
                        // rolled past so use the one from pricing call
                        yield = fillerYield;
                    } else {
                        throw ModelException(method, "Past value missing for '"+
                                             id +"' on " + simDate.toString());
                    }
                } else {
                    double   pv = discount->pv(startDate, endDate);
                    double   ndays = endDate.daysDiff(startDate);
                    yield = (double)mmBasis / ndays * (1.0 / pv - 1.0);
                    
                    if (doSetFiller) {
                        fillerYield = yield;
                        fillerSet = true;
                    }
                }
                
                FRIfaces::IRValueSP rValue(new FR::RConstDouble(yield));
                frCtrl->setRValue(iSimDate, this, rValue);
            }
        }
        return 0; // indicates that we've saved the value ourselves
    }

private:
    // fields
    const string               id;      // transient - set upon construction
    const string               mmYieldName;
    const string               ccyName;
    const int                  startDateOffset; // # sim dates
    const int                  endDateOffset;   // # sim dates
    const int                  mmBasis;         // 360 or 365
    const CashFlowArrayConstSP pastValues;

    // internal
    mutable double       fillerYield; // $unregistered
    mutable bool         fillerSet; // $unregistered

    FRMMYieldVariable(): LValueBase(TYPE), startDateOffset(0), endDateOffset(0),
                         mmBasis(360), fillerYield(0), fillerSet(false) {}

    static IObject* defaultFRMMYieldVariable(){
        return new FRMMYieldVariable();
    }
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FRMMYieldVariable, clazz);
        SUPERCLASS(FR::LValueBase);
        EMPTY_SHELL_METHOD(defaultFRMMYieldVariable);
        FIELD(id, "");
        FIELD_MAKE_TRANSIENT(id);
        FIELD(mmYieldName, "mm yield name");
        FIELD(ccyName, "currency Name");
        FIELD(startDateOffset, "startDateOffset");
        FIELD(endDateOffset, "endDateOffset");
        FIELD(mmBasis, "360/365");
        FIELD(pastValues, "past yield values");
        FIELD_MAKE_OPTIONAL(pastValues);
        clazz->setPublic(); // make visible to EAS/spreadsheet
        Addin::registerConstructor("FR_VAR_MM_YIELD",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a FR Money Market Yield Variable",
                                   TYPE);
    }
    
};

CClassConstSP const FRMMYieldVariable::TYPE = CClass::registerClassLoadMethod(
    "FRMMYieldVariable", typeid(FRMMYieldVariable), load);
    
///////////////////////////////////////////////////////
// FRMMRateVariable. A different flavour of the above, with support
// for on-line interface Note that currently this will fail if a value
// is requested in the past.  We are therefore requiring the rules to
// be amended with each MMRate reference in the past.
///////////////////////////////////////////////////////
class FRMMRateVariable: public FR::LValueBase,
                        virtual public IGetMarket {
public:
    static CClassConstSP const TYPE;

    /** Historically (before parser) shallow copies were employed 
        to avoid an exponential explosion in the memory demands of 'operation'
        classes. This is not really relevant with the parser, but the shallow
        clone remains. Note that for this to be correct all the fields need
        to be const.
        However, upon the introduction of some contained market data (YieldCurveWrapper)
        we no longer have all fields being 'const' (getMarket updates the value) and
        tweaking relies on clone() to return a new instance. So here (and in a few
        other instances) we need to restore the deep copy. 
        An alternative would be to remove the clone override on the parent, but this
        is a more conservative approach. */
    IObject* clone() const{
        return CObject::clone();
    }

    /** Returns the type that this variable represents */
    FRIfaces::VarType getType() const{
        return FRIfaces::doubleType;
    }

    /** returns a string that can be used to identify this object */
    virtual const string& getID() const{
        return id;
    }

    virtual void validatePop2Object(){
        // initialise our id
        const_cast<string&>(id) = irName; // validatePop2Object is privileged
        dcc = DayCountConventionSP(
            DayCountConventionFactory::make(dayCountConv));
        rateMat = MaturityPeriodSP(new MaturityPeriod(maturityPeriod));
        if (!compounding.empty() &&
            !CString::equalsIgnoreCase(compounding, 
                                       CompoundBasis::toString(CompoundBasis::SIMPLE)) &&
            !CString::equalsIgnoreCase(compounding, 
                                       CompoundBasis::toString(CompoundBasis::ANNUAL))) {
            throw ModelException("FRMMRateVariable::validatePop2Object",
                                 "Compounding must be empty, " + 
                                 CompoundBasis::toString(CompoundBasis::SIMPLE) + 
                                 " or " +
                                 CompoundBasis::toString(CompoundBasis::ANNUAL));
        }
    }

    virtual FRIfaces::ILValue* getLValue(FRController* frCtrl) const{
        throw ModelException("FRMMRateVariable::getLValue", "Variable "+
                             id+" is read-only");
    }

    /** returns the input to define the variable in IMS */
    virtual StringArray getIMSInput() const {
        StringArray result;
        
        result.push_back("MMRate");
        result.push_back(id);
        string data(ccyName.getName()+","+maturityPeriod+","+
                    dayCountConv+","+compounding);
        if (pastFixings.get() &&
            !pastFixings->empty()) {
            data = data + ",Fixings";
            result.push_back(data); 
            StringArray fixings(CashFlow::toStringArray(*pastFixings.get()));
            for(int i=0; i<fixings.size(); i++) {
                result.push_back(fixings[i]);
            }
        } else {
            result.push_back(data); 
        }

        return result;
    }

    // IGetMarket
    virtual void getMarket(const IModel* model, const MarketData* market) {
        ccyName.getData(model, market);
    }

    /** Class for holding expected discount factors along each point*/
    class MMRateSV: public FR::RValueDouble, 
                    public virtual FRIfaces::ILValueDouble,
                    public virtual FRController::IStateVarClient{
        struct MyRT{
            TGetValue*                func;
            TSetValue*                setFunc;
            SVExpectedDiscFactorSP discFactorSV;
            double                    minusReciprocalYearFrac; // = -1/yearFrac
            bool                      isSimpleCompounding;

            explicit MyRT(double yearFrac,
                          bool   isSimpleCompounding): 
                func(&getValue), setFunc(&setValue), 
                minusReciprocalYearFrac(-1.0/yearFrac),
                isSimpleCompounding(isSimpleCompounding){}
            
            static void setValue(void* structure, double value){
                throw ModelException("MMRateSV::setValue", "Cannot assign a"
                                     " value to a const variable");
            }
            static double getValue(void* structure){
                MyRT* rt = (MyRT*)structure;
                // need to review/test this carefully to make sure it's safe
                // to ignore begin()/end()
                // Get [log of] discount factor between dates
                double logOfPV = rt->discFactorSV->firstDF();
                // and then turn into a rate.
                if (rt->isSimpleCompounding) {
                    double simple = (1.0 - exp(-logOfPV)) * rt->minusReciprocalYearFrac;
                    return simple; // using simple convention
                    // simple: coupon = r*yearFrac
                    // = exp(-logOfPV) - 1.0
                }
                // NB pow(pv,b) = pv^b = exp(b*ln pv)
                double rate = exp(logOfPV * rt->minusReciprocalYearFrac) - 1.0;
                return rate; // using annual convention
                // annual: coupon = (1+r)^yearFrac - 1.0
                // = exp(logOfPV * -1/yearFrac)^yearFrac - 1.0
                // = exp(logOfPV * -1) - 1.0 = 1/pv -1.0 (more useful?)
            }
            
            void* operator new(size_t size){
                return FR::MemMgr::alloc(size);
            }
            void operator delete(void *ptr){
                FR::MemMgr::dealloc(ptr);
            }
        };
        SVGenExpectedDiscFactorSP    discFactorGen;      //!< Generator for spot
        MyRT*                     rt;
    public:
        /** Populates the collector with the state variables required by the
            various assets */
        virtual void collectStateVars(IStateVariableCollectorSP svCollector) const{
            svCollector->append(discFactorGen.get()); // ask for spot SV
        }

        /** To be called when the path generator changes (currently before doing
            the past, and before doing the future) */
        virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen){
            rt->discFactorSV = discFactorGen->getSVExpectedDiscFactor(newPathGen);
        }
        
        // simple constructor
        MMRateSV(double                 yearFrac,
                 bool                   isSimpleCompounding,
                 SVGenExpectedDiscFactorSP discFactorGen): 
            discFactorGen(discFactorGen), 
            rt(new MyRT(yearFrac, isSimpleCompounding)){}
        
        ~MMRateSV(){
            delete rt;
        }
        // get the variable expressed as a double
        virtual double getValue() {
            return MyRT::getValue(rt);
        }
        virtual FRIfaces::IRValueDouble::RT* getRT(){
            return (FRIfaces::IRValueDouble::RT*)rt;
        }
        // 'set' the variable using a double - this throws an exception
        virtual void setValue(double value){
            rt->setFunc(rt, value);
        }
        void set(const IObjectConstSP& object){
            rt->setFunc(rt, 0.0); // this throws an exception
        }
        virtual void setReset(char* reset){} // Path object handles this
    };

protected:
    /* Suppose that we can look up FRController::IStateVarClient based upon 
       this. Then look up each time we get here. If null, create one.
       It needs to have */

    /** creates FlexRule::IRValue which represents the value
        of this expression at a specific timepoint */
    virtual FRIfaces::IRValue* createRValue(int           index,
                                            FRController* frCtrl) const{
        static const string method = "FRMMRateVariable::createRValue";

        const FRIfaces::IProductView* productView = frCtrl->productView();
        int                numDates = productView->getSimDates()->size();
        const DateTime&    today = productView->getValueDate();
        DateTimeArrayConstSP simDates(productView->getSimDates());
        bool isSimpleCompounding = CString::equalsIgnoreCase(compounding, 
                                                             CompoundBasis::toString(CompoundBasis::SIMPLE));
        YieldCurveConstSP  yc = ccyName.getSP();
        // sanity check
        if (!yc.get()) {
            throw ModelException(method, 
                                 "No market data available for YieldCurve " + ccyName.getName());
        }
        if (frCtrl->stateVarUsed()){
            const DateTime& calcDate = (*simDates)[index];
            if (calcDate > today) {
                //// can just create variable for this specific time point
                // XXX! SRM currently requires same dates for rateStart & calc in SVGenDiscFactor
                // DateTime rateStartDate(yc->settles(calcDate));
                DateTime rateStartDate(calcDate);
                DateTime rateMatDate(yc->rateMaturity(rateStartDate,
                                                      rateMat.get(),
                                                      0 /* bad day */));
                SVGenExpectedDiscFactorSP dfGen(
                    new SVGenExpectedDiscFactor(calcDate, // when to compute
                                             rateStartDate, // pv from here
                                             yc,
                                             DateTimeArray(1, rateMatDate),
                                             true)); // compute log
                double yearFrac = dcc->years(rateStartDate, rateMatDate);
                // All other values mean annual compounding (the default)
                MMRateSV* rValue = new MMRateSV(yearFrac, isSimpleCompounding, dfGen);
                /* register with controller */
                frCtrl->addToStateVars(rValue, false /* do not free - we are
                                                        returning this object */);
                return rValue;
            } else {
                // Expect an entry in pastFixings for this ...
                double rate;
                bool isMissing = false;
                try {
                    if (!pastFixings.get()) {
                        isMissing = true;
                    } else {
                        rate = CashFlow::interpolate(*pastFixings.get(),
                                                     calcDate);
                    }
                } catch (exception& ) {
                    isMissing = true;
                }
                if (isMissing) {
                    // Creating this FR::LValueDouble results in the value
                    // throwing an exception if it is used (and not set)
                    return new FR::LValueDouble(id.c_str());
                }
                return new FR::RConstDouble(rate);
            }
        }
        for (int iSimDate = 0; iSimDate < numDates; iSimDate++) {
            // avoid putting values in more than once - this can happen if
            // this routine is called again after it has failed (eg missing
            // sample). A calling method might ignore the error in the [vain]
            // hope that the particular expression will not be evaluated at
            // run time (eg consider 'a = i>3? value[i-3]: 0.0')
            if (!frCtrl->getRValue(this, iSimDate)){
                const DateTime& fixDate = (*simDates)[iSimDate];
                if (fixDate > today) {
                    double irFwd = yc->fwd(fixDate,
                                           rateMat.get(),
                                           0, // use YC bad day conv
                                           dcc.get(),
                                           isSimpleCompounding? 
                                           CompoundBasis::SIMPLE:
                                           CompoundBasis::ANNUAL);
                    FRIfaces::IRValueSP rValue(new FR::RConstDouble(irFwd));
                    frCtrl->setRValue(iSimDate, this, rValue);
                } else {
                    // Expect an entry in pastFixings for this ...
                    double rate;
                    bool isMissing = false;
                    try {
                        if (!pastFixings.get()) {
                            isMissing = true;
                        } else {
                            rate = CashFlow::interpolate(*pastFixings.get(), 
                                                         fixDate);
                        }
                    } catch (exception& ) {
                        isMissing = true;
                    }
                    FRIfaces::IRValueSP rValue;
                    if (isMissing) {
                        // Creating this FR::LValueDouble results in the value
                        // throwing an exception if it is used (and not set)
                        rValue = FRIfaces::IRValueSP(new FR::LValueDouble(id.c_str()));
                    } else {
                        rValue = FRIfaces::IRValueSP(new FR::RConstDouble(rate));
                    }
                    frCtrl->setRValue(iSimDate, this, rValue);
                }
            }
        }
        return 0; // indicates that we've saved the value ourselves
    }

private:
    // fields
    const string         id;      // transient - set upon construction
    const string         irName;          //e.g. "R"
    YieldCurveWrapper    ccyName;         //e.g. "currency-GBP-ERB". Not const, note.

    const string         maturityPeriod;  //e.g. "3M"
    const string         dayCountConv;    //e.g. "Act/360"
    const string         compounding;     //[optional] e.g. "Annual", "Simple"
    CashFlowArrayConstSP pastFixings;     //[optional] 

    //transient
    MaturityPeriodSP     rateMat;
    DayCountConventionSP dcc;

    FRMMRateVariable(): LValueBase(TYPE) {}

    FRMMRateVariable(const string               irName,
                     const string               ccyName,
                     const string               maturityPeriod, //e.g. "3M"
                     const string               dayCountConv,  //e.g. "Act/360"
                     const string               compounding,
                     CashFlowArrayConstSP       pastFixings): 
        LValueBase(TYPE), irName(irName), ccyName(ccyName),
        maturityPeriod(maturityPeriod), dayCountConv(dayCountConv),
        compounding(compounding), pastFixings(pastFixings) {
        // finish the job - notably assigning to 'id' also
        validatePop2Object();
    }
    
    static IObject* defaultFRMMRateVariable(){
        return new FRMMRateVariable();
    }

    static FR::LValueBase* create(const string& varName, const StringArray& value) {
        static const string routine = "FRMMRateVariable::create";
        if (value[0].empty()) {
            throw ModelException(routine, "Value should contain ccy name, "
                                 "DCC, maturity period and [optionally] compounding basis."
                                 "May also give lookup id for fixings (this requires compounding basis to be explicit)");
        }
        string params[5];
        // [0]=ccyName, [1]=maturityPeriod, 
        // [2]=dayCountConv , [3]=compounding basis
        // in that order and comma separated
        // An optional [4]th param is the lookup id for 
        // past fixings. That will trigger value.size()>1
        // containing past fixing info. We can ignore the
        // param here.
        const char* pos = value[0].c_str();
        int paramId = 0;
        while (*pos) {
            char   c;
            // skip whitespace
            while ((c = *pos) == ' ' || c == '\t'){
                pos++;
            }
            if (c == ','){
                // skip over our chosen separator and continue
                pos++;
                // count params
                paramId++;
                if (paramId>4) {
                    throw ModelException(routine, "Too many commas with "
                                         "variable " + varName + 
                                         " and value of " + value[0]);
                }
            } else {
                // next contiguous block is the text we want
                const char* varBegin = pos;
                do{
                    pos++; /* Get another character. */
                    c = *pos;
                } while (c != '\0' && c != ',' && c != ' ' && c != '\t');
                const char* varEnd = pos;
                // done with this one
                params[paramId] = string(varBegin, varEnd - varBegin);
            } 
        }

        CashFlowArraySP pastFixings;
        StringArray fixings;
        switch (paramId) {
        case 2:
            // provide a default for optional 'compounding'
            params[3] = CompoundBasis::toString(CompoundBasis::ANNUAL);
            // no fixings can be specified in this case
            break;
        case 3:
            // no fixings specified in this case
            break;
        case 4:
            // fixings indicated...but we ignore the param (used earlier)
            fixings = value;
            if (fixings.size()<2) {
                throw ModelException(routine, "No fixings indicated!? Internal error!");
            }
            fixings.erase(fixings.begin(), fixings.begin()+1);
            pastFixings = CashFlowArraySP(new CashFlowArray(CashFlow::fromStringArray(fixings)));
            break;
        default:
            throw ModelException(routine, "For variable " + varName + 
                                 " require ccyName, maturityPeriod, "
                                 "dayCountConv[, compounding basis][,fixings lookup id] "
                                 "in that order and comma "
                                 "separated but given " + value[0]);
        }

        // A simple "return new ..." gives gcc compilation warning 
        FR::LValueBase* mmr = new FRMMRateVariable(varName,
                                                   params[0],
                                                   params[1],
                                                   params[2],
                                                   params[3],
                                                   pastFixings);
        return mmr;
    }

    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FRMMRateVariable, clazz);
        SUPERCLASS(FR::LValueBase);
        IMPLEMENTS(IGetMarket);
        EMPTY_SHELL_METHOD(defaultFRMMRateVariable);
        FIELD(id, "");
        FIELD_MAKE_TRANSIENT(id);
        FIELD(irName, "variable's name");
        FIELD(ccyName, "e.g. currency-GBP-ERB");
        FIELD(maturityPeriod, "e.g. 3M");
        FIELD(dayCountConv, "e.g. Act/360");
        FIELD(compounding, "Empty, Annual or Simple");
        FIELD_MAKE_OPTIONAL(compounding);
        FIELD(rateMat, "rate maturity");
        FIELD_MAKE_TRANSIENT(rateMat);
        FIELD(dcc, "dcc");
        FIELD_MAKE_TRANSIENT(dcc);
        FIELD(pastFixings, "Historical rate fixings");
        FIELD_MAKE_OPTIONAL(pastFixings);
        clazz->setPublic(); // make visible to EAS/spreadsheet
        Addin::registerConstructor("FR_VAR_MM_RATE",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a FR Money Market Rate Variable",
                                   TYPE);
        FRVarFactory::registerCreateMethod(FRMMRateVariable::create, "MMRate");
    }
};

CClassConstSP const FRMMRateVariable::TYPE = CClass::registerClassLoadMethod(
    "FRMMRateVariable", typeid(FRMMRateVariable), load);
    
///////////////////////////////////////////////////////
// FRSwapRateVariable. Similar to the above but it's a swap rate!
///////////////////////////////////////////////////////
class FRSwapRateVariable: public FR::LValueBase,
                          virtual public IGetMarket{
public:
    static CClassConstSP const TYPE;
    static const string FRONT_STUB;
    static const string BACK_STUB;

    /** See comment in FRMMRateVariable */
    IObject* clone() const{
        return CObject::clone();
    }

    /** Returns the type that this variable represents */
    FRIfaces::VarType getType() const{
        return FRIfaces::doubleType;
    }

    /** returns a string that can be used to identify this object */
    virtual const string& getID() const{
        return id;
    }

    virtual void validatePop2Object(){
        static const string method("FRSwapRateVariable::validatePop2Object");
        // initialise our id
        const_cast<string&>(id) = name; // validatePop2Object is privileged
        dcc = DayCountConventionSP(
            DayCountConventionFactory::make(dayCountConv));
        try {
            swapMat = MaturityPeriodSP(new MaturityPeriod(maturityPeriod));
        } catch(exception& e) {
            throw ModelException(e, method, 
                                 "Failed to identify swap maturity period : " +
                                 maturityPeriod);
        }
        try {
            couponFreq = MaturityPeriodSP(new MaturityPeriod(1, couponPeriod));
        } catch(exception& e) {
            throw ModelException(e, method, 
                                 "Failed to identify coupon period : " + 
                                 couponPeriod + 
                                 ". Note we require NO number, just a letter.");
        }
        if (CString::equalsIgnoreCase(stubLocation, FRONT_STUB)) {
            stubAtEnd = false;
        } else if (CString::equalsIgnoreCase(stubLocation, BACK_STUB)) {
            stubAtEnd = true; 
        } else {
            throw ModelException(method, 
                                 "stubLocation (" + stubLocation + 
                                 ") should be " + FRONT_STUB +
                                 " or " + BACK_STUB);
        }
    }

    virtual FRIfaces::ILValue* getLValue(FRController* frCtrl) const{
        throw ModelException("FRSwapRateVariable::getLValue", "Variable "+
                             id+" is read-only");
    }

    /** returns the input to define the variable in IMS */
    virtual StringArray getIMSInput() const {
        StringArray result;
        
        result.push_back("SwapRate");
        result.push_back(id);
        string data(ccyName.getName()+","+maturityPeriod+","+
            couponPeriod+","+stubType+","+
            stubLocation+","+dayCountConv);
        if (pastFixings.get() &&
            !pastFixings->empty()) {
            data = data + ",Fixings";
            result.push_back(data); 
            StringArray fixings(CashFlow::toStringArray(*pastFixings.get()));
            for(int i=0; i<fixings.size(); i++) {
                result.push_back(fixings[i]);
            }
        } else {
            result.push_back(data); 
        }

        return result;
    }

    // IGetMarket
    virtual void getMarket(const IModel* model, const MarketData* market) {
        ccyName.getData(model, market);
    }

    /** Class for holding swap rates along each point*/
    class SwapRateSV: public FR::RValueDouble, 
                      public virtual FRIfaces::ILValueDouble,
                      public virtual FRController::IStateVarClient{
        struct MyRT{
            TGetValue*                func;
            TSetValue*                setFunc;
            SVGenIRSwap::IStateVarSP     irSwapSV;

            explicit MyRT(): 
                func(&getValue), setFunc(&setValue){}
            
            static void setValue(void* structure, double value){
                throw ModelException("SwapRateSV::setValue", "Cannot assign a"
                                     " value to a const variable");
            }
            static double getValue(void* structure){
                static const string method = "SwapRateSV::getValue";
                MyRT* rt = (MyRT*)structure;

                double parYield;
                double annuity;
                rt->irSwapSV->parYield(parYield, annuity);
                return parYield;
            }
            
            void* operator new(size_t size){
                return FR::MemMgr::alloc(size);
            }
            void operator delete(void *ptr){
                FR::MemMgr::dealloc(ptr);
            }
        };
        SVGenIRSwapSP    irSwapGen;
        MyRT*         rt;
    public:
        /** Populates the collector with the state variables required by the
            various assets */
        virtual void collectStateVars(IStateVariableCollectorSP svCollector) const{
            svCollector->append(irSwapGen.get()); 
        }

        /** To be called when the path generator changes (currently before doing
            the past, and before doing the future) */
        virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen){
            rt->irSwapSV = irSwapGen->getIRSwapSV(newPathGen);
        }
        
        // simple constructor
        SwapRateSV(SVGenIRSwapSP irSwapGen): 
            irSwapGen(irSwapGen), rt(new MyRT()){}
        
        ~SwapRateSV(){
            delete rt;
        }
        // get the variable expressed as a double
        virtual double getValue() {
            return MyRT::getValue(rt);
        }
        virtual FRIfaces::IRValueDouble::RT* getRT(){
            return (FRIfaces::IRValueDouble::RT*)rt;
        }
        // 'set' the variable using a double - this throws an exception
        virtual void setValue(double value){
            rt->setFunc(rt, value);
        }
        void set(const IObjectConstSP& object){
            rt->setFunc(rt, 0.0); // this throws an exception
        }
        virtual void setReset(char* reset){} // Path object handles this
    };

protected:
    /* Suppose that we can look up FRController::IStateVarClient based upon 
       this. Then look up each time we get here. If null, create one.
       It needs to have */

    /** creates FlexRule::IRValue which represents the value
        of this expression at a specific timepoint */
    virtual FRIfaces::IRValue* createRValue(int           index,
                                            FRController* frCtrl) const{
        static const string method = "FRSwapRateVariable::createRValue";

        const FRIfaces::IProductView* productView = frCtrl->productView();
        int                numDates = productView->getSimDates()->size();
        const DateTime&    today = productView->getValueDate();
        DateTimeArrayConstSP simDates(productView->getSimDates());
        YieldCurveConstSP  yc = ccyName.getSP();
        // sanity check
        if (!yc.get()) {
            throw ModelException(method, 
                                 "No market data available for YieldCurve " + ccyName.getName());
        }

        if (frCtrl->stateVarUsed()){
            const DateTime& calcDate = (*simDates)[index];
            if (calcDate > today) {
                //// can just create variable for this specific time point
                // XXX! SRM currently requires same dates for rateStart & calc in SVGenDiscFactor
                // DateTime rateStartDate(yc->settles(calcDate));
                DateTime rateStartDate(calcDate);
                DateTime swapMatDate(yc->rateMaturity(rateStartDate,
                                                      swapMat.get(),
                                                      0 /* bad day */));
                SVGenIRSwapSP irSwapGen(new SVGenIRSwap(yc, // coupon curve
                                                  yc, // => use general formula
                                                  rateStartDate,
                                                  swapMatDate,
                                                  couponPeriod, // fixedPayInterval
                                                  dayCountConv,
                                                  stubType,
                                                  stubAtEnd,
                                                  "None", //accrueBadDayConv
                                                  "None", //payBadDayConv
                                                  HolidaySP(Holiday::noHolidays()),
                                                  true/* NOT USED! isCashSettled*/));

                SwapRateSV* rValue = new SwapRateSV(irSwapGen);
                /* register with controller */
                frCtrl->addToStateVars(rValue, false /* do not free - we are
                                                      returning this object */);
                return rValue;
            } else {
                // Expect an entry in pastFixings for this ...
                double rate;
                bool isMissing = false;
                try {
                    if (!pastFixings.get()) {
                        isMissing = true;
                    } else {
                        rate = CashFlow::interpolate(*pastFixings.get(),
                                                     calcDate);
                    }
                } catch (exception& ) {
                    isMissing = true;
                }
                if (isMissing) {
                    // Creating this FR::LValueDouble results in the value
                    // throwing an exception if it is used (and not set)
                    return new FR::LValueDouble(id.c_str());
                }
                return new FR::RConstDouble(rate);
            }
        }
        for (int iSimDate = 0; iSimDate < numDates; iSimDate++) {
            // avoid putting values in more than once - this can happen if
            // this routine is called again after it has failed (eg missing
            // sample). A calling method might ignore the error in the [vain]
            // hope that the particular expression will not be evaluated at
            // run time (eg consider 'a = i>3? value[i-3]: 0.0')
            if (!frCtrl->getRValue(this, iSimDate)){
                const DateTime& fixDate = (*simDates)[iSimDate];
                if (fixDate > today) {
                    // don't handle stubType here. What to do?
                    DateTime rateStartDate(yc->settles(fixDate));
                    DateTime swapMatDate(yc->rateMaturity(rateStartDate,
                                                          swapMat.get(),
                                                          0 /* bad day */));
                    double swapRate = yc->couponRate(rateStartDate,
                                                     swapMatDate,
                                                     *couponFreq.get(),
                                                     stubAtEnd,
                                                     dcc.get());
                    FRIfaces::IRValueSP rValue(new FR::RConstDouble(swapRate));
                    frCtrl->setRValue(iSimDate, this, rValue);
                } else {
                    // Expect an entry in pastFixings for this ...
                    double rate;
                    bool isMissing = false;
                    try {
                        if (!pastFixings.get()) {
                            isMissing = true;
                        } else {
                            rate = CashFlow::interpolate(*pastFixings.get(), 
                                                         fixDate);
                        }
                    } catch (exception& ) {
                        isMissing = true;
                    }
                    FRIfaces::IRValueSP rValue;
                    if (isMissing) {
                        // Creating this FR::LValueDouble results in the value
                        // throwing an exception if it is used (and not set)
                        rValue = FRIfaces::IRValueSP(new FR::LValueDouble(id.c_str()));
                    } else {
                        rValue = FRIfaces::IRValueSP(new FR::RConstDouble(rate));
                    }
                    frCtrl->setRValue(iSimDate, this, rValue);
                }
            }
        }
        return 0; // indicates that we've saved the value ourselves
    }

private:
    // fields
    const string         id;      // transient - set upon construction
    const string         name;            //e.g. "R"
    YieldCurveWrapper    ccyName;         //e.g. "currency-GBP-ERB". Not const, note.
    const string         maturityPeriod;  //e.g. "2Y"
    const string         couponPeriod;    //of coupons, "M", "Q", "S", "A"
    const string         stubType;        //e.g. "None", "Bond", "Swap"
    const string         stubLocation;    //e.g. "frontStub", "backStub"
    const string         dayCountConv;    //e.g. "Act/360"
    CashFlowArrayConstSP pastFixings;     //[optional] 

    //transient
    MaturityPeriodSP     swapMat;    // from maturityPeriod
    MaturityPeriodSP     couponFreq; // from frequency
    bool                 stubAtEnd;  // from stubLocation
    DayCountConventionSP dcc;

    FRSwapRateVariable(): LValueBase(TYPE) {}

    FRSwapRateVariable(string               name,
                       string               ccyName,
                       string               maturityPeriod,
                       string               couponPeriod,
                       string               stubType,
                       string               stubLocation,
                       string               dayCountConv,
                       CashFlowArrayConstSP pastFixings):
        LValueBase(TYPE), name(name), ccyName(ccyName),
        maturityPeriod(maturityPeriod), couponPeriod(couponPeriod),
        stubType(stubType),
        stubLocation(stubLocation), dayCountConv(dayCountConv), 
        pastFixings(pastFixings) {
        // finish the job - notably assigning to 'id' also
        validatePop2Object();
    }
    
    static IObject* defaultFRSwapRateVariable(){
        return new FRSwapRateVariable();
    }

    static FR::LValueBase* create(const string& varName, const StringArray& value) {
        static const string routine = "FRSwapRateVariable::create";
        string params[7];
        int paramId = 0;
        if (!value[0].empty()) {
            // [0]=ccyName, 
            // [1]=maturityPeriod, 
            // [2]=couponPeriod, 
            // [3]=stubType,
            // [4]=stubLocation,
            // [5]=dayCountConv 
            // [6]=fixings id
            // in that order and comma separated
            const char* pos = value[0].c_str();
            while (*pos) {
                char   c;
                // skip whitespace
                while ((c = *pos) == ' ' || c == '\t'){
                    pos++;
                }
                if (c == ','){
                    // skip over our chosen separator and continue
                    pos++;
                    // count params
                    paramId++;
                    if (paramId>6) {
                        throw ModelException(routine, "Too many commas with "
                                             "variable " + varName + 
                                             " and value of " + value[0]);
                    }
                } else {
                    // next contiguous block is the text we want
                    const char* varBegin = pos;
                    do{
                        pos++; /* Get another character. */
                        c = *pos;
                    } while (c != '\0' && c != ',' && c != ' ' && c != '\t');
                    const char* varEnd = pos;
                    // done with this one
                    params[paramId] = string(varBegin, varEnd - varBegin);
                } 
            }
        }
        CashFlowArraySP pastFixings;
        StringArray fixings;
        if (paramId == 5) {
            // no pastFixings but otherwise ok
        } else if (paramId == 6) {
            // fixings indicated...but we ignore the param (used earlier)
            fixings = value;
            if (fixings.size()<2) {
                throw ModelException(routine, "No fixings indicated!? Internal error!");
            }
            fixings.erase(fixings.begin(), fixings.begin()+1);
            pastFixings = CashFlowArraySP(new CashFlowArray(CashFlow::fromStringArray(fixings)));
        } else {
            // here if value.empty() too
            throw ModelException(routine, "For variable " + varName + 
                                 " require ccyName, maturityPeriod, "
                                 "couponPeriod, stubType, stubLocation, "
                                 "dayCountConv[,fixings id] in that order and comma "
                                 "separated but given " + 
                                 (value[0].empty()?"an empty string":value[0]));
        }

        // A simple "return new ..." gives gcc compilation warning 
        FR::LValueBase* var = new FRSwapRateVariable(varName,
                                                     params[0], 
                                                     params[1], 
                                                     params[2],
                                                     params[3],
                                                     params[4],
                                                     params[5],
                                                     pastFixings);
        return var;
    }

    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FRSwapRateVariable, clazz);
        SUPERCLASS(FR::LValueBase);
        IMPLEMENTS(IGetMarket);
        EMPTY_SHELL_METHOD(defaultFRSwapRateVariable);
        FIELD(id, "");
        FIELD_MAKE_TRANSIENT(id);
        FIELD(name, "variable's name");
        FIELD(ccyName, "e.g. currency-GBP-ERB");
        FIELD(maturityPeriod, "Maturity of the CMS, e.g. 2Y");
        FIELD(couponPeriod, "M, Q, S, A (no number)");
        FIELD(stubType, "None, Bond or Swap");
        FIELD(stubLocation, "frontStub or backStub");
        FIELD(dayCountConv, "e.g. Act/360");
        FIELD_NO_DESC(swapMat);
        FIELD_MAKE_TRANSIENT(swapMat);
        FIELD_NO_DESC(couponFreq);
        FIELD_MAKE_TRANSIENT(couponFreq);
        FIELD_NO_DESC(stubAtEnd);
        FIELD_MAKE_TRANSIENT(stubAtEnd);
        FIELD_NO_DESC(dcc);
        FIELD_MAKE_TRANSIENT(dcc);
        FIELD(pastFixings, "Historical rate fixings");
        FIELD_MAKE_OPTIONAL(pastFixings);
        clazz->setPublic(); // make visible to EAS/spreadsheet
        Addin::registerConstructor("FR_VAR_SWAP_RATE",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a Flex Swap Rate Variable",
                                   TYPE);
        FRVarFactory::registerCreateMethod(FRSwapRateVariable::create, "SwapRate");
    }
};

const string FRSwapRateVariable::FRONT_STUB = "frontStub";
const string FRSwapRateVariable::BACK_STUB = "backStub";

CClassConstSP const FRSwapRateVariable::TYPE = CClass::registerClassLoadMethod(
    "FRSwapRateVariable", typeid(FRSwapRateVariable), load);
    
///////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////


bool loadFRBondVariable() {
    return (FRBondVariable::TYPE != 0);
}

DRLIB_END_NAMESPACE

