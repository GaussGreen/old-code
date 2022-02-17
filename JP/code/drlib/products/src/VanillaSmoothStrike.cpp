//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VanillaSmoothStrike.cpp
//
//   Description : Smooth the strike for big vanlla positions
//                 It's really just the weighted sum of vanillas with
//                 different strikes
//
//   Author      : Andrew J Swain
//
//   Date        : 18 August 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Generic1Factor.hpp"
#include "edginc/PassThroughModel.hpp"
#include "edginc/UtilFuncs.hpp"
#include "edginc/Vanilla.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"

DRLIB_BEGIN_NAMESPACE

class VanillaSmoothStrike: public Generic1Factor,
                           public PassThroughModel::IIntoProduct,
                           public LastSensDate {
public:
    static CClassConstSP const TYPE; 

    virtual void validatePop2Object(){
        static const string method = "VanillaSmoothStrike::validatePop2Object";
        try { 
            if (weights.empty() || strikes.empty()) {
                throw ModelException(method, "must have at least one weight "
                                     "or strike");
            }
                
            if (weights.size() != strikes.size()) {
                throw ModelException(method, "weights and strikes must be same "
                                     "length");
            }

        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    virtual void Validate() {
        static const string method("VanillaSmoothStrike::Validate");
        try { 
            AssetUtil::assetCrossValidate(asset.get(),
                                          fwdStarting,
                                          startDate,
                                          valueDate,
                                          discount,
                                          this);
            if (fwdStarting && oneContract) {
                throw ModelException(method, "Can't be forward starting and "
                                     "one contract");
            }

            if (!american && isExercised && valueDate < maturity) {
                throw ModelException(method,
                                     "option is european but is already "
                                     "exercised - today ("+valueDate.toString()+
                                     ") is before maturity (" + 
                                     maturity.toString() + ")");
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }    
    }

    /** Indicates whether VEGA_MATRIX is sensible for this instrument.*/
    bool avoidVegaMatrix(const IModel*model){
        return false;  // it's a vanilla
    }

    /** Returns all strikes the imnt is sensitive to  */
    DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                      const IModel*      model){
        static const string method("VanillaSmoothStrike::getSensitiveStrikes");
        try {
            DoubleArraySP sensStrikes(new DoubleArray(0));
            if (avoidVegaMatrix(model)) {
                throw ModelException(method, 
                                     "VEGA_MATRIX is not valid for this "
                                     "instrument");
            }

            if (maturity.isLess(valueDate) || isExercised) {
                // if it's dead do nothing
            }
            else {
                // add in each strike
                DateTime imntStart = fwdStarting ? startDate : valueDate;
                for (int i = 0; i < strikes.size(); i++) {
                    LinearStrikeVolRequestSP volRequest(new LinearStrikeVolRequest(strikes[i], 
                                                                                   imntStart, 
                                                                                   maturity,
                                                                                   fwdStarting));

                    SensitiveStrikeDescriptor sensStrikeDesc;
                    sensStrikeDesc.forwardOnly = false;
                    asset->getSensitiveStrikes(volRequest.get(), outputName, 
                                               sensStrikeDesc, sensStrikes);
                }
            }

            return sensStrikes;
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        clazz->setDescription("vanilla smooth strike");
        REGISTER(VanillaSmoothStrike, clazz);
        SUPERCLASS(Generic1Factor);
        IMPLEMENTS(PassThroughModel::IIntoProduct);
        IMPLEMENTS(LastSensDate);
        EMPTY_SHELL_METHOD(defaultVanillaSmoothStrike);
        FIELD(isCall, "isCall");
        FIELD(american, "american");
        FIELD(maturity, "maturity");
        FIELD(strike, "strike");
        FIELD(weights, "weights");
        FIELD(strikes, "strikes");
        FIELD(isExercised, "isExercised");
        FIELD_MAKE_OPTIONAL(isExercised);
        FIELD(finalLevel, "spot at maturity");
        FIELD_MAKE_OPTIONAL(finalLevel);
    }

    static IObject* defaultVanillaSmoothStrike(){
        return new VanillaSmoothStrike();
    }
    
private:
    friend class VanillaSmoothStrikeClosedForm;
    friend class VanillaSmoothStrikePassThrough;

    VanillaSmoothStrike():Generic1Factor(TYPE) {}; 
    VanillaSmoothStrike(const VanillaSmoothStrike& rhs);
    VanillaSmoothStrike& operator=(const VanillaSmoothStrike& rhs);

    void requests(Control* control, CResults* results) const {
        static const string method = "VanillaSmoothStrike::requests";
        try {
            DateTime settles = instSettle->settles(maturity, asset.get());

            OutputRequest* request =
                control->requestsOutput(OutputRequest::PAYMENT_DATES);
            if (request) {
                DateTimeArray dates(1, settles);
                OutputRequestUtil::recordPaymentDates(control,results,&dates); 
            }
        
            request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
            if (request && (valueDate.isGreaterOrEqual(maturity) || isExercised)) {
                CashFlowArray payments;
                double value = GetIntrinsic(finalLevel, strike, isCall, true);

                payments.push_back(CashFlow(settles, value));
                
                OutputRequestUtil::recordKnownCashflows(control,
                                                        results,
                                                        discount->getCcy(),
                                                        &payments); 
            }

            InstrumentUtil::recordFwdAtMat(control,
                                           results,
                                           maturity,
                                           valueDate,
                                           asset.get());
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    // pricing is done by the instrument as it's a wrapper around the vanilla
    // model throws control of pricing to here from the 'product'
    void price(PassThroughModel* model, Control* control, Results* results) const{
        static const string method = "VanillaSmoothStrike::price";        
        try  {
            if (priceDeadInstrument(control, results)) {
                return; // all done;
            }

            double value = 0.0;
            for (int i = 0; i < weights.size(); i++) {
                // build up each mini-vanilla in turn
                DateTimeArray exerDate(1, maturity);
                DoubleArray   exerValue(1, strikes[i]);
                
                Schedule exer(exerDate, exerValue, Schedule::INTERP_LINEAR);

                CVanillaSP miniv(CVanilla::make(valueDate,
                                                isCall,
                                                american,
                                                &exer,
                                                fwdStarting,
                                                startDate,
                                                oneContract,
                                                notional,
                                                initialSpot,
                                                asset.get(),
                                                ccyTreatment,
                                                discount.get(),
                                                instSettle.get()));

                Results miniRes;
                model->calculate(miniv.get(), control, &miniRes);
                double fv = miniRes.retrievePrice();

                value += weights[i] * fv;
            }

            results->storePrice(value, discount->getCcy());    
            requests(control, results);
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    DateTime endDate(const Sensitivity* sensControl) const {
        DateTime instEnd  = instSettle->settles(maturity, asset.get());
        DateTime assetEnd = asset->settleDate(maturity);
        DateTime end = instEnd.isGreater(assetEnd) ? instEnd : assetEnd;
        return end;       
    }

    virtual bool priceDeadInstrument(CControl* control, CResults* results) const {       
        if (maturity.isGreater(valueDate) && !isExercised) {
            return false;
        }

        double   value;
        DateTime settles = instSettle->settles(maturity, asset.get());

        if (valueDate >= settles) {
            // settled already
            value = 0.0;
        }
        else if (isExercised || valueDate >= maturity) {
            value = GetIntrinsic(finalLevel, strike, isCall, true);
            value *= discount->pv(valueDate, settles);
            value *= InstrumentUtil::scalePremium(oneContract,
                                                  false,
                                                  notional,
                                                  0.0,
                                                  initialSpot);
        } 
        else {
            value = 0.0;
        }            

        requests(control, results);
        results->storePrice(value, discount->getCcy());
        return true;   
    }   

    /** Implementation of PassThroughModel::IntoProduct interface */
    PassThroughModel::IProduct* createProduct(
        PassThroughModel* model) const;

private:
    bool        isCall;
    bool        american;
    DateTime    maturity;
    double      strike;  // the real one
    DoubleArray weights; // weight of each mini-vanilla
    DoubleArray strikes; // strike of each mini-vanilla

    bool        isExercised;
    double      finalLevel;
};

CClassConstSP const VanillaSmoothStrike::TYPE = CClass::registerClassLoadMethod(
    "VanillaSmoothStrike", typeid(VanillaSmoothStrike), VanillaSmoothStrike::load);


/** private class */
class VanillaSmoothStrikePassThrough: public PassThroughModel::IProduct{
private:
    const VanillaSmoothStrike* vss; // a reference

public:
    VanillaSmoothStrikePassThrough(const VanillaSmoothStrike* vss): vss(vss){}

    void price(PassThroughModel* model,
               Control*          control, 
               Results*          results)const{
        vss->price(model, control, results);
    }
};

/** Implementation of PassThroughModel::IntoProduct interface */
PassThroughModel::IProduct* VanillaSmoothStrike::createProduct(
    PassThroughModel* model) const {
    return new VanillaSmoothStrikePassThrough(this);
}

// for class loading 
bool VanillaSmoothStrikeLoad() {
    return (VanillaSmoothStrike::TYPE != 0);
}

DRLIB_END_NAMESPACE
