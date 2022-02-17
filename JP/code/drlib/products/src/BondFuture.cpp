//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : BondFuture.cpp
//
//   Description : Bond future instrument
//
//   Author      : Ian Stares
//
//   Date        : 25 August 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ClosedForm.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/Bond.hpp"
#include "edginc/InstrumentSettlement.hpp"
#include "edginc/RootFinder.hpp"
#include "edginc/RhoParallel.hpp"
#include "edginc/FutureAsAsset.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/PropertyTweakHypothesis.hpp"
#include "edginc/Control.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/Results.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/RateParallel.hpp"

DRLIB_BEGIN_NAMESPACE

class BondFuture: public CInstrument, 
                  public virtual ClosedForm::IIntoProduct,
                  public virtual LastSensDate,
                  public virtual Theta::Shift,
                  public virtual Func1D::NoDeriv,
                  public virtual IFutureAsAsset {
public:
    static CClassConstSP const TYPE;
    
    /** copy market data relevant to the instrument */
    virtual void GetMarket(const IModel* model, const CMarketDataSP market) {
        market->GetReferenceDate(valueDate);
        discount.getData(model, market);
        ctdBond->getMarket(model, market.get());
        futuresSettle->getMarket(model, market.get());
    }

    virtual void Validate() {
        static const string method = "BondFuture::Validate";
        try {
            if ( !Maths::isPositive(ctdConversionFactor) ) {
                throw ModelException(method, "Invalid CTD conversion factor " + Format::toString(ctdConversionFactor));
            }
            if ( !ctdBond->getMaturityDate().isGreater(futuresExpiryDate) ) {
                throw ModelException(method, "Bond maturity is on or before futures expiry");
            }
            if ( !Maths::isPositive(futuresPrice) ) {
                throw ModelException(method, "Invalid futures price " + Format::toString(futuresPrice));
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** Implementation of ClosedForm::IntoProduct interface */
    virtual ClosedForm::IProduct* createProduct(ClosedForm* model) const;

    virtual DateTime getValueDate() const {
        return valueDate;
    }

    /** when to stop tweaking */
    virtual DateTime endDate(const Sensitivity* sensControl) const {
        return ctdBond->getMaturityDate();
    }

    virtual bool sensShift(Theta* shift) {
        try {
            valueDate = shift->rollDate(valueDate);
        }
        catch (exception& e) {
            throw ModelException(e, "BondFuture::sensShift (theta)");
        }    
        return true; // our components have theta type sensitivity
    }

    // to implement the interface needed for Brent root finding
    virtual double operator()(double zSpreadGuess) const{
        return priceFromZSpread(zSpreadGuess) - futuresPrice;
    }

    /** Returns the name of the instrument's discount currency. */
    virtual string discountYieldCurveName() const {
        return discount.getName();
    }

    //----------- to implement the interface IFutureAsAsset -----------------------//
  
    /** Returns a date after which the instrument can no longer be used as
        an asset, here the future settlement date  */
    virtual DateTime maturityDate() const {
        //return futuresSettle->settles(futuresExpiryDate, 0);
        return futuresExpiryDate;
    }

    /** Returns the yield curve used for discounting */
    virtual YieldCurveConstSP getDiscount() const {
        return discount.getSP();
    }

    //------------ end of IFutureAsAsset -------------------------------------------//

private:
    friend class BondFutureHelper;
    friend class BondFutureClosedForm;

 // for reflection
    BondFuture(): CInstrument(TYPE), ctdConversionFactor(0.0), futuresPrice(0.0), 
                  zSpread(0.0), haveZSpread(false) {};
    BondFuture(const BondFuture& rhs);
    BondFuture& operator=(const BondFuture& rhs);

    // deal with output requests
    void requests(Control* control, CResults* results) const {
        static const string method = "BondFuture::requests";
        try {
            OutputRequest* request =
                    control->requestsOutput(OutputRequest::BOND_FUTURE_Z_SPREAD);
            if (request) {
                results->storeRequestResult(request, IObjectSP(CDouble::create(zSpread)));
            }

            request = control->requestsOutput(OutputRequest::BOND_DURATION);
            if (request) {
                try {
                    // BOND_DURATION reported is the Macaulay duration of the CTD bond from maturity
                    // using the discount curve with the z-spread applied (so consistent with bond future price)
                    results->storeRequestResult(request, ctdBond->duration(valueDate,discount.getSP()));

                    YieldCurveSP copyCurve(copy(discount.get()));
                    PropertyTweakHypothesis<RateParallel>(zSpread, OutputName::SP(discount->getCcy())).applyTo(copyCurve);
            
                    double duration = ctdBond->duration(futuresExpiryDate, copyCurve);
                    // Scale by conversion factor
                    duration /= ctdConversionFactor;
                    results->storeRequestResult(request, duration);

                } catch (exception& e) {
                    UntweakableSP untweakableSpread(new Untweakable(e));
                    results->storeRequestResult(request,untweakableSpread);
                }
            }

            request = control->requestsOutput(OutputRequest::PVBP);
            if (request) {
                try {
                    // Difference in Qlib vs. Concorde methodology for bond future PVBP calculation:
                    // Concorde defines it to be the spot PVBP of the CTD bond * conversion factor.
                    // We use the PVBP of the forward bond price to exclude risk wrt bond cashflows payable before maturity.

                    // Get forward price of bond implied by market price of future 
                    // (assume future is quoted as %par)
                    double ctdPriceAtMat = futuresPrice * ctdConversionFactor;

                    // Scale the bond price for non-100 face value bonds
                    ctdPriceAtMat *= ctdBond->getFaceValue()/100;
                    
                    double pvbp = ctdBond->pvbp(ctdPriceAtMat, 
                                                false,      // dirty
                                                futuresExpiryDate);

                    // Scale back to %par units and by conversion factor
					pvbp = (pvbp/ctdConversionFactor) * (100/ctdBond->getFaceValue());

                    // Finally make absolute and scale by 100
                    results->storeRequestResult(request, 100*fabs(pvbp));
                } catch (exception& e) {
                    UntweakableSP untweakable(new Untweakable(e));
                    results->storeRequestResult(request,untweakable);
                }
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    // Shift the yield curve inputs by the Z-spread and calculate the futures price by pv'ing the CTD
    // and scaling by the conversion factor - remember to reset the yield curve once we're done
    double priceFromZSpread(double zSpread) const {
        YieldCurveSP copyCurve(copy(discount.get()));
        PropertyTweakHypothesis<RateParallel>(
            zSpread, OutputName::SP(discount->getCcy())).applyTo(copyCurve);

        DateTime settleDate = futuresSettle->settles(futuresExpiryDate, 0);
        double ctdDirtyFwd = ctdBond->presentValue(settleDate, copyCurve);

        // future is quoted as %par, so divide by face value (which may not be 100)
        // and scale by 100%
        ctdDirtyFwd *= 100.0/ctdBond->getFaceValue();

        return ctdDirtyFwd / ctdConversionFactor;
    }

    // calculate Z-spread from given futures price. 
    // i.e. find the spread so that priceFromZSpread matches the input futures price
    void calculateZSpread() const {
        static const double TOLERANCE = 1.0e-8; 
      
        try {
            // use Brent to solve for the Z-spread
            ZBrent solver(TOLERANCE);
            try { // first go at bracketing
                zSpread = solver.solve(*this, -0.1, 0.1);
            }
            catch (exception& ) { // expand to silly range
                zSpread = solver.solve(*this, -0.5, 0.5);
            }
            // now update our object if successful
            haveZSpread = true;
        }
        catch (exception& e) {
            haveZSpread = false;
            zSpread = 0.0;

            // just for error messaging
            YieldCurveConstSP yc(discount.get());
            double approxValue = ctdBond->presentValue(futuresExpiryDate, yc);
            double approxPct   = 100.0*approxValue/ctdBond->getFaceValue();

            throw ModelException(e, "BondFuture::calculateZSpread",
                                 "Couldn't solve for futures price of " + 
                                 Format::toString(futuresPrice) + 
                                 " for bond with value of " +                                  
                                 Format::toString(approxPct) + "%");
        }
    }


    void price(Control* control, CResults* results) const {
        static const string method = "BondFuture::price";
        try {
            double value = 0.0;
        
            if (valueDate.isLess(futuresExpiryDate)) {
                if (!haveZSpread) {
                    // calculate Z-spread for later use
                    calculateZSpread();
                }
                if (!haveZSpread) {
                // the Z-spread has not been calculated. This is a problem
                    throw ModelException(method, "Z-spread has not previously been calculated");
                }
                value = priceFromZSpread(zSpread);
            }
        
            results->storePrice(value, discount->getCcy());

            if (control && control->isPricing() ) {
                requests(control, results);
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }
    

    BondSP                  ctdBond;
    YieldCurveWrapper       discount;
    DateTime                valueDate;
    const DateTime          futuresExpiryDate;
    InstrumentSettlementSP  futuresSettle;
    double                  ctdConversionFactor;
    double                  futuresPrice;
    mutable double          zSpread;
    mutable bool            haveZSpread;
};


/** private class */
class BondFutureClosedForm: public ClosedForm::IProduct{
private:
   const BondFuture* bf; // a reference

public:
    BondFutureClosedForm(const BondFuture* bf): bf(bf){}

    void price(ClosedForm* model,
               Control*    control, 
               CResults*   results) const{
        bf->price(control, results);
    }
};
    
/** Implementation of ClosedForm::IntoProduct interface */
ClosedForm::IProduct* BondFuture::createProduct(
    ClosedForm* model) const{
    return new BondFutureClosedForm(this);
}

class BondFutureHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        clazz->setDescription("Bond future");
        REGISTER(BondFuture, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(ClosedForm::IIntoProduct);
        IMPLEMENTS(LastSensDate);
        IMPLEMENTS(Theta::Shift);
        IMPLEMENTS(IFutureAsAsset);
        EMPTY_SHELL_METHOD(defaultBondFuture);
        FIELD(ctdBond, "The cheapest to deliver bond");
        FIELD(discount, "identifies discount curve");
        FIELD(valueDate, "valuation date");
        FIELD_MAKE_OPTIONAL(valueDate);
        FIELD(futuresExpiryDate, "expiry date of futures");
        FIELD(futuresSettle, "Settlement details for futures");
        FIELD(ctdConversionFactor, "CTD conversion factor");
        FIELD(futuresPrice, "the market futures price");
        FIELD(zSpread, "the Z-spread");
        FIELD_MAKE_TRANSIENT(zSpread);
        FIELD(haveZSpread, "has Z-spread been calculated");
        FIELD_MAKE_TRANSIENT(haveZSpread);
    }

    static IObject* defaultBondFuture(){
        return new BondFuture();
    }
};

CClassConstSP const BondFuture::TYPE = CClass::registerClassLoadMethod(
    "BondFuture", typeid(BondFuture), BondFutureHelper::load);

/* for class loading */
bool BondFutureLoad() {
    return (BondFuture::TYPE != 0);
}


DRLIB_END_NAMESPACE
