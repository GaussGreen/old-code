//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : TaxWrapper.cpp
//
//   Description : Linear system accounting for tax liability of an option
//                 at various MTM dates.
//                 Started from Andrew's TrailFee code rev 1.5
//
//   Date        : April 2003
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/GenericNAssetAndImnt.hpp" // to get into Pyramid - no actual assets needed
#include "edginc/ClosedFormLN.hpp"
#include "edginc/Maths.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/SensitiveStrikes.hpp"
#include "edginc/OutputRequestUtil.hpp"
#include "edginc/ITaxableInst.hpp"  // to extract data from underlying instrument
#include <algorithm>

DRLIB_BEGIN_NAMESPACE

class TaxWrapper: public GenericNAssetAndImnt,
                  virtual public CClosedFormLN::IIntoProduct,
                  virtual public ISensitiveStrikes,
                  virtual public LastSensDate,
                  virtual public Theta::Shift {
public:
    static CClassConstSP const TYPE; 

    virtual void validatePop2Object() {
        static const string routine = "TaxWrapper::validatePop2Object";
        // apply any explicit ordering info (required to avoid IMS limitation)
        if (!taxPaymentDatesOrder.empty()) {
            if (taxPaymentDatesOrder.size() != taxPaymentDates.size()) {
                throw ModelException(routine,
                                     "When provided, tax payment dates order array len=" +
                                     Format::toString(taxPaymentDatesOrder.size()) + 
                                     " must have length equal to the number of tax payment dates #=" +
                                     Format::toString(taxPaymentDates.size()));
            }
            // We require each number mentioned once and start from 1. 
            // This will (implicitly) catch out-of-bounds too
            IntArray checkOrder = taxPaymentDatesOrder;
            sort(checkOrder.begin(), checkOrder.end());
            for(int i=0; i<checkOrder.size(); i++) {
                if (checkOrder[i] != i+1) {
                    throw ModelException(routine,
                                         "taxPaymentDatesOrder has a value " +
                                         Format::toString(checkOrder[i]) + 
                                         " when expected " + Format::toString(i+1));
                }
            }
            orderedTaxPaymentDates = DateTimeArray(taxPaymentDates.size());
            for(int j=0; j<taxPaymentDates.size(); j++) {
                orderedTaxPaymentDates[taxPaymentDatesOrder[j]-1] = taxPaymentDates[j];
            }
        } else {
            orderedTaxPaymentDates = taxPaymentDates;
        }
    }

    virtual void Validate() {
        static const string method = "TaxWrapper::Validate";
        try {
            int i;
            if (taxMTMs.size()<1) {
                throw ModelException(method, "Require at least one tax mtm date!");                                                     
            }
            // 'n' tax MTM dates & rates
            // 'n' tax payment dates
            if (taxMTMs.size() != orderedTaxPaymentDates.size()) {
                throw ModelException(method, "Require equal numbers of tax payment and mtm dates but given " +
                                     Format::toString(orderedTaxPaymentDates.size()) + " tax payments and " +
                                     Format::toString(taxMTMs.size()) + " MTMs");
            }
            if (mtmRefDate>getValueDate()) {
                throw ModelException(method, 
                                     "MTM Reference date " + 
                                     mtmRefDate.toString() +
                                     " must not be after today!" + 
                                     getValueDate().toString());
            }
            if (mtmRefDate>=taxMTMs[0].date) {
                throw ModelException(method, 
                                     "MTM Reference date " + 
                                     mtmRefDate.toString() +
                                     " must be before first tax MTM date " + 
                                     taxMTMs[0].date.toString());
            }

            /* May have "deferred" case where early MTM dates produce tax
               liability only at deal maturity, so later payment date than
               later MTM dates (i.e. overlapping), so do NOT validate :-
            DateTime::ensureIncreasing(taxPaymentDates,
                                       "Tax payment dates",
                                       true);
            **/
            CashFlow::ensureDatesIncreasing(taxMTMs,
                                            "Tax MTM dates",
                                            true/*fail if empty*/);

            for(i=0; i<taxMTMs.size(); i++) {
                if (orderedTaxPaymentDates[i]<taxMTMs[i].date) {
                    throw ModelException(method, 
                                         "Tax payment date [" + Format::toString(i+1) +
                                         "] " + orderedTaxPaymentDates[i].toString() +
                                         " must not be before corresponding MTM date " + 
                                         taxMTMs[i].date.toString());
                }
            }
            for(i=0; i<pastMTMValues.size(); i++) {
                if (i>0 &&
                    pastMTMValues[i].date<=pastMTMValues[i-1].date) {
                    throw ModelException(method, 
                                         "Past MTM Value dates must be strictly increasing but [" 
                                         + Format::toString(i+1) + "] " + 
                                         pastMTMValues[i].date.toString() +
                                         " <= [" + Format::toString(i) + "] " + 
                                         pastMTMValues[i-1].date.toString());
                }
                // While a 0 MTM value is actually legitimate, disallowing 0 makes it possible
                // to check for missing past values which is more important. This requires
                // that a genuinely past MTM of 0 has to be entered as slightly non-0.
                if (pastMTMValues[i].date <= getValueDate() &&
                    Maths::isZero(pastMTMValues[i].amount)) {
                    throw ModelException(method, 
                                         "The past MTM Value for date " +
                                         pastMTMValues[i].date.toString() + 
                    " is zero, indicating a missing level. Please provide a non-zero value.");
                }
            }
        }    
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** indicates whether VEGA_MATRIX is sensible for this instrument */
    virtual bool avoidVegaMatrix(const IModel* model) {
        ISensitiveStrikes* vegaMatrixImnt = 
            dynamic_cast<ISensitiveStrikes*>(imnt.get());

        if (vegaMatrixImnt) {
            return vegaMatrixImnt->avoidVegaMatrix(useOverrideModel?
                                                   model:imntModel.get());
        }
        else {
            return true;  // can't do it
        }
    }
        
    /** returns all strikes on the vol surface to which 
        this instrument is sensitive */
    virtual DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                              const IModel*      model) {
        if (avoidVegaMatrix(model)) {
            throw ModelException("TaxWrapper::getSensitiveStrikes", 
                                 "VEGA_MATRIX is not valid for this instrument");
        }

        ISensitiveStrikes& strikes = dynamic_cast<ISensitiveStrikes&>(*imnt);

        return strikes.getSensitiveStrikes(outputName, 
                                           useOverrideModel?
                                           model:imntModel.get());
    }

    /** roll through time (setting historic values) */
    virtual bool sensShift(Theta* theta) {
        try{
            DateTime preValueDate = valueDate;
            valueDate = theta->rollDate(valueDate);
            // now do something with pastMTMValues
            bool pastRollDate = false;
            for (int i = 0; !pastRollDate && i < pastMTMValues.size(); i++){
                if ((pastMTMValues[i].date.isGreater(preValueDate) && 
                     !pastMTMValues[i].date.isGreater(valueDate)) ||
                    (pastMTMValues[i].date.equals(valueDate)    
                     && Maths::isZero(pastMTMValues[i].amount))) {
                    pastMTMValues[i].amount = mtmForRoll;
                } else if (pastMTMValues[i].date.isGreater(valueDate) ) {
                    pastRollDate = true;
                }
            }
        } catch (exception& e){
            throw ModelException(e, "TaxWrapper::sensShift(Theta) failed.");
        }
        return true; // continue to tweak components which implement Theta
    }


    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(TaxWrapper, clazz);
        SUPERCLASS(GenericNAssetAndImnt);
        IMPLEMENTS(CClosedFormLN::IIntoProduct);
        IMPLEMENTS(ISensitiveStrikes);
        IMPLEMENTS(LastSensDate);
        EMPTY_SHELL_METHOD(defaultTaxWrapper);
        FIELD(mtmRefDate, "mtmRefDate");
        FIELD(taxMTMs, "taxMTMs");
        FIELD(taxPaymentDates, "taxPaymentDates");
        FIELD(taxPaymentDatesOrder, "Order of taxPaymentDates");
        FIELD_MAKE_OPTIONAL(taxPaymentDatesOrder);
        FIELD_NO_DESC(orderedTaxPaymentDates);
        FIELD_MAKE_TRANSIENT(orderedTaxPaymentDates);
        FIELD(pastMTMValues, "pastMTMValues");
        FIELD(useOverrideModel, "else price using instruments true model");
        FIELD(debugOn,          "true=>produce debug output");
        FIELD_MAKE_OPTIONAL(debugOn);
        FIELD(mtmForRoll, "mtmForRoll");
        FIELD_MAKE_TRANSIENT(mtmForRoll);
    }

    static IObject* defaultTaxWrapper(){
        return new TaxWrapper();
    }
    
private:
    friend class TaxWrapperClosedForm; 

    TaxWrapper():GenericNAssetAndImnt(TYPE), 
        useOverrideModel(false), 
        debugOn(false), 
        taxPaymentDatesOrder(0), orderedTaxPaymentDates(0),
        mtmForRoll(0.0) {}; 
    TaxWrapper(const TaxWrapper& rhs);
    TaxWrapper& operator=(const TaxWrapper& rhs);

    DateTime endDate(const Sensitivity* sensitivity) const {
        DateTime      lastDate;
        LastSensDate* lsd = dynamic_cast<LastSensDate*>(imnt.get());
        if (lsd) {
            SensitivityConstSP shift(sensitivity);
            // if we're not using the TaxWrapper's model , we need to use the 
            // one we're pricing with to get the last sens date
            if (!useOverrideModel) {
                shift = SensitivityConstSP(sensitivity->spawn(imntModel.get()));
            }

            lastDate = lsd->endDate(shift.get());
        } else {
            throw ModelException("TaxWrapper::endDate",
                                 "underlying imnt (" + 
                                 imnt->getClass()->getName() + 
                                 ") does not implement LastSensDate");
        }
        // order doesn't matter here
        for (int i=0; i<taxPaymentDates.size(); i++) {
            if (lastDate<taxPaymentDates[i]) {
                lastDate = taxPaymentDates[i];
            }
        }
        return lastDate;
    }

    /** Implementation of ClosedFormLN::IntoProduct interface */
    CClosedFormLN::IProduct* createProduct(CClosedFormLN* model) const;

private:
    DateTime          mtmRefDate;
    CashFlowArray     taxMTMs;         //.date is MTM date, .amount is tax rate
    DateTimeArray     taxPaymentDates;
    CashFlowArray     pastMTMValues;   // can I use IPastValues?
    bool              useOverrideModel;
    bool              debugOn;

    // optional
    IntArray          taxPaymentDatesOrder; // can determine order of taxPaymentDates
    // transient
    DateTimeArray     orderedTaxPaymentDates; // after possible application of ...Order

    // for possible use with theta rolling
    mutable double    mtmForRoll;
};

CClassConstSP const TaxWrapper::TYPE = CClass::registerClassLoadMethod(
    "TaxWrapper", typeid(TaxWrapper), TaxWrapper::load);


/** private class */
class TaxWrapperClosedForm: public CClosedFormLN::IProduct{
private:
    const TaxWrapper* inst; // a reference
    // tax MTM dates customised to this instrument instance
    DateTimeArray     taxMTMDates;
    DateTimeArray     taxPaymentDates;  // after any ordering done in 'inst'

    // ("name") as per doc
    int                numTaxPeriods;        // convenient
    DoubleMatrix       In;                   // unit matrix nxn
    DoubleMatrix       In1;                  // unit matrix (n+1)x(n+1)
    DoubleMatrix       taxRates;             // ("X") used to give tax on diff in mtm values each interval
    DoubleMatrix       taxValues;            // ("V") used to give tax values of payments made each interval
    DoubleMatrix       fvToMTMDates;         // ("Y") fv factors from today to mtm dates
    DoubleMatrix       pvPayToMTMDates;      // ("U") pv factors from pay to mtm dates
    mutable DoubleArray optPrices;           // ("B") MTM values (past then future)
    int                iFirstFutureMTM;
    double             pvFromFirstFutureMTM; // may be final und opt pay date if that is next

    // XXX would be nice to be able to do "a+b"
    static DoubleArray sum(const DoubleArray& a, const DoubleArray& b) {
        if (a.size() != b.size()) {
            throw ModelException("DoubleArray sum",
                                 "Arrays not of equal lengths( " + Format::toString(a.size()) +
                                 " and " + Format::toString(b.size()) + ")" );
        }
        DoubleArray s = a;
        for(int i=0; i<a.size(); i++) {
            s[i] += b[i];
        }
        return s;
    }

public:
    TaxWrapperClosedForm(const TaxWrapper* inst): 
        inst(inst),
        taxPaymentDates(inst->orderedTaxPaymentDates),
        numTaxPeriods(taxPaymentDates.size()) {
        static const string method = "TaxWrapperClosedForm::TaxWrapperClosedForm";
        int i, j;

        // We need to know the final settlement date of the underlying option 
        const ITaxableInst::Basic* instBasic = dynamic_cast<ITaxableInst::Basic*>(inst->imnt.get());
        if (!instBasic) {
            throw ModelException(method, 
                                 "Underlying instrument of type (" + inst->imnt->getClass()->getName() +
                                 ") does not support tax treatment!");
        }
        const DateTime& undOptMat = instBasic->getFinalPaymentDate();

        // Can't cope with pricing after underlying option final settlement date
        // Because ... the valuation of tax MTMs in final period is done at undOptMat.
        if (inst->valueDate > undOptMat) {
            throw ModelException(method, 
                                 "Cannot price after underlying instrument final pay date " + undOptMat.toString());
        }

        // For greater convenience in calculations, since can treat mtmRefDate & taxMTMs uniformly
        taxMTMDates = DateTimeArray(inst->taxMTMs.size()+1);
        taxMTMDates[0] = inst->mtmRefDate;
        for(i=0; i<inst->taxMTMs.size(); i++) {
            taxMTMDates[i+1] = inst->taxMTMs[i].date;
        }
        // While we're here we may as well do the check that
        // und opt mat falls inside last MTM period.
        int iLastMTM = taxMTMDates.size()-1;
        if (undOptMat>taxMTMDates[iLastMTM] ||
            undOptMat<=taxMTMDates[iLastMTM-1]) {
            throw ModelException(method, 
                                 "Underlying instrument final pay date " + undOptMat.toString() +
                                 " must fall in final MTM period (" + 
                                 taxMTMDates[iLastMTM-1].toString() +
                                 " , " + taxMTMDates[iLastMTM].toString() + "]");
        }

        /** Note alloc is DoubleMatrix(numCols, numRows), 
            and to reference row i, col j use M[j][i] */
        In = DoubleMatrix(numTaxPeriods,numTaxPeriods); // all 0.0s, except for ...
        for(i=0;i<numTaxPeriods;i++) {
            In[i][i] = 1.0;
        }
        In1 = DoubleMatrix(numTaxPeriods+1,numTaxPeriods+1); // all 0.0s, except for ...
        for(i=0;i<numTaxPeriods+1;i++) {
            In1[i][i] = 1.0;
        }
        taxRates = DoubleMatrix(numTaxPeriods+1,numTaxPeriods); // all 0.0s, except for ...
        for(i=0;i<numTaxPeriods;i++) {
            taxRates[i][i] = -inst->taxMTMs[i].amount;
            taxRates[i+1][i] = inst->taxMTMs[i].amount;
        }
        taxValues = DoubleMatrix(numTaxPeriods,numTaxPeriods); // all 0.0s, except for ...
        for(i=0;i<numTaxPeriods;i++) {
            for(j=0;j<numTaxPeriods;j++) {
                if (taxPaymentDates[j] < taxMTMDates[i+1] &&
                    taxPaymentDates[j] >=  taxMTMDates[i]) {
                    taxValues[j][i] = inst->taxMTMs[i].amount;
                }
            }
        }
        fvToMTMDates = DoubleMatrix(numTaxPeriods+1, numTaxPeriods+1);  // no pv unless ...
        for(i=0;i<numTaxPeriods+1;i++) {
            for(j=0;j<numTaxPeriods+1;j++) {
                if (i==j) {
                    const DateTime* t = &taxMTMDates[i];
                    if (undOptMat < *t) {
                        t = &undOptMat;
                    }
                    if (*t>inst->valueDate) {
                        fvToMTMDates[j][i] = 1.0 / inst->discount->pv(inst->valueDate, *t);
                    } else {
                        fvToMTMDates[j][i] = 1.0;
                    }
                } else {
                    fvToMTMDates[j][i] = 0.0;
                }
            }
        }

        pvPayToMTMDates = DoubleMatrix(numTaxPeriods,numTaxPeriods+1); // all 0.0s, except for ...
        for(i=0;i<numTaxPeriods+1;i++) {
            for(j=0;j<numTaxPeriods;j++) {
                // For the "final" tax year we actually need to stop at undOptMat
                const DateTime* t = &taxMTMDates[i];
                if (undOptMat < *t) {
                    t = &undOptMat;
                }
                if (*t > inst->valueDate && 
                    taxPaymentDates[j] >= *t) {
                    pvPayToMTMDates[j][i] = inst->discount->pv(*t, 
                                                 taxPaymentDates[j]);
                }
            }
        }

        iFirstFutureMTM = numTaxPeriods+1;
        optPrices = DoubleArray(numTaxPeriods+1,0.0);
        for(i=0, j=0; i<numTaxPeriods+1; i++) {
            // For the "final" tax year we actually need to stop at undOptMat
            const DateTime* t = &taxMTMDates[i];
            if (undOptMat < *t) {
                t = &undOptMat;
            }
            if (*t > inst->valueDate) {
                // undOptPrice done later
                // capture the first only
                if (iFirstFutureMTM>numTaxPeriods) {
                    iFirstFutureMTM = i;
                    pvFromFirstFutureMTM = inst->discount->pv(*t);
                }
            } else {
                while(j<inst->pastMTMValues.size() && 
                      inst->pastMTMValues[j].date < *t) {
                    j++;
                }
                if (j>=inst->pastMTMValues.size() ||
                    inst->pastMTMValues[j].date>*t) {
                    throw ModelException(method,
                                         "Past MTM value missing for " + t->toString());
                }
                optPrices[i] = inst->pastMTMValues[j].amount;
            }
        }
        if (iFirstFutureMTM>numTaxPeriods) {
            throw ModelException(method,
                                 "All MTM dates are past!");
        }


    }

    void price(CClosedFormLN* model,
               Control*       control, 
               CResults*      results)const{
        static const string method = "TaxWrapperClosedForm::price";

        try  {
            // Possibly know more about the underlying instrument - namely coupons it may be
            // paying, which allows more complete tax treatment...
            // Note this is done within the pricing (rather than initialisation) so that 
            // if the TaxableInst needs to react to tweaks we get correct sensitivities.
            const ITaxableInst::WithCoupons*  instWithCoupons  = dynamic_cast<ITaxableInst::WithCoupons*>(inst->imnt.get());
            CashFlowArrayConstSP coupons; // any discrete payments made by underlying instrument (may be empty)
            DoubleArray          taxOnCoupons;         // ("TC") tax paid on coupons in each MTM period
            DoubleArray          sumFutureCoupons;     // ("FC") sum of coupons paid after each MTM date
            if (!instWithCoupons) {
                // only know FV so is correct only for non-coupon paying instruments
                // We supply an empty set of coupons ...
                coupons = CashFlowArrayConstSP(new CashFlowArray(0));
            } else {
                coupons = instWithCoupons->getCoupons();
                // Check none fall after "maturity"...
                if (coupons->back().date > taxMTMDates.back()) {
                    throw ModelException(method,
                                         "Coupon (" + coupons->back().date.toString() +
                                         ") cannot be after final tax MTM date (" +
                                         taxMTMDates.back().toString() + ")");
                }
            }
            taxOnCoupons = DoubleArray(numTaxPeriods, 0.0);
            int i,j;
            for(i=0, j=0; i<numTaxPeriods; i++) {
                taxOnCoupons[i] = 0.0;
                // A single pass through for taxOnCoupons
                while (j<coupons->size() &&
                       ((*coupons)[j].date >= taxMTMDates[i]) &&
                       ((*coupons)[j].date < taxMTMDates[i+1])) {
                    taxOnCoupons[i] += (*coupons)[j].amount;
                    j++;
                }
                taxOnCoupons[i] *= inst->taxMTMs[i].amount;
            }
            sumFutureCoupons = DoubleArray(numTaxPeriods+1, 0.0);
            // Cumulative sum means working backwards is most efficient
            // Validation above on coupons means we know 
            // sumFutureCoupons[numTaxPeriods]=0
            for(i=numTaxPeriods-1, j=coupons->size()-1; i>=iFirstFutureMTM; i--) {
                // taxMTMDates are in increasing order, so can start the next
                // one off where the later one finished.
                // Also, iFirstFutureMTM >= 1 (by definition/validation)
                sumFutureCoupons[i] = sumFutureCoupons[i+1] * 
                        inst->discount->pv(taxMTMDates[i], taxMTMDates[i+1]);
                while (j>=0 &&
                       (*coupons)[j].date > taxMTMDates[i]) {
                    sumFutureCoupons[i] += (*coupons)[j].amount * 
                        inst->discount->pv(taxMTMDates[i], (*coupons)[j].date);
                    j--;
                }
            }

            // useful
            DoubleMatrix V1 = (In - taxValues).computeInverse();
            DoubleMatrix G1 = pvPayToMTMDates * V1;
            // Ref doc eq# 7
            DoubleMatrix G = (In1 - G1 * taxRates).computeInverse();

            // At the moment the formulation treats 2 distinct cases - either
            // no coupons (and have a fair value) or all coupons (so 0 fair value)
            // The natural full solution is to know the fair value at future dates
            // (which is what we're simulating with "option+coupons") but we're a way off that.
            // Note that by making the condition here base off the instrument type implementing
            // ITaxableInst we are stopping the ability to support an instrument which can be
            // treated as "option+coupon". 


            // Need the underlying option price to complete 'optPrices' and possibly for theta
            // and since it's (currently) cheap for "coupons" cases we don't mind ....
            CControlSP ctrl(Control::makeFromFlags("", 0.0));
            IModel*    myModel = inst->useOverrideModel ? model : inst->imntModel.get();
            
            ResultsSP baseResults(myModel->Run(inst->imnt.get(), ctrl.get()));
            double undOptPrice = baseResults->retrievePrice();
            if (!instWithCoupons) {
                for(i=iFirstFutureMTM; i<optPrices.size(); i++) {
                    optPrices[i] = undOptPrice;
                }
            } // else 0.0s

            // Option mtm values including tax. (ref doc eq# 10)
            DoubleArray M = G.mult(sum(sum(fvToMTMDates.mult(optPrices), sumFutureCoupons), G1.mult(taxOnCoupons)));

            // Actual tax payments (ref doc eq# 8.1) 
            DoubleArray P = V1.mult(sum(taxRates.mult(M), taxOnCoupons));

            // The price of tax is simply the value of the future tax payments
            double payoff = 0.0;
            for(int p=0;p<taxPaymentDates.size(); p++) {
                if (taxPaymentDates[p]>inst->valueDate) {
                    payoff += P[p] * 
                        inst->discount->pv(taxPaymentDates[p]);
                }
            }
            results->storePrice(payoff, 
                                inst->discount->getCcy());   

            if (control->isPricing()){

                try {
                    // PAYMENT_DATES is a list of all dates on which payments may occur
                    // including past and potential future dates.
                    OutputRequest* request =
                        control->requestsOutput(OutputRequest::PAYMENT_DATES);
                    if (request && !request->getHasFinished()) {
                        OutputRequestUtil::recordPaymentDates(control,results,
                                                              &taxPaymentDates);
                    }

                    CashFlowArray knownCashFlows(0);
                    
                    // the index is 1-based since we need consider here only
                    // mtm dates which have corresponding payment dates, and the
                    // first mtm date is a reference only. Note there is validation
                    // to ensure mtmRefDate is past.
                    for(int i=1; i<iFirstFutureMTM; i++) {
                        CashFlow cf(taxPaymentDates[i-1], P[i-1]);
                        knownCashFlows.push_back(cf);
                    }

                    // KNOWN_CASHFLOWS should have dates a subset of PAYMENT_DATES
                    // and be supplied for all past cash flows, and any future ones
                    // that are determined.
                    // Here this means flows for each known tax payment.
                    request =control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
                    if (request && !request->getHasFinished() && 
                        knownCashFlows.size()>0) {
                        OutputRequestUtil::recordKnownCashflows(control,
                                                                results,
                                                                inst->discount->getCcy(),
                                                                &knownCashFlows); 
                    }

                    // record price and und opt price for possible use for theta  
                    // when rolling across MTM date
                    inst->mtmForRoll = undOptPrice + payoff;
                }
                catch (exception& e) {
                    throw ModelException(e, method);
                }
                
                {
                    // for information only
                    DoubleArray Apast = optPrices;
                    DoubleArray Afuture = optPrices;
                    for(int i=0; i<optPrices.size(); i++) {
                        if (i<iFirstFutureMTM) {
                            Afuture[i] = 0.0;
                        } else {
                            Apast[i] = 0.0;
                        }
                    }

                    OutputRequest* request = control->requestsOutput(OutputRequest::TAX_GROSS_UP);
                    if (request) {
                        if (Maths::isZero(undOptPrice)) {
                            results->storeRequestResult(request, 
                                                        IObjectSP(new NotApplicable()));
                        } else {
                            DoubleArray forGrossUp = G.mult(fvToMTMDates.mult(Afuture));
                            results->storeRequestResult(request, 
                                                        CDoubleSP(CDouble::create(forGrossUp[iFirstFutureMTM] / undOptPrice * 
                                                                                  pvFromFirstFutureMTM)));
                        }
                    }
                    request = control->requestsOutput(OutputRequest::TAX_NEXT_CASHFLOW);
                    if (request) {
                        if (Maths::isZero(undOptPrice)) {
                            results->storeRequestResult(request, 
                                                        IObjectSP(new NotApplicable()));
                        } else {
                            DoubleArray forLoan = G.mult(Apast);
                            results->storeRequestResult(request, 
                                                        CDoubleSP(CDouble::create(forLoan[iFirstFutureMTM] / undOptPrice *
                                                                                  pvFromFirstFutureMTM)));
                        }
                    }
                }
                if (inst->debugOn) {
                    OutputNameConstSP Name_M(new OutputName("TAX_M"));
                    results->storeCcylessResult(IObjectSP(new DoubleArray(M)),
                                                Results::DEBUG_PACKET,
                                                Name_M);
                    OutputNameConstSP Name_P(new OutputName("TAX_P"));
                    results->storeCcylessResult(IObjectSP(new DoubleArray(P)),
                                                Results::DEBUG_PACKET,
                                                Name_P);
                    OutputNameConstSP Name_G(new OutputName("TAX_G"));
                    results->storeCcylessResult(IObjectSP(new DoubleMatrix(G)),
                                                Results::DEBUG_PACKET,
                                                Name_G);
                    OutputNameConstSP Name_X(new OutputName("TAX_X"));
                    results->storeCcylessResult(IObjectSP(new DoubleMatrix(taxRates)),
                                                Results::DEBUG_PACKET,
                                                Name_X);
                    OutputNameConstSP Name_ALPHA(new OutputName("TAX_ALPHA"));
                    results->storeCcylessResult(IObjectSP(CDouble::create(undOptPrice)),
                                                Results::DEBUG_PACKET,
                                                Name_ALPHA);
                    OutputNameConstSP Name_V(new OutputName("TAX_V"));
                    results->storeCcylessResult(IObjectSP(new DoubleMatrix(taxValues)),
                                                Results::DEBUG_PACKET,
                                                Name_V);
                    OutputNameConstSP Name_Z(new OutputName("TAX_Z"));
                    results->storeCcylessResult(IObjectSP(new DoubleMatrix(pvPayToMTMDates)),
                                                Results::DEBUG_PACKET,
                                                Name_Z);
                    OutputNameConstSP Name_Y(new OutputName("TAX_Y"));
                    results->storeCcylessResult(IObjectSP(new DoubleMatrix(fvToMTMDates)),
                                                Results::DEBUG_PACKET,
                                                Name_Y);
                    OutputNameConstSP Name_FC(new OutputName("TAX_FC"));
                    results->storeCcylessResult(IObjectSP(new DoubleArray(sumFutureCoupons)),
                                                Results::DEBUG_PACKET,
                                                Name_FC);
                    OutputNameConstSP Name_TC(new OutputName("TAX_TC"));
                    results->storeCcylessResult(IObjectSP(new DoubleArray(taxOnCoupons)),
                                                Results::DEBUG_PACKET,
                                                Name_TC);
                }
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }
};

/** Implementation of ClosedForm::IntoProduct interface */
CClosedFormLN::IProduct* TaxWrapper::createProduct(CClosedFormLN* model) const
{
    return new TaxWrapperClosedForm(this);
}

// for class loading 
bool TaxWrapperLoad() {
    return (TaxWrapper::TYPE != 0);
}

DRLIB_END_NAMESPACE
