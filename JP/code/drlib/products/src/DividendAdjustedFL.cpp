//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DividendAdjustedFL.cpp
//
//   Description : Dividend adjusted flow through option instrument
//
//   Date        : 09/27/2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/DividendAdjusted.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/SimpleEquity.hpp"
#include "edginc/StruckEquity.hpp"
#include "edginc/Actual365F.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/NumericalIntegrationLN.hpp"
#include "edginc/ImpliedSample.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/Black.hpp"



#pragma optimize( "g", off )

DRLIB_BEGIN_NAMESPACE

///////// local helpers /////////////////////
/* calculate expected delta for fwdDate and maturty matDate */
static double fwdDelta( bool    isCall, 
                        double  fwd,
                        double  strike,
                        double  variance,
                        double  fwdDateOverMatDate)
{
    if (Maths::isNegative(strike))
    {    // average option already guaranteed to be in the money
        return 1.0;
    }
    else
    {
        double sqrtVar = sqrt(variance);
              
        double d1 = (log(fwd/strike) + (0.5 * variance))/sqrtVar;
        double d11 = d1 - sqrtVar*fwdDateOverMatDate;

        double delta = (isCall ? N1(d11) : -N1(-d11));
    
        return delta;
    }
}

/** returns true if ex date is found in which case corresponding pay date is assigned. 
            false otherwise, with pay date set to 0 */
static bool findPayDate(const DateTime& exDate, const CashFlowArray& exDateFlow, 
                        const CashFlowArray& payDateFlow, DateTime& foundPayDate)
{
    if (exDateFlow.size() != payDateFlow.size())
        throw ModelException("findPayDate", "ex date array and pay date array have different sizes");

    for (int i = 0; i < exDateFlow.size(); i++)
    {
        if (exDateFlow[i].date.getDate() == exDate.getDate())
        {    
            foundPayDate = payDateFlow[i].date;
            return true;
        }
            
    }

    foundPayDate = DateTime(0,0);
    return false;
}

///////// end local helpers /////////////////////

////////////////////// DividendAdjustedFL class //////////////////

class DividendAdjustedFL: public DividendAdjBase,
                          public CClosedFormLN::IIntoProduct
{
public:
    static CClassConstSP const TYPE;

    /** Implementation of ClosedForm::IntoProduct interface */
    virtual CClosedFormLN::IProduct* createProduct(CClosedFormLN* model) const;

    /** rolls value date and sets initial spot for Theta,
        return true if sub objects need to be tweaked */
    virtual bool sensShift(Theta* shift);

    /** instrument validation override */
    virtual void Validate();

    /** substitute equity div's by assumed div's */
    void substituteEqDiv(AssetSP eq, bool omitSpecialDivDate = false, 
                                     const DateTime& specialDivDate = DateTime(0,0)) const ;


    DoubleArraySP getPastPeriodPayments() const;

    /** return all dates which a payment may be made */
    DateTimeArraySP getPaymentDates();

    void calcFlowData( const DoubleArraySP paidSoFar, // unadjusted (same unit as assumed div)
         const DoubleArray& weightDiv,
         const DoubleArray& weightAssumed,
         CashFlowArraySP    flowThroughCashFlows) const;

    bool dateInHist(const DateTime& aDate, const CashFlowArray& cfArr) const;

private:
    friend class DividendAdjustedFLHelper;
    friend class DividendAdjustedFLClosedForm;

    DividendAdjustedFL();
    DividendAdjustedFL(const DividendAdjustedFL& rhs);
    DividendAdjustedFL& operator=(const DividendAdjustedFL& rhs);
    static void load(CClassSP& clazz);
    
protected:

    static const string NOTIONAL;
    static const string DELTA_ASSUMED;
    static const string DELTA_EQ_DIV;
    static const string ZERO_MU;

    // input in addition to those in DividendAdjBase already
    string                  flowThruType;
    CashFlowArray           histDelta;
    CashFlowArray           histFlow; // historic recorded flows before fx/delta weighting

    // transient fields to allow theta shift to use it
    mutable double          deltaToday; // this is delta for cash flow, not necessarily equal to delta of the option
    mutable double          divDiffToday; // difference in div amount
    DividendArray   eqOrigDivArr; // needed for value change in zero mu pricing call.

   DividendAdjustedFL(CClassConstSP const type);  // for derived types
};

typedef smartPtr<DividendAdjustedFL> DividendAdjustedFLSP;

const string DividendAdjustedFL::NOTIONAL = "NOTIONAL";
const string DividendAdjustedFL::DELTA_ASSUMED = "DELTA_ASSUMED";
const string DividendAdjustedFL::DELTA_EQ_DIV = "DELTA_EQ_DIV";
const string DividendAdjustedFL::ZERO_MU = "ZERO_MU";

////////////// class methods ////////////////////////////
// validation
void DividendAdjustedFL::Validate()
{
    const string method = "DividendAdjustedFL::Validate";
    DividendAdjBase::Validate();

    bool isNotional = (flowThruType == DividendAdjustedFL::NOTIONAL);
    bool isZeroMu = (flowThruType == DividendAdjustedFL::ZERO_MU);
    bool isStruck = (ccyTreatment == CAsset::CCY_TREATMENT_STRUCK);
    bool isProtected = (ccyTreatment == CAsset::CCY_TREATMENT_PROTECTED);
    bool isDelta = !(isNotional || isZeroMu);
    
    if (isDelta && (isStruck || isProtected))
    {
        throw ModelException(method, "Currency struck or protected options are currently not supported "
            "for flow through types DELTA_ASSUMED and DELTA_EQ_DIV.");
    }

    if (fwdAdj && !isZeroMu)
    {
        throw ModelException(method, "Fwd adj is only currently supported "
            "for flow through type ZERO_MU.");
    }

    // need this for vega matrix strikes
    avgInst = buildAvgSpot(avgOut, strike);

    if(!isNotional && !isZeroMu && hasYield())
        throw ModelException(method, flowThruType + " doesn't support div yields"); 

    // validation,  can further verify that dates match !
    if (!isZeroMu)
    {
        if(isStruck && histFlow.size() != histFXRates.size())
            throw ModelException(method, "num of past fx rates " + Format::toString(histFXRates.size()) +
                        " is not equal to num of known payments " + Format::toString(histFlow.size()));

        if(!isNotional && histFlow.size() != histDelta.size())
            throw ModelException(method, "num of past delta's " + Format::toString(histDelta.size()) +
                        " is not equal to num of known payments " + Format::toString(histFlow.size()));
    }

    if(assumedDivs.size() == 0)
        throw ModelException(method, "no input for assumed divs.");

    if(fwdStarting && assumedDivs[0].date < startDate)
        throw ModelException(method, "assumed div cannot be before fwd start date.");

    if(isZeroMu)
    {
        for (int i =0; i<numExDivDates.size(); i++)
        {
            if ((int)numExDivDates[i].amount != 1)
                throw ModelException(method, "for flowThruType=ZERO_MU, each div period must contain exactly one "
                    " dividend. Period starting " + numExDivDates[i].date.toString() + " contains " +
                    Format::toString((int)numExDivDates[i].amount) + " dividends.");
        }
        if(eqDivArr.size() != assumedDivs.size())
            throw ModelException(method, "for flowThruType=ZERO_MU, assumed div array must match num of divs in the period");

        // for fwd starting, cannot have divs between value date and fwd start date
        if(fwdStarting)
        { 
            // extract divs between value date and start date 
            DividendListSP divs = AssetUtil::getAllDivsBetweenDates(asset.get(),
                                                                    valueDate,
                                                                    startDate);

            if (divs->getDivAmounts()->size() > 0)
                throw ModelException(method, "ZERO_MU type cannot have dividends between value date and fwd start date.");
        }

        /* to do: remove this comment
           
        if (divadjType != "NONE")
            throw ModelException(method, "ZERO_MU type cannot have UP or DOWN for divadjType.");

        */
    }
}

/** private class */
class DividendAdjustedFLClosedForm: public CClosedFormLN::IProduct{
private:
    IModel*                    model;
    DividendAdjustedFL*        divadjFL; // a reference, this is not const

public:
    DividendAdjustedFLClosedForm(IModel* model, const DividendAdjustedFL* inst): 
        model(model) {
        divadjFL = const_cast<DividendAdjustedFL*>(inst); }

    void price(CClosedFormLN*   model,
               Control*        control, 
               CResults*       results) const;

   
    double computeDelta(DoubleArray& weightDiv,
                        DoubleArray& weightAssumed) const;

    // calc the value difference between spot + assumed and spot + eq_div across an ex-div date for ZERO_MU type
    double valueChange(double price, Control* control) const;

    // same kind of interp as for average spot, also set strike used for interp
    CVolRequestLN* volInterp(double& strikeUsed) const;

    // for zero mu, it is required to massage cf array to properly handle UP/DOWN.
    void adjustForUpDown(CashFlowArraySP &payments, const string& divAdjType) const;
};

// price call entry point
void DividendAdjustedFLClosedForm::price(CClosedFormLN*   model,
                                            Control*        control, 
                                            CResults*       results) const
{
    static const string method = "DividendAdjustedFLClosedForm::price";
    try {

        // keep pricing unless a date after mat date, for cash flow output
        if (divadjFL->valueDate.getDate() > divadjFL->matDate.getDate())
        {// if settled already, return 0. 
            if (divadjFL->valueDate >= divadjFL->instSettle->settles(divadjFL->matDate, divadjFL->asset.get()))
            {   
                results->storePrice(0.0, divadjFL->discount->getCcy());
                return;
            }
        }

        int i;
        double premium = 0.0;
        double cashToday = 0.0;
        DateTime valDate = divadjFL->valueDate;

        // get orig eq div.  This is used only for computing value change for zero mu
        if (!!control && control->isPricing())
        {
            divadjFL->eqOrigDivArr = divadjFL->getEqDivs();
        }
        
        divadjFL->eqDivArr = divadjFL->getEqDivs();

      
        // for debugging
#if 0
        {
            static i = 0;
            string dawFile = "c:\\tmp\\daw-dbgD" + Format::toString(i) + ".xml";
            XMLWriter xml(dawFile.c_str());
            divadjFL->write("OBJECT", &xml);
            string dawFile2 = "c:\\tmp\\daw-asset-dbgD" + Format::toString(i) + ".xml";
            XMLWriter xml2(dawFile2.c_str());
            divadjFL->asset.get()->write("OBJECT", &xml2);
            i++;
        }
#endif

        CashFlowArray flowAmounts; // with ex date
        CashFlowArray flowAmountsOnPayDates; // with pay date
        // only for notional case do we need callput sign (since delta already has the correct sign).
        DoubleArray    weightDiv(divadjFL->eqDivArr.size(), divadjFL->isCall?1.0:-1.0);
        DoubleArray    weightAssumed(divadjFL->assumedDivs.size(), divadjFL->isCall?1.0:-1.0);

        bool isStruck = (divadjFL->ccyTreatment == CAsset::CCY_TREATMENT_STRUCK);
        bool isZeroMu = (divadjFL->flowThruType == DividendAdjustedFL::ZERO_MU);
        bool isDeltaThru = (divadjFL->flowThruType == DividendAdjustedFL::DELTA_ASSUMED 
                    || divadjFL->flowThruType == DividendAdjustedFL::DELTA_EQ_DIV);
        // compute scaling factor 
        double scalingFactor = InstrumentUtil::scalePremium(
                                divadjFL->oneContract,
                                divadjFL->fwdStarting,
                                divadjFL->notional,
                                divadjFL->fwdStarting ? divadjFL->asset->fwdValue(divadjFL->startDate) : 1.0,
                                divadjFL->initialSpot);

        // compute hist flow, note that for zero_mu case, dates are pay dates, otherwise dates are ex-dates
        // also amount is in the same unit as assumed div except for zero_mu case, where it is cashflow
        if (isZeroMu) // treat zero_mu
        {// price using assumed divs 
            flowAmounts = divadjFL->histFlow;
            flowAmountsOnPayDates = flowAmounts;
            divadjFL->substituteEqDiv(divadjFL->asset.getSP()); // note that only amount replaced in this case, dates are not
        }
        else
        {
            // get hist flow first
            for (i=0; i<divadjFL->histFlow.size(); i++)
            {
                if (divadjFL->histFlow[i].date >= valDate)
                    break; // spreadsheet tend to have future dates filled

                flowAmounts.push_back(divadjFL->histFlow[i]);
                if (isStruck)
                    flowAmounts[i].amount *= divadjFL->histFXRates[i].amount;
                if (isDeltaThru)
                    flowAmounts[i].amount *= divadjFL->histDelta[i].amount;
            }
            // calc weights for delta types
            if(isDeltaThru)
            {// set weights to delta, keep today's delta for possible output
                double delta = computeDelta(weightDiv, weightAssumed);
                if (control && control->isPricing()) // for recording and theta use
                    divadjFL->deltaToday = delta;
            }

            CashFlowArraySP futFlow(new CashFlowArray);
            divadjFL->calcFlowData(divadjFL->getPastPeriodPayments(), weightDiv, weightAssumed, futFlow);


            // combine past and future cash flow
            flowAmounts.insert(flowAmounts.end(), futFlow->begin(), futFlow->end());

            

        }
        // multiply by scale factor, look up pay date (and store in flowAmountsOnPayDate), and calc sum pv'ed
        flowAmountsOnPayDates.clear();
        CashFlowArraySP knownPayments(new CashFlowArray(0)); // to capture known cashflows.
        for (i=0; i<flowAmounts.size(); i++)
        {
            DateTime payDate = flowAmounts[i].date; // for zero mu case, flowAmounts are on pay date.
            if (!isZeroMu) // need to reset payDate (since flowAmounts contain ex dates)
            {
                // scale factor needs to be applied to all amounts (past data as well,
                // since the units are based on one contract)
                flowAmounts[i].amount *= scalingFactor;

                // process flowAmounts with divadj type (UP, DOWN, BOTH).
                if (divadjFL->divadjType == "UP") flowAmounts[i].amount = Maths::max(flowAmounts[i].amount,0.0);
                if (divadjFL->divadjType == "DOWN") flowAmounts[i].amount = Maths::min(flowAmounts[i].amount,0.0);


                // lookup pay date: if not in list, use ex instead (for assumed)
                const DateTime& exDate = flowAmounts[i].date;
                if (! divadjFL->findDivPayDate(exDate, divadjFL->eqDivArr, payDate))
                {
                    payDate = exDate;
                }
            }
            // now add to price 
            const double& amount = flowAmounts[i].amount;
            flowAmountsOnPayDates.push_back(CashFlow(payDate, amount));

            // output all (even past) payments
            knownPayments->push_back(CashFlow(payDate, amount));

            if(payDate > valDate) 
            {
                premium += amount * 
                           divadjFL->discount->pv(valDate, payDate);
            }
        }
        // build an average spot inst and calc price
        model->Price(divadjFL->buildAvgSpot(divadjFL->avgOut, divadjFL->strike).get(), control, results);

        double avgPrice = results->retrievePrice();
        // now add price to premium
        premium += avgPrice;

        // see if there is an extra amount for zero mu today, if the instrument has not been updated yet
         double zeroMuCashToday = 0.0;

         if(isZeroMu && !divadjFL->fwdAdj && divadjFL->recordHist()) // is today an exdate? 
         {
             DateTime divPayDate;
             if (!divadjFL->findDivPayDate(valDate, divadjFL->eqDivArr, divPayDate))
             {
                 divPayDate = valDate; // should never get here; to do: for cleaner code, recordHist() could return pay date.
             }
             
             if (!divadjFL->dateInHist(divPayDate, divadjFL->histFlow))
             {      
                 zeroMuCashToday = valueChange(avgPrice, control); // imnt is not updated yet, need to compute amount 
             }
             // else it's already 0.0 
         }
         
         // and store the price plus any possible additional amount from zero mu generated cf's
         results->storePrice(premium + zeroMuCashToday, divadjFL->discount->getCcy());
         
         // record price after adding to zeroMuCashToday, to be computed here:
        

        if ( control && control->isPricing() ) {
            OutputRequest* request = NULL;

            // output ex-date cashflows. (for debug)
            OutputNameConstSP flOut(new OutputName("Divadj-CashFlows-ExDates"));
            CashFlowArraySP out = CashFlowArraySP(new CashFlowArray(flowAmounts));
            results->storeGreek(out, Results::DEBUG_PACKET, flOut);

            // output pay-date cashflows. (for debug)
            OutputNameConstSP flOut2(new OutputName("Divadj-CashFlows-PayDates"));
            CashFlowArraySP out2 = CashFlowArraySP(new CashFlowArray(flowAmountsOnPayDates));
            results->storeGreek(out2, Results::DEBUG_PACKET, flOut2);

            // ******** output to help debug *********
            if (isDeltaThru)
            {
                // output future weights (for debug)
                OutputNameConstSP delta_div_name(new OutputName("Divadj-futureWeight-eqDiv"));
                DoubleArraySP delta_div = DoubleArraySP(new DoubleArray(weightDiv));
                results->storeGreek(delta_div, Results::DEBUG_PACKET, delta_div_name);
                OutputNameConstSP delta_asm_name(new OutputName("Divadj-futureWeight-eqAsm"));
                DoubleArraySP delta_asm = DoubleArraySP(new DoubleArray(weightAssumed));
                results->storeGreek(delta_asm, Results::DEBUG_PACKET, delta_asm_name);
            }

            // ******** if any hist data should be recorded *********
            bool recordHist = divadjFL->recordHist();
            // below allow recording at matDate for fwdAdj
            
            // initialize payDate to matDate.  
            // to do: write getCurrPayDate to set payDate right here.
            DateTime payDate = divadjFL->instSettle->settles(divadjFL->matDate, divadjFL->asset.get());

            DateTime divPayDate;
            if (!divadjFL->findDivPayDate(valDate, divadjFL->eqDivArr, divPayDate))
                divPayDate = valDate;

            if(divadjFL->fwdAdj && flowAmounts.size() > 0 && !isZeroMu)
            {
                cashToday = flowAmounts[0].amount; // just one cash flow in fwdAdj case, (except in zero mu case)
            }

            divadjFL->divDiffToday = cashToday;

            if (recordHist)
            {
                if(divadjFL->fwdAdj)
                {
                    double fwdZeroRate = DividendAdjBase::getFwdZeroRate(divPayDate, payDate, divadjFL);

                    DividendAdjBase::processOutputRequest(OutputRequest::HIST_DISC_RATE_TO_MAT,
                        fwdZeroRate,
                        control, 
                        results);
                   
                    DividendAdjBase::processOutputRequest(OutputRequest::DRO_HIST_DISC_RATE_TO_MAT,
                                         fwdZeroRate,
                                         control, 
                                         results);
                }
                if (isZeroMu)
                {
                    if (!divadjFL->fwdAdj)
                        payDate = divPayDate;

                    double df = divadjFL->discount->pv(valDate, payDate);
                    cashToday = valueChange(avgPrice, control) / df; // amount to be paid on pay date
                    divadjFL->divDiffToday = cashToday; // just for recording
                }
                else
                {
                    if(!divadjFL->fwdAdj)
                    {
                        divadjFL->divDiffToday = divadjFL->getHist(flowAmounts, valDate);
                        cashToday = divadjFL->divDiffToday;
                        if(!findPayDate(valDate, flowAmounts, flowAmountsOnPayDates, payDate))
                            throw ModelException(method, "cannot find cash flow pay date");
                    }
                    if (StruckEquity::TYPE->isInstance(divadjFL->asset.get()))
                    {
                        double fx = DividendAdjBase::getFXValue(divPayDate, divadjFL->asset.get());
                    
                    DividendAdjBase::processOutputRequest(OutputRequest::HIST_SPOT_FX,
                        fx,
                        control, 
                        results);
                   
                    DividendAdjBase::processOutputRequest(OutputRequest::DRO_HIST_SPOT_FX,
                                         fx,
                                         control, 
                                         results);
   
                        
                    divadjFL->divDiffToday /= fx; // convention here is to record in underlying ccy
                    }
                    // if delta to be recorded
                    if (isDeltaThru)
                    {
                        DividendAdjBase::processOutputRequest( OutputRequest::HIST_DELTA,
                                                                divadjFL->deltaToday,
                                                                control,
                                                                results );

                        DividendAdjBase::processOutputRequest( OutputRequest::DRO_HIST_DELTA,
                                                                divadjFL->deltaToday,
                                                                control,
                                                                results );

                        if (Maths::isZero(divadjFL->deltaToday))
                        {
                            throw ModelException("DividendAdjustedFLClosedForm::price",
                                "Delta cannot be 0.00000000000");
                        }
                        divadjFL->divDiffToday /= divadjFL->deltaToday; // convention is amount diff in div unit
                    }
                }
                
                if (isZeroMu)
                {
                    // now that we have computed payDate, add it to knownPayments
                    CashFlow flow(payDate, cashToday);
                    knownPayments->insert(knownPayments->begin(), flow); // insert so it will be in order

                    adjustForUpDown(knownPayments, divadjFL->divadjType);
                }

            }


            /** HIST_DIV_DIFF */
            
            // first, determine whether or not we should record it.
            bool shouldRecordHistDivDiff = false;
            if (recordHist)
            {    // on ex dates, always record for zero mu, and only for non fwd adj for the other types.
                if(isZeroMu || !divadjFL->fwdAdj)
                {
                    shouldRecordHistDivDiff = true;
                }
            }
            
            if (!isZeroMu && divadjFL->fwdAdj && divadjFL->matDate.getDate() == valDate.getDate())
            {    // always record on maturity for non zero mu, fwdAdj case
                shouldRecordHistDivDiff = true;
            }

            // now, record it if required            
            if (shouldRecordHistDivDiff)
            {

                CashFlowSP flow(new CashFlow(isZeroMu? payDate : valDate, divadjFL->divDiffToday));
                // this is for histFlow recording
                 if ((request = control->requestsOutput(OutputRequest::HIST_DIV_DIFF))){
                    // put amount in HIST_DIV_DIFF: to do -- store whole cf when pyr can
                    results->storeRequestResult(request, flow->amount);
                }   
                 if ((request = control->requestsOutput(OutputRequest::HIST_DIV_DIFF_ANSI_DATE))){
                    // put amount in HIST_DIV_DIFF: to do -- store whole cf when pyr can
                    DateTime::MonthDayYear mdy = flow->date.toMDY();
                    double ansiDate = (double) (mdy.day + 100 * mdy.month + 10000 * mdy.year);

                    results->storeRequestResult(request, ansiDate);
                }   

                if ((request = control->requestsOutput(OutputRequest::DRO_HIST_DIV_DIFF))){
                    results->storeRequestResult(request, flow);
                }   
            }


            /** KNOWN_CASHFLOWS and PAYMENT_DATES */
            // always output 
            if((request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS))){
                OutputRequestUtil::recordKnownCashflows(control, 
                                                        results, 
                                                        divadjFL->discount->getCcy(), 
                                                        knownPayments.get()); 

            
            if((request = control->requestsOutput(OutputRequest::PAYMENT_DATES))){
                OutputRequestUtil::recordPaymentDates(control,results, 
                    divadjFL->getPaymentDates().get());



            }

            }
            
        }                    

    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

void DividendAdjustedFLClosedForm::adjustForUpDown(CashFlowArraySP &payments, const string& divadjType) const
{
    int i;
    for (i=0; i<payments->size(); i++)
    {
        // delete neg cf for UP, pos cf for DOWN
        if (  
             (   divadjType == "UP"    &&    Maths::isNegative( (*payments)[i].amount )   ) ||
             (   divadjType == "DOWN"  &&   !Maths::isNegative( (*payments)[i].amount )   )
           )
        {
            payments->erase(payments->begin() + i);
            i--;
        }
    }
}
        

// calc the value difference between spot + assumed and spot + eq_div across an ex-div date for ZERO_MU type
double DividendAdjustedFLClosedForm::valueChange(double price, Control* control) const
{
    // copy inst as we'll change it
    DividendAdjustedFLSP inst = DividendAdjustedFLSP(copy(divadjFL));
    // find today's div
    int idxFound = -1;
    for (int j = 0;j < inst->eqOrigDivArr.size();j++)
    {
        if(inst->eqOrigDivArr[j].getExDate().equals(inst->valueDate,false))
        {
           idxFound = j;
           break;
        }
    }
    if (idxFound == -1) // should never happen really
        throw ModelException("DividendAdjustedFLClosedForm::valueDifference", 
                "value Date is not an ex div date for recording cash flow");

    /** calculate price diff betwen spot+assumed div and spot+eq div */
    CResults results;
    // pv factor
    double df = inst->discount->pv(inst->eqOrigDivArr[idxFound].getExDate(), inst->eqOrigDivArr[idxFound].getPayDate());
    double spot = inst->asset->getSpot();
    double assumedSpot = spot + df*(inst->eqOrigDivArr[idxFound].getDivAmount()
                                    - inst->assumedDivs[idxFound].amount);

    SpotLevel assumedDivSpot(assumedSpot);

    // We only want to override the equity spot (don't want to also override any FXAsset spot!)
    OutputNameSP equityName(new OutputName(inst->asset->getTrueName()));
    assumedDivSpot.findAndShift(inst, equityName);
    model->Price(inst->buildAvgSpot(inst->avgOut, inst->strike).get(), control, &results);
    double assumedPrice = results.retrievePrice();

    return (assumedPrice - price);
}

// calcualte delta's for weighting div pass through. Delta's are pure model delta in pay off ccy for one option
// returns today's delta, to be recorded if needed
double DividendAdjustedFLClosedForm::computeDelta(DoubleArray& weightDiv,
                                                  DoubleArray& weightAssumed) const
{
    static const string method = "DividendAdjustedFLClosedForm::computeDelta";

    double absStrike;
    
    CVolRequestLNSP volRequest(volInterp(absStrike));
    // interpolate the vol using our LN request
    CVolProcessedBSSP  volBS(divadjFL->asset->getProcessedVol(volRequest.get()));

    DateTime imntStartDate = divadjFL->valueDate;
    /* if it's forward starting, convert percentage strikes to absolute values */
    if (divadjFL->fwdStarting)
    {
        imntStartDate = divadjFL->startDate;
        absStrike *= divadjFL->asset->fwdValue(divadjFL->startDate);
    }

    // variance 
    double variance = divadjFL->avgOut->averageVariance(volBS.get(),
                                                imntStartDate,
                                                false);

    AssetSP asset = AssetSP( copy(divadjFL->asset.get())); // we may replace divs, so copy.

    // if we are using assumed div list for delta, need to reset div list in asset's equity
    if (divadjFL->flowThruType == DividendAdjustedFL::DELTA_ASSUMED) 
        divadjFL->substituteEqDiv(asset);
     
    // compute weighted fwd values (average fwd)
    double fwdForDelta = divadjFL->avgOut->futureSampleSum(asset.get(), divadjFL->valueDate );

    Actual365F dayCount; // should be good enough, or use volBS->GetTimeMetric()->yearFrac()
    
    const DateTime& lastOutDate = divadjFL->avgOut->getLastDate();
    double T = dayCount.years(imntStartDate, lastOutDate);
    double futureWeight = divadjFL->avgOut->futureWeight(divadjFL->valueDate);
    
    int futAssmd = 0;

    // get on the fly num ex div dates to properly handle mu special
    CashFlowArraySP numExDivDatesInUseSP = divadjFL->getNumDivDates(divadjFL->eqDivArr);
    CashFlowArray& numExDivDatesInUse = *numExDivDatesInUseSP;

    // calc future delta's for assumed
    for (int i =0; i<divadjFL->assumedDivs.size(); i++)
    {
        DateTime lastDayPd = (i == divadjFL->assumedDivs.size() - 1) ? lastOutDate.rollDate(-1):
                                    divadjFL->assumedDivs[i + 1].date.rollDate(-1); 
    
        if (lastDayPd.getDate() >= divadjFL->valueDate.getDate())
        {
            // need calc only if period is empty, so that last date of assumed can be a pay date
            if ((int)numExDivDatesInUse[i].amount == 0)
            {
                double t = dayCount.years(imntStartDate, lastDayPd);
                //fwd delta
                weightAssumed[i] = fwdDelta(divadjFL->isCall, fwdForDelta, absStrike , variance, t/T);
                // average down delta's
                weightAssumed[i] *= divadjFL->avgOut->futureWeight(lastDayPd)/futureWeight;

            }
        }
        else // use history
        {
            futAssmd = i+1; // keep track of future starting pt 
            if (futAssmd == divadjFL->assumedDivs.size())
                futAssmd --; // use the last one, bust should never some here ?

        }
    }

    // calc delta's for equity divs
    for (int j = 0; j < divadjFL->eqDivArr.size(); j++)
    {
        DateTime exDate = divadjFL->eqDivArr[j].getExDate();
        if (exDate >= divadjFL->assumedDivs[futAssmd].date)
        {
            if (exDate.getDate() < divadjFL->valueDate.getDate())
            {// use history
                weightDiv[j] = divadjFL->getHist(divadjFL->histDelta, exDate);
            }
            else
            {
                double t = dayCount.years(imntStartDate, exDate);
                //fwd delta
                weightDiv[j] = fwdDelta(divadjFL->isCall, fwdForDelta, absStrike, variance, t/T);
                // average down delta's
                weightAssumed[j] *= divadjFL->avgOut->futureWeight(exDate)/futureWeight;
            }
        }
    }
    // return today's delta for adjustment
    return fwdDelta(divadjFL->isCall, fwdForDelta, absStrike, variance, 0.0); 
}

// same kind of vol request as average spot, also set strike used for interp
CVolRequestLN* DividendAdjustedFLClosedForm::volInterp(double& adjStrike) const {
    static const string method = "DividendAdjustedClosedForm::::volInterp";
    try {
        CVolRequestLN* interp = 0;
        
        double sumSoFar = divadjFL->avgOut->sumToDate(divadjFL->valueDate);

        if (divadjFL->avgOut->countFutureSamples(divadjFL->valueDate) > 0) {
            double futureWeight = divadjFL->avgOut->futureWeight(divadjFL->valueDate);
        
            if (Maths::isZero(futureWeight)) {
                throw ModelException(method, 
                                     "total of future avg out weights "
                                     "must be > 0.0");
            }

            adjStrike = (divadjFL->strike - sumSoFar)/futureWeight;
        }
        else {
            adjStrike = divadjFL->strike - sumSoFar;
        }
            
        // if effective strike is -ve, use ATM vol 
        if (Maths::isPositive(adjStrike)) {
            DateTime effectiveEnd;
            if (divadjFL->fwdStarting) {
                effectiveEnd = divadjFL->avgOut->expectedEndDate(divadjFL->valueDate);
            }

            interp = new LinearStrikeVolRequest(adjStrike, 
                                                divadjFL->fwdStarting ? divadjFL->startDate : 
                                                divadjFL->valueDate, 
                                                divadjFL->fwdStarting ? effectiveEnd :
                                                divadjFL->avgOut->getLastDate(),
                                                divadjFL->fwdStarting);
        }
        else {
            interp = new ATMVolRequest();
        }            
        return interp;
    }
    catch (exception& e) {
        throw ModelException(&e, method);
    }
}

// calculate weights (delta) and cashflows 
// this method is similar to calcStrike() in DividendAdjusted
void DividendAdjustedFL::calcFlowData(
         const DoubleArraySP paidSoFar, // needed for flow through, unadjusted (same unit as assumed div)
         const DoubleArray& weightDiv,
         const DoubleArray& weightAssumed,
         CashFlowArraySP    flowThroughCashFlows) const
{
    const string method = "DividendAdjustedFL::calcFlowData";

    if (paidSoFar->size() != assumedDivs.size())
    {
        throw ModelException(method, "when computing cashflows (as in the flow through case, paidSoFar size"
            " (equals " + Format::toString(paidSoFar->size()) + " but should be equal to the size of "
            "the assumed divs array (equals " + Format::toString(assumedDivs.size()) + ").");
    }

    int i,j;

    DoubleArray divDollars(eqDivArr.size());
    DoubleArray fxDiv(eqDivArr.size());
    DoubleArray fxAsmd(assumedDivs.size());
    DoubleArray gfDiv(eqDivArr.size());
    DoubleArray gfAsmd(assumedDivs.size());

    // get dates lists first
    DateTimeArray eqExDates(eqDivArr.size());
    for (i=0; i<eqExDates.size(); i++)
        eqExDates[i] = eqDivArr[i].getExDate();

    DateTime currAssumedStartDate = valueDate;
    DateTimeArray assumedDates(assumedDivs.size());
    for (i=0; i<assumedDates.size(); i++)
    {
        assumedDates[i] = assumedDivs[i].date;
        if (assumedDates[i] <= valueDate)
            currAssumedStartDate = assumedDates[i];
    }

    // last avg out date
    const DateTime& lastOutDate = avgOut->getLastDate();

    // get FX and fwd discount for the actual dividend dates
    calcFXandDFdiv(fxDiv,gfDiv);

    double pastDivDollars = 0.0;
    // calc dollar amount
    for (j = 0;j < eqDivArr.size();j++)
    {
        // if we are computing (at hit) cashflows, don't include past divs in the sum
        if (fwdAdj || eqExDates[j].getDate() >= valueDate.getDate())
        {
            divDollars[j] = div2dollar(eqDivArr[j], asset.get());
        }
        else if(eqExDates[j].getDate() >= currAssumedStartDate.getDate())
        {
            pastDivDollars += div2dollar(eqDivArr[j], asset.get());;
        }
    }

    //calc the total value of the actual div in the assumed periods
    DoubleArray periodSum(assumedDates.size());
    subTotal(assumedDates, lastOutDate, true, eqExDates, divDollars, periodSum); 

    // get FX and fwd discount for the end of the assumed periods
    calcFXandDFassumed(periodSum, fxAsmd, gfAsmd);

    // weights for the propotional strike adjustment, on eq div intervals
    DoubleArray prop(eqDivArr.size(), 1.0); // init to 1.0 if no adjustments needed
    if (avgWeightAdjusted)
        propAdjust(eqExDates, avgOut.get()->getDates(), avgOut.get()->getWeights(), prop); 

    // weights for the propotional strike adjustment, for assumed intervals
    DoubleArray propexp(assumedDivs.size(), 1.0); // init to 1.0 if no adjustments needed
    if (avgWeightAdjusted)
        propAdjust(assumedDates, avgOut.get()->getDates(), avgOut.get()->getWeights(), propexp); 
    
    int k = 0;
    double sum1 = 0.0;
    double sumexp = 0.0; 
    double m;
    //main loop 
    double fwdFlowThru = 0.0;
    for (i = 0; i<assumedDivs.size(); i++)
    { 
        bool isLast = (i==assumedDivs.size()-1);
        const DateTime& endDate = (isLast ? lastOutDate : assumedDates[i+1]);
        // compute pastDivAdjust, and paidSoFarThisPd, required for correct adjustment to flow through payment
        double pastDivAdjust = 
            ( assumedDates[i].getDate() < valueDate.getDate() &&
              (i == assumedDates.size() - 1 || valueDate.getDate() <= assumedDates[i+1].getDate()) 
            ) ?
            pastDivDollars : 0.0;
        double paidSoFarThisPd = !!paidSoFar? (*paidSoFar)[i] : 0.0;

        // get on the fly num ex div dates to properly handle mu special
        CashFlowArraySP numExDivDatesInUseSP = getNumDivDates(eqDivArr);
        CashFlowArray& numExDivDatesInUse = *numExDivDatesInUseSP;

        // assumed period contains actual dividends
        if (!Maths::isZero(numExDivDatesInUse[i].amount))
        {
            // if there are only 0 divs in the qtr, then equally weight payments by setting periodSum to 1 / numDivs.
            double denom = Maths::isZero(periodSum[i]) ? 1.0 / numExDivDatesInUse[i].amount : periodSum[i];

            for (j = k; j<eqDivArr.size(); j++)
            {   
                double numer = 1.0;
                double x = 1.0;
                if (!Maths::isZero(periodSum[i]))
                {
                    numer = Maths::isZero(periodSum[i]) ? 1.0 : divDollars[j];
                    x = 1.0;
                }
                else
                {   // div sum is 0: instead of  {D - AD/sum div} get {- A / [num div]}
                    numer = 1.0;
                    x = 0.0;
                }

                if(
                    eqExDates[j] >= assumedDates[i] &&
                    ((!isLast && eqExDates[j] < endDate) ||
                    (isLast && eqExDates[j] <= endDate)) // last date is inclusive
                    
                    &&
                    
                    !(!fwdAdj && eqExDates[j].getDate() < valueDate.getDate()) // no flow on past divs
                 ) 
                {
                    if(eqDivArr[j].getDivType() == Dividend::AMOUNT)
                        m = 1.0; // not yield case
                    else
                        m = 0.0; // yield case

                    // adjustment
                    sum1 += prop[j]*fxDiv[j]*numer*(m-assumedDivs[i].amount/denom)*gfDiv[j];
                
                    double expectedStrikeAdjust = prop[j]*fxDiv[j]*numer*
                                    (x-(assumedDivs[i].amount + paidSoFarThisPd - pastDivAdjust)/denom);
                    sumexp += expectedStrikeAdjust * gfDiv[j];
                    if(!fwdAdj)
                    {
                        flowThroughCashFlows->push_back(CashFlow(eqDivArr[j].getExDate(),
                                                    expectedStrikeAdjust*weightDiv[j]));
                    }
                    else
                        fwdFlowThru += expectedStrikeAdjust*weightDiv[j] * gfDiv[j];
                }
                if(eqExDates[j] >= endDate)
                {
                    k = j;
                    break;
                }
            }
        }
        else
        {// no actual dividends in the assumed period
            double expectedStrikeAdjust = -propexp[i] * fxAsmd[i] * 
                                (assumedDivs[i].amount + paidSoFarThisPd - pastDivAdjust);
            sum1 +=   expectedStrikeAdjust * gfAsmd[i];
            sumexp += expectedStrikeAdjust * gfAsmd[i];
            if(!fwdAdj)
            {
                flowThroughCashFlows->push_back(CashFlow(isLast? endDate:endDate.rollDate(-1),
                                 expectedStrikeAdjust*(weightAssumed[i])));
            }
            else
                fwdFlowThru += expectedStrikeAdjust*weightAssumed[i]* gfAsmd[i];
        }
    }

    // only one cash flow for fwdAdj
    if(fwdAdj)
    {
        const DateTime&  matPayDate = avgOut->getLastDate(); // should be final settlementDate ?
        flowThroughCashFlows->push_back(CashFlow(matPayDate, fwdFlowThru));
    }
}

/** substitute equity div's by assumed div's */
// for ZERO_MU flowType, eqDiv exDates and pay dates are kept
void DividendAdjustedFL::substituteEqDiv(AssetSP eq, bool omitSpecialDivDate,
                     const DateTime& specialDivDate) const

{
    int i = 0;
    // create new div list, using assumed divs.
    DividendArray divarr(assumedDivs.size());
    for (i = 0; i < assumedDivs.size(); i ++)
    {
        if(flowThruType == ZERO_MU){

            divarr[i] = Dividend(eqDivArr[i].getExDate(),
                                 eqDivArr[i].getPayDate(),
                                 Dividend::AMOUNT,
                                 assumedDivs[i].amount);
        }
        else{
            divarr[i] = Dividend(assumedDivs[i].date,
                                 assumedDivs[i].date,
                                 Dividend::AMOUNT,
                                 assumedDivs[i].amount);
        }
    }

    // Delete special dividend
    if (omitSpecialDivDate)
    {
        for ( i = 0; i < divarr.size(); i ++ )
        {
            if (divarr[i].getExDate().equals(specialDivDate, false)) // ignore time
            {
                divarr.erase(divarr.begin() + i);
            }
        }
    }
    
    
    // see if asset can return its equity.  If so, replace the divs.
    IHaveEquity* haveEquity = dynamic_cast<IHaveEquity*>(eq.get());
    if (haveEquity)
    {
        DividendListSP newDivList = DividendListSP(new DividendList(divarr));
        haveEquity->getEquity()->setDivList(newDivList.get());
    }
    else
    {
        throw ModelException("DividendAdjustedFL::substituteEqDiv",
                            "couldn't cast asset to IHaveEquity");
    }
}

/** for each period, returns the total payments made, in growth ccy */
DoubleArraySP DividendAdjustedFL::getPastPeriodPayments() const
{
    int i;
    DoubleArraySP pastPeriodPayments = DoubleArraySP(new DoubleArray(assumedDivs.size()));

    if (fwdAdj)
        return pastPeriodPayments;

    //calc the sum of all payments made in each period
    DateTimeArray eqExDates(eqDivArr.size());
    DoubleArray      exDatePayments(eqDivArr.size());

    for (i=0; i<eqExDates.size(); i++)
    {
        eqExDates[i] = eqDivArr[i].getExDate();
        
        if (!fwdAdj && eqExDates[i].getDate() < valueDate.getDate())
        {// use history
            exDatePayments[i] = getHist(histFlow, eqDivArr[i].getExDate());
        }
        else
        {
            exDatePayments[i] = 0.0;
        }
    }

    DateTimeArray assumedDates(assumedDivs.size());
    for (i=0; i<assumedDates.size(); i++)
        assumedDates[i] = assumedDivs[i].date;
    
    subTotal(assumedDates, avgOut->getLastDate(), true, eqExDates, exDatePayments, *pastPeriodPayments); 

    return pastPeriodPayments;
}

DateTimeArraySP DividendAdjustedFL::getPaymentDates()
{
    DateTimeArraySP allPayDates = DateTimeArraySP(new DateTimeArray(0));
    DividendArray divArray = getEqDivs();
    if (!fwdAdj)
    {
        for (int i=0; i<divArray.size(); i++)
        {
            allPayDates->push_back(divArray[i].getPayDate());
        }
    }

    // add maturity date
    allPayDates->push_back(instSettle->settles(matDate, asset.get()));
    
    return allPayDates;
}

bool DividendAdjustedFL::dateInHist(const DateTime& aDate, const CashFlowArray& cfArr) const
{
    bool found = false;
    for (int i = 0; i < cfArr.size(); i++)
    {
        if (cfArr[i].date.getDate() == aDate.getDate())
        {
            found = true;
            break;
        }
    }
    
    return found;
}
/** Rolls the value date and sets initial spot if rolling over start date */
bool DividendAdjustedFL::sensShift(Theta* shift)
{    
    DateTime newDate = shift->rollDate(valueDate);
    // if today is an ex-date, compute pay date.  thetaHelper will use it for computing fwd rates.
    DateTime divPayDate;
    if (! findDivPayDate(valueDate, getEqDivs(), divPayDate))
    {
        divPayDate = valueDate;
    }

    DateTime payDate = fwdAdj ?
                       instSettle->settles(matDate, asset.get())
                       : divPayDate;
    
    thetaHelper(shift, 
                newDate, 
                true,       // useFwd
                divPayDate);   // dividend pay date

    if (newDate == valueDate) 
        return true; 

    if(recordHist())
    {
        if (flowThruType == DividendAdjustedFL::ZERO_MU)
        {
            histFlow.push_back(CashFlow(payDate, divDiffToday));
        }
        else if(!fwdAdj)
        {
            histFlow.push_back(CashFlow(valueDate, divDiffToday));
        }

        if (flowThruType == DividendAdjustedFL::DELTA_ASSUMED
            || flowThruType == DividendAdjustedFL::DELTA_EQ_DIV)
            histDelta.push_back(CashFlow(valueDate, deltaToday));

    }
    
    // for all ex dates between valueDate and newDate, push back histFlow.
    DividendArray eqDvs = getEqDivs();
    
    // currently, don't support hopping over ex dates for types other than NOTIONAL or ZERO_MU
    // also, don't support hopping for cases with more than 1 ex div date / qtr

    bool grOneDv = false;
    for (int j = 0; j < numExDivDates.size(); j++) {
        if (numExDivDates[j].amount >= 2.0) {
            grOneDv = true;
        }
    }

    double sgn = (flowThruType == DividendAdjustedFL::NOTIONAL && !isCall)? -1.0 : 1.0;
    bool pastNewDate = false;
    int  i = 0;
    while ( !pastNewDate && i < eqDvs.size() ) {
        const DateTime &exDate = eqDvs[i].getExDate();
        if ( exDate.getDate() > valueDate.getDate() && exDate.getDate() < newDate.getDate() ) {
            if (grOneDv) {
                throw ModelException("DividendAdjustedFL::sensShift", 
                    "Functionality not yet supported!  \n"
                    "Cannot roll over ex-dates after tomorrow for theta when sched has > 1 div / qtr.");
            }
                
            if (flowThruType != DividendAdjustedFL::NOTIONAL && flowThruType != DividendAdjustedFL::ZERO_MU) {
                throw ModelException("DividendAdjustedFL::sensShift", 
                    "Functionality not yet supported!  \n"
                    "Cannot roll over ex-dates after tomorrow for theta for DELTA types.");
            }

            if (fwdAdj) {
                throw ModelException("DividendAdjustedFL::sensShift", 
                    "Functionality not yet supported!  \n"
                    "Cannot roll over ex-dates after tomorrow for theta when is Fwd Adj = TRUE");
            }

            // get the div's pay date for use below
            DateTime thisDivPayDate;
            if (! findDivPayDate(exDate, getEqDivs(), thisDivPayDate))
            {
                thisDivPayDate = exDate;
            }
            
            // currently do not support unpaid zero mu payments (to do: can approximate with forward delta.)
            if( flowThruType == DividendAdjustedFL::ZERO_MU &&  !newDate.isGreater(thisDivPayDate) ) {
                throw ModelException("DividendAdjustedFL::sensShift", 
                    "Functionality not yet supported!  \n"
                    "Cannot roll over ex-dates that haven't been paid yet for ZER0_MU case.");
            }

            double divDiff = sgn*computeDivDiff(eqDvs[i]);
            if (divadjType == "UP")        { divDiff = Maths::max(divDiff, 0.0); }
            else if (divadjType == "DOWN") { divDiff = Maths::min(divDiff, 0.0); }

            // update histFlow
            histFlow.push_back(CashFlow(exDate, divDiff));

            // update FX Rates
            if (StruckEquity::TYPE->isInstance(asset.get())) {
               
                DateTime thisPayDate = fwdAdj ?
                    instSettle->settles(matDate, asset.get())
                    : thisDivPayDate;
                            
                    histFXRates.push_back( CashFlow(exDate, 
                                getFXValue(thisPayDate, asset.get())) );

            }

        } else if ( exDate.isGreater(newDate) ) {
            pastNewDate = true; // don't need to look anymore, as this and remaining ex dates > new date.
        }
        ++i;
    }

    // roll today 
    valueDate = newDate;
    
    return true;
}

// for reflection
DividendAdjustedFL::DividendAdjustedFL(): DividendAdjBase(TYPE),
                                                    flowThruType(DELTA_ASSUMED),
                                                    deltaToday(0.0),
                                                    divDiffToday(0.0)
                                                    {usePayDate = true;}

class DividendAdjustedFLHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(DividendAdjustedFL, clazz);
        SUPERCLASS(DividendAdjBase);
        IMPLEMENTS(CClosedFormLN::IIntoProduct);
        EMPTY_SHELL_METHOD(defaultDividendAdjustedFL);
        FIELD(isCall,           "Is it a call option");
        FIELD(divadjType,       "Type of the adjustment");
        FIELD(fwdAdj,           "Use forward adjustment");
        FIELD(avgOut,                  "average-out samples list");
        FIELD(flowThruType,      "model type: DELTA_ASSUMED or NOTIONAL");
        FIELD(histDelta,        "Historic delta values for computing div delta flow through amonuts");
        FIELD_MAKE_OPTIONAL(histDelta);
        FIELD(histFlow,          "historic recorded flows before fx/delta weighting");
        FIELD_MAKE_OPTIONAL(histFlow);

        // transient
        FIELD(deltaToday, "delta");
        FIELD_MAKE_TRANSIENT(deltaToday);
        FIELD(divDiffToday, "dividend difference recorded on ex-date, same as cash today for zero_mu type");
        FIELD_MAKE_TRANSIENT(divDiffToday);
        FIELD(eqDivArr, "");
        FIELD_MAKE_TRANSIENT(eqDivArr);
        FIELD(eqOrigDivArr, "");
        FIELD_MAKE_TRANSIENT(eqOrigDivArr);


    }

    static IObject* defaultDividendAdjustedFL(){
        return new DividendAdjustedFL();
    }
};

CClassConstSP const DividendAdjustedFL::TYPE = CClass::registerClassLoadMethod(
    "DividendAdjustedFL", typeid(DividendAdjustedFL), DividendAdjustedFLHelper::load);
   
/** Implementation of ClosedForm::IntoProduct interface */
CClosedFormLN::IProduct* DividendAdjustedFL::createProduct(
    CClosedFormLN* model) const{

    return new DividendAdjustedFLClosedForm(model, this);
}


// for derived class
DividendAdjustedFL::DividendAdjustedFL(CClassConstSP const type): 
    DividendAdjBase(type),
    flowThruType(ZERO_MU),
    deltaToday(0.0),
    divDiffToday(0.0) {usePayDate = true;}

////////////////////// DFWSpecial class //////////////////

class DFWSpecial: public DividendAdjustedFL, 
                  virtual public NumericalIntegrationLN::IIntoProduct
{
public:
    static CClassConstSP const TYPE;

    friend class AdjustedStrikeNumerical;
    /** Implementation of ClosedForm::IntoProduct interface */
    virtual NumericalIntegrationLN::IProduct* createProduct(NumericalIntegrationLN* model) const;

    // used to explicitly call parent's GetMarket with ClosedFormLN model to avoid spline
    void GetMarket(const IModel*          model, 
                   const CMarketDataSP    market);

private:
    // utility methods
    double specialDivAmount() const;
    double specialAssumedAmount() const;
    DateTime specialDivPayDate() const;
    void zeroEqDivs(AssetSP eq) const;
    void setVol(CAsset* asset, double vol) const;
    CAssetSP makeFlatVolAsset(CAssetSP asset) const;
        
    // new fields
    double strikeShiftForVolInterp;
    DateTime specialDivDate;
    bool useFwdVariance; // if true properly compute fwd variance at day before div date.
                         // otherwise, compute vol from today 
    double regularDivSize; // the regular small div that is used in computing yield 
                           // in computing forward at day before div date
    bool useExpStrikeForPrice;

    // transient
    mutable double cachedNewStrike; 

    friend class DFWSpecialHelper;
    friend class DFWSpecialClosedForm;

    DFWSpecial();
    DFWSpecial(const DFWSpecial& rhs);
    DFWSpecial& operator=(const DFWSpecial& rhs);
    static void load(CClassSP& clazz);
    
protected:

    // input in addition to those in DividendAdjBase already
};

typedef smartPtr<DFWSpecial> DFWSpecialSP;

// for reflection
DFWSpecial::DFWSpecial(): DividendAdjustedFL(TYPE), 
                          strikeShiftForVolInterp(3.0),
                          specialDivDate(DateTime(0,0)),
                          useFwdVariance(false),
                          regularDivSize(.16),
                          useExpStrikeForPrice(false),
                          cachedNewStrike(0.0){};


////////////// class methods ////////////////////////////


////////////////////// DFWSpecialHelper class //////////////////

class DFWSpecialHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(DFWSpecial, clazz);
        SUPERCLASS(DividendAdjustedFL);
        IMPLEMENTS(NumericalIntegrationLN::IIntoProduct);
        EMPTY_SHELL_METHOD(defaultDFWSpecial);
        FIELD(strikeShiftForVolInterp,           "Volatility for fair value calculation is interpolated at K - strikeShiftForVolInterp");
        FIELD(specialDivDate,       "Date of special dividend");
        FIELD(useFwdVariance,       "True, compute fwd variance when computing BS at div date");
        FIELD(regularDivSize,       "Size of regular contract div used to to compute yield for fwd computation at day before div date");
        FIELD_MAKE_OPTIONAL(regularDivSize);
        FIELD(useExpStrikeForPrice, "True, use Black Scholes value with expected strike for price.");
        FIELD_MAKE_OPTIONAL(useExpStrikeForPrice);
        FIELD(cachedNewStrike, "saved new strike for theta on day before ex date");
        FIELD_MAKE_TRANSIENT(cachedNewStrike);
    }

    static IObject* defaultDFWSpecial(){
        return new DFWSpecial();
    }
};

CClassConstSP const DFWSpecial::TYPE = CClass::registerClassLoadMethod(
    "DFWSpecial", typeid(DFWSpecial), DFWSpecialHelper::load);
   
bool DividendAdjustedFLLoad()
{
    return DFWSpecial::TYPE && true; // && true to shut the compiler up
}

void DFWSpecial::GetMarket(const IModel*          model, 
                                 const CMarketDataSP    market)
{
    // this goes against the philosophy of the library as it is
    // supposed to be up to the model to choose the market data
    IModelSP   cfln(new CClosedFormLN("VolPreferred"));

    // explicitly call parent with ClosedFormLN to avoid spline
    DividendAdjBase::GetMarket(cfln.get(), market);
}

/////////////////////////////// class AdjustedStrikeNumerical //////////////
class AdjustedStrikeNumerical: virtual public NumericalIntegrationLN::IProduct
{

    typedef enum{ADJUSTED_STRIKE, PRICE, SKEW} TComputeMode;

private:
    const DFWSpecial* dfwSpecial;
    static double solvePrecision;
    static int maxIter;

    TComputeMode computeMode;

    DateTime dayBeforeSpecialDiv;
    double specialDivSize;
    double specialAssumedSize;

    CAssetSP assetNew; // asset with assumed dividends.  
    CAssetSP assetFlatVol; // asset with flat vol, interpolated at original strike - strike shift at maturity



    double fwdNoDivsDayBeforeSpecialDiv;
    double fwdNoDivsMat;

    double fwdDayBeforeSpecialDiv; 
    double fwdMat; 

public:
       AdjustedStrikeNumerical(DFWSpecial* dfwSpecial) : 
       dfwSpecial(dfwSpecial) 
       {
           HolidayConstSP hols(AssetUtil::getHoliday(dfwSpecial->asset.get()));
           dayBeforeSpecialDiv = hols->addBusinessDays(dfwSpecial->specialDivDate, -1);
           specialDivSize = dfwSpecial->specialDivAmount(); 
           specialAssumedSize = dfwSpecial->specialAssumedAmount(); 

               
           assetNew = CAssetSP (copy(dfwSpecial->asset.get()));
           assetFlatVol = dfwSpecial->makeFlatVolAsset(assetNew);
           
           fwdNoDivsDayBeforeSpecialDiv = 1.0 / dfwSpecial->discount->pv(dayBeforeSpecialDiv);
           fwdNoDivsMat = 1.0 / dfwSpecial->discount->pv(dfwSpecial->avgOut->getLastDate());

           dfwSpecial->substituteEqDiv(assetNew, true, dfwSpecial->specialDivDate);
           
           if (!isSpecial())
           {
               fwdDayBeforeSpecialDiv = 1.0;
               fwdMat = 1.0;
           }
           else
           {
               fwdDayBeforeSpecialDiv = assetNew->fwdValue(dayBeforeSpecialDiv);
               fwdMat = assetNew->fwdValue(dfwSpecial->avgOut->getLastDate());
           }
       }
       
       bool isSpecial()
       {
           return !dfwSpecial->valueDate.isGreater(dayBeforeSpecialDiv);
       }

       void price(NumericalIntegrationLN* model,
           Control*                control, 
           CResults*               results) 
        {
            // calculate the discount factor back to today
            double discFactor = dfwSpecial->instSettle->pv(dfwSpecial->valueDate,
                dfwSpecial->avgOut->getLastDate(),
                dfwSpecial->discount.get(), 
                dfwSpecial->asset.get());
            
            // compute scale factor
           double scalingFactor = InstrumentUtil::scalePremium(dfwSpecial->oneContract,
               false,
               dfwSpecial->notional,
               1.0,
               dfwSpecial->initialSpot);
 

           if (!isSpecial())
           {

                DateTime maturity = dfwSpecial->avgOut->getLastDate();
               
                double fwd = assetNew->fwdValue(maturity);
                double strikeToUse = dfwSpecial->strike;

                if (dfwSpecial->valueDate.equals(dfwSpecial->specialDivDate,false) &&
                    model->doingTimeShift())
                {   // if we are rolling to ex-date, reduce fwd by fwd value of large div
                    strikeToUse = dfwSpecial->cachedNewStrike;
                    fwd -= (dfwSpecial->specialDivAmount() 
                        / dfwSpecial->discount->pv(dfwSpecial->specialDivPayDate(),maturity));
                      // as our parent doesn't use restoreable shift
                }

                LinearStrikeTSVolRequest volRequest(
                   strikeToUse,
                   dfwSpecial->valueDate,
                   maturity,
                   false);  
            
                CVolProcessedBSSP volBS(assetNew->getProcessedVol(&volRequest));
                double var = volBS->CalcVar(dfwSpecial->valueDate, maturity);            
                
                double zero_mu_price = Black::price(
                    dfwSpecial->isCall, fwd, strikeToUse, 1.0, var);

                results->storePrice(
                    zero_mu_price * scalingFactor * discFactor , dfwSpecial->discount->getCcy());

                return;
           }
           else if(dfwSpecial->valueDate.equals(dayBeforeSpecialDiv, false))
           {
                // simply compute new strike, and then call bsExpectedStrike for premium:
               computeMode = ADJUSTED_STRIKE;
               double new_strike = payoff(assetNew->getSpot());
               double price_on_day_before_special_div = bsExpectedStrike(new_strike);
                results->storePrice(
                    price_on_day_before_special_div * scalingFactor * discFactor , dfwSpecial->discount->getCcy());

               
                // other outputs
                if (control && control->isPricing() ) 
                {
                    // adjStrike again
                    OutputNameConstSP adjStrikeO(new OutputName("debug_exp_strike"));
                    results->storeGreek(CDoubleSP(CDouble::create(new_strike)), Results::DEBUG_PACKET, 
                        adjStrikeO);
                    
                    OutputRequest* request = NULL;
                    
                    if ((request = control->requestsOutput(OutputRequest::ADJUSTED_STRIKE))){
                        results->storeRequestResult(request, new_strike);
                    }
                    
                    dfwSpecial->cachedNewStrike = new_strike;  // store for theta
                }
                
                return;
            
           }
           // computed adjusted price
           // adjStrike
           double adjStrike = expectedAdjustedStrike(model);
           
           // option price
           double premium = 0.0;
           
           if (dfwSpecial->useExpStrikeForPrice)
           {
               premium = bsExpectedStrike(adjStrike);
           }
           else
           {
               premium = optionPrice(model); 
           }
               
           results->storePrice(premium * scalingFactor * discFactor , dfwSpecial->discount->getCcy());
           
           // try catch todo!!!!
           
           // other outputs
           if (control && control->isPricing() ) 
           {
               OutputRequest* request = NULL;
               
               if ((request = control->requestsOutput(OutputRequest::ADJUSTED_STRIKE))){
                   results->storeRequestResult(request, adjStrike);
               }

                           
                 // adjStrike again
               OutputNameConstSP adjStrikeO(new OutputName("debug_exp_strike"));
               results->storeGreek(CDoubleSP(CDouble::create(adjStrike)), Results::DEBUG_PACKET, 
                   adjStrikeO);
               
               // f1
               OutputNameConstSP debug_fwd1(new OutputName("debug_fwd1"));
               results->storeGreek(CDoubleSP(CDouble::create(fwdDayBeforeSpecialDiv)), Results::DEBUG_PACKET, 
                   debug_fwd1);
               
               // f2
               OutputNameConstSP debug_fwd2(new OutputName("debug_fwd2"));
               results->storeGreek(CDoubleSP(CDouble::create(fwdMat)), Results::DEBUG_PACKET, 
                   debug_fwd2);
               
               // pv factor
               OutputNameConstSP debug_pv(new OutputName("debug_pv"));
               results->storeGreek(CDoubleSP(CDouble::create(discFactor)), Results::DEBUG_PACKET, 
                   debug_pv);

               // skew
               
               // first combpute discount factor to day before ex-date
            double discFactorToDayBeforeEx = dfwSpecial->instSettle->pv(dfwSpecial->valueDate,
                dayBeforeSpecialDiv,
                dfwSpecial->discount.get(), 
                dfwSpecial->asset.get());

               OutputNameConstSP debug_skew(new OutputName("debug_skew"));
               results->storeGreek(CDoubleSP(CDouble::create(discFactorToDayBeforeEx * premiumFromSkew(model))), Results::DEBUG_PACKET, 
                   debug_skew);
               
               // zero mu price.  this + skew should equal price (actually we are fudging: this would be true
               // if we used the same vol to integrate at div - 1.  However, we are using vol from smile to 
               // do this integration; so in fact the price returned is actually better than zero mu price
           
               double origPremNoPv = origPremium();
               OutputNameConstSP debug_origPremium(new OutputName("debug_origPremium"));
               results->storeGreek(CDoubleSP(CDouble::create(origPremNoPv * scalingFactor * discFactor)), Results::DEBUG_PACKET, 
                   debug_origPremium);
               
               
           }
           
        
        }
       
        double expectedAdjustedStrike(NumericalIntegrationLN* ni)
        {
            computeMode = ADJUSTED_STRIKE;
            return  ni->integrate(this);
        }

        double optionPrice(NumericalIntegrationLN* ni)
        {
            computeMode = PRICE;
            return  ni->integrate(this);
        }
        double premiumFromSkew(NumericalIntegrationLN* ni)
        {
            computeMode = SKEW;
            return  ni->integrate(this);
        }
        // returns undiscounted zero_mu price with the strike shift in vol interp.
        double origPremium()    
        {
            DateTime maturity = dfwSpecial->avgOut->getLastDate();
            
            LinearStrikeTSVolRequest volRequest(
                dfwSpecial->strike - dfwSpecial->strikeShiftForVolInterp, // no shift
                dfwSpecial->valueDate,
                maturity,
                false);  
            
            CVolProcessedBSSP volBS(assetNew->getProcessedVol(&volRequest));
            
            double origVariance = volBS->CalcVar(dfwSpecial->valueDate, maturity);                                   
            
            return Black::price(dfwSpecial->isCall, fwdMat, dfwSpecial->strike, 1.0, origVariance ); 
        }

        // bsExpectedStrike
        double bsExpectedStrike(double k)    
        {
            DateTime maturity = dfwSpecial->avgOut->getLastDate();
            
            LinearStrikeTSVolRequest volRequest(
                k,
                dfwSpecial->valueDate,
                maturity,
                false);  
            
            CVolProcessedBSSP volBS(assetNew->getProcessedVol(&volRequest));
            
            double var = volBS->CalcVar(dfwSpecial->valueDate, maturity);            
            
            double fwdToMatWithSpecialDiv = fwdMat - specialDivSize
                                                *fwdNoDivsMat/fwdNoDivsDayBeforeSpecialDiv;

            return Black::price(dfwSpecial->isCall, fwdToMatWithSpecialDiv, k, 1.0, var); 
        }

        /** following methods drive the numerical integration */
       /** gives the distribution of the fwd */
        PDFCalculator* pdfCalculator() const 
        {       
            
            static const string method = "AdjustedStrikeNumerical::pdfCalculator";
            try {

                // use assetFlatVol in creating pdfCalculator
                LinearStrikeTSVolRequest volRequest(
                    0.0, // interp level of 0.0 as a dummy.  
                    dfwSpecial->valueDate,
                    dayBeforeSpecialDiv,
                    false);  
                
                PDFRequestLNStrike pdfRequest(&volRequest);
                
                return assetFlatVol->pdfCalculator(&pdfRequest);
            }
            catch (exception& e) {
                throw ModelException(e, method);
            }
        }

        /** given a spot level, return the payoff */
        double payoff(double spot) const 
        {
            //DateTime dayBeforeSpecialDiv = dfwSpecial->specialDivDate.rollDate(-1); // should do this in constructor
            
            // todo: set to EOD???
            //double specialDivSize = 3.0; // todo!!!!!
            DateTime maturity = dfwSpecial->avgOut->getLastDate();
            
            if (spot < specialDivSize)
            {
                return dfwSpecial->strike;          
            }

            double timeToMat = (maturity.getDate() - dayBeforeSpecialDiv.getDate())/365.0;
            double regularDivDiscount = exp(-timeToMat * dfwSpecial->regularDivSize / spot);

            // contractual def of forward used for K computation    
            double growthFactor = regularDivDiscount * fwdNoDivsMat / fwdNoDivsDayBeforeSpecialDiv;

            double origFwd = spot * growthFactor;
            double origVariance = fwdVarianceToMat(dfwSpecial->strike, dayBeforeSpecialDiv, maturity);

            double origPrice = Black::price(dfwSpecial->isCall, origFwd, dfwSpecial->strike, 1.0, origVariance); 
            
            // now set up solver
            // for purposes of computing the strike, newFwd reduces the spot by the difference of the special div and the assumed div on the 
            // special div date.
            double newFwd = (spot - ( specialDivSize - specialAssumedSize )) * growthFactor;  
            double hiBound = dfwSpecial->strike;
            double loBound = 0.0;
    
            double iter = 0;

            while (hiBound - loBound > AdjustedStrikeNumerical::solvePrecision && 
                                        iter < AdjustedStrikeNumerical::maxIter)
            {
                double tempStrike = (hiBound + loBound) / 2.0;
                //double newVariance = fwdVarianceToMat(tempStrike, dayBeforeSpecialDiv, maturity);
                double newPrice = Black::price(dfwSpecial->isCall, newFwd, tempStrike, 1.0, origVariance); 
                if (newPrice > origPrice) 
                {
                    loBound = tempStrike; 
                }
                else
                {
                    hiBound = tempStrike;
                }
                iter ++;
            }

            double newStrike = hiBound;
            if (computeMode == ADJUSTED_STRIKE)
            {
                return newStrike;
            }
            else // price or skew
            {
                double newVariance = fwdVarianceToMat(newStrike + dfwSpecial->strikeShiftForVolInterp, dayBeforeSpecialDiv, maturity);
                double realNewFwd = (spot - specialDivSize) * fwdMat / fwdDayBeforeSpecialDiv;  // includes future assumed divs, after special date
                // need to compute other stuff
                double newPrice = Black::price(dfwSpecial->isCall, realNewFwd, newStrike, 1.0, newVariance);

                double origPriceWithCorrectFwd = Black::price(dfwSpecial->isCall, realNewFwd, newStrike, 1.0, origVariance);
                if(computeMode == SKEW)
                {
                    return newPrice - origPriceWithCorrectFwd;
                }
                else
                {
                    return newPrice;
                }
            }
        }

       /** centre point of distribution */
        double centre(const DateTime& date) const 
        {
            return assetNew->fwdValue(dayBeforeSpecialDiv); // fwd at special div date - 1
        }
        
        /** variance at T for given K.  Note that this does not shift by strikeShiftForVolInterp since it is computed before 
            the special ex-div date */
        double variance(double strike, const DateTime& date) const {
            static const string method = "AdjustedStrikeNumerical::variance";
            try {
                LinearStrikeTSVolRequest volRequest(
                    strike, // no shift
                    dfwSpecial->valueDate,
                    date,
                    false);  

                CVolProcessedBSSP volBS(assetFlatVol->getProcessedVol(&volRequest));
                
                return volBS->CalcVar(dfwSpecial->valueDate, date);                                   
            }
            catch (exception& e) {
                throw ModelException(e, method);
            }   
        }
        
       /** when are we getting the distrubution for */
        DateTime time() const 
        {
            return dayBeforeSpecialDiv;
        }

        // destructor?
        
        // utility function        
        double fwdVarianceToMat(double strike, const DateTime& fwdDate, const DateTime& maturity) const
        {
            
            LinearStrikeTSVolRequest volRequest(
                strike  - dfwSpecial->strikeShiftForVolInterp,
                dfwSpecial->valueDate,
                maturity,
                false);  
            
            CVolProcessedBSSP volBS(assetNew->getProcessedVol(&volRequest));
            
            double varianceToMat = volBS->CalcVar(dfwSpecial->valueDate, maturity);

            if (dfwSpecial->useFwdVariance)
            {
                double spot = dfwSpecial->asset->getSpot();
                double scale = spot / (spot - specialDivSize) ;
                scale *= scale;
                
                LinearStrikeTSVolRequest volRequestFront(
                    strike,  // no shift for front strike
                    dfwSpecial->valueDate,
                    fwdDate,
                    false);  
                
                CVolProcessedBSSP volBSFront(assetNew->getProcessedVol(&volRequestFront));
                
                double varianceToFwdDate = scale * volBSFront->CalcVar(dfwSpecial->valueDate, fwdDate);
                return varianceToMat - varianceToFwdDate;
            }
            else
            {
                double t1 = volBS->calcTradingTime(fwdDate, maturity);
                double t2 = volBS->calcTradingTime(dfwSpecial->valueDate, maturity);

                // note:  we are interpolating vol at maturity.  (Since the 2 maturities
                // straddling the ex-div - 1 date may have different strike conventions).
                // Hence we have to convert variance to maturity to fwd variance:

                return varianceToMat * t1 / t2;
            }
        }

};   

/** Implementation of ClosedForm::IntoProduct interface */
NumericalIntegrationLN::IProduct* DFWSpecial::createProduct(
    NumericalIntegrationLN* model) const{

//    NumericalIntegrationLN::IProduct* divadjFL = DividendAdjustedFL::createProduct(model);
//    return new DFWSpecialClosedForm(model, this, divadjFL);
    
    return new AdjustedStrikeNumerical(const_cast<DFWSpecial* >(this));    
}


double DFWSpecial::specialAssumedAmount() const

{
    double specialAssumedAmount = 0.0;
    bool found = false;
    for (int i = 0; i < assumedDivs.size(); i ++)
    {
            if (assumedDivs[i].date.equals(specialDivDate, false)) // ignore time
            {
                specialAssumedAmount = assumedDivs[i].amount;
                found = true;
                break;
            }
    }

    if (found)
    {
        return specialAssumedAmount;
    }
    else
    {
        throw ModelException("DFWSpecial::specialAssumedAmount", 
            "Special Div Date " + specialDivDate.toString() + " cannot be found in the assumed divs!");
    }
}

double DFWSpecial::specialDivAmount() const

{
    double specialAmount = 0.0;
    bool found = false;
    DividendArray divarr = getEqDivs();
    for (int i = 0; i < eqDivArr.size(); i ++)
    {
            if (divarr[i].getExDate().equals(specialDivDate, false)) // ignore time
            {
                specialAmount = divarr[i].getDivAmount();
                found = true;
                break;
            }
    }

    if (found)
    {
        return specialAmount;
    }
    else
    {
        throw ModelException("DFWSpecial::specialDivAmount", 
            "Special Div Date " + specialDivDate.toString() + " cannot be found in the equity's divs!");
    }
}

DateTime DFWSpecial::specialDivPayDate() const

{
    DateTime specialPayDate(0,0);
    bool found = false;
    DividendArray divarr = getEqDivs();
    for (int i = 0; i < eqDivArr.size(); i ++)
    {
            if (divarr[i].getExDate().equals(specialDivDate, false)) // ignore time
            {
                specialPayDate = divarr[i].getPayDate();
                found = true;
                break;
            }
    }

    if (found)
    {
        return specialPayDate;
    }
    else
    {
        throw ModelException("DFWSpecial::specialDivAmount", 
            "Special Div Date " + specialDivDate.toString() + " cannot be found in the equity's divs!");
    }
}


//** set vol  */
void DFWSpecial::setVol(CAsset* asset, double vol) const
{
    EquityBase* eq = dynamic_cast<EquityBase*>(asset);
    VolLevel shift(vol);
    VolLevel::Shift* tweak = dynamic_cast<VolLevel::Shift*>(const_cast<CVolBase*>(eq->getVol().get()));
    tweak->sensShift(&shift);
}

CAssetSP DFWSpecial::makeFlatVolAsset(CAssetSP asset) const
{
    // first compute vol interpolated at maturity at original strike:
    DateTime maturity = avgOut->getLastDate();
    
    LinearStrikeTSVolRequest volRequestMat(
        strike - strikeShiftForVolInterp,
        valueDate,
        maturity,
        false);  
    
    CVolProcessedBSSP volBS(asset->getProcessedVol(&volRequestMat));
    
    double vol = volBS->CalcVol(valueDate, maturity);            
    
    // now create an asset with flat vol, using vol interpolated at maturity
    CAssetSP assetFlatVol((copy(asset.get())));                
    setVol(assetFlatVol.get(), vol);

    return assetFlatVol;
}

int AdjustedStrikeNumerical::maxIter = 200;
double AdjustedStrikeNumerical::solvePrecision = .00001;

/************************************************/

class DividendAdjustedStrike : public DividendAdjustedFL
{
private:
    class Product : public CClosedFormLN::IProduct
    {
    public:
        Product( const DividendAdjustedStrike & inst ) :
          inst( const_cast< DividendAdjustedStrike & >( inst ) ) {}

    private:
        DividendAdjustedStrike & inst;

        /** implementation of CClosedFormLN::IProduct interface */
        void price(
            CClosedFormLN* model,
            Control*       control, 
            CResults*      results) const
        {
            static const string method = "DividendAdjustedStrike::Product::price";
            try
            {
                // keep pricing unless a date after mat date, for cash flow output
                if( inst.valueDate.getDate() > inst.matDate.getDate() )
                {
                    // if settled already, return 0.
                    if( inst.valueDate >= inst.instSettle->settles(inst.matDate, inst.asset.get()) )
                    {   
                        results->storePrice( 0., inst.discount->getCcy() );
                        return;
                    }
                }

                HolidayConstSP hols = AssetUtil::getHoliday( inst.asset.get() );
                DateTime valueDate = inst.valueDate;
                DateTime exDivDate = hols->addBusinessDays( valueDate, inst.adjDateOffset );

                // try to adjust the strike
                if( control && control->isPricing() )
                {
                    int size = inst.eqOrigDivArr.size();
                    for( int i = 0; i < size; ++i )
                    {
                        // whether to do strike adjustment
                        if( ! inst.eqOrigDivArr[ i ].getExDate().equals( exDivDate, false ) )
                            continue;

                        // if adjusted strike is not recorded yet
                        if( ! inst.dateInHist( exDivDate, inst.histStrikes ) )
                        {
                            // find adjusted strike
                            double adjStrike = inst.findAdjustedStrike( results, hols, i ); 
                            if( ! Maths::equals( inst.strike, adjStrike ) )
                            {
                                // record strike for sens calculations
                                if( ! inst.histStrikes.size() || inst.histStrikes.back().date.getDate() < exDivDate.getDate() )
                                    inst.histStrikes.push_back( CashFlow( exDivDate, adjStrike ) );

                                // set the strike from history if needed
                                inst.setHistStrike();

                                // record adjusted strike
                                OutputRequest * request;
                                if( request = control->requestsOutput( OutputRequest::ADJUSTED_STRIKE ) )
                                    results->storeRequestResult( request, adjStrike );

                                // record strike history
                                if( request = control->requestsOutput( OutputRequest::DRO_HIST_DIV_DIFF ) )
                                    results->storeRequestResult( request, CashFlowArraySP::attachToRef( &inst.histStrikes ) );
                                if( request = control->requestsOutput( OutputRequest::HIST_DIV_DIFF ) )
                                    results->storeRequestResult( request, adjStrike );
                                if( request = control->requestsOutput( OutputRequest::HIST_DIV_DIFF_ANSI_DATE ) )
                                {
                                    DateTime::MonthDayYear mdy = exDivDate.toMDY();
                                    double ansiDate = (double)( mdy.day + 100 * mdy.month + 10000 * mdy.year );
                                    results->storeRequestResult( request, ansiDate );
                                }   
                            }
                        }
                        break;
                    }
                }

                // build an average spot inst and calc price
                model->Price( inst.buildAvgSpot( inst.avgOut, inst.strike ).get(), control, results );
                double price = results->retrievePrice();

                // store the price
                results->storePrice( price, inst.discount->getCcy() );
            }
            catch( exception & e )
            {
                throw ModelException( e, method );
            }
        }
    };

    static const int maxIter;
    static const double strikePrecision;
    static const double pricePrecision;

    double findAdjustedStrike(
        CResults * results,
        const HolidayConstSP & hols,
        int divIndex ) const
    {
        // check if need moving to ex div date
        if( adjDateOffset && adjOnExDate )
        {
            // copy instrument
            smartPtr< DividendAdjustedStrike > inst( copy( this ) );

            // move to ex div date
            Theta theta( adjDateOffset, HolidaySP::constCast( hols ) );
            theta.applyScenario( inst );

            // set adjOnExDate to false as instrument is already moved to ex div date
            inst->adjOnExDate = false;

            // call itself again
            return inst->findAdjustedStrike( results, hols, divIndex );
        }

        CResults res;

        // get equity name
        OutputNameSP equityName( new OutputName( asset->getTrueName() ) );

        // discount factor from pay date to ex-div date
        double df = discount->pv( valueDate, eqOrigDivArr[ divIndex ].getPayDate() );

        // current spot
        double spot = asset->getSpot();

        // amount of actual dividend alredy reflected in the spot
        double actualDivPaid = adjDateOffset > 0 ? 0. : eqOrigDivArr[ divIndex ].getDivAmount();

        // build AverageSpot instrument
        AverageSP avgInst( buildAvgSpot( avgOut, strike ) );
        const IAverageSpot & avgSpot = dynamic_cast< const IAverageSpot & >( *avgInst );

        // record variance at original strike
        double variance = avgSpot.calcAvgVariance(
            fixVolInterpAtStrike ? strike - strikeShiftForVolInterp : strike );

        // set equity spot to ( spot + pv( actual div paid - assumed div ) )
        SpotLevel assumedLevel( spot + df * ( actualDivPaid - assumedDivs[ divIndex ].amount ) );
        assumedLevel.findAndShift( avgInst, equityName );

        // calculate assumed price
        avgSpot.priceLN( 0, &res, strike, variance );
        double assumedPrice = res.retrievePrice();
        results->storeGreek(
            CDoubleSP( CDouble::create( assumedPrice ) ),
            Results::DEBUG_PACKET, 
            OutputNameConstSP( new OutputName("DBG_ASSUMED_PRICE_AFTER_DIV") ) );

        // set equity spot to ( spot + pv( actual div paid - actual div ) )
        SpotLevel actualLevel( spot + df * ( actualDivPaid - eqOrigDivArr[ divIndex ].getDivAmount() ) );
        actualLevel.findAndShift( avgInst, equityName );

        // calculate actual price
        avgSpot.priceLN( 0, &res, strike, variance );
        double actualPrice = res.retrievePrice();
        results->storeGreek(
            CDoubleSP( CDouble::create( actualPrice ) ),
            Results::DEBUG_PACKET, 
            OutputNameConstSP( new OutputName("DBG_ACTUAL_PRICE_AFTER_DIV") ) );

        // check for 'UP' / 'DOWN' only adjustment
        if( actualPrice > assumedPrice && divadjType == "UP" ||
            actualPrice < assumedPrice && divadjType == "DOWN" )
            return strike;

        // pick lo and hi levels for strike
        double diff = ::fabs( assumedDivs[ divIndex ].amount - eqOrigDivArr[ divIndex ].getDivAmount() );
        double loStrike = Maths::max( strike - 100. * diff, 0. );
        results->storeGreek(
            CDoubleSP( CDouble::create( loStrike ) ),
            Results::DEBUG_PACKET, 
            OutputNameConstSP( new OutputName("DBG_LO_STRIKE") ) );
        double hiStrike = strike + 100. * diff;
        results->storeGreek(
            CDoubleSP( CDouble::create( hiStrike ) ),
            Results::DEBUG_PACKET, 
            OutputNameConstSP( new OutputName("DBG_HI_STRIKE") ) );

        // find sign of correlation between strike and price
        double newStrike = hiStrike;
        avgSpot.priceLN( 0, &res, newStrike,
            fixVolInterpAtStrike ? variance : avgSpot.calcAvgVariance( newStrike ) );
        double price = res.retrievePrice();
        bool posCorr = price > actualPrice;

        // solve for strike using bi-section method
        int iter = 0;
        while( ++iter <= maxIter &&
            ( hiStrike - loStrike > strikePrecision || ::fabs( price - assumedPrice ) > pricePrecision ) )
        {
            newStrike = ( hiStrike + loStrike ) / 2.;
            avgSpot.priceLN( 0, &res, newStrike,
                fixVolInterpAtStrike ? variance : avgSpot.calcAvgVariance( newStrike ) );
            price = res.retrievePrice();

            if( posCorr ^ ( price > assumedPrice ) )
                loStrike = newStrike;
            else
                hiStrike = newStrike;
        }

        results->storeGreek(
            CIntSP( CInt::create( iter ) ),
            Results::DEBUG_PACKET, 
            OutputNameConstSP( new OutputName("DBG_ITER") ) );

        if( iter > maxIter )
        {
            throw ModelException( "DividendAdjustedStrike::Product::FindAdjustedStrike",
                "Solution cannot be found - exceeded maximum number of iterations " + Format::toString( maxIter ) );
        }

        results->storeGreek(
            CDoubleSP( CDouble::create( newStrike ) ),
            Results::DEBUG_PACKET, 
            OutputNameConstSP( new OutputName("DBG_ADJUSTED_STRIKE") ) );

        return newStrike;
    }

    void setHistStrike()
    {
        // find historical strike
        int size = histStrikes.size();
        for( int i = 0; i < size; ++i )
        {
            if( histStrikes[i].date.getDate() <= valueDate.getDate() )
                strike = histStrikes[i].amount;
        }
    }

public:
    static CClassConstSP const TYPE;

    /** instrument validation override */
    virtual void Validate()
    {
        // set ZERO_MU mode
        flowThruType = ZERO_MU;

        if( ! numExDivDates.size() )
        {
            // set 1 div per period
            int size = assumedDivs.size();
            numExDivDates.resize( size );
            for( int i = 0; i < size; ++i )
            {
                numExDivDates[i].date = assumedDivs[i].date;
                numExDivDates[i].amount = 1.;
            }
        }

        // call parent for validation
        DividendAdjustedFL::Validate();
    }

    /** rolls value date and sets initial spot for Theta,
        return true if sub objects need to be tweaked */
    virtual bool sensShift( Theta * shift )
    {
        static const string method = "DividendAdjustedStrike::sensShift";
        try
        {
            // calculate the new date
            DateTime newDate = shift->rollDate( valueDate );

            // If fwd start date falls between value date and theta date 
            // must multiply the percentage strike by spot (theta) or 
            // fwd at start date (theta forward spot) to convert to an 
            // absolute strike, as will no longer be forward starting.
            // Also have to set initial spot to same level.

            if( fwdStarting && valueDate <= startDate && startDate <= newDate )
            {
                fwdStarting = false;
                initialSpot = asset->getThetaSpotOnDate(shift, startDate);
                strike *= initialSpot;
            }  

            // fill in any samples
            avgOut->roll( asset.get(), valueDate, newDate, ! shift->useAssetFwds() );

            // set new value date 
            valueDate = newDate;

            // store actual divs if not stored yet
            if( ! eqOrigDivArr.size() )
                eqOrigDivArr = getEqDivs();
            if( ! eqDivArr.size() )
                eqDivArr = eqOrigDivArr;

            int size = eqOrigDivArr.size();

            // validate number of assumed divs
            if( assumedDivs.size() != size )
                throw ModelException( method, "Number of assumed divs must match number of actual divs in the period" );

            int lastHistStrikeDate = histStrikes.size() ? histStrikes.back().date.getDate() : 0;

            // set past assumed divs equal to actual divs
            for( int i = 0; i < size; ++i )
            {
                // validate assumed div dates
                if( ! assumedDivs[ i ].date.equals( eqOrigDivArr[ i ].getExDate(), false ) )
                    throw ModelException( method, "Assumed div dates must match actual ex div dates" );

                if( assumedDivs[ i ].date.getDate() > newDate.getDate() )
                    break;

                // replace assumed div amount only if strike has been adjusted already
                if( assumedDivs[ i ].date.getDate() <= lastHistStrikeDate )
                    assumedDivs[ i ].amount = eqOrigDivArr[ i ].getDivAmount();
            }

            // substitute actual divs with assumed divs
            substituteEqDiv( asset.getSP() );

            // set the strike from history if needed
            setHistStrike();
        }
        catch( exception & e )
        {
            throw ModelException( e, method );
        }    
        return true; // our components have theta type sensitivity
    }

    /** implementation of ClosedForm::IntoProduct interface */
    virtual CClosedFormLN::IProduct* createProduct(CClosedFormLN* model) const
    {
        return new Product( *this );
    }

protected:
    CashFlowArray histStrikes;
    int adjDateOffset;
    bool adjOnExDate;
    bool fixVolInterpAtStrike;
    double strikeShiftForVolInterp;

    DividendAdjustedStrike( CClassConstSP const type = TYPE ) :
        DividendAdjustedFL( type ),
        adjDateOffset( 1 ),
        adjOnExDate( true ),
        fixVolInterpAtStrike( true ),
        strikeShiftForVolInterp( 0. )
    {}

public:
    // use DividendAdjustedBase as superclass to avoid some mandatory fields declared in DividendAdjustedFL
    static void load(CClassSP& clazz)
    {
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(DividendAdjustedStrike, clazz);
        EMPTY_SHELL_METHOD(create);

        SUPERCLASS(DividendAdjBase);
        IMPLEMENTS(CClosedFormLN::IIntoProduct);

        // base fields
        FIELD(isCall, "Is it a call option");
        FIELD(avgOut, "Average-Out samples list");
        FIELD(divadjType, "Adjust strike only if price to be moved UP/DOWN on exDivDate, default NONE");
        FIELD_MAKE_OPTIONAL(divadjType);

        // own filed
        FIELD(histStrikes, "Historical strikes, last one is the current strike");
        FIELD_MAKE_OPTIONAL(histStrikes);
        FIELD(adjDateOffset, "Perform strike adjustment on (exDivDate - adjDateOffset), default 1");
        FIELD_MAKE_OPTIONAL(adjDateOffset);
        FIELD(adjOnExDate, "Specifies whether to perform strike adjustment as if on exDivDate even if adjDateOffset != 0, default TRUE");
        FIELD_MAKE_OPTIONAL(adjOnExDate);
        FIELD(fixVolInterpAtStrike, "Whether to fix volatility interpolated at (strike - strikeShiftForVolInterp) for fair value calculation, default TRUE");
        FIELD_MAKE_OPTIONAL(fixVolInterpAtStrike);
        FIELD(strikeShiftForVolInterp, "Interpolate volatility for fair value calculation at (strike - strikeShiftForVolInterp), default 0.0");
        FIELD_MAKE_OPTIONAL(strikeShiftForVolInterp);

        // transient
        FIELD_NO_DESC(eqOrigDivArr);
        FIELD_MAKE_TRANSIENT(eqOrigDivArr);
    }

    static IObject* create() { return new DividendAdjustedStrike(); }
};

const int DividendAdjustedStrike::maxIter = 200;
const double DividendAdjustedStrike::strikePrecision = .00001;
const double DividendAdjustedStrike::pricePrecision = .0000001;

CClassConstSP const DividendAdjustedStrike::TYPE = CClass::registerClassLoadMethod(
    "DividendAdjustedStrike", typeid(DividendAdjustedStrike), DividendAdjustedStrike::load);
   
bool DividendAdjFlowThroughLoad()
{
    return DividendAdjustedFL::TYPE && DividendAdjustedStrike::TYPE && true; // && true to shut the compiler up
}

DRLIB_END_NAMESPACE
