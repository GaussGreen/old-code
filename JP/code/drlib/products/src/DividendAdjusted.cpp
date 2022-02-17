//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DividendAdjusted.cpp
//
//   Description : Dividend adjusted option instrument
//
//   Date        : 11/20/2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/DividendAdjusted.hpp"
#include "edginc/Maths.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/ProtEquity.hpp"
#include "edginc/StruckEquity.hpp"
#include "edginc/Actual365F.hpp"

#if defined(_MSC_VER)
#pragma optimize( "g", off )
#endif


DRLIB_BEGIN_NAMESPACE

///////// local helpers /////////////////////

// count number of divs between dates
static int countDivs(const DividendArray& eqDiv, const DateTime& date1, 
                     const DateTime& date2, bool inclDate2)
{
    int count = 0;
    for (int j = 0; j<eqDiv.size(); j++)
    {
        if((eqDiv[j].getExDate() >= date1) &&
           (eqDiv[j].getExDate() < date2) ||
           (inclDate2 && eqDiv[j].getExDate() == date2))

            count++;
    }
    return count;
}

///////// local helpers end /////////////////////

////////////////// more helpers ///////////////

// write a requested double to results
    // override base
void DividendAdjBase::validatePop2Object()
{
    if (assumedDivs.size() == 0)
    {
        throw ModelException("DividendAdjBase::Validate", "Cannot have empty assumed div schedule!");
    }

}



void DividendAdjBase::processOutputRequest(const string&   outputRequest,
                                     const double    value,      
                                     Control*        control, 
                                     CResults*       results)
{
    OutputRequest* request = NULL;
    if ((request = 
            control->requestsOutput(outputRequest))){
                results->storeRequestResult(request, value);
    }   


}


// returns the difference between a div and an assumed div, according to our convention of multiple divs / qtr.
double DividendAdjBase::computeDivDiff(
                             int    numDvPd,   // number of divs in the period
                             double amtDvPd,   // total dollar amount of divs in the period
                             double amtDv,     // dollar amount of the div
                             double amtAssm)   // dollar amount of the assumed divs in the period
{
    if (Maths::isZero(amtDvPd))
    {
        return Maths::isZero(amtAssm)? 0.0 : -amtAssm/(double)numDvPd;
    }
    else
    {
        return (amtDv)*(1.0 - amtAssm/amtDvPd);
    }
}

    // returns the difference between a div and an assumed div, according to our convention of multiple divs / qtr.
double DividendAdjBase::computeDivDiff(const Dividend& exDiv)
{
    const DateTime& exDate = exDiv.getExDate();
    DividendArray eqDivArr = getEqDivs();
    
    int iPd = 0;
    while ( iPd < assumedDivs.size() - 1 && assumedDivs[iPd+1].date.getDate() < exDate.getDate() ) { 
        iPd++;
    }
    // note: if exDate is beyond all assumed div dates, we will correctly get the last period.
    
    double amtDvPd = 0.0;
    int lastPd = assumedDivs.size() - 1;
    int iEx = 0;
    for ( iEx = 0; iEx < eqDivArr.size(); iEx++ )
    {
        if ( (iPd == lastPd && eqDivArr[iEx].getExDate().getDate() >= assumedDivs[iPd].date.getDate()) ||
             (iPd != lastPd && eqDivArr[iEx].getExDate().getDate() >= assumedDivs[iPd].date.getDate() &&    
                               eqDivArr[iEx].getExDate().getDate() < assumedDivs[iPd+1].date.getDate()) )
        {
            amtDvPd += div2dollar(eqDivArr[iEx], asset.get());
        }
    }

    return computeDivDiff((int)numExDivDates[iPd].amount, amtDvPd, div2dollar(exDiv, asset.get()), assumedDivs[iPd].amount);
}



// util function to calc proportional strike weight adjustments.
// first sum to each sub date and the 1-sum.
void DividendAdjBase::propAdjust(const DateTimeArray& periodDates, const DateTimeArray& subDates,
                       const DoubleArray& subValues, DoubleArray& prop)
{
    int i,j;
    double sum;

    prop.resize(periodDates.size());
    for (i = 0;i < periodDates.size();i++)
    {
        sum = 0.0;
        for (j = 0;j < subDates.size();j++)
        {  
            if (subDates[j] < periodDates[i])
                sum += subValues[j];
        }
        prop[i] = 1.0 - sum;
    }

}
//get historic data
double DividendAdjBase::getHist(const CashFlowArray& histValues, DateTime date)
{
    for (int j = 0;j < histValues.size();j++)
    {
        if (histValues[j].date.equals(date, false))
        {
            return histValues[j].amount;
        }
    }
    
    throw ModelException("DividendAdjBase::getHist",
                         "Date " + date.toString() + 
                         " not found in history data");

    return 0.0; // never reach here
}

// get growth factor from date1 to date2
double DividendAdjBase::growthFactor(const DateTime& date1, const DateTime& date2, const DividendAdjBase* inst)
{
    // for protected get stock growth factor
    if (ProtEquity::TYPE->isInstance(inst->asset.get()))
    {
        const ProtEquity* eq = static_cast<const ProtEquity*>(inst->asset.get());
        return (1.0/eq->getYC()->pv(date1, date2));    
    }
    else   
        return (1.0/inst->discount->pv(date1, date2));
}

// get zero rate from valueDate to date
double DividendAdjBase::getZeroRate(const DateTime& date, const DividendAdjBase* inst)
{
    if (ProtEquity::TYPE->isInstance(inst->asset.get()))
    {
        const ProtEquity* eq = static_cast<const ProtEquity*>(inst->asset.get());
        return eq->getYC()->zero(date);    
    }
    else   
        return inst->discount->zero(date);
}

// get fwd zero rate from date1 to date2
double DividendAdjBase::getFwdZeroRate(const DateTime& date1, 
                                    const DateTime& date2, 
                                    const DividendAdjBase* inst)
{
    DayCountConventionSP dcc(DayCountConventionSP(new Actual365F()));

    if (ProtEquity::TYPE->isInstance(inst->asset.get()))
    {
        const ProtEquity* eq = static_cast<const ProtEquity*>(inst->asset.get());
        return eq->getYC()->fwd(date1, date2, dcc.get(), CompoundBasis::CONTINUOUS);    
    }
    else   
        return inst->discount->fwd(date1, date2, dcc.get(), CompoundBasis::CONTINUOUS);   
}

// get spot fx rate for the date
double DividendAdjBase::getFXValue(const DateTime& date, const CAsset* struckAsset)
{

    if (StruckEquity::TYPE->isInstance(struckAsset))
    {
        return static_cast<const StruckEquity*>(struckAsset)->getFX().get()->fwdValue(date);
    }
    else 
        return 1.0;
}

double DividendAdjBase::div2dollar(const Dividend& div, const CAsset* asset)
{
    const string method = "DividendAdjBase::div2dollar";

    // see if asset can return its equity.  
    const IHaveEquity* hasEquity = dynamic_cast<const IHaveEquity*>(asset);
    if (!hasEquity)
    {
        throw ModelException(method, "couldn't cast asset to IHaveEquity");
    }

    if(div.getDivType() == Dividend::AMOUNT)
    {
        return div.getDivAmount();
    }
    else
    {
        DateTime beforeExDivDate(div.getExDate().getDate(),DateTime::BEFORE_EX_DIV_TIME);
        return div.getDivAmount() * hasEquity->getEquity()->fwdValue(beforeExDivDate); 
    }    
}

// calculate sub total using data from subDates & subValues, summed according to periodDates.
// last date of periodDates is inputted separately to allow inclusive/exclusive counting
void DividendAdjBase::subTotal(const DateTimeArray& periodDates, const DateTime& lastDate, bool lastDateInclusive, 
                     const DateTimeArray& subDates, const DoubleArray& subValues, DoubleArray& total)
{
    int i,j;
    int curr = 0;
    double sum;

    total.resize(periodDates.size());
    for (i = 0; i<periodDates.size()-1; i++)
    {
        sum = 0.0;
        for (j = curr; j<subDates.size(); j++)
        {  
            if (subDates[j] >= periodDates[i] && subDates[j] < periodDates[i+1])
                sum += subValues[j];
            if (subDates[j] >= periodDates[i+1])
            {
                curr = j;
                break;
            }
        }
        total[i] = sum;
    }

    // do last date separately to allow inclusive
    sum = 0.0;
    for (j = 0; j<subDates.size(); j++)
    {  
        if ((subDates[j] >= periodDates[periodDates.size()-1] && subDates[j]<lastDate)
            || (subDates[j] == lastDate && lastDateInclusive))
            sum += subValues[j];
    }
    total[periodDates.size()-1] = sum;
}

/** finds corresponding pay dates for the supplied ex dates, using the supplied divArray.  
    If the pay date cannot be found, uses ex date. */
void DividendAdjBase::lookUpPayDates(const DateTimeArray& exDates, const DividendArray& divArray, DateTimeArray& foundPayDates)
{
    foundPayDates.resize(exDates.size());
    for (int i = 0; i < exDates.size(); i++)
    {
        if (! findDivPayDate(exDates[i], divArray, foundPayDates[i]) )
        {    // couldn't find the pay date
            foundPayDates[i] = exDates[i];
        }
            
    }
}

/** returns true if ex date is in the div array, whence foundPayDate is assigned to exDate's pay date. 
            false otherwise, with foundPayDate set to 0 */
bool DividendAdjBase::findDivPayDate(const DateTime& exDate, const DividendArray& divArray, DateTime& foundPayDate)
{
    for (int i = 0; i < divArray.size(); i++)
    {
        if (divArray[i].getExDate().getDate() == exDate.getDate())
        {    
            foundPayDate = divArray[i].getPayDate();
            return true;
        }
            
    }

    foundPayDate = DateTime(0,0);
    return false;
}

////////////////////////////////////////////////////////////
///////////////////// DividendAdjusted class ////////////////
////////////////////////////////////////////////////////////
class DividendAdjusted: public DividendAdjBase,
                        public CClosedFormLN::IIntoProduct
{
public:

    friend class DividendAdjustedClosedForm;

    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);

    static IObject* defaultDividendAdjusted(){
        return new DividendAdjusted();
    }

    // override base
    virtual void validatePop2Object();

    /** instrument validation override */
    virtual void Validate();

    /** Implementation of ClosedForm::IntoProduct interface */
    virtual CClosedFormLN::IProduct* createProduct(CClosedFormLN* model) const;

    // calculate avg out sample dates and weights
    SampleListSP calcOutSample() const;

    // calc the adjusted strikes
    void calcStrikes(double& avgStrike, double& absStrike, double& expectedStrike,
                     double& strikeAdjustmentToDate, double& currAdjustment) const;


protected:

    DividendAdjusted() :DividendAdjBase(TYPE),
                        adjustmentType(0),
                        useFwdValOfDivs(0),
                        useEquivYields(0) {};

    DividendAdjusted(const DividendAdjusted& rhs);
    DividendAdjusted& operator=(const DividendAdjusted& rhs);

    // below are for old interface
    int             isCall;
    int             adjustmentType;
    int             useFwdValOfDivs;
    int             useEquivYields;

    CashFlowArray   aveSchedOrig;
    DoubleArray     aveSchedOrigWeights;
};


///////////////////// DividendAdjBase class methods////////////////
void DividendAdjBase::Validate()
{
    int i, numOfDiv;

    static const string method = "DividendAdjBase::Validate";
    // just check the things that aren't/cannot be checked in 
    // validatePop2Object
    if (!asset){
        throw ModelException(method, "Asset is null");
    }
    if (!discount){
        throw ModelException(method, "Discount YC is null");
    }

    AssetUtil::assetCrossValidate(asset.get(),
                                  fwdStarting,
                                  startDate,
                                  valueDate,
                                  discount,
                                  this);
    
    if (FXAsset::TYPE->isInstance(asset.get())) {
        throw ModelException(method,
                             "Options on FX assets are not allowed yet");
    }

    int numavg = avgOut->getDates().size(); 
    if (numavg == 0)
    {
        throw ModelException(method, "Number of average out dates can not be zero");
    }

    // if fwd starting, can't be one contract
    if ( fwdStarting && oneContract)
    {
        throw ModelException(method, "Can't be forward starting and "
                             "one contract");
    }


    if (assumedDivs.size() == 0)
        throw ModelException(method,"number of the assumed dividends must be > 0.");

    
    if (assumedDivs.size() != numExDivDates.size())
        throw ModelException(method,"size of the assumed dividend periods must be the same as size of numExDivDates.");

    for (i = 0;i < assumedDivs.size(); i++)
    {
        // force time of div going to ex to when it should be
        numExDivDates[i].date = DateTime(numExDivDates[i].date.getDate(), 
                                       Dividend::DIVIDEND_EXDIV_TIME);
        // force time of div going to ex to when it should be
        assumedDivs[i].date = DateTime(assumedDivs[i].date.getDate(), 
                                      Dividend::DIVIDEND_EXDIV_TIME);
        if (assumedDivs[i].date >= avgOut->getLastDate())
            throw ModelException(method, "assumed dividend date "+assumedDivs[i].date.toString()+
                                 " can not be >= maturity date" + avgOut->getLastDate().toString());

        // check inputed number of dividends in the assumed period
        if (!(numExDivDates[i].date == assumedDivs[i].date)) {
            throw ModelException(method,
                                 Format::toString(i) + 
                                 "th period start date for assumed "
                                 "dividends does not match start date "
                                 "for numExDivDates input.");    
        }
    }

    // get divs, this is called again in each pricing call 
    eqDivArr = getEqDivs();

    // check number of dividends in the assumed periods matches market data
    for (i = 0;i < assumedDivs.size();i++)
    {
        const DateTime& start = assumedDivs[i].date;
        const DateTime& end = (i ==assumedDivs.size() - 1? avgOut->getLastDate() : assumedDivs[i+1].date);
        numOfDiv = countDivs(eqDivArr, start, end, i==assumedDivs.size() - 1);

        if (numOfDiv != (int)numExDivDates[i].amount) {
            throw ModelException(method,
                             "In period " + Format::toString(i) + " (which starts on  " +
                             start.toString() + " and ends on " + end.toString() + ", the contract assumed " + 
                             Format::toString((int)numExDivDates[i].amount) + 
                             " ex-div dates, but found " + Format::toString(numOfDiv) +
                              " in the equity.");    
        }
    }

    if (matDate.getDate() < avgOut->getLastDate().getDate())
        throw ModelException(method, "maturity date cannot be before last average out date.");
        
    if (!(ccyTreatment == CAsset::CCY_TREATMENT_STRUCK) && assumedDivsInStruckCcy)
        throw ModelException(method, "assumedDivsInStruckCcy can only be true for ccy struck case.");
    
    if (assumedDivsInStruckCcy)
    {
        int i;
        for (i =0; i<numExDivDates.size(); i++)
        {
            if ((int)numExDivDates[i].amount != 1)
                throw ModelException(method, "for assumedDivsInStruckCcy case, each div period must contain exactly one "
                " dividend. Period " + Format::toString(i) + " starting " + numExDivDates[i].date.toString() + " contains " +
                Format::toString((int)numExDivDates[i].amount) + " dividends.");
        }

        if(eqDivArr.size() != assumedDivs.size())
            throw ModelException(method, "for flowThruType=ZERO_MU, assumed dividend array must match"
                                 " the number of dividends in the period");

        for (i =0; i<eqDivArr.size(); i++)
        {
            if (eqDivArr[i].getExDate().getDate() != assumedDivs[i].date.getDate())
                throw ModelException(method, "for assumedDivsInStruckCcy case, equity dividend date " + 
                                     eqDivArr[i].getExDate().toString() +  " must match assumed dividend date "+
                                     assumedDivs[i].date.toString());
        }
    }
}

/** get equity dividends */
DividendArray DividendAdjBase::getEqDivs() const
{
    // extract all divs 
    DividendListSP divs = AssetUtil::getAllDivsBetweenDates(asset.get(),
                                                        assumedDivs[0].date.rollDate(-1),
                                                        avgOut->getLastDate());

    return divs.get()->getArray();
}

/** gets num div dates on the fly.  This allows proper treatment for MuSpecial */
CashFlowArraySP DividendAdjBase::getNumDivDates(const DividendArray& currentEqDivs) const
{
    CashFlowArraySP numDivDates(new CashFlowArray(0));
    
    for (int i = 0;i < assumedDivs.size();i++)
    {
        const DateTime& start = assumedDivs[i].date;
        const DateTime& end = (i ==assumedDivs.size() - 1? avgOut->getLastDate() : assumedDivs[i+1].date);
        int numDiv = countDivs(currentEqDivs, start, end, i==assumedDivs.size() - 1);
        numDivDates->push_back(CashFlow(numExDivDates[i].date, numDiv));
    }

    return numDivDates;
}

bool DividendAdjBase::hasYield() const
{
    const string method = "DividendAdjBase::hasYield";

    bool hasYield = false;


    for (int i = 0; i < eqDivArr.size(); i++) {
        if ((!(eqDivArr[i].getDivType() == Dividend::AMOUNT)) && 
            (!(eqDivArr[i].getExDate() > valueDate))) {
            throw ModelException(method, "historic dividend (" + 
                                 eqDivArr[i].getExDate().toString() + 
                                 ") must be dollar amount (today is " +
                                 valueDate.toString() + ")");
        }
        if(eqDivArr[i].getDivType() == Dividend::PERCENT) 
            hasYield = true;

        if(eqDivArr[i].getDivType() == Dividend::CONTINUOUS)
            throw ModelException(method, "CONTINUOUS dividends not supported");
    } 

    return hasYield;
}

// build a average spot instrument
AverageSP DividendAdjBase::buildAvgSpot(SampleListConstSP outSample, double strikeAvg) const
{
    AverageSP avg(Average::makeAvgSpot(isCall,
                                       matDate,
                                       strikeAvg,
                                       outSample.get(),
                                       instSettle.get(),
                                       premiumSettle.get(),
                                       asset.get(),
                                       ccyTreatment,
                                       discount.get(),
                                       valueDate,
                                       fwdStarting,
                                       startDate,
                                       oneContract,
                                       notional,
                                       initialSpot));
    
    return avg;
}

void DividendAdjBase::GetMarket(const IModel*          model, 
                                 const CMarketDataSP    market)
{
    const string method = "DividendAdjBase::GetMarket";

    market->GetReferenceDate(valueDate);

    if (asset.usingCache())
    {
        CAsset::getAssetMarketData(model, market.get(), ccyTreatment, 
                                   discount, asset);
    }

    discount.getData(model, market);

    instSettle->getMarket(model, market.get());
    if( !instSettle )
    {
        throw ModelException(method, "Instrument settlement is NULL");
    }

    if (premiumSettle.get())
    {
        premiumSettle->getMarket(model, market.get());
    }
}

DateTime DividendAdjBase::getValueDate() const
{
    return valueDate;
}

/** when to stop tweaking */
DateTime DividendAdjBase::endDate(const Sensitivity* sensControl) const {
    
    DateTime instEnd  = instSettle->settles(matDate, asset.get());
    DateTime assetEnd = asset->settleDate(matDate);
    DateTime end = instEnd.isGreater(assetEnd) ? instEnd : assetEnd;
    return end;
}

// calculate FX spot and growth factor for eq div dates
void DividendAdjBase::calcFXandDFdiv(DoubleArray& fx, DoubleArray& gf) const
{
    Actual365F dayCount;

    DateTime  matPayDate = avgOut->getLastDate(); // to do: should be matDate + settle ?
    if (usePayDate) // true only for new FL, to do old div adj should use this as well
        matPayDate = instSettle->settles(matDate, asset.get());

    for (int j = 0; j < eqDivArr.size(); j++)
    {
        gf[j] = 1.0;
        DateTime exDate = eqDivArr[j].getExDate();
        if (exDate.getDate() < valueDate.getDate()) 
        {// use history
            if (StruckEquity::TYPE->isInstance(asset.get()))
                fx[j] = getHist(histFXRates, exDate);
            else
                fx[j] = 1.0;
            // for growth factor
            if (fwdAdj)
            {
                double zRate = getHist(histDiscRates, exDate);
                // to do : should use pay date with fwd zero rate ! but this will change numbers
                const DateTime& cashPayDate = (usePayDate ? eqDivArr[j].getPayDate() : exDate);
                gf[j] = exp(log(zRate + 1.0)*dayCount.years(cashPayDate, matPayDate)); 
            }
        }
        else
        {
            // to do : should use pay date with fwd zero rate and fwd fx rate
            const DateTime& cashPayDate = (usePayDate ? eqDivArr[j].getPayDate() : exDate);
            fx[j] = getFXValue(cashPayDate, asset.get());
            if (fwdAdj)
                gf[j] = growthFactor(cashPayDate, matPayDate, this);
        }
    }
}

// calculate FX spot and growth factor for assumed dates
void DividendAdjBase::calcFXandDFassumed(const DoubleArray& periodSum,
                                       DoubleArray& fx,
                                       DoubleArray& gf) const
{
    Actual365F dayCount;

    // get on the fly num ex div dates to properly handle mu special
    CashFlowArraySP numExDivDatesInUseSP = getNumDivDates(getEqDivs());
    CashFlowArray& numExDivDatesInUse = *numExDivDatesInUseSP;

    const DateTime& lastOutDate = avgOut->getLastDate();
    DateTime  matPayDate = avgOut->getLastDate(); // to do: should be matDate + settle ?
    if (usePayDate) // true only for new FL, to do old div adj should use this as well
        matPayDate = instSettle->settles(matDate, asset.get());

    // process assumed div dates
    for (int j = 0; j < assumedDivs.size(); j++)
    {
        const DateTime& cashPayDate = (j==assumedDivs.size()-1 ? lastOutDate :
                                        assumedDivs[j+1].date.rollDate(-1));

        gf[j] = 1.0;
        if ( Maths::isZero(numExDivDatesInUse[j].amount) && cashPayDate.getDate() < valueDate.getDate())
        {//if period is empty and need history (note we count the divs, rather than compute sum)
            if (StruckEquity::TYPE->isInstance(asset.get()))
                fx[j] = getHist(histFXRates, cashPayDate);
            else
                fx[j] = 1.0;
            // for growth factor
            if (fwdAdj)
            {
                double zRate = getHist(histDiscRates, cashPayDate);
                gf[j] = exp(log(zRate + 1.0)*dayCount.years(cashPayDate, matPayDate));
            }
        }
        else
        {
            fx[j] = getFXValue(cashPayDate, asset.get());
            if (fwdAdj)
                gf[j] = growthFactor(cashPayDate, matPayDate, this);
        }
    }
}

//////////////////////////////////////////////////
/** private product class */
//////////////////////////////////////////
class DividendAdjustedClosedForm: public CClosedFormLN::IProduct{
private:
    IModel*                  model;
    const DividendAdjusted*   divadjInst; // a reference

public:
    DividendAdjustedClosedForm(IModel* model, const DividendAdjusted*  inst): 
        model(model), divadjInst(inst){}

    void price(CClosedFormLN*   model,
               Control*        control, 
               CResults*       results) const;
};

// price entry point
void DividendAdjustedClosedForm::price(CClosedFormLN*   model,
                                       Control*        control, 
                                       CResults*       results) const
{
    static const string method = "DividendAdjustedClosedForm::price";
    try {

        // if assumedDivsInStruckCcy, adjust assumed divs.
        if (divadjInst->assumedDivsInStruckCcy)
        {
            // divided by fx fwd at ex date to convert assumed divs to base ccy
            const_cast<DividendAdjusted*>(divadjInst)->assumedDivs2BaseCcy(); 

            if (control && control->isPricing()) // write assumed divs in base ccy to debug packet
            {    
                // output ex-date cashflows. (for debug)
                OutputNameConstSP assumedDivsBaseOut(new OutputName("assumedDivsBase"));
                CashFlowArraySP out = CashFlowArraySP(new CashFlowArray(divadjInst->assumedDivs));
                results->storeGreek(out, Results::DEBUG_PACKET, assumedDivsBaseOut);
    
            }
        }

        DateTime settlementDate = divadjInst->instSettle->settles(divadjInst->matDate, divadjInst->asset.get());
        const DateTime&  matPayDate = divadjInst->avgOut->getLastDate(); // should be settlementDate ?

        if (divadjInst->valueDate >= settlementDate)
        {   // settled already
            results->storePrice(0.0, divadjInst->discount->getCcy());
            return;
        }

        // get divs
        const_cast<DividendAdjusted*>(divadjInst)->eqDivArr = divadjInst->getEqDivs();

        // for debugging
#if 0
        {
            static i = 0;
            string dawFile = "c:\\tmp\\d-dbgD" + Format::toString(i) + ".xml";
            XMLWriter xml(dawFile.c_str());
            divadjInst->write("OBJECT", &xml);
            string dawFile2 = "c:\\tmp\\d-asset-dbgD" + Format::toString(i) + ".xml";
            XMLWriter xml2(dawFile2.c_str());
            divadjInst->asset.get()->write("OBJECT", &xml2);
            i++;
        }
#endif

        //debug
        double  strike; // strike to be passed to average model
        double  absStrike; // absolute strike  (=strike if not fwd starting)
        double  initStrike; // keeps track of initial strike
        double  expectedStrike; // expected absolute strike
        double  strikeAdjustmentToDate; // strike adjustment up to date
        double  currAdjustment; // adjustment for current period
        double  newScale = 1.0;  // scale factor for notional
        SampleListSP outSample; // averaging out sample

        // calc strike
        divadjInst->calcStrikes(strike, absStrike, expectedStrike, strikeAdjustmentToDate, currAdjustment);

        if (divadjInst->scaleType != "NONE")
        {
            // calc out samples
            outSample = divadjInst->avgOut;
            initStrike = divadjInst->strike;
            newScale = initStrike / strike;
        }

        else
        {
            // calc out samples
            outSample = divadjInst->calcOutSample(); 
        }

        // build an average spot inst and calc price
        model->Price(divadjInst->buildAvgSpot(outSample, strike).get(), control, results);
            
        // Do Notional adjustment if scale A
        double avgPrice = results->retrievePrice();
        results->storePrice(avgPrice * newScale, divadjInst->discount->getCcy());

        if (control && control->isPricing() ) {
            OutputRequest* request = NULL;

            if ((request = control->requestsOutput(OutputRequest::ADJUSTED_STRIKE))){
                results->storeRequestResult(request, absStrike);
            }

            if ((request = control->requestsOutput(OutputRequest::EXPECTED_STRIKE))){
                results->storeRequestResult(request, expectedStrike);
            }

            if ((request = control->requestsOutput(OutputRequest::EXPECTED_SCALE_FACTOR))){
                results->storeRequestResult(request, newScale);
            }

            if ((request = control->requestsOutput(OutputRequest::SCALE_FACTOR_UPTO_DATE))){
                results->storeRequestResult(request, strikeAdjustmentToDate);
            }

            // always output these useful outputs
            if ((request = control->requestsOutput(OutputRequest::STRIKE_ADJUSTMENT_UPTO_DATE))){
                if (divadjInst->scaleType == "STRIKE_DIVIDED") {
                    results->storeRequestResult(request, divadjInst->strike / strikeAdjustmentToDate);
                }
                else if (divadjInst->scaleType == "STRIKE_MULTIPLIED") {
                    results->storeRequestResult(request, divadjInst->strike * strikeAdjustmentToDate);
                }
                else {
                    results->storeRequestResult(request, divadjInst->strike - strikeAdjustmentToDate);
                }
            }
            
            if ((request = control->requestsOutput(OutputRequest::CURRENT_STRIKE_ADJUSTMENT))){
                results->storeRequestResult(request, currAdjustment);
            } 

            // record hist data if needed
            if( divadjInst->recordHist())
            {
                if (divadjInst->fwdAdj)
                {
//                     if ((request = 
//                         control->requestsOutput(OutputRequest::HIST_DISC_RATE_TO_MAT))){
//                        results->storeRequestResult(request, 
//                                            DividendAdjBase::getZeroRate(matPayDate, divadjInst));
//                    }   
                    
                    double zeroRate = DividendAdjBase::getZeroRate(matPayDate, divadjInst);

                    DividendAdjusted::processOutputRequest(OutputRequest::HIST_DISC_RATE_TO_MAT,
                                                           zeroRate,
                                                           control, 
                                                           results);
                   
                    DividendAdjusted::processOutputRequest(OutputRequest::DRO_HIST_DISC_RATE_TO_MAT,
                                                           zeroRate,
                                                           control, 
                                                           results);
                }
                
                if (StruckEquity::TYPE->isInstance(divadjInst->asset.get()))
                {
                    if ((request = 
                         control->requestsOutput(OutputRequest::HIST_SPOT_FX))){
                        results->storeRequestResult(request, DividendAdjBase::getFXValue(divadjInst->valueDate, 
                                                                                         divadjInst->asset.get()));
                    }   
                }
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


/** indicates whether VEGA_MATRIX is sensible for this instrument */
bool DividendAdjBase::avoidVegaMatrix(const IModel* model)
{
    /* this should possibly be false for local vol models etc. */
    return false;
}

/** returns all strikes on the vol surface to which 
    this instrument is sensitive */
DoubleArraySP DividendAdjBase::getSensitiveStrikes(OutputNameConstSP outputName,
                                                    const IModel*      model)
{
    return Average::getSensitiveStrikes(avgInst.get(), outputName, model);
}

/** Rolls the value date and sets initial spot if rolling over start date */
bool DividendAdjBase::sensShift(Theta* shift)
{
    DateTime newDate = shift->rollDate(valueDate);
    thetaHelper(shift, newDate);

    // roll today 
    valueDate = newDate;

    return true;
}

// if useFwd is true, record forward to fwdDate in setting fx and df */
void DividendAdjBase::thetaHelper(Theta*          shift, 
                                  const DateTime& newDate, 
                                  bool            isNewFL, 
                                  const DateTime& divPayDate)
{    
    DateTime  matPayDate = avgOut->getLastDate(); // this is old div adj way
    if (isNewFL) // this will be for new flow thru type
        matPayDate = instSettle->settles(matDate, asset.get());

    if (fwdStarting && newDate.isGreaterOrEqual(startDate) &&
        startDate.isGreaterOrEqual(valueDate))
    {
        fwdStarting = false;
        initialSpot = asset->getThetaSpotOnDate(shift, startDate);
        strike *= initialSpot;
    }

    avgOut->roll(asset.get(), valueDate, newDate, !shift->useAssetFwds());
    
    if (newDate == valueDate)
        return; 

    if (recordHist())
    {
        if (fwdAdj)
        {
            double zr = (isNewFL) ? getFwdZeroRate(divPayDate, matPayDate, this) : getZeroRate(matPayDate, this);
            histDiscRates.push_back(CashFlow(valueDate, zr));
        }

        if (StruckEquity::TYPE->isInstance(asset.get()))
            histFXRates.push_back( CashFlow(valueDate, getFXValue(isNewFL ? divPayDate : valueDate, asset.get())) );

    }

    // update for ex-dates after today
    DividendArray eqDvs = getEqDivs();
    bool pastNewDate = false;
    int  i = 0;
    while ( !pastNewDate && i < eqDvs.size() ) {
        const DateTime &exDate = eqDvs[i].getExDate();
        if ( exDate.getDate() > valueDate.getDate() && exDate.getDate() < newDate.getDate() ) {
            // update FX Rates
            if (ccyTreatment == CAsset::CCY_TREATMENT_STRUCK) {
                 histFXRates.push_back( CashFlow(exDate, 
                                getFXValue(exDate, asset.get())) );
            }
            
            // update Disc Rates
            if (fwdAdj) {
                double zr = getFwdZeroRate(exDate, matPayDate, this); // to do: for FL, need forward from payDate.
                histDiscRates.push_back(CashFlow(exDate, zr));
            }
        } else if ( exDate.isGreater(newDate) ) {
            pastNewDate = true; // don't need to look anymore, as this and remaining ex dates > new date.
        }
        ++i;
    }

    
    // special case: no divs in a period
    for( i=0; i<numExDivDates.size(); i++ ) {
        if (Maths::isZero(numExDivDates[i].amount)) {
            bool isLast = (i==numExDivDates.size()-1);
            const DateTime& endDate = (isLast ? avgOut->getLastDate() : assumedDivs[i+1].date);
            if (valueDate.getDate() < endDate.getDate() && newDate.getDate() >= endDate.getDate()) {
                // need to update if we roll over a future end of no div period
                const DateTime& cashPayDate = isLast ? avgOut->getLastDate() : assumedDivs[i+1].date.rollDate(-1);
                // updated FX Rates
                if (ccyTreatment == CAsset::CCY_TREATMENT_STRUCK) {
                     histFXRates.push_back( CashFlow(cashPayDate, getFXValue( cashPayDate, asset.get() ) ) );
                }
                // update Disc Rates
                if (fwdAdj) {
                    double zr = getFwdZeroRate(cashPayDate, matPayDate, this); // to do: for FL, need forward from payDate.
                    histDiscRates.push_back(CashFlow(cashPayDate, zr));
                }
            }
        }
    }
    
    // if notional adjusted type
    if (scaleType != "NONE") {

        double spot    = 0.0;
        bool   useSpot = !shift->useAssetFwds();

        if (useSpot) {
            spot = asset->getSpot();
        } 

        // fill spot historical sample point if needed
        bool pastRollDate = false;
        int  j = 0;
        while (!pastRollDate && j < spotBeforeExDates.size() ) {
            if ((spotBeforeExDates[j].date.isGreater(valueDate) && 
                 !spotBeforeExDates[j].date.isGreater(newDate)) ||
                (spotBeforeExDates[j].date.equals(valueDate)    && 
                 Maths::isZero(spotBeforeExDates[j].amount))) {
                if (useSpot) {
                    spotBeforeExDates[j].amount = spot;
                }
                else {
                    spotBeforeExDates[j].amount = asset->fwdValue(spotBeforeExDates[j].date);
                }
            } 
            else if (spotBeforeExDates[j].date.isGreater(newDate)) {
                pastRollDate = true;
            }
            ++j;
        }

        // fill spot assumed sample point if needed
        bool pastRollDateAsmd = false;
        int  k = 0;
        while (!pastRollDateAsmd && k < spotBeforeAssumedDates.size() ) {
            if ((spotBeforeAssumedDates[k].date.isGreater(valueDate) && 
                 !spotBeforeAssumedDates[k].date.isGreater(newDate)) ||
                (spotBeforeAssumedDates[k].date.equals(valueDate)    && 
                 Maths::isZero(spotBeforeAssumedDates[k].amount))) {
                if (useSpot) {
                    spotBeforeAssumedDates[k].amount = spot;
                }
                else {
                    spotBeforeAssumedDates[k].amount = asset->fwdValue(spotBeforeAssumedDates[k].date);
                }
            } 
            else if (spotBeforeAssumedDates[k].date.isGreater(newDate)) {
                pastRollDateAsmd = true;
            }
            ++k;
        }

    }

    return;
};

CSensControl* DividendAdjBase::AlterControl(const IModel*       modelParams,
                                             const CSensControl* sensControl) const
{
    SensControlPerName* alteredControl = NULL;
    if (Delta::TYPE->isInstance(sensControl)         &&
        CClosedFormLN::TYPE->isInstance(modelParams) )
    {
        const Delta* delta = 
            dynamic_cast<const Delta*>((IObject*)sensControl);
        ShiftSizeCollectorSP shiftSizeVisitor(new ShiftSizeCollector(
            delta,
            strike,
            fwdStarting?
            ShiftSizeCollector::FWD_START_ADJUSTMENT:
            ShiftSizeCollector::SPOT_START_ADJUSTMENT));
        asset->accept(shiftSizeVisitor.get());

        if ( Maths::isPositive(shiftSizeVisitor->getShiftSize()) )
        {
            alteredControl = new Delta(shiftSizeVisitor->getShiftSize());
            alteredControl->
                setMarketDataName(sensControl->getMarketDataName());
        }
    }

    return alteredControl;
}

/////////////// dividend adjusted class methods ////////
// calculate avg out sample dates and weights
///// additonal data mapping first, data mapping is just to maintain old interface
void DividendAdjusted::validatePop2Object()
{
    int numavg = aveSchedOrig.size();
    DateTimeArray dates(numavg);
    DoubleArray   values(numavg);

    for(int i = 0;i < numavg;i++)
    {
        dates[i] = aveSchedOrig[i].date;
        values[i] = aveSchedOrig[i].amount;
    }

    avgOut = SampleListSP(new SampleList(dates,values,aveSchedOrigWeights)); 

    DividendAdjBase::isCall = (isCall == 1);

    if (adjustmentType == 0)
        divadjType = "NONE";
    else if (adjustmentType == 1)
        divadjType = "UP";
    else
        divadjType = "DOWN";     

    fwdAdj = (useFwdValOfDivs == 1);

    // now call the parent's validatePop2Object
    DividendAdjBase::validatePop2Object();

}

//////////////////////
void DividendAdjusted::Validate()
{
    DividendAdjBase::Validate();

    // need this for vega matrix strikes
    double strikeAvg, absStrike, expectedStrike, dummy;
    calcStrikes(strikeAvg, absStrike, expectedStrike, dummy, dummy);
    SampleListSP outSample;
    if (scaleType != "NONE") {
        outSample = avgOut; 
    }
    else {
        outSample = calcOutSample(); 
    }
    avgInst = buildAvgSpot(outSample, strikeAvg);

    if((!(divadjType == "NONE")) && (useEquivYields == 1 || hasYield()))
        throw ModelException("DividendAdjBase::Validate()", "Dividend adjustment type must be NONE for div yields");
}

SampleListSP DividendAdjusted::calcOutSample() const
{
    const string method = "DividendAdjusted::calcOutSample";

    int k = 0;
    int j;

    // calc the number of dividend yields
    for (j = 0;j < eqDivArr.size();j++)
    {
        if ((eqDivArr[j].getExDate() <= valueDate) && (eqDivArr[j].getDivType() != Dividend::AMOUNT))
            throw ModelException(method,"past dividend must be amount (not yield)"); 
         
        if ((useEquivYields == 1 && eqDivArr[j].getExDate()>valueDate)
            || eqDivArr[j].getDivType() != Dividend::AMOUNT) 
        {
            k++; 
        }
    }

    // copy of average out schedule
    DateTimeArray sDates = avgOut.get()->getDates();
    DoubleArray sValues = avgOut.get()->getValues();
    DoubleArray sWeights = avgOut.get()->getWeights();

    if (k > 0)
    {// there are yield div to deal with
        DateTimeArray dates(k + 1);
        DoubleArray values(k + 1);
        DoubleArray weights(k + 1);

        DoubleArray fxDiv(eqDivArr.size());
        DoubleArray gfDiv(eqDivArr.size());
        calcFXandDFdiv(fxDiv,gfDiv); //get FX and growth factor

        k = 0;
        //adjusting weights
        const ProtEquity* protEq = dynamic_cast<const ProtEquity*>(asset.get()); 
        for (j = 0; j < eqDivArr.size(); j++)
        {
            if(eqDivArr[j].getDivType() != Dividend::AMOUNT)
            {
                dates[k] = eqDivArr[j].getExDate();
                values[k] = 0.0;
                weights[k] = gfDiv[j] * eqDivArr[j].getDivAmount()/(1.0 - eqDivArr[j].getDivAmount());
                k++;
            }
            else if (useEquivYields == 1 && eqDivArr[j].getExDate() > valueDate)
            {
                dates[k] = eqDivArr[j].getExDate();
                values[k] = 0.0;
                DateTime exDate(eqDivArr[j].getExDate().getDate(), DateTime::BEFORE_EX_DIV_TIME);
                double fwd = (protEq ? protEq->unadjustedFwdValue(exDate) : asset->fwdValue(exDate));

                DateTime timetmp(eqDivArr[j].getExDate().getDate(),DateTime::BEFORE_EX_DIV_TIME);
                weights[k] = gfDiv[j] * eqDivArr[j].getDivAmount()/(fwd/fxDiv[j] - eqDivArr[j].getDivAmount()); 
                k++;
            } 
        }
        // insert or combine weight adjustments with average sample array
        int n = 0;
        int numInsert = 0;
        for (j=0; j<avgOut.get()->getDates().size(); j++)
        {
            while (n<k && dates[n] < sDates[j+numInsert])
            {
                sDates.insert(sDates.begin()+j+numInsert, dates[n]);
                sValues.insert(sValues.begin()+j+numInsert, values[n]);
                sWeights.insert(sWeights.begin()+j+numInsert, weights[n]);
                n ++;
                numInsert ++;
            }
            if (n<k && dates[n] == sDates[j+numInsert])
            {
                sWeights[j+numInsert] += weights[n];
                n ++;
            }
        }
    }

    //new sample list
    return SampleListSP(new SampleList(sDates,sValues,sWeights));
}


// calc the adjusted strikes
void DividendAdjusted::calcStrikes(double& avgStrike, double& absStrike, double& expectedStrike,
    double& strikeAdjustmentToDate, double& currAdjustment) const
{
    // since assumedDivsInStruckCcy needs a new assumedDivs, must regenerate it for each price call
    // from original schedule
    CashFlowArraySP assumedDivsInUseSP;
    if (assumedDivsInStruckCcy)
    {
        assumedDivsInUseSP = assumedDivs2BaseCcy();
    }
    else
    {
        assumedDivsInUseSP = CashFlowArraySP(copy(&assumedDivs));
    }

    CashFlowArray& assumedDivsInUse = *assumedDivsInUseSP;

    const string method = "DividendAdjusted::calcStrikes";

    int i,j;

    DoubleArray divDollars(eqDivArr.size());
    DoubleArray fxDiv(eqDivArr.size());
    DoubleArray fxAsmd(assumedDivsInUse.size());
    DoubleArray gfDiv(eqDivArr.size());
    DoubleArray gfAsmd(assumedDivsInUse.size());

    // get dates lists first
    DateTimeArray eqExDates(eqDivArr.size());
    for (i=0; i<eqExDates.size(); i++)
        eqExDates[i] = eqDivArr[i].getExDate();

    DateTimeArray assumedDates(assumedDivsInUse.size());
    for (i=0; i<assumedDates.size(); i++)
        assumedDates[i] = assumedDivsInUse[i].date;

    // last avg out date
    const DateTime& lastOutDate = avgOut->getLastDate();

    // init values
    avgStrike = strike;

    double fwdAtStart = 1.0;
    // if it's forward starting, convert percentage strikes 
      // to absolute values based on spot at start date 
    if (fwdStarting)
    {
        fwdAtStart = asset->fwdValue(startDate);
        avgStrike *= fwdAtStart;
    }

    // init expected strike
    expectedStrike = avgStrike;

    // get FX and fwd discount for the actual dividend dates
    calcFXandDFdiv(fxDiv,gfDiv);

    // calc dollar amount
    for (j = 0;j < eqDivArr.size();j++)
    {
        divDollars[j] = div2dollar(eqDivArr[j], asset.get());
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
    DoubleArray propexp(assumedDivsInUse.size(), 1.0); // init to 1.0 if no adjustments needed
    if (avgWeightAdjusted)
        propAdjust(assumedDates, avgOut.get()->getDates(), avgOut.get()->getWeights(), propexp); 
    
    // initialize strikeAdjustmentToDate, to be updated in main loop below.
    if (scaleType != "NONE") {
        strikeAdjustmentToDate = 1.0;
        currAdjustment = 1.0;
    }
    else {
        strikeAdjustmentToDate = 0.0;
        currAdjustment = 0.0;
    }

    int k = 0;
    double adjust = strikeAdjustmentToDate;
    double adjustExp = strikeAdjustmentToDate;
    double m;
    double spotAtEx;
    double plusOrMinus = 1.0;

    if (scaleType == "STRIKE_MULTIPLIED") {
        plusOrMinus = -1.0;
    }

    // get on the fly num ex div dates to properly handle mu special
    CashFlowArraySP numExDivDatesInUseSP = getNumDivDates(eqDivArr);
    CashFlowArray& numExDivDatesInUse = *numExDivDatesInUseSP;

    //main loop 
    for (i = 0; i<assumedDivsInUse.size(); i++)
    { 
        bool isLast = (i==assumedDivsInUse.size()-1);
        const DateTime& endDate = (isLast ? lastOutDate : assumedDates[i+1]);

        // if there are only 0 divs in the qtr, then equally weight payments by setting periodSum to 1 / numDivs.
        double denom = Maths::isZero(periodSum[i]) ? 1.0 / numExDivDatesInUse[i].amount : periodSum[i];

        // assumed period contains actual dividends
        if (!Maths::isZero(numExDivDatesInUse[i].amount))
        {
            for (j = k; j<eqDivArr.size(); j++)
            {
                double numer = 1.0;
                if (!Maths::isZero(periodSum[i]))
                {
                    numer = divDollars[j];
                }
                else
                {   
                    numer = 1.0;
                }

                if(eqExDates[j] >= assumedDates[i] &&
                    ((!isLast && eqExDates[j] < endDate) ||
                    (isLast && eqExDates[j] <= endDate))) // last date is inclusive
                {
                    
                    if (scaleType != "NONE") 
                    {
                        DateTime beforeExDate((AssetUtil::getHoliday(asset.get())->addBusinessDays(eqExDates[j],-1)).getDate(), DateTime::END_OF_DAY_TIME);
                        if (eqExDates[j] <= valueDate)
                        {
                            if(spotBeforeExDates[j].date.getDate() == beforeExDate.getDate()) 
                            {
                                if (spotBeforeExDates[i].amount > 0.0)
                                {
                                    spotAtEx = spotBeforeExDates[i].amount;
                                }
                                else
                                {
                                    throw ModelException(method, "Past Historic Samples must be strictly superior to 0.");
                                }
                            }
                            else 
                            {
                                throw ModelException(method, "Past Historic Samples must be at Ex-date - 1 Business day");
                            }
                        }
                        else
                        {
                            spotAtEx = asset->fwdValue(beforeExDate);
                        }

                        // adjustment for the strike
                        adjust *= 1.0 + plusOrMinus * prop[j] * fxDiv[j] * (divDollars[j] - (assumedDivsInUse[i].amount) * numer/denom) / spotAtEx;
                    }

                    else 
                    {    
                        if(eqDivArr[j].getDivType() == Dividend::AMOUNT
                            && !(useEquivYields == 1 && eqExDates[j] > valueDate)) // not yield case
                            m = 1.0;
                        else
                            m = 0.0;
        
                        // adjustment for the strike
                        adjust += prop[j]*fxDiv[j]*gfDiv[j]*(divDollars[j]*m - assumedDivsInUse[i].amount*numer/denom);
                        adjustExp += prop[j]*fxDiv[j]*gfDiv[j]*(divDollars[j] - assumedDivsInUse[i].amount*numer/denom);
                    }
                    
                    // for past ex-dates, add to strikeAdjustmentToDate 
                    if (eqExDates[j] <= valueDate)
                    {
                        strikeAdjustmentToDate = adjust;
                    }
               
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
            double thisAsmAdj;

            if (scaleType != "NONE") 
            {
                DateTime beforeExDate((AssetUtil::getHoliday(asset.get())->addBusinessDays(assumedDates[i],-1)).getDate(), DateTime::END_OF_DAY_TIME);
                if (assumedDates[i].getDate() <= valueDate.getDate())
                {
                    if(spotBeforeAssumedDates[i].date.getDate() == beforeExDate.getDate()) 
                    {
                        if (spotBeforeAssumedDates[i].amount > 0.0)
                        {
                            spotAtEx = spotBeforeAssumedDates[i].amount;
                        }
                        else
                        {
                            throw ModelException(method, "Past Historic Samples must be strictly superior to 0.");
                        }
                    }
                    else 
                    {
                        throw ModelException(method, "Past Historic Samples must be at Assumed div date - 1 Business day");
                    }
                }
                else
                {
                    spotAtEx = asset->fwdValue(beforeExDate);                
                }

                thisAsmAdj = 1.0 - plusOrMinus * propexp[i] * fxAsmd[i] * (assumedDivsInUse[i].amount) / spotAtEx;
                adjust   *= thisAsmAdj;

                // for past ex-dates, add to strikeAdjustmentToDate 
                if (assumedDates[i].getDate() <= valueDate.getDate())
                {
                    strikeAdjustmentToDate = adjust;
                }
            }

            else 
            {
                thisAsmAdj = -propexp[i] * fxAsmd[i] * (assumedDivsInUse[i].amount) * gfAsmd[i];
                adjust   += thisAsmAdj;
                adjustExp += thisAsmAdj;

                // for past ex-dates, add to strikeAdjustmentToDate 
                if (endDate.getDate() - 1 <= valueDate.getDate())
                {
                    strikeAdjustmentToDate = adjust;
                }
            }

            if (endDate.getDate()  == valueDate.getDate())
            {
                currAdjustment = thisAsmAdj;
            }
        }
    }

    //dividend adjustment type for the adjusted strike
    if (scaleType == "STRIKE_DIVIDED") 
    {
        if (divadjType == "UP"){
            avgStrike /= Maths::min(adjust,1.0);
            strikeAdjustmentToDate = Maths::min(strikeAdjustmentToDate,1.0);
        }
        else if (divadjType == "DOWN"){
            avgStrike /= Maths::max(adjust,1.0);
            strikeAdjustmentToDate = Maths::max(strikeAdjustmentToDate,1.0);
        }
        else{
            avgStrike /= adjust;
        }
    }
    else if (scaleType == "STRIKE_MULTIPLIED") 
    {
        if (divadjType == "UP"){
            avgStrike *= Maths::max(adjust,1.0);
            strikeAdjustmentToDate = Maths::max(strikeAdjustmentToDate,1.0);
        }
        else if (divadjType == "DOWN"){
            avgStrike *= Maths::min(adjust,1.0);
            strikeAdjustmentToDate = Maths::min(strikeAdjustmentToDate,1.0);
        }
        else{
            avgStrike *= adjust;
        }
    }
    else
    {
        if (divadjType == "UP"){
            avgStrike -= Maths::min(adjust,0.0);
            expectedStrike -= Maths::min(adjustExp,0.0);
            strikeAdjustmentToDate = Maths::min(strikeAdjustmentToDate,0.0);
        }
        else if (divadjType == "DOWN"){
            avgStrike -= Maths::max(adjust,0.0);
            expectedStrike -= Maths::max(adjustExp,0.0);
            strikeAdjustmentToDate = Maths::max(strikeAdjustmentToDate,0.0);
        }
        else{
            avgStrike -= adjust;
            expectedStrike -= adjustExp;
        }
    }

    // keep abslute strike for printing
    absStrike = avgStrike;
    // scale strike for calling average model if fwd starting
    avgStrike /=fwdAtStart;
}


/* true = this is a ex-div date or last date-1 of assumed period with no div */
bool DividendAdjBase::recordHist() const
{
    // get on the fly num ex div dates to properly handle mu special
    CashFlowArraySP numExDivDatesInUseSP = getNumDivDates(getEqDivs());
    CashFlowArray& numExDivDatesInUse = *numExDivDatesInUseSP;
    
    bool record = false;

    for (int j = 0;j < eqDivArr.size();j++)
    {
        if(eqDivArr[j].getExDate().equals(valueDate,false))
        {
           record = true;
           break;
        }
    }
    for (int i = 0;i < assumedDivs.size() - 1;i++)
    {
        // note: the last day of this period is assumedDivs[i+1] - 1.
        if((assumedDivs[i+1].date.getDate()-1) == valueDate.getDate() && (int)numExDivDatesInUse[i].amount == 0)
        {
            record = true;
            break;
        }
    }

    // if we're at maturity with no divs in last period, record as well   
    if (matDate.equals(valueDate, false) && (int)numExDivDatesInUse.back().amount == 0)
    {
        record = true;
    }
    return record;
}


/** for assumedDivsInStruckCcy case, we convert assumed divs into base ccy before pricing */
CashFlowArraySP DividendAdjBase::assumedDivs2BaseCcy() const
{
    CashFlowArraySP newAssumedDivs = CashFlowArraySP(new CashFlowArray(assumedDivs));
    for (int i=0; i < (*newAssumedDivs).size(); i++)
    {
        double fx = 1.0;
        DateTime& exDate = (*newAssumedDivs)[i].date;
        if (exDate.getDate() < valueDate.getDate()) 
        {    // use history
            fx= getHist(histFXRates, exDate);
        }
        else
        {
            fx = getFXValue(exDate, asset.get());
        }
        if (Maths::isZero(fx))
            throw ModelException("DividendAdjBase::assumedDivs2BaseCcy", "Zero fx rate for " + (*newAssumedDivs)[i].date.toString());

        (*newAssumedDivs)[i].amount /= fx;
    }

    return newAssumedDivs; 
}

////////////////// registrations ////////////////
// for derived class
DividendAdjBase::DividendAdjBase(CClassConstSP const type): 
    Generic1Factor(type),
    strike(0.0),
    avgWeightAdjusted(false),
    usePayDate(false),
    scaleType("NONE"),
    isCall(0),
    divadjType("NONE"),
    fwdAdj(true),
    assumedDivsInStruckCcy(false)
{}

    class DividendAdjBaseHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(DividendAdjBase, clazz);
        SUPERCLASS(Generic1Factor);
        IMPLEMENTS(ISensitiveStrikes);
        IMPLEMENTS(Theta::Shift);
        IMPLEMENTS(LastSensDate);
        FIELD(strike,           "Strike");
        FIELD(matDate,          "maturity");
        FIELD(assumedDivs,      "assumed dividend list");
        FIELD(histDiscRates,    "Historic discount rates");
        FIELD_MAKE_OPTIONAL(histDiscRates);
        FIELD(histFXRates,      "Historic Fx rates");
        FIELD_MAKE_OPTIONAL(histFXRates);
        FIELD(numExDivDates,     "Number of Eq Divs in period");
        FIELD_MAKE_OPTIONAL(numExDivDates);
        FIELD(scaleType,  "Notional and strike adjustment type");
        FIELD_MAKE_OPTIONAL(scaleType);
        FIELD(spotBeforeExDates,  "past samples at actual ex-dates for notional and strike scaling");
        FIELD_MAKE_OPTIONAL(spotBeforeExDates);
        FIELD(spotBeforeAssumedDates,  "past samples at assumed div dates for notional and strike scaling");
        FIELD_MAKE_OPTIONAL(spotBeforeAssumedDates);

        FIELD(avgWeightAdjusted,     "Are dividend adjustments/payments proportional to avg weight left.  NOT AVAILABLE");
        FIELD_MAKE_OPTIONAL(avgWeightAdjusted);

        FIELD(assumedDivsInStruckCcy, "If True, then for ccy struck, assumed divs" 
                     "are specified in the struck ccy.");
        FIELD_MAKE_OPTIONAL(assumedDivsInStruckCcy);

        // for clone. these fields are registered in child class if need to appear at interface
        // to allow use of old parmeters
        FIELD(isCall, "Is it a call option");
        FIELD_MAKE_TRANSIENT(isCall);
        FIELD(divadjType,       "0: adjust both up and down.  1: up only.  2: down only");
        FIELD_MAKE_TRANSIENT(divadjType);
        FIELD(fwdAdj,           "Use forward growth factor adjustment for divs");
        FIELD_MAKE_TRANSIENT(fwdAdj);
        FIELD(avgOut,                  "average-out samples list");
        FIELD_MAKE_TRANSIENT(avgOut);

        // transient, used lots of places
        FIELD(eqDivArr,   "equity dividend array");
        FIELD_MAKE_TRANSIENT(eqDivArr);
        // for calling Average::getSensitiveStrikes()
        FIELD(avgInst,                  "avg instrument, used for getSensitiveStrike");
        FIELD_MAKE_TRANSIENT(avgInst);

        // to be removed once correctly treating disc rates etc ?
        FIELD(usePayDate, "false = use ex-dividend date,true = div pay-date");
        FIELD_MAKE_TRANSIENT(usePayDate);

    }
};

CClassConstSP const DividendAdjBase::TYPE = CClass::registerClassLoadMethod(
    "DividendAdjBase", typeid(DividendAdjBase), DividendAdjBaseHelper::load);
   

/** Implementation of ClosedForm::IntoProduct interface */
CClosedFormLN::IProduct* DividendAdjusted::createProduct(
    CClosedFormLN* model) const{

    return new DividendAdjustedClosedForm(model, this);
}

///////////////////// DividendAdjusted class members ////////////////
CClassConstSP const DividendAdjusted::TYPE = CClass::registerClassLoadMethod(
    "DividendAdjusted", typeid(DividendAdjusted), load);

void DividendAdjusted::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(DividendAdjusted, clazz);
    SUPERCLASS(DividendAdjBase);
    IMPLEMENTS(CClosedFormLN::IIntoProduct);
    EMPTY_SHELL_METHOD(defaultDividendAdjusted);
    FIELD(isCall,  "Is it a call option"); // old is an integer, new is bool in base
    FIELD(adjustmentType,   "Type of the strike adjustment"); // old integer, new is a string in base
    FIELD(useFwdValOfDivs,  "Use forward adjustment");
    FIELD(useEquivYields,   "All dividends as yields");
    FIELD(aveSchedOrig,      "average-out samples");
    FIELD(aveSchedOrigWeights, "average-out weights");

}

//// loader
extern bool DividendAdjFlowThroughLoad();

bool DividendAdjustedLoad()
{
    bool status = DividendAdjBase::TYPE &&
                  DividendAdjusted::TYPE &&
                  DividendAdjFlowThroughLoad();

    return status;
}


DRLIB_END_NAMESPACE

