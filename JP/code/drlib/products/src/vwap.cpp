//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : vwap.cpp
//
//   Description : volume weighted avg price instrument
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : September 26, 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/vwap.hpp"
#include "edginc/Maths.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/Tree1fLN.hpp"
#include "edginc/BadDayFollowing.hpp"
DRLIB_BEGIN_NAMESPACE

void VWAP::Validate()
{
    static const string method = "VWAP::Validate";
    // just check the things that aren't/cannot be checked in 
    // validatePop2Object
    if (!asset){
        throw ModelException(method, "Asset is null");
    }
    if (!discount){
        throw ModelException(method, "Discount YC is null");
    }

    DateTime matDate = lastDate;
    AssetUtil::assetCrossValidate(asset.get(),
                          false, //fwdStarting
                          startDate,
                          valueDate,
                          discount,
                          this);  
    
    if ( FXAsset::TYPE->isInstance(asset.get()) )
    {
        throw ModelException(method,
                             "Options on FX assets are not allowed yet");
    }

    // do some initialization.
    nDays = 0;
    BadDayFollowing bdc;
    
    int i;
    int daysDiff;
    DateTime tmpDate = valueDate.rollDate(-1);
    DateTimeArray sampleTimes = vwapSchedule->getDates();
    while ((tmpDate = bdc.adjust(tmpDate.rollDate(1), hols.get())).isLess(lastDate)) {
        daysDiff = tmpDate.getDate() - valueDate.getDate();
        nDays++;

        for (i = 0; i < vwapSchedule->length(); i++) {
            sampleDates.push_back(sampleTimes[i].rollDate(daysDiff));
        }
    }

    if (hols->businessDaysDiff(startDate, lastDate)+1 != minProfitPerShare->size()){
        throw ModelException(method, "minProfitPerShare array is not the same length as number of days "
            "between start and last dates" + Format::toString(hols->businessDaysDiff(startDate, lastDate)+1) + 
            "!=" +Format::toString(minProfitPerShare->size()));
    }

    currentDateOffset = hols->businessDaysDiff(startDate, valueDate);
}

void VWAP::validatePop2Object()
{
    static const string method("VWAP::validatePop2Object");
    int                 numDates;

    // can't get exercise schedule from Market - fail if it is NULL
    if( !vwapSchedule )
    {
        throw ModelException(method, "Exercise schedule is NULL");
    }

    // can't get instrument settlement from Market - fail if it is NULL
    if( !instSettle )
    {
        throw ModelException(method, "Instrument settlement is NULL");
    }

    // check that we have at least one entry in the exercise schedule
    numDates = vwapSchedule->length();
    if ( numDates < 1)
    {
        throw ModelException(method, "Exercise schedule is empty");
    }

}

void VWAP::GetMarket(const IModel*          model, 
                         const CMarketDataSP    market)
{
    static const string method("VWAP::GetMarket");
    market->GetReferenceDate(valueDate);

    CAsset::getAssetMarketData(model, market.get(), ccyTreatment, 
                               discount, asset);

    discount.getData(model, market);
    instSettle->getMarket(model, market.get());


    hols.getData(model, market);
    
}

DateTime VWAP::getValueDate() const
{
  return valueDate;
}

/** when to stop tweaking */
DateTime VWAP::endDate(const Sensitivity* sensControl) const {
    DateTime matDate = lastDate;
    DateTime instEnd  = instSettle->settles(matDate, asset.get());
    DateTime assetEnd = asset->settleDate(matDate);
    DateTime end = instEnd.isGreater(assetEnd) ? instEnd : assetEnd;
    return end;
}

void VWAP::addOutputRequests(Control* control,
                                 Results* results,
                                 const double& fairValue,
                                 const double& indVol) const
{
    DateTime matDate = lastDate;

    // take care of additional outputs
    if ( control && control->isPricing() )
    {
        // IND_VOL
        OutputRequest* request = NULL;
        if (control->requestsOutput(OutputRequest::IND_VOL, request)) 
        {
            results->storeRequestResult(request, indVol); 
        }
        
        // FWD_AT_MAT
        InstrumentUtil::recordFwdAtMat(control,
                                       results,
                                       matDate,
                                       valueDate,
                                       asset.get());

    }
}


/** indicates whether VEGA_MATRIX is sensible for this instrument */
bool VWAP::avoidVegaMatrix(const IModel* model)
{
    /* this should possibly be false for local vol models etc. */
    return false;
}

/** returns all strikes on the vol surface to which 
    this instrument is sensitive */
DoubleArraySP VWAP::getSensitiveStrikes(OutputNameConstSP outputName,
                                            const IModel*      model)
{
    DoubleArraySP sensStrikes = DoubleArraySP(new DoubleArray(0));

    if (avoidVegaMatrix(model)) {
        throw ModelException("VWAP::getSensitiveStrikes", 
                             "VEGA_MATRIX is not valid for this instrument");
    }

    double   strike       = asset->getSpot();

    // create a vol request object to be passed on
    LinearStrikeVolRequestSP volRequest(new LinearStrikeVolRequest(strike, 
                                                                   valueDate, 
                                                                   lastDate,
                                                                   false));

    SensitiveStrikeDescriptor sensStrikeDesc;
    sensStrikeDesc.forwardOnly = false;

    asset->getSensitiveStrikes(volRequest.get(), outputName, 
                               sensStrikeDesc, sensStrikes);

    return sensStrikes;
}

/** Rolls the value date and sets initial spot if rolling over start date */
bool VWAP::sensShift(Theta* shift)
{    
    DateTime newDate = shift->rollDate(valueDate);

    // no theta for vwap. 

    if (newDate.getDate() != valueDate.getDate()) {
        throw ModelException("VWAP::sensShift", 
                             "No Theta for VWAP");
    }

    // could maybe set a sample that occurs at the same time as valueDate if it is yet to be set as
    // this routine is called with a 0 shift before pricing.
    
    return true;
};


// for reflection
VWAP::VWAP(): 
CInstrument(TYPE){}

// Monte Carlo
class VWAP_MC : public IMCProduct, virtual public IMCProductLN {
private:
    const VWAP*     inst;
    double              strike;
    double              mult;
    bool                useRatio;

    double histSharesBoughtToday;
    double histVWAP;
    double histAmountForSharesToday;

    DoubleArray weights;
    double totalWeight;


public:
    

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen)const{
        // empty
    }

    // equivalent to InstIntoMCProduct 
    VWAP_MC(const VWAP*           inst,
              const IRefLevelConstSP&   refLevel,    // how to 'avg in'
              const SimSeriesConstSP&   simSeries,   // simulation dates
              const IPastValuesConstSP& mcPastValues, // historic values
              const DateTime&           matDate):    // maturity date
    IMCProduct(inst->asset.get(),
              inst->valueDate,
              inst->discount.get(),
              refLevel,
              simSeries,
              mcPastValues, // historic values
              inst->instSettle.get(),
              matDate),
    inst(inst) {
        strike = inst->vwapSchedule->lastValue();

        useRatio = false;
        mult = 1.0;
        
        weights = inst->vwapSchedule->getValues();
        int i;
        totalWeight = 0.;
        for (i=0; i < inst->vwapSchedule->length(); i++) {
             totalWeight += weights[i];
        }

        histSharesBoughtToday = 0.;
        histVWAP = 0.;
        histAmountForSharesToday = 0.;

    }

    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&                prices); 

    
    int getMinMaxIndex(double stockPrice) const
    {
        int i;
        int minMaxIndex = 0;
        for (i=0; i<inst->shareMinMax->size(); i++) {
            if (stockPrice < (*inst->shareMinMax)[i]) {
                break;
            } else {
                minMaxIndex = i;
            }
        }
        return minMaxIndex;
    }
    
    double getMinShares(double stockPrice) const
    {
        return (*inst->minPD)[getMinMaxIndex(stockPrice)];
    }

    // for the LogNormal path generator
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGenerator,
                                    int                     iAsset) const {

        CVolRequestLNArray   reqarr(1); // one interp level/path per asset here
        
        reqarr[0] = CVolRequestLNSP(new LinearStrikeVolRequest(
            inst->asset->getSpot(), // atm  
            inst->valueDate, 
            inst->lastDate,
            false));  // it's not forward starting
        
        return reqarr;
    }

    /** extra output requests */
    virtual void recordExtraOutput(CControl*     control,
                                   Results*      results,
                                   const IMCPrices& prices) const;
    
    double calcProfit(double sharePrice, 
                      double sharesBoughtToday, 
                      double amountForSharesToday, 
                      double vwap,
                      int j,
                      int *nDaysAbove, 
                      double *sharesTraded,
                      double *amountPaidByCompany) const;

};

/** extra output requests */
void VWAP_MC::recordExtraOutput(CControl*     control,
                                Results*      results,
                                const IMCPrices& prices) const
{

    double minShares;
    double sharesToBuy;
    int i;
    int nextSampleIdx;
    double sharesTraded = inst->sharesSoFar; // excluding today
    double amountPaidByCompany = inst->dollarsSoFar; // excluding today


    nextSampleIdx = -1;
    for (i = 0; i < inst->vwapSchedule->length(); i++) {
        if (inst->sampleDates[i] > inst->valueDate) {
            nextSampleIdx = i;
            break;
        }
    }
    
    if (nextSampleIdx != -1) {
        minShares = getMinShares(inst->asset->getSpot());
        if (inst->hasThreeDayRule == true) {
            if (inst->nDaysSoFar < 3 && inst->asset->getSpot() >= inst->threeDayLevel) {
                minShares = 0.; // can't buy
            }
        }
        minShares = Maths::min(minShares, inst->maxShares - sharesTraded );
        minShares = Maths::min(minShares, (inst->maxDollars - amountPaidByCompany )/inst->asset->getSpot());
        sharesToBuy = minShares*weights[nextSampleIdx]/totalWeight;
    } 
    else
    {
        double sharesTraded = inst->sharesSoFar;
        double amountPaidByCompany = inst->dollarsSoFar;
        int nDaysAbove = inst->nDaysSoFar;

        calcProfit(inst->asset->getSpot(),
                   histSharesBoughtToday, 
                   histAmountForSharesToday,
                   histVWAP,
                   0,
                   &nDaysAbove,
                   &sharesTraded,
                   &amountPaidByCompany);

        sharesToBuy = sharesTraded - inst->sharesSoFar - histSharesBoughtToday;
    }
    
    OutputNameConstSP notionalOutput(new OutputName("SHARES_TO_BUY"));
    results->storeScalarGreek(sharesToBuy, Results::DEBUG_PACKET, notionalOutput);

}


void VWAP_MC::payoff(const IPathGenerator*  pathGen,
                     IMCPrices&                prices)
{
    double price = 0.;
    
    if (pathGen->doingPast() == true) {
        int    firstStep = pathGen->begin(0); // only one asset
        int    lastStep = pathGen->end(0); // only one asset
        int i;

        double sharesTraded = inst->sharesSoFar; // excluding today
        double amountPaidByCompany = inst->dollarsSoFar; // excluding today

        double minShares;
        double sharesToBuyNow;

        for (i=firstStep; i<lastStep; i++) {
            minShares = getMinShares(pathGen->Path(0,0)[i]);
            if (inst->hasThreeDayRule == true) {
                if (inst->nDaysSoFar < 3 && pathGen->Path(0,0)[i] >= inst->threeDayLevel) {
                    minShares = 0.; // can't buy
                }
            }
            minShares = Maths::min(minShares, inst->maxShares - sharesTraded );
            minShares = Maths::min(minShares, (inst->maxDollars - amountPaidByCompany )/pathGen->Path(0,0)[i]);
            sharesToBuyNow = minShares*weights[i]/totalWeight;
            histSharesBoughtToday += sharesToBuyNow;
            histVWAP +=  pathGen->Path(0,0)[i]*weights[i]/totalWeight;
            histAmountForSharesToday += pathGen->Path(0,0)[i]*sharesToBuyNow;
        }

    } else {
        
        int i,j;
        int iStart;
        double profit;
        int npd = inst->vwapSchedule->length();
        
        double sharesTraded = inst->sharesSoFar;

        double minShares;
        double sharesBoughtToday;
        double sharesToBuyNow;
        double amountForSharesToday;
        double vwap;
        double amountPaidByCompany = inst->dollarsSoFar;


        int nDaysAbove;

        
        for (j=0; j<inst->nDays; j++) {
            sharesBoughtToday = 0.;
            amountForSharesToday = 0.;
            vwap = 0;
            
            // I doubt Theta would work for this
            if (j == 0) {
                iStart = pathGen->begin(0);
                sharesBoughtToday = histSharesBoughtToday;
                vwap = histVWAP;
                amountForSharesToday = histAmountForSharesToday;
                nDaysAbove = inst->nDaysSoFar;
            } else {
                iStart = 0;
            }

            if (sharesTraded < inst->maxShares && amountPaidByCompany < inst->maxDollars) {
               
                for (i = iStart; i<npd; i++) {
                    minShares = getMinShares(pathGen->Path(0,0)[j*npd+i]);
                    if (inst->hasThreeDayRule == true) {
                        if (nDaysAbove < 3 && pathGen->Path(0,0)[j*npd+i] >= inst->threeDayLevel) {
                            minShares = 0.; // can't buy
                        }
                    }
                    minShares = Maths::min(minShares, inst->maxShares - sharesTraded);
                    minShares = Maths::min(minShares, (inst->maxDollars-amountPaidByCompany)/pathGen->Path(0,0)[j*npd+i]);
                    sharesToBuyNow = minShares*weights[i]/totalWeight;
                    sharesBoughtToday += sharesToBuyNow;
                    vwap +=  pathGen->Path(0,0)[j*npd+i]*weights[i]/totalWeight;
                    amountForSharesToday += pathGen->Path(0,0)[j*npd+i]*sharesToBuyNow;
                }
                
                profit = calcProfit(pathGen->Path(0,0)[(j+1)*npd-1],
                    sharesBoughtToday, 
                    amountForSharesToday,
                    vwap,
                    j,
                    &nDaysAbove,
                    &sharesTraded,
                    &amountPaidByCompany);

                price += profit*inst->discount->pv(inst->valueDate, inst->valueDate.rollDate(j));
            } 
        } 
    }
    prices.add(price * mult);
}


double VWAP_MC::calcProfit(double sharePrice, 
                           double sharesBoughtToday, 
                           double amountForSharesToday,
                           double vwap,
                           int j,
                           int *nDaysAbove, 
                           double *sharesTraded,
                           double *amountPaidByCompany) const
{
    int minMaxIndex;
    bool noSharesToday;
    double maxSharesToday;
    double minSharesToday;
    double profit;

    if (sharePrice >= inst->threeDayLevel) {
        *nDaysAbove = *nDaysAbove + 1;
        if (*nDaysAbove <= 3) {
            noSharesToday = true;
        } 
        else {
            noSharesToday = false;
        }
    } 
    else {
        *nDaysAbove = 0;
        noSharesToday = false;
    }
    
    if (noSharesToday == true) {
        // must sell off any shares bought today
        profit = sharesBoughtToday*sharePrice - amountForSharesToday;

    } else {

        minMaxIndex = getMinMaxIndex(sharePrice);
    
    
        if (*amountPaidByCompany + amountForSharesToday > inst->maxDollars) {
            maxSharesToday = (inst->maxDollars - *amountPaidByCompany)/vwap;
        } else {
            maxSharesToday = (*inst->maxPD)[minMaxIndex];
        }

        if (sharesBoughtToday + *sharesTraded > inst->maxShares) {
            maxSharesToday = Maths::min(maxSharesToday, inst->maxShares - *sharesTraded);
        }

        minSharesToday = Maths::min((*inst->minPD)[minMaxIndex], maxSharesToday);

        if (minSharesToday == maxSharesToday) {
            // must deliver maxSharesToday
            double nSharesForMarket = sharesBoughtToday - maxSharesToday;

            *sharesTraded += maxSharesToday*1.000000001; // add a small amount to make sure the 
            *amountPaidByCompany += maxSharesToday*vwap*1.000000001;
        
            profit = maxSharesToday*vwap + nSharesForMarket*sharePrice 
                - amountForSharesToday;
        } 
        else // can buy or sell
        {
            double tmpSharesToTrade;

            double profitMaxShares;
            double profitMinShares;
            double profitDoNothing;
            bool canDoNothing;
        
            profitMaxShares = (sharesBoughtToday - maxSharesToday)*sharePrice +
                maxSharesToday*vwap - amountForSharesToday;

            profitMinShares = (sharesBoughtToday - minSharesToday)*sharePrice +
                minSharesToday*vwap - amountForSharesToday;

            if (sharesBoughtToday >= minSharesToday && sharesBoughtToday <= maxSharesToday) {
                profitDoNothing = sharesBoughtToday*vwap - amountForSharesToday;
                canDoNothing = true;
            } else {
                profitDoNothing = -99999999999.;
                canDoNothing = false;
            }


            if (profitMinShares > profitMaxShares) { // it's linear so we don't need to check profitDoNothing
                // always go to the min number if it makes the most profit
                profit = profitMinShares;
                tmpSharesToTrade = minSharesToday;
            } 
            else if (profitMaxShares > (*inst->minProfitPerShare)[j+inst->currentDateOffset]*maxSharesToday) {
                // buy to the max if you make at least minProfitPerShare
                profit = profitMaxShares;
                tmpSharesToTrade = maxSharesToday;   
            } 
            else if (canDoNothing == false) {
                if (sharesBoughtToday >= maxSharesToday) {
                    if (profitMinShares > 0.) {
                        // Could it be worthwhile to take a loss? I can't imagine anyone being able to 
                        // cross that psychological barrier to do so.
                   
                        // sell down to min. 
                        profit = profitMinShares;
                        tmpSharesToTrade = minSharesToday;
                    }
                    else {
                        // sell to max
                        profit = profitMaxShares;
                        tmpSharesToTrade = maxSharesToday;   
                    }
                
                }
                else
                {
                    // buy to min
                    profit = profitMinShares;
                    tmpSharesToTrade = minSharesToday;
                }
            }
            else {
                // do nothing
                profit = profitDoNothing;
                tmpSharesToTrade = sharesBoughtToday;
            }
        
            *sharesTraded += tmpSharesToTrade;
            *amountPaidByCompany += tmpSharesToTrade*vwap;
        }   
    }

    return profit;
}

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* VWAP::createProduct(const MonteCarlo* model) const {
    static const string routine("VWAP::createProduct");
    int i;


    
    // v simple simSeries
    SimSeriesSP simSeries(new SimSeries(1));
    simSeries->addDates(sampleDates);

    // v simple RefLevel
    IRefLevelSP refLevel(IRefLevel::Util::makeFwdStart(startDate));
    
    // need to add start date and maturity to past values
    DateTimeArray  allDates(1, startDate);
    DoubleArray    allValues(1, 999999.);

    DateTimeArray sampleTimes = vwapSchedule->getDates();

    for (i = 0; i < vwapSchedule->length(); i++) {
        if (sampleTimes[i] <= valueDate) {
            allDates.push_back(sampleTimes[i]);
            allValues.push_back((*histSharePrice)[i]);
        }
    }

    DoubleMatrix allValuesAsMatrix(allValues);

    IPastValuesSP pastValues(
        IPastValues::Util::makeSimple(allDates, allValuesAsMatrix));
    
    IMCProduct*  prod = new VWAP_MC(this, refLevel, simSeries, 
                                     pastValues, vwapSchedule->lastDate());

    return prod;
}


/** Returns the name of the instrument's discount currency */
string VWAP::discountYieldCurveName() const {
    return discount.getName();
}



class VWAPHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VWAP, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(IMCIntoProduct);
        IMPLEMENTS(ISensitiveStrikes);
        IMPLEMENTS(Theta::Shift);
        IMPLEMENTS(LastSensDate);
        EMPTY_SHELL_METHOD(defaultVanilla);
        FIELD(valueDate,        "valuation Date");
        FIELD(startDate,        "Option start date");
        FIELD(vwapSchedule,            "vwap sample schedule");
        FIELD(asset,            "Underlying of option");
        FIELD(discount,         "Discount curve");
        FIELD(ccyTreatment,     "Currency Treatment");
        FIELD(instSettle, "Instrument settlement at maturity");
        FIELD(lastDate,         "last date purchases can be made");
        FIELD(nDays,"Private");
        FIELD_MAKE_TRANSIENT(nDays);
        FIELD(hols,"Holidays");
        FIELD(sampleDates,"Private");
        FIELD_MAKE_TRANSIENT(sampleDates);
        FIELD(maxShares,"maximum number of shares to trade");
        FIELD(sharesSoFar,"shares traded already");
        FIELD(shareMinMax,"share price for min/max per day");
        FIELD(minPD,"minimum number of shares to trade per day");
        FIELD(maxPD,"maximum number of shares to trade per day");
        FIELD(minProfitPerShare, "minimum profit per extra share traded")
        FIELD(maxDollars,"maximum number of dollars company will spend");
        FIELD(dollarsSoFar,"dollars company has spent already");
        FIELD(histSharePrice, "historical prices on value date");
        FIELD(hasThreeDayRule, "has a three-day rule");
        FIELD(threeDayLevel,"share price must close above this level for three days for shares to be delivered if closing price is above this level");
        FIELD_MAKE_OPTIONAL(threeDayLevel);
        FIELD(nDaysSoFar,"number of days consecutively above three day limit so far");
        FIELD_MAKE_OPTIONAL(nDaysSoFar);
        FIELD(currentDateOffset,"Private");
        FIELD_MAKE_TRANSIENT(currentDateOffset);
    }

    static IObject* defaultVanilla(){
        return new VWAP();
    }
};

CClassConstSP const VWAP::TYPE = CClass::registerClassLoadMethod(
    "VWAP", typeid(VWAP), VWAPHelper::load);
bool  VWAPLoad() {
    return (VWAP::TYPE != 0);
   }


DRLIB_END_NAMESPACE

