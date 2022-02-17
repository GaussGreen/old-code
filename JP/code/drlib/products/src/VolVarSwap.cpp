//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolVarSwap.cpp
//
//   Description : VolVarSwap contract
//
//   Author      : Ning shen
//
//   Date        : 12 Feb 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolVarSwap.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/Black.hpp"
#include "edginc/CashSettlePeriod.hpp"
#include "edginc/MarketDataFetcherLNSpline.hpp"
#include "edginc/VolSV.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/imslerror.hpp"
#include "edginc/NonPricingModel.hpp"
#include "edginc/BadDayConventionFactory.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/StubFactory.hpp"
#include "edginc/SwapTool.hpp"
#include "edginc/ProtEquity.hpp"
#include "edginc/ProtAsset.hpp"
#include "edginc/XCB.hpp"
#include "edginc/Fund.hpp"
#include "edginc/FourierProcessSV.hpp"
#include "edginc/FourierProcessSVJ.hpp"
#include "edginc/FourierProcessSVCJ.hpp"
#include "edginc/Scenario.hpp"
#include "edginc/VolRequestTime.hpp"
#include "edginc/ClientRunnable.hpp"
#include "edginc/ObservationBuilder.hpp"
#include "edginc/MDFUtil.hpp"
#include "edginc/PastSamplesEvent.hpp"
#include "edginc/VegaMatrixLite.hpp"
#include "edginc/SimpleEquity.hpp"

DRLIB_BEGIN_NAMESPACE

#define xxxRG

/** instrument validation */
void VolVarShell::Validate() {
    static const string method = "VolVarShell::Validate";

    // check that settlement is cash
    if (instSettle->isPhysical()) {
        throw ModelException(method,
                             "Only cash settlement is allowed");
    }

    if (ccyTreatment == CAsset::CCY_TREATMENT_STRUCK) {
        throw ModelException(method, "swap can't be ccy struck");
    }

    if (FXAsset::TYPE->isInstance(asset.get())){
        throw ModelException(method, "FX underlying not allowed");
    }

    if (samples.size() < 2 && !subtractMeanVol){
        throw ModelException(method,
                             "Need at least 2 samples in the vol sample");
    }
    if (samples.size() < 3 && subtractMeanVol){
        throw ModelException(method,
                             "Need at least 3 samples in the vol sample");
    }

    int i;
    for (i = 1; i < samples.size(); i++) {
        /* fail if sample dates are not increasing */
        if (samples[i].date <= samples[i-1].date)
        {
            throw ModelException(method,
                "vol sample point " +
                Format::toString(i) +
                " has date (" +
                samples[i].date.toString() +
                ") which is before (or on) previous "
                "sample date (" +
                samples[i-1].date.toString() + ")");
        }
    }

    for (i = 0;
    i < samples.size() && valueDate.isGreater(samples[i].date);
    i++) {
        if (!Maths::isPositive(samples[i].amount)) {
            throw ModelException(method,
                "historic sample on (" +
                samples[i].date.toString() +
                ") for " + asset->getTrueName() + " is missing");
        }
    }

    if (Maths::isNegative(strikeVol)) {
        throw ModelException(method,
            "strike is not >= 0 (" +
            Format::toString(strikeVol) + ")");
    }

    if (Maths::isZero(strikeVol) && !dontScaleByStrike) {
        throw ModelException(method, "can't scale by a zero strike");
    }

    AssetUtil::assetCrossValidate(asset.get(),
                                  false,
                                  valueDate,
                                  valueDate,
                                  discount,
                                  this);

    if (!isVanilla) {
        // Populate total (expected) number of returns
        numTotalReturns = samples.size()-1;
    }

    if(observationsPerYear == 0) {
        throw ModelException(method, "PPY (observationsPerYear) must not be zero");
    }
}

/** indicates whether VEGA_MATRIX is sensible for this instrument */
bool VolVarShell::avoidVegaMatrix(const IModel* model) {
    if (ImpliedIntegration::TYPE->isInstance(model)) {
        return false;
    } else if (ClosedFormIntegrateLN::TYPE->isInstance(model)) {
        return false;
    } else {
        return true;
    }
}

/** returns all strikes on the vol surface to which this instrument is sensitive */
DoubleArraySP VolVarShell::getSensitiveStrikes(OutputNameConstSP outputName,
                                          const IModel*      model) {
    static const string method("VolVarShell::getSensitiveStrikes");
    try {
        if (avoidVegaMatrix(model)) {
            throw ModelException("VEGA_MATRIX is not valid for this instrument.");
        }
        return getSensitiveStrikesHelper(outputName, valueDate, samples.back().date, asset.get(), model);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

DoubleArraySP VolVarShell::getSensitiveStrikesHelper(OutputNameConstSP outputName,
                                               const DateTime& valueDate,
                                               const DateTime& maturityDate,
                                               const CAsset*   asset,
                                               const IModel*         model) {
    static const string method("VolVarShell::getSensitiveStrikesHelper");
    try {
        DoubleArraySP sensStrikes;
        if (ImpliedIntegration::TYPE->isInstance(model)) {
            const int tailWidth = 20;   // arbitrary
            const int midWidth  = 15;   // arbitrary
            ImpliedIntegration& ii = dynamic_cast<ImpliedIntegration&>(*const_cast<IModel*>(model));
            sensStrikes = ii.sensitiveStrikes(valueDate,maturityDate,asset,midWidth,tailWidth);
        } else if (ClosedFormIntegrateLN::TYPE->isInstance(model)) {
            ClosedFormIntegrateLN& cfi = dynamic_cast<ClosedFormIntegrateLN&>(*const_cast<IModel*>(model));
            sensStrikes = cfi.sensitiveStrikes(outputName, asset,valueDate,maturityDate);
        }
        return sensStrikes;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }

}

/** Returns rolls value date and sets initial spot for Theta,
    return true if sub objects need to be tweaked */
bool VolVarShell::sensShift(Theta* shift) {
    DateTime newDate = shift->rollDate(valueDate);

    // fill vol sample point if needed
    double spot    = 0.0;
    bool   useSpot = !shift->useAssetFwds();

    if (useSpot) {
        spot = asset->getSpot();
    }

    bool pastRollDate = false;
    int  i = 0;
    while (!pastRollDate && i < samples.size() ) {
        if ((samples[i].date.isGreater(valueDate) &&
             !samples[i].date.isGreater(newDate)) ||
            (samples[i].date.equals(valueDate)    &&
             Maths::isZero(samples[i].amount))) {
            if (useSpot) {
                samples[i].amount = spot;
            }
            else {
                samples[i].amount = asset->fwdValue(samples[i].date);
            }
        }
        else if (samples[i].date.isGreater(newDate)) {
            pastRollDate = true;
        }
        ++i;
    }

    return Generic1Factor::sensShift(shift);
}

DateTime VolVarShell::endDate(const Sensitivity* sensControl) const {
    const DateTime& matDate = samples[samples.size()-1].date;
    const DateTime& instEnd  = instSettle->settles(matDate, asset.get());
    const DateTime& assetEnd = asset->settleDate(matDate);
    const DateTime& end = instEnd.isGreater(assetEnd) ? instEnd : assetEnd;
    return end;
}

// determine static hedge
DoubleArray VolVarShell::staticHedge(
    CMarketDataSP      market,
    IModel*            model,
    const DoubleArray& strikes,
    const BoolArray&   isCall,
    double             loShares,
    double             hiShares) {
    static const string method("VolVarShell::staticHedge");
    try {
        if (Maths::isNegative(loShares)) {
            throw ModelException(method, "low shares (" +
                Format::toString(loShares) +
                ") is -ve");
        }
        if (Maths::isNegative(hiShares)) {
            throw ModelException(method, "high shares (" +
                Format::toString(hiShares) +
                ") is -ve");
        }

        if (strikes.size() < 2) {
            throw ModelException(method, "need at least two strikes");
        }

        int i;
        for (i = 1; i < strikes.size(); i++) {
            if (strikes[i] <= strikes[i-1]) {
                throw ModelException(method,
                    "strikes must be in ascending order"
                    "- strike " + Format::toString(i) +
                    " ("+Format::toString(strikes[i]) +
                    ") <= strike "+Format::toString(i-1)+
                    " ("+Format::toString(strikes[i-1]));
            }
        }

        if (!Maths::isPositive(strikes[0])) {
            throw ModelException(method, "strikes must be > 0");
        }

        if (isCall.size() != strikes.size()) {
            throw ModelException(method,
                "strikes and call/put must be same length");
        }

        // grab that market data
        model->getInstrumentAndModelMarket(market.get(), this);


        DoubleArray hedges(strikes.size()+1);
        DoubleArray k(strikes.size()+1);
        DoubleArray payoff(strikes.size()+1);

        DateTime maturity = samples[samples.size()-1].date;
        double   fwd      = asset->fwdValue(maturity);

        HolidayConstSP hols(AssetUtil::getHoliday(asset.get()));

        int    daysInSwap = hols->businessDaysDiff(valueDate, maturity);
        double time       = (double)observationsPerYear/(double)daysInSwap;

        k[0] = strikes[0]/2;
        for (i = 0; i < strikes.size(); i++) {
            k[i+1] = strikes[i];
        }

        for (i = 0; i < k.size(); i++) {
            payoff[i] = -2.0*(log(k[i]/fwd)-k[i]/fwd+1)*time;
        }

        // given bounds
        hedges[0] = loShares;
        hedges[hedges.size()-1] = hiShares;

        // given i strikes & call/puts, return i+1 hedges
        DoubleArray callDelta(k.size());
        DoubleArray putDelta(k.size());

        for (i = 0; i < k.size()-1; i++) {
            callDelta[i] = isCall[i-1] ? (payoff[i+1]-payoff[i])/(k[i+1]-k[i]) : 0.0;
        }
        for (i = 1; i < k.size(); i++) {
            putDelta[i]  = !isCall[i-1] ? -(payoff[i]-payoff[i-1])/(k[i]-k[i-1]) : 0.0;
        }

        for (i = 1; i < hedges.size()-1; i++) {
            if (isCall[i-1]) {
                hedges[i] = callDelta[i]-callDelta[i-1];
            }
            else {
                hedges[i] = putDelta[i]-putDelta[i+1];
            }
        }

        return hedges;
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

/** calculate historical vol */
double VolVarShell::historicalVol(int& numReturns) const {
    return historicalVol(numReturns,
                         valueDate,
                         asset->getSpot());
}

/************************************************************************************/
/** Div Adjuster to remove continuous dividends for dividend adjustment computation */
/************************************************************************************/
VolVarShell::DivAdjuster::DivAdjuster() {}


void VolVarShell::DivAdjuster::adjustDividend(Dividend& div){
    if (div.getDivType() == Dividend::CONTINUOUS)
    {
        div.convertToDollar(0.0);
    }
}


VolVarShell::DivAdjuster::~DivAdjuster(){}


/** calculate historical vol, given a value date */
double VolVarShell::historicalVol(int&            numReturns,
                                  const DateTime& valueDate,
                                  double          spot) const {
    static const string method("VolVarShell::historicalVol");
    try {
        double vol;

        if (samples.size() <= 1) {
            throw ModelException(method,
                "number of samples (" +
                Format::toString(samples.size()) +
                ") <= 1");
        }

        // see if there's no history
        if (samples[0].date >= valueDate) {
            numReturns = 0;
            return 0.0;
        }

        // extract all past divs in (first sample date, value date].
        VolVarShell::DivAdjuster divAdjuster;

        //   // Then initialise our DividendCollector
        DividendCollector collector
          (asset.get(),
          &divAdjuster, // optional
          valueDate,
          samples.front().date,
          valueDate);

        DividendListSP divs = collector.getDividends();

        const DividendArray& divArray = divs->getArray();


        numReturns = 0;

        int    numPastSamples = 1;
        double logRelative;
        double sumSqr = 0.0;
        double sum    = 0.0;

        // Set up an index to match div dates with Sample Dates.
        int iDiv = 0;
        int iDivStart = 0;

        // Include all samples if pricing at time of expiry.
        // Note if expiry is 23:59:58 and valueDate is EOD then expiry<valueDate so last return is included
        bool onExpiry = samples.back().date.equals(valueDate);

        for (int i = 1; i < samples.size() &&
            (valueDate.isGreater(samples[i].date) || onExpiry);
            i++)
        {
            double divAmount = 0.0;
            if (!noDivAdj)
            {
                // for each sample, need to compute sum of divs in
                // the range (previous sample date, this sample date]
                bool overThisSample = false;
                divAmount = 0.0; // reset to 0.0 for each sample date
                iDiv = iDivStart;

                while (iDiv < divArray.size() && !overThisSample)
                {
                    const int &thisSample = samples[i].date.getDate();
                    const int &previousSample = samples[i-1].date.getDate(); // note that i starts at 1
                    const int &divDate =  divArray[iDiv].getExDate().getDate();

                    if (previousSample < divDate && divDate <= thisSample)
                    {
                        // accumulate div
                        divAmount += divArray[iDiv].getDivAmount();
                    }

                    if (divDate > thisSample)
                    {
                        overThisSample = true;
                        iDivStart = iDiv;
                    }
                    iDiv++;
                }
            }

            numPastSamples++;

			logRelative = divAdjOnExDate? log((samples[i].amount + divAmount)/samples[i-1].amount):
                                          log(samples[i].amount/(samples[i-1].amount - divAmount));
            sumSqr      += logRelative * logRelative;
            sum         += logRelative;
            numReturns++;
        }

        // ideally we want to calculate historic vol from the historic
        // samples and the future vol from the vol surface between value
        // date and the last sample date.  However, this leaves out the
        // vol between the last historic sample and value date.  To cover
        // this we set the sample after value date to the current spot
        // and calculate the historic vol up the sample after value date.
        // This essentially captures the vol between last historic sample
        // sample and value date. Set sample on or after value date to spot

        // if adjusting hist variance for divs, must adjust current return
		double divToday = 0.0;
		if (!noDivAdj)
		{
			bool overToday = false;
			int iDiv = 0;
			while (iDiv < divArray.size() && !overToday)
			{
               int divDate =  divArray[iDiv].getExDate().getDate();
			   if (divDate == valueDate.getDate())
			   {
				   divToday += divArray[iDiv].getDivAmount();
			   }

			   if (divDate > valueDate.getDate())
			   {
				   overToday = true;
			   }

               iDiv++;
			}
		}

        if (numPastSamples < samples.size()) {
			logRelative = divAdjOnExDate? log( (spot+divToday) /samples[numPastSamples-1].amount):
		                                  log(spot/(samples[numPastSamples-1].amount - divToday));

            sumSqr      += logRelative * logRelative;
            sum         += logRelative;
            numReturns++;
        }

        if (isVanilla) {
            numReturns = numPastReturns;
        } else {
            VolVarShell* dummyPtr = const_cast<VolVarShell*>(this);
            dummyPtr->numPastReturns = numReturns;
        }

        // need at least 2 returns to get mean
        if (!subtractMeanVol || numReturns < 2) {
            vol = sqrt(sumSqr/numReturns);
        }
        else {
            double numTotalReturnsD = (double)(numTotalReturns);
            double numReturnsD = (double)(numReturns);
            vol = sqrt(numTotalReturnsD / (numTotalReturnsD - 1.0) / numReturnsD * (sumSqr - sum*sum / numReturnsD));
        }
        return vol * sqrt((double)observationsPerYear);
     }
     catch (exception& e) {
         throw ModelException(e, method);
     }
}




class VolVarShellHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(VolVarShell, clazz);
        SUPERCLASS(Generic1Factor);
        IMPLEMENTS(LastSensDate);
        IMPLEMENTS(ISensitiveStrikes);
        FIELD(samples, "vol samples");
        FIELD(strikeVol, "contract strike in vol unit");
        FIELD(observationsPerYear, "number of sample points per year");
        FIELD(subtractMeanVol, "if mean vol is subtracted");
        FIELD(payoffType, "FORWARD=vol/var swap (default), CALL = call option, PUT = put option, FORWARD_CAP = cap portion only of a capped var swap");
        FIELD_MAKE_OPTIONAL(payoffType);
        FIELD(dontScaleByStrike, "TRUE = don't divide by 2K");
        FIELD_MAKE_OPTIONAL(dontScaleByStrike);
        FIELD(cap, "Cap relative to strike");
        FIELD_MAKE_OPTIONAL(cap);
        FIELD(noDivAdj, "Addition additional effect: E[sum(log(1-D)^2)] due to discrete dividends to variance.");
        FIELD_MAKE_OPTIONAL(noDivAdj);
		FIELD(divAdjOnExDate, "true - add to ex-date sample; false - subtract from previous sample");
        FIELD_MAKE_OPTIONAL(divAdjOnExDate);
        FIELD(isVanilla, "whether a vanilla variance swap");
        FIELD_MAKE_TRANSIENT(isVanilla);
        FIELD(numPastReturns, "number of past returns");
        FIELD_MAKE_TRANSIENT(numPastReturns);
        FIELD(numTotalReturns, "total (expected) number of returns");
        FIELD_MAKE_TRANSIENT(numTotalReturns);
    }
};

CClassConstSP const VolVarShell::TYPE = CClass::registerClassLoadMethod(
    "VolVarShell", typeid(VolVarShell), VolVarShellHelper::load);

typedef smartPtr<VolVarShell> VolVarShellSP;

class StaticHedge : public CObject {
    static CClassConstSP const TYPE;

    VolVarShellSP imnt;
    CMarketDataSP market;
    IModelSP      model;
    DoubleArray   strikes;
    BoolArray     isCall;
    double        loShares;
    double        hiShares;

    StaticHedge() : CObject(TYPE) {}

    // addin for hedges
    static IObjectSP hedgeAddin(StaticHedge* params) {
        DoubleArray hedges = params->imnt->staticHedge(params->market,
                                                       params->model.get(),
                                                       params->strikes,
                                                       params->isCall,
                                                       params->loShares,
                                                       params->hiShares);

        DoubleArraySP output(new DoubleArray(hedges.size()));
        *output = hedges;
        return IObjectSP(output);
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic();
        clazz->setDescription("compute static hedge for vol and variance swaps");
        REGISTER(StaticHedge, clazz);
        SUPERCLASS(CObject);
        FIELD(imnt, "imnt");
        FIELD(market, "market");
        FIELD(model, "model");
        FIELD(strikes, "strikes");
        FIELD(isCall, "call/put per strike");
        FIELD(loShares, "low shares");
        FIELD(hiShares, "high shares");
        EMPTY_SHELL_METHOD(defaultStaticHedge);

        Addin::registerClassObjectMethod("VOLVAR_SWAP_STATIC_HEDGE",
                                         Addin::RISK,
                                         "static hedges for vol/variance swap",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)hedgeAddin);
    }

    static IObject* defaultStaticHedge(){
        return new StaticHedge();
    }
};

CClassConstSP const StaticHedge::TYPE = CClass::registerClassLoadMethod(
    "StaticHedge", typeid(StaticHedge), load);


// Type registration
CClassConstSP const VarCapCVModel::TYPE =
    CClass::registerClassLoadMethod("VarCapCVModel", typeid(VarCapCVModel),
    VarCapCVModel::load);

const string VarCapCVModel::DEFAULT = "default";
const string VarCapCVModel::CONTROL_VARIATE = "control_variate";
const string VarCapCVModel::MULTIPLICATIVE_SCALING = "multiplicative_scaling";
const string VarCapCVModel::ADJUST_MEAN_VOL = "adjust_mean_vol";

// Product
class VarCapCVModel::Product {
public:
    Product(VarianceSwapSP  varCap,
            CAssetSP        marketAsset):
    varCap(varCap),
    marketAsset(marketAsset) {}

    virtual ~Product() {}

    VarianceSwapSP  varCap;          // Variance cap - supporting multiple integrands
    CAssetSP        marketAsset;     // Asset with VolSurface
};

// variance swap

void VarianceSwap::avoidVegaMatrixLite(const IModel* model) {
    static const string method = "VarianceSwap::avoidVegaMatrixLite";
    
    if (avoidVegaMatrix(model)) {
        // Check basic VEGA_MATRIX
        throw ModelException(method, "Instrument does not support VEGA_MATRIX and hence not VEGA_MATRIX_LITE");
    } else if (!ClosedFormIntegrateLN::TYPE->isInstance(model)) {
        // Make sure we have ClosedFormIntegrateLN
        throw ModelException (method, "VEGA_MATRIX_LITE is only supported for model ClosedFormIntegrateLN");
    } else if (!SimpleEquity::TYPE->isInstance(asset.get())) {
        // Allow LITE only for SimpleEquity
        throw ModelException (method, "Only SimpleEquity underlyings supported");
    } else if(!CString::equalsIgnoreCase(payoffType, "FORWARD")) {
        // Disallow Caps or Capped Swaps
        throw ModelException (method, "Only Swap style instruments supported. Payoff type is " + payoffType);
    }
}



void VarianceSwap::validatePop2Object(){
    //always set it to false
    //this can only be overriden by createProduct
    useCV = false;

    //create copy of asset wrapper
    smartPtr<CAssetWrapper> tmpPtr(copy(&asset));
    marketAsset = *tmpPtr;
}

/** retrieve market data */
void VarianceSwap::GetMarket(const IModel*       model,
                             const CMarketDataSP market) {
    static const string routine("Variance::GetMarket");

    const VarCapCVModel* cvModel = dynamic_cast<const VarCapCVModel*>(model);
    //VarCapCVModel* cvModel = dynamic_cast<VarCapCVModel*>(const_cast<IModel*>(model));
    if(cvModel) {
        GetCVMarket(cvModel, market);
        return;
    }

    market->GetReferenceDate(valueDate);
    CAsset::getAssetMarketData(model, market.get(), ccyTreatment,
                               discount, asset);
    discount.getData(model, market);
    instSettle->getMarket(model, market.get());

    VarianceSwapUtil::validateBasis(model, market.get(), asset);
}

// GetMarket method for VarCapCVModel
void VarianceSwap::GetCVMarket(const VarCapCVModel*   cvModel,
                 const CMarketDataSP    market) {
    market->GetReferenceDate(valueDate);

    // Use both vols from the market
    CAsset::getAssetMarketData(cvModel->volModel.get(), market.get(), ccyTreatment,
                               discount, asset);

    CAsset::getAssetMarketData(cvModel->fwdModel.get(), market.get(), ccyTreatment,
                               discount, marketAsset);

    discount.getData(cvModel->fwdModel.get(), market);
    instSettle->getMarket(cvModel->fwdModel.get(), market.get());

    VarianceSwapUtil::validateBasis(cvModel->fwdModel.get(), market.get(), asset);
}

/** price given the future vol squared */
/** NOTE: we allow Model to be NULL if we don't need it */                 
void VarianceSwap::price(double varFuture, 
                         IModel* model,
                         Control* control, 
                         Results* results, 
                         VarSwapBasis::VarSwapBasisProcSP basisProc,
                         VanillaContractsRecorderSP recorder) const {
    static const string method = "VarianceSwap::price";

    if (!recorder.get() )
    {
        recorder = VanillaContractsRecorder::createVanillaOptionRecorder(control);
    }
   
    try {
        if (payoffType != "FORWARD") {
            // must be a forward
            throw ModelException(method, "This model only prices payoffType FORWARD." );
        }

        int numReturns=0;
        DateTime maturity   = samples[samples.size()-1].date;
        double   volPast    = historicalVol(numReturns);
        double   pastWeight = 0;

        pastWeight = numReturns/(double)(numTotalReturns);

        // if zero trading time, then past weight must equal 1.0
        ATMVolRequest volReq;
        CVolProcessedBSSP vol(asset->getProcessedVol(&volReq));
        double t = vol->calcTradingTime(isFwdStarting() ? samples.front().date : valueDate, samples.back().date);
        if (!Maths::isPositive(t)) {
            pastWeight = 1.0;
        }

        // now add effect due to dividends.  Note: discreteDivsAdjustment computes the effect of the future dividend
        // on the total variance, to avoid confusion with both past / future weights, and switching from a continuous
        // integration to a discrete summand.  Thus we have to normalize to determine modification to
        // varFuture
        double effectFromDivs = VarSwapUtilities::futureDiscreteDivsAdjustment(asset.get(), valueDate, samples.front().date,
            samples.back().date, observationsPerYear, numTotalReturns + 1);

        if (noDivAdj && Maths::isPositive(1.0 - pastWeight)) {
            // normalization in discreteDivsAdjustment is always observationsPerYear/numTotalReturns
            // therefore we do not have to distinguish between subtractMeanVol = true/false
            varFuture += effectFromDivs / (double)(observationsPerYear) * (double)(numTotalReturns);
        }

        double totalVolSquared = 0.0;
        double optionScalingFactor = 1.0;
        if (!Maths::isPositive(t)) {
            totalVolSquared = volPast*volPast;
        } else {
            double contributionPast = volPast*volPast*pastWeight;

            double contributionFuture = 0.0;
            if (subtractMeanVol) {
                optionScalingFactor = (double)(observationsPerYear) / (double)(numTotalReturns - 1);
                contributionFuture = varFuture * (double)(observationsPerYear) / (double)(numTotalReturns - 1);
            } else {
                optionScalingFactor = (double)(observationsPerYear) / (double)(numTotalReturns);
                contributionFuture = varFuture * (double)(observationsPerYear) / (double)(numTotalReturns);
            }

            double corrTermTotal = 0.0;
            if (subtractMeanVol && (numReturns>1)) {
                if (numReturns>0) {
                    double spotAtSamplesStart = samples[0].amount;
                    double currentSpot = asset->getSpot();
                    double fwdAtMaturity = asset->fwdValue(maturity);
                    double atmVol = vol->CalcVol(valueDate, maturity);
                    double timeToMaturity = vol->calcTradingTime(valueDate, samples.back().date);

                    double corrTerm1 = log(currentSpot/spotAtSamplesStart) * log(currentSpot/spotAtSamplesStart) / (double)(numReturns);
                    double corrTerm2 = log(fwdAtMaturity/spotAtSamplesStart) * log(fwdAtMaturity/spotAtSamplesStart)
                        - log(fwdAtMaturity/spotAtSamplesStart) * atmVol * atmVol * timeToMaturity
                        + atmVol*atmVol*timeToMaturity + pow(atmVol,4.0) * timeToMaturity*timeToMaturity / 4.0;
                    corrTerm2 /= (double)(numTotalReturns);
                    corrTermTotal = (corrTerm1 - corrTerm2) * (double)(observationsPerYear) / (double)(numTotalReturns - 1) ;
                }
            }

            totalVolSquared = contributionPast + contributionFuture + corrTermTotal;
        }

        bool dead = valueDate.isGreater(instSettle->settles(maturity, asset.get()));
        double pv = 1.0;
        if(!dead) {
            // Alive instrument
            pv = instSettle->pv(valueDate,maturity,discount.get(),asset.get());

        }

        double scaleFactor = 100.0;
        double result = pv * VarianceSwapUtil::priceVarSwapSimple(totalVolSquared, strikeVol, notional,
                scaleFactor, dontScaleByStrike);

        // Scale option portfolio
        optionScalingFactor *= 100.0 * notional * pv;
        if(!dontScaleByStrike) {
            optionScalingFactor /= 2.0 * strikeVol;
        }
        recorder->scaleNbContracts(optionScalingFactor);

        if(dead) {
            results->storePrice(0.0, discount->getCcy());
        } else {
            results->storePrice(result, discount->getCcy());
        }

        // take care of additional outputs
        double volFuture = 0.0;
        if(numTotalReturns - numPastReturns) {
            double ratio = (double)(observationsPerYear)/(double)(numTotalReturns - numPastReturns);
            if (subtractMeanVol) {
                ratio *= (double)(numTotalReturns)/(double)(numTotalReturns -1);
            }
            volFuture = sqrt( ratio * varFuture );
        }

        //Vega Matrix Lite
        if (control)
        {
            SensitivitySP sens(control->getCurrentSensitivity());
            VegaMatrixLiteSP vml(dynamic_cast<VegaMatrixLite *>(sens.get()));
            if (vml.get()) {
                CModelLN *lnModel = dynamic_cast<CModelLN *>(model);
                QLIB_VERIFY(model != NULL, "'model' should not be NULL when computing a VegaMatrixLite sensitivity.");
                QLIB_VERIFY(lnModel != NULL, "Model must derive from CModelLN to be used with VegaMatrixLite.");
                
                VanillaInfo::storeVegaMatrix(
                    vml,
                    const_cast<VarianceSwap*>(this),
                    lnModel,
                    recorder,
                    CAssetSP::constCast(asset.getSP()),
                    valueDate,
                    endDate(vml.get()),
                    results);
                
                return;
            }
        }
        
        if (control && control->isPricing()) {
            OutputRequest* request =
                control->requestsOutput(OutputRequest::VOL_IN_FUTURE);
            if (request) {
                results->storeRequestResult(request, volFuture);
            }
            request = control->requestsOutput(OutputRequest::IND_VOL);
            if (request) {
                results->storeRequestResult(request, volFuture);
            }
            request = control->requestsOutput(OutputRequest::VOL_IN_PAST);
            if (request) {
                results->storeRequestResult(request, volPast);
            }
            request = control->requestsOutput(OutputRequest::PAST_WEIGHT);
            if (request) {
                results->storeRequestResult(request, pastWeight);
            }
            request = control->requestsOutput(OutputRequest::DISCOUNT_FACTOR);
            if (request && !dead) {
                results->storeRequestResult(request, pv);
            }
            request = control->requestsOutput(OutputRequest::TOTAL_VOL);
            if (request) {
                results->storeRequestResult(request, sqrt(totalVolSquared));
            }
            request = control->requestsOutput(OutputRequest::STRIKE_VOL);
            if (request) {
                results->storeRequestResult(request, strikeVol);
            }
            request = control->requestsOutput(OutputRequest::EXPECTED_N);
            if (request) {
                results->storeRequestResult(request, numTotalReturns);
            }
            request = control->requestsOutput(OutputRequest::IMPLIED_N);
            // The number of returns 'implied' by samples
            if (request) {
                results->storeRequestResult(request, samples.size()-1);
            }
            request = control->requestsOutput(OutputRequest::PAYMENT_DATES);
            if (request) {
                DateTime paymentDate = instSettle->settles(maturity, asset.get());
                DateTimeArray date(1, paymentDate);
                OutputRequestUtil::recordPaymentDates(control,results,&date);
            }
            request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
            if (request) {
                if (valueDate.isGreaterOrEqual(maturity)) {
                    DateTime paymentDate = instSettle->settles(maturity, asset.get());
                    CashFlow cf(paymentDate, result/pv);
                    CashFlowArray cfl(1, cf);
                    OutputRequestUtil::recordKnownCashflows(control,
                        results,
                        discount->getCcy(),
                        &cfl);
                }
            }

            // FSA regulatory numbers
            request = control->requestsOutput(OutputRequest::FSA_VALUE);
            if (request) {
                // Fair value
                // We don't know that FV is appropriate for caps but it seems reasonable
                results->storeRequestResult(request, result);
            }

            if (payoffType == "FORWARD") {
                // We don't know what to report for caps
                request = control->requestsOutput(OutputRequest::FSA_PRR);
                if (request) {
                    // Position Risk Requirement
                    double reserve = 0.0;
                    double totalVol = sqrt(totalVolSquared);
                    if (totalVol < 0.5) {
                        reserve = 0.08;
                    } else if (totalVol < 0.75) {
                        reserve = 0.12;
                    } else /*if (totalVol < 1.0)*/ {
                        reserve = ceil(2.0*totalVol)*0.08;   // (16% if TOTAL_VOL in (75%, 100%), 24% if in (100%, 150%) etc.)
                    }

                    reserve *= scaleFactor * notional;
                    // If scaling by strike the result will be multiplied by vega notional not variance notional.
                    // So need to scale down by 2K
                    if (!dontScaleByStrike) {
                        reserve /= 2.0 * strikeVol;
                    }

                    results->storeRequestResult(request, reserve);
                }
            }

            // Record basis information at maturity i.e. ignore forward starting VarSwaps for here
            if(basisProc.get()) {
                double tradYear = vol->calcTradingTime(valueDate, samples.back().date);
                basisProc->recordRequests(tradYear, control, results);
            }

            // DELAY_PRICE
            InstrumentUtil::delayPriceHelper(control,
                results,
                result,
                valueDate,
                discount.get(),
                asset.get(),
                premiumSettle.get());
            // FWD_AT_MAT
            InstrumentUtil::recordFwdAtMat(control,
                results,
                maturity,
                valueDate,
                asset.get());

            // output div effect to debug packet
            OutputNameConstSP divEffOut(new OutputName("effectFromDivs"));
            results->storeGreek(CDoubleSP(CDouble::create(effectFromDivs)), Results::DEBUG_PACKET, divEffOut);
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** price a dead instrument until settlement - exercised, expired, knocked out etc.
    returns true if it is dead (and priced), false if it is not dead */
bool VarianceSwap::priceDeadInstrument(CControl* control, CResults* results) const {
    static const string method = "VarianceSwap::priceDeadInstrument";

    DateTime maturity = samples[samples.size()-1].date;
    DateTime settles  = instSettle->settles(maturity, asset.get());

    if (maturity.isGreater(valueDate)) {
        return false;
    }
#if 1
    else if (valueDate.isGreaterOrEqual(settles)) {
#else
    else if (valueDate.isGreater(settles)) {
#endif
        if (payoffType == "FORWARD") {
            // Call price method even for settled instruments so
            // we can get the OutputRequests
            price(0.0, NULL, control, results);
        }
        results->storePrice(0.0, discount->getCcy());
        return true;
    }
    else {
        if (payoffType == "FORWARD") {
            price (0.0, NULL, control, results);
        } else if (payoffType == "FORWARD_CAP") {
            int numReturns;
            double volPast = historicalVol(numReturns);
            double price = Maths::max(volPast*volPast - cap*strikeVol*strikeVol,0.);
            double pv = instSettle->pv(valueDate,
                                       maturity,
                                       discount.get(),
                                       asset.get());
            price *= pv*notional*100.0;
            if (!dontScaleByStrike) {
                price /= (2.0*strikeVol);
            }
            results->storePrice(price, discount->getCcy());
            if (control && control->isPricing()) {
                OutputRequest* request = control->requestsOutput(OutputRequest::VOL_IN_PAST);
                if (request) {
                    results->storeRequestResult(request, volPast);
                }
                request = control->requestsOutput(OutputRequest::PAST_WEIGHT);
                if (request) {
                    results->storeRequestResult(request, 1.);
                }
                request = control->requestsOutput(OutputRequest::VOL_IN_FUTURE);
                if (request) {
                    results->storeRequestResult(request, 0.);
                }
                request = control->requestsOutput(OutputRequest::TOTAL_VOL);
                if (request) {
                    results->storeRequestResult(request, volPast);
                }
                request = control->requestsOutput(OutputRequest::DISCOUNT_FACTOR);
                if (request) {
                    results->storeRequestResult(request, pv);
                }

                // KNOWN_CASHFLOWS
                request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
                if (request) {
                    DateTime paymentDate = instSettle->settles(maturity, asset.get());
                    CashFlow cf(paymentDate, price / pv);
                    CashFlowArray cfl(1, cf);
                    OutputRequestUtil::recordKnownCashflows(control,
                        results,
                        discount->getCcy(),
                        &cfl);
                }

                // PAYMENT_DATES
                request = control->requestsOutput(OutputRequest::PAYMENT_DATES);
                if (request) {
                    DateTime paymentDate = instSettle->settles(maturity, asset.get());
                    DateTimeArray date(1, paymentDate);
                    OutputRequestUtil::recordPaymentDates(control,results,&date);
                }
            }
        } else {
            throw ModelException(method, "Unknown payoff type " + payoffType);
        }

        return true;
    }
}


/** Weight for call / put options portfolio */
class VarianceSwap::VarSwapIntegrandWeight: public Function1DDouble {
public:
    /** Full constructor */
    VarSwapIntegrandWeight(const Range& integrationDomain):
    Function1DDouble(integrationDomain) {}

    /** Implements 1/k weights */
    virtual double operator()(double relativeStrike) const {
        double relativeStrikeWeight = 1.0 / Maths::square(relativeStrike);
        return relativeStrikeWeight;
    }
};


// Constructor used by VanillaVarSwap
VarianceSwap::VarianceSwap(const DateTime& valueDate,
             const DateTime& startDate,
             const bool fwdStarting,
             const double initialSpot,
             const string& ccyTreatment,
             InstrumentSettlementSP instSettle,
             InstrumentSettlementSP premiumSettle,
             CAssetWrapper asset,
             YieldCurveWrapper discount,
             CashFlowArray samples,
             const double strikeVol,
             const int observationsPerYear,
             const bool subtractMeanVol,
             const string& payoffType,
             const bool dontScaleByStrike,
             const double cap,
             const bool noDivAdj,
             const bool divAdjOnExDate,
             const bool isVanilla,
             const int numPastReturns,
             const int numTotalReturns) : VolVarShell(TYPE) {

    // Generic1Factor
    this->valueDate = valueDate;
    this->startDate = startDate;
    this->fwdStarting = fwdStarting;
    this->oneContract = true;
    this->notional = 1.0;
    this->initialSpot = initialSpot;
    this->ccyTreatment = ccyTreatment;
    this->instSettle = instSettle;
    this->premiumSettle = premiumSettle;
    this->asset = asset;
    this->discount = discount;

    // VolVarShell
    this->samples = samples;
    this->strikeVol = strikeVol;
    this->observationsPerYear = observationsPerYear;
    this->subtractMeanVol = subtractMeanVol;
    this->payoffType = payoffType;
    this->dontScaleByStrike = dontScaleByStrike;
    this->cap = cap;
    this->noDivAdj = noDivAdj;
    this->divAdjOnExDate = divAdjOnExDate;
    this->isVanilla = isVanilla;
    this->numPastReturns = numPastReturns;
    this->numTotalReturns = numTotalReturns;
}

double VarianceSwap::varFromCalls(double                            lowStrike,
                                  double                            highStrike,
                                  int                               nbSteps,
                                  const DateTime&                   maturity,
                                  double                            fwd,
                                  string                            integrationMethod,
                                  const bool                        allowNegFwdVar,
                                  double                            absPrecisionAdjusted,
                                  double                            relPrecisionAdjusted,
                                  bool                              isDoingVegaMatrix,
                                  VarSwapBasis::VarSwapBasisProcSP  basisProc,
                                  VanillaContractsRecorderSP        recorder) const {
    static const string method = "VarianceSwap::varFromCalls";
    try {
        double varFuture = 0.0;
        LinearStrikeVolRequestSP volRequest(new LinearStrikeVolRequest(0.0,
                                                                       valueDate,
                                                                       maturity,
                                                                       false));
        // by default, allowNegFwdVar is set to false, so we have to update
        volRequest->allowNegativeFwdVar(allowNegFwdVar);

        if (isDoingVegaMatrix) {
            integrationMethod = ClosedFormIntegrateLN::DISCRETE_INTEGRATION;
        }

        // Determine price cutoff: assume single factor - this would have already failed for NFactors

        double cutoff = 0.0;
        double tradYear = 0.0;
        if(basisProc.get()) {
            ATMVolRequest volReq;
            CVolProcessedBSSP vol(asset->getProcessedVol(&volReq));
            tradYear = vol->GetTimeMetric()->yearFrac(valueDate,maturity);
            cutoff = basisProc->interpPriceCutoff(tradYear);
        }

        // Update lowStrike
        lowStrike = Maths::max(lowStrike, cutoff);

        if(CString::equalsIgnoreCase(integrationMethod, ClosedFormIntegrateLN::DISCRETE_INTEGRATION)) {
            if(nbSteps<2) {
                nbSteps = 10; // arbitrary
            }
            if(lowStrike<=0.0) {
                throw ModelException("The low Strike must be strictly positive!");
            }

            double dK = (highStrike - lowStrike) / ((double) (nbSteps - 1));

            double integratedValue = 0.0;
            double curLongPosition = 1.0/lowStrike;
            double callPaySoFar = 0.0;

            // Create VolInterp once for efficiency
            // choose how to interpolate the vol - go for traditional route for now

            double scalingFactor = - recorder->getScalingFactor() * 2.0;
            recorder->setScalingFactor(scalingFactor);

            for (int i=0; i<nbSteps-1; i++) {
                // Set the strike
                double curStrike = lowStrike + dK*(double)i;
                volRequest->setStrike(curStrike);

                // interpolate the vol for each strike
                CVolProcessedBSSP volBS(asset->getProcessedVol(volRequest.get()));
                // double variance = volBS->CalcVar(valueDate, maturity);
                double vol = volBS->CalcVol(valueDate, maturity);
                double yearFrac = volBS->calcTradingTime(valueDate, maturity);
                double variance = Maths::square(vol) * yearFrac;
                if (Maths::isNegative(variance)) {
                    throw ModelException("negative variance");
                }

                // Nb of calls
                double nCalls = (log((curStrike+dK)/lowStrike) - curLongPosition*(curStrike + dK - lowStrike) + callPaySoFar)/dK;
                callPaySoFar = callPaySoFar + nCalls*(curStrike - lowStrike);
                curLongPosition += nCalls;
                
                // finally call Black model
#if 0
                double curCall = Black::price(true, fwd, curStrike, 1., variance);
                integratedValue += curCall * nCalls;
#endif
                double thisPrice = BlackPrice(true, fwd, curStrike, 1.0, vol, yearFrac, 
                    nCalls, maturity, recorder);

                integratedValue += thisPrice;
            }

            varFuture = - 2.0*(log(lowStrike/fwd)-1 + fwd/lowStrike + integratedValue);
        } else {
            // create an integrator by name and pass on precision in terms of integrator
            Integrator1DSP myIntegrator = Integrator1D::createIntegrator(integrationMethod,
                                                                         absPrecisionAdjusted,
                                                                         relPrecisionAdjusted);
            // Create range
            refCountPtr<Range> interval;
            if (IntFuncInf::TYPE->isInstance(myIntegrator)) {
                // for IMSL, create open zero infinity interval
                interval = refCountPtr<Range>(new Range(OpenBoundary(cutoff/fwd), Infinity(Infinity::Plus)));
            } else {
                // otherwise, create closed lowstrike/fwd highstrike/fwd interval
                // worst case: the integrator fails, if it expects an open boundary
                interval = refCountPtr<Range>(new Range(ClosedBoundary(lowStrike/fwd), ClosedBoundary(highStrike/fwd)));
            }
            // create integrand and integrate
            double scalingFactor = recorder->getScalingFactor() * 2.0;
            recorder->setScalingFactor(scalingFactor);
            VarianceSwap::VarSwapIntegrandWeightConstSP optionWeights(new 
                VarianceSwap::VarSwapIntegrandWeight(*interval));
            StaticReplicationIntegrandSP integrand(new StaticReplicationIntegrand(
                optionWeights, recorder, asset.getSP(), valueDate, maturity, allowNegFwdVar));
            varFuture = 2 * myIntegrator->integrate(*integrand);
            // VarSwapIntegrand myIntegrand(*interval, fwd, valueDate, maturity, asset, volRequest);
            // varFuture = 2 * myIntegrator->integrate(myIntegrand);
        }

        if(basisProc.get() && Maths::isPositive(tradYear)) {
            // Deal with VarSwapBasis now
            double volBasis = basisProc->interpVolBasis(tradYear);

            // Compute number of future returns
            // DO NOT USE numTotalReturns here: the instrument does not get updated when it is fwdStarting !!!
            // Compute number of future returns using holidays
            TimeMetricConstSP metric = basisProc->GetTimeMetric();
            HolidayWrapper hols = metric->getHolidays();
            int numFutureReturns = hols->businessDaysDiff(valueDate, maturity);

            // Convert basis from trading time to PPY time
            double ppyTime = 0.0;
            if (subtractMeanVol && numFutureReturns > 0) {
                ppyTime = (double)(numFutureReturns - 1) / 252.0;
            } else {
                ppyTime = (double)(numFutureReturns) / 252.0;
            }

            // Compute correction due to basis
            double sqrtVarBasis = volBasis * sqrt(ppyTime);
            double basisCorrection = sqrtVarBasis * sqrtVarBasis + 2.0 * sqrtVarBasis * sqrt(varFuture);
            varFuture += basisCorrection;
        }

        return varFuture;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


// ** using Fourier method for SVJ **
template <class Process, class Product>
class VarOptIntegrand:public Function1DDouble {
public:
    VarOptIntegrand(const  Process& process,
                    const  Product& product,
                    const  DateTime& maturity,
					bool   isVolSwap,
                    double strike):
        Function1DDouble(Range(OpenBoundary(0.0), Infinity(Infinity::Plus))),
        process(process), product(product), matdate(maturity), isVolSwap(isVolSwap), strike(strike) {}

    double operator()(double  x) const {
		if(isVolSwap) {
			// based on Matytsin 2000 Merrill Lynch presentation for VolatilitySwap
			return (1.0 - exp(process.cumulant(product, -x, matdate) + x * strike).real()) / (x * sqrt(x));
		}
		else {
#ifdef RG
			static double omega = +0.05;    // XXX make this an input
			// based on similar derivation to that in 'Single Integrand Approach' (Section 2.1.1.)
			// Venardos (2002) 'The Theory Behind the Fourier Engine'
			Complex z(omega, x);
			return (exp(process.cumulant(product, z, matdate) + z * strike) / (z * z)).real();
#else
			// based on Matytsin 2000 Merrill Lynch presentation
			Complex z(0.0, x);
			return ( (1.0 - exp(process.cumulant(product, z, matdate)-z*strike))/(x * x) ).real();
#endif
		}
    }

private:

    const  Process& process;
    const  Product& product;
    DateTime        matdate;
	bool            isVolSwap;
    double          strike;
};


// Fourier payoff

// equivalent to InstIntoFourierProduct
VarOptFP::VarOptFP(const VarianceSwap* inst,
		           bool                isVolSwap,
 				   bool			   	   useMatytsinOnly,
                   const DateTime&     matDate):    // maturity date
    FourierProduct(inst->asset.get(),
                   inst->valueDate,
                   inst->discount.get(),
                   inst->instSettle.get()),
    inst(inst),
    maturity(matDate),
	isVolSwap(isVolSwap),
    useMatytsinOnly(useMatytsinOnly) {
    static const string method = "VarOptFP::VarOptFP";

	if ((inst->payoffType != "FORWARD_CAP") && (!isVolSwap)) {
		// forward starting not allowed
		throw ModelException(method, "Fourier models only allowed for variance caps." );
	}

}

void VarOptFP::validate(FourierProcess* process) const{
    const static string method("VarOptFP::validate");
    int i = 0;
    for(; i < mAsset->NbAssets(); i++) {
        string CcyTtreatment(mAsset->assetGetCcyTreatment(i));
        if(CcyTtreatment == "S") {
            throw ModelException(method,
                                 "Currency Struck assets not supported.");
        }
    }
    // give the process the chance to do some extra validation
    // and to initialize its transient fields
    process->validate(this);
}

	void VarOptFP::price(const FourierEngine* model,
                   Control*             control,
                   Results*             results){

    // valueDate >= matDate is taken care of here
    if(inst->priceDeadInstrument(control, results)){
        return; // dead instrument priced
    }

    FourierProduct::price(model,
                          control,
                          results);
}

/** Constructs integrands for var option */
Function1DDoubleArrayConstSP VarOptFP::Integrand(const FourierEngine * model,
                                       const Integrator1D*  integrator) {
    static const string method = "VarOptFP::Integrand";

    try{
        const FourierProcess& process = model->getProcess();

        if (inst->subtractMeanVol) {
            totalYears = (double)(inst->numTotalReturns - 1)/(double)inst->observationsPerYear;
        } else {
             totalYears = (double)inst->numTotalReturns/(double)inst->observationsPerYear;
        }

        int      numReturns;
        volPast    = inst->historicalVol(numReturns);
        pastWeight = (double)numReturns/(double)(inst->numTotalReturns);

        double histVar = volPast * volPast * pastWeight * totalYears;
        strike = inst->cap * inst->strikeVol * inst->strikeVol * totalYears - histVar;

        if (inst->noDivAdj && Maths::isPositive(1.0 - pastWeight)) {

            // now add effect due to dividends.  Note: discreteDivsAdjustment computes the effect of the future dividend
            // on the total variance, to avoid confusion with both past / future weights, and switching from a continuous
            // integration to a discrete summand.  Thus we have to normalize to determine modification to
            // varFuture
            double effectFromDivs = VarSwapUtilities::futureDiscreteDivsAdjustment(inst->asset.get(), inst->valueDate,
                inst->samples.front().date, inst->samples.back().date, inst->observationsPerYear, inst->numTotalReturns + 1);

            // normalization in discreteDivsAdjustment is always observationsPerYear/numTotalReturns
            // therefore we do not have to distinguish between subtractMeanVol = true/false
            strike -= effectFromDivs / (double)(inst->observationsPerYear) * (double)(inst->numTotalReturns);
        }

        Function1DDoubleArraySP functions(new Function1DDoubleArray(1));

        if (inst->isFwdStarting()){
            years = process.getTimeMetric().yearFrac(inst->samples.front().date, maturity);
            const FwdStFourierProductQuadVar& thisProd = *this;
            const FwdStFourierProcessQuadVar* thisProc = dynamic_cast<const FwdStFourierProcessQuadVar*>(&process);

            if(!thisProc) {
                throw ModelException(method, "Process does not support FwdStFourierProcessQuadVar interface." );
            }

            (*functions)[0] = Function1DDoubleSP(new VarOptIntegrand<FwdStFourierProcessQuadVar, FwdStFourierProductQuadVar>
                                                      (*thisProc,
                                                       thisProd,
                                                       maturity,
												       isVolSwap,
                                                       strike));
        } else {
            years = process.getTimeMetric().yearFrac(inst->valueDate, maturity);
            const StFourierProductQuadVar& thisProd = *this;
            const StFourierProcessQuadVar* thisProc = dynamic_cast<const StFourierProcessQuadVar*>(&process);

            if(!thisProc) {
                throw ModelException(method, "Process does not support StFourierProcessQuadVar interface." );
            }

            (*functions)[0] = Function1DDoubleSP(new VarOptIntegrand<StFourierProcessQuadVar, StFourierProductQuadVar>
                                                      (*thisProc,
                                                       thisProd,
                                                       maturity,
												       isVolSwap,
                                                       strike));
        }

        return functions;
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

/** Post process method for VolSwapFP */
void VarOptFP::postResults(const FourierEngine* model,
                           const Integrator1D*  integrator,
                           const FourierProductIntegrator1D::IntegralArray& integrals,
                           CControl*            control,
                           CResults*            results) {
    static const string method = "VarOptFP::postResults";
    double price = 0.0;

    try {
        const FourierProcess& process = model->getProcess();
        double varHM;
        
        double divVar = 0.0;
        if (inst->isFwdStarting()){
            const FwdStFourierProductQuadVar& thisProd = *this;
            const FwdStFourierProcessQuadVar* thisProc = dynamic_cast<const FwdStFourierProcessQuadVar*>(&process);
            if(!thisProc) {
                throw ModelException(method, "Process does not support FwdStFourierProcessQuadVar interface" );
            }
            varHM = thisProc->expectation(thisProd, maturity);
        } else {
            const StFourierProductQuadVar& thisProd = *this;
            const StFourierProcessQuadVar* thisProc = dynamic_cast<const StFourierProcessQuadVar*>(&process);
            if(!thisProc) {
                throw ModelException(method, "Process does not support StFourierProcessQuadVar interface" );
            }
            varHM = thisProc->expectation(thisProd, maturity);
        }
        //future realized variance or quadratic variation (NORMALIZED BY TIME)
        double futVar = varHM * years;
		if (isVolSwap) {
			price = integrals[0] / (2 * sqrt(Maths::PI));
			if (!useMatytsinOnly) {
				//price -= sqrt(varHM * years - strike);
                price -= sqrt(futVar - strike); 
			}
			price /= sqrt(totalYears);
		}
		else {
			price = integrals[0] / Maths::PI;
			//price = price +.5*(varHM * years - strike);
    		price = price +.5*(futVar - strike); 
            price /= totalYears;
		}

        // See if option price is positive
        if((price <= 0.0) && (!isVolSwap)) {
            if(price > -0.01) {
                price = 0.0;
            } else {
                throw ModelException("(Scaleless) Price of variance  cap is negative " +
                                     Format::toString(price) +
                                     ". Check integrator parameters e.g. tolerance, nbSteps.");
            }
        }

        // Keep some info aside before we apply notions, scalings and pv
        // E ( AvgTotalVar - Strike^2 * cap )
        double unscaledPrice = price;
        // E (AvgTotalVar)

        // Do dividend adjustments
        if(inst->noDivAdj) {
            // No dividend adjustment = use full return = add the effect of divs to vanillas
            divVar = VarSwapUtilities::futureDivsVariance(
               inst->asset.get(), inst->valueDate, inst->valueDate, maturity);
        }
        
        //pastVar is realized Var accumulated in the past (volPast is normalized, hence multiply by time)
        //pastWeight * totalYears = discrete fraction of time from first obsDate to valueDate
        double pastVar = volPast * volPast * pastWeight * totalYears;
        
        //futVar is future integrated variance or future quadratic variation
        //obtained using expectation from the relevant stochVol model (expectation is normalized, hence multiply by time)  
        //years = continuous fraction of time from first obsDate to maturity
        /*
        double futVar = varHM * years;
        */

        //divVar is future variance due to future dividends (specific to quadVar)
        //(no need to multiply by time, sum log(Retun)^2 has the dimension of a variance)
        
        //total realized + future var (or quadVar) normalized by fraction of years
        //double unscaledFwd   = ((varHM + divVar) * years + volPast * volPast * pastWeight * totalYears) / totalYears;
        double unscaledFwd   = (pastVar + futVar + divVar) / totalYears;

        double pv = inst->instSettle->pv(inst->valueDate,
                                         maturity,
                                         inst->discount.get(),
                                         inst->asset.get());

        price *= pv*inst->notional*100.0;

        if (!inst->dontScaleByStrike) {
            price /= (2.0*inst->strikeVol);
        }

        results->storePrice(price, discount->getCcy());

        // take care of additional outputs
        if (control && control->isPricing()) {
            // Past Vol
            OutputRequest* request = control->requestsOutput(OutputRequest::VOL_IN_PAST);
            if (request) {
                results->storeRequestResult(request, volPast);
            }

            // Past Weight
            request = control->requestsOutput(OutputRequest::PAST_WEIGHT);
            if (request) {
                results->storeRequestResult(request, pastWeight);
            }

            /* Compute Past and Future Vols */

            //contributionPast is pastVar normalized by ExpNPast / PPY ie volPast * volPast * pastWeight
            double contributionPast = volPast*volPast*pastWeight;
            //contributionFuture
            double contributionFuture = 0.0;
            if(inst->numTotalReturns - inst->numPastReturns) {
                // Has future
                // Note if pricing on maturity date we have numPastReturns = expectedN
                // The way time is divided means that VOL_IN_FUTURE = 0 and TOTAL_VOL = VOL_IN_PAST in this case
                // varHM * years is future realized variance or quadratic variation / divVar is future variance due to dividends
                if (inst->subtractMeanVol) {
                    contributionFuture = (futVar+ divVar) * 
                        (double)inst->observationsPerYear/(double)(inst->numTotalReturns - 1);
                } else {
                    contributionFuture = (futVar + divVar) * 
                        (double)inst->observationsPerYear/(double)(inst->numTotalReturns);
                }
            }

            // Future Vol
            request = control->requestsOutput(OutputRequest::VOL_IN_FUTURE);
            if (request) {
                    
                double volFuture = 0.0;
                if(inst->numTotalReturns - inst->numPastReturns) {
                    // Has future
                    // expN / futureN
                    double adj = (double)inst->numTotalReturns/(double)(inst->numTotalReturns - inst->numPastReturns);
                    volFuture = sqrt(contributionFuture*adj);
                }
                results->storeRequestResult(request, volFuture);
            }

			// IND_VOL: not support such that IND_VOL can be reported for (composite) capped variance swap.

            double totalVolSquared = contributionPast + contributionFuture;

            // Total Vol
            request = control->requestsOutput(OutputRequest::TOTAL_VOL);
            if (request) {
                results->storeRequestResult(request, sqrt(totalVolSquared));
            }

            // Discount
            request = control->requestsOutput(OutputRequest::DISCOUNT_FACTOR);
            if (request) {
                results->storeRequestResult(request, pv);
            }

            // Expected N
            request = control->requestsOutput(OutputRequest::EXPECTED_N);
            if (request) {
                results->storeRequestResult(request, inst->numTotalReturns);
            }

            // IMPLIED_N. The number of returns 'implied' by samples
            request = control->requestsOutput(OutputRequest::IMPLIED_N);
            if (request) {
                results->storeRequestResult(request, inst->samples.size()-1);
            }

/*
            OLD CODE, kept for the moment
        
            // Option Implied Vol
            request = control->requestsOutput(OutputRequest::VAR_OPTION_IND_VOL);
            if (!isVolSwap && request) {
                if(Maths::isPositive(years)) {
                    double variance = 0.0;
                    try {
                        bool success = Black::impliedVariance(true,                                            // option call or put
                                                              unscaledFwd,                                     // forward price
                                                              inst->cap * inst->strikeVol * inst->strikeVol,   // strike
                                                              1.0,                                             // pv
                                                              0.3 * 0.3 * years,                               // initial var guess
                                                              unscaledPrice,                                   // option price
                                                              2.0 * 0.3 * 1.0e-5 * years,                      // var accuracy
                                                              variance);

                        if(success) {
                            double impliedVol = sqrt(variance / years);
                            results->storeRequestResult(request, impliedVol);
                        }
                    } catch(exception&) {
                        // Don't bother
                    }
                }
            }

*/

            //GAD, correction for IND_VOL
            // Option Consistent Implied Vol (based on future realized variance / past realized variance accounted for by modifying strike)
            request = control->requestsOutput(OutputRequest::VAR_OPTION_IND_VOL);
            if (request) {
                if(Maths::isPositive(years)) {
                    double variance = 0.0;
                    try{
                        //consistent values
                        //consistent fwd is based on future realized variance
                        double consFwdPrice = (unscaledFwd - pastVar / totalYears) / (1.0 - pastWeight) ; //past var has been taken out and put in the strike
                        //consistent strike accounts for past realized variance
                        double consStrike = (inst->cap * inst->strikeVol * inst->strikeVol - pastVar / totalYears) / (1.0 - pastWeight) ;   
                        //consistent price is scaled by pastN
                        double consUnscaledPrice = unscaledPrice / (1.0 - pastWeight);

                        //test on strike, if negative option is worthless and we give zero implied vol
                        if(Maths::isPositive(consStrike)){
                            bool success = Black::impliedVariance(true,                                  // option call or put
                                                                consFwdPrice,//changed                           // forward price
                                                                consStrike, //changed                           // strike
                                                                1.0,                                             // pv
                                                                0.3 * 0.3 * years,                               // initial var guess
                                                                consUnscaledPrice, //changed                     // option price
                                                                2.0 * 0.3 * 1.0e-5 * years,                      // var accuracy
                                                                variance);
                            if(success) {
                                double impliedVol = sqrt(variance / years);
                                results->storeRequestResult(request, impliedVol);
                            }
                        }
                    }
                    catch(exception&) {
                        // Don't bother
                    }
                }
            }
            //end of GAD

			// payment_dates
			request = control->requestsOutput(OutputRequest::PAYMENT_DATES);
			if (request) {
				DateTime maturity   = inst->samples[inst->samples.size()-1].date;
				DateTime paymentDate = inst->instSettle->settles(maturity, inst->asset.get());
				DateTimeArray date(1, paymentDate);
				OutputRequestUtil::recordPaymentDates(control,results,&date);
			}

            // DELAY_PRICE
            InstrumentUtil::delayPriceHelper(control,
                                             results,
                                             price,
                                             inst->valueDate,
                                             inst->discount.get(),
                                             inst->asset.get(),
                                             inst->premiumSettle.get());

            try {
                // FWD_AT_MAT
                InstrumentUtil::recordFwdAtMat(control,
                                               results,
                                               maturity,
                                               inst->valueDate,
                                               inst->asset.get());
            } catch(exception& e) {
                // Don't rethrow
                ModelException ee = ModelException::addTextToException(e, method + ": Computation of forward price failed");

                IObjectSP fwdObj(new Untweakable(ee));

                results->storeRequestResult(control->requestsOutput(OutputRequest::FWD_AT_MAT),
                                            fwdObj,
                                            OutputNameSP(new OutputName(inst->asset.get()->getName())));
            }
        }
    }

    catch (exception& e){
        throw ModelException(e, method);
    }
}

// ** using Fourier method to price variance cap square payoff, E[((V-K)+)^2]
template <class Process, class Product>
class VarOptSqrIntegrand:public Function1DDouble {
public:
    VarOptSqrIntegrand(const Process& process,
                       const Product& product,
                       const DateTime& maturity,
                       const double& strike,
                       const double mean):
        Function1DDouble(Range(OpenBoundary(0.0), Infinity(Infinity::Plus))),
        process(process), product(product), matdate(maturity), strike(strike), mean(mean) {}

    double operator()(double  x) const {
#ifdef RG
        static double omega = +0.05;    // XXX make this an input
        // based on similar derivation to that in 'Single Integrand Approach' (Section 2.1.1.)
        // Venardos (2002) 'The Theory Behind the Fourier Engine'. Not tested
        Complex z(omega, x);
        //return (exp(process.cumulant(product, z, matdate) - z * strike) / (z*(z-1.)*(z-2.))).real();
#else
        Complex z(0.0, x);
        return ( (2.*mean-strike)*x - (exp(process.cumulant(product, z, matdate))*(1.+exp(-z*strike))).imag() )/(x * x * x);
#endif
    }

private:

    const Process& process;
    const Product& product;
    DateTime       matdate;
    double         strike;
    double         mean;
};

// Fourier payoff for CV adjustment
class VarOptCVFP: public FourierProduct,
                  public FourierProductIntegrator1D,
                  public StFourierProductQuadVar,
                  public FwdStFourierProductQuadVar {
private:
    const VarianceSwap* inst;
    DateTime            maturity;
    double              mult;
    double              years;
    double              totalYears;
    double              strike;
    double              volPast;
    double              pastWeight;

public:
    static const string VAR_SWAP;
    static const string VAR_SWAP_SQR;
    static const string VAR_CAP;
    static const string VAR_CAP_SQR;
    static const string MODEL_SWAP;
    static const string STRIKE;

    // equivalent to InstIntoFourierProduct
    VarOptCVFP(const VarianceSwap*       inst,
               const DateTime&           matDate):    // maturity date
        FourierProduct(inst->asset.get(),
                       inst->valueDate,
                       inst->discount.get(),
                       inst->instSettle.get()),
        inst(inst),
        maturity(matDate){
        static const string method = "VarOptCVFP::VarOptCVFP";

        if (inst->payoffType != "FORWARD_CAP") {
            // forward starting not allowed
            throw ModelException(method, "Fourier models only allowed for variance caps." );
        }

    }

    virtual void validate(FourierProcess* process) const {
        const static string method("VarOptCVFP::validate");
        int i = 0;
        for(; i < mAsset->NbAssets(); i++) {
            string CcyTtreatment(mAsset->assetGetCcyTreatment(i));
            if(CcyTtreatment == "S") {
                throw ModelException(method,
                                     "Currency Struck assets not supported.");
            }
        }
        // give the process the chance to do some extra validation
        // and to initialize its transient fields
        process->validate(this);
    }

    virtual void price(const FourierEngine* model,
                       Control*             control,
                       Results*             results){

        // valueDate >= matDate is taken care of here
        if(inst->priceDeadInstrument(control, results)){
            return; // dead instrument priced
        }

        FourierProduct::price(model,
                              control,
                              results);
    }

    /** Constructs integrands for var option */
    Function1DDoubleArrayConstSP Integrand(const FourierEngine * model,
                                           const Integrator1D*  integrator) {
        static const string method = "VarOptCVFP::Integrand";

        try{
            const FourierProcess& process = model->getProcess();

            totalYears = (double)inst->numTotalReturns/(double)inst->observationsPerYear;

            int      numReturns;
            volPast    = inst->historicalVol(numReturns);
            pastWeight = (double)numReturns/(double)(inst->numTotalReturns);

            double histVar = volPast * volPast * pastWeight * totalYears;
            strike = inst->cap * inst->strikeVol * inst->strikeVol * totalYears - histVar;

            if (inst->noDivAdj && Maths::isPositive(1.0 - pastWeight)) {

                // now add effect due to dividends.  Note: discreteDivsAdjustment computes the effect of the future dividend
                // on the total variance, to avoid confusion with both past / future weights, and switching from a continuous
                // integration to a discrete summand.  Thus we have to normalize to determine modification to
                // varFuture
                double effectFromDivs = VarSwapUtilities::futureDiscreteDivsAdjustment(inst->asset.get(), inst->valueDate,
                    inst->samples.front().date, inst->samples.back().date, inst->observationsPerYear, inst->numTotalReturns + 1);

                // normalization in discreteDivsAdjustment is always observationsPerYear/numTotalReturns
                // therefore we do not have to distinguish between subtractMeanVol = true/false
                strike -= effectFromDivs / (double)(inst->observationsPerYear) * (double)(inst->numTotalReturns);
            }

            Function1DDoubleArraySP functions(new Function1DDoubleArray(3));

            if (inst->isFwdStarting()){
                const FwdStFourierProductQuadVar& thisProd = *this;
                const FwdStFourierProcessQuadVar* thisProc = dynamic_cast<const FwdStFourierProcessQuadVar*>(&process);
                if(!thisProc) {
                    throw ModelException(method, "Process does not support FwdStFourierProcessQuadVar interface." );
                }
                years = process.getTimeMetric().yearFrac(inst->samples.front().date, maturity);
                double  mean = thisProc->expectation(thisProd, maturity) * years;

                //variance cap payoff: E[(V-K)+]
                (*functions)[0] = Function1DDoubleSP(new VarOptIntegrand<FwdStFourierProcessQuadVar, FwdStFourierProductQuadVar>
                                                          (*thisProc,
                                                           thisProd,
                                                           maturity,
													       false,
                                                           strike));

                //variance cap square payoff: E[((V-K)+)^2]
                (*functions)[1] = Function1DDoubleSP(new VarOptSqrIntegrand<FwdStFourierProcessQuadVar, FwdStFourierProductQuadVar>
                                                          (*thisProc,
                                                           thisProd,
                                                           maturity,
                                                           strike,
                                                           mean));

                //variance square payoff: E[V^2]
                (*functions)[2] = Function1DDoubleSP(new VarOptSqrIntegrand<FwdStFourierProcessQuadVar, FwdStFourierProductQuadVar>
                                                          (*thisProc,
                                                           thisProd,
                                                           maturity,
                                                           0.0,
                                                           mean));
            } else {
                const StFourierProductQuadVar& thisProd = *this;
                const StFourierProcessQuadVar* thisProc = dynamic_cast<const StFourierProcessQuadVar*>(&process);
                if(!thisProc) {
                    throw ModelException(method, "Process does not support StFourierProcessQuadVar interface." );
                }
                years = process.getTimeMetric().yearFrac(inst->valueDate, maturity);
                double  mean = thisProc->expectation(thisProd, maturity) * years;

                //variance cap payoff: E[(V-K)+]
                (*functions)[0] = Function1DDoubleSP(new VarOptIntegrand<StFourierProcessQuadVar, StFourierProductQuadVar>
                                                          (*thisProc,
                                                           thisProd,
                                                           maturity,
													       false,
                                                           strike));

                //variance cap square payoff: E[((V-K)+)^2]
                (*functions)[1] = Function1DDoubleSP(new VarOptSqrIntegrand<StFourierProcessQuadVar, StFourierProductQuadVar>
                                                          (*thisProc,
                                                           thisProd,
                                                           maturity,
                                                           strike,
                                                           mean));

                //variance square payoff: E[V^2]
                (*functions)[2] = Function1DDoubleSP(new VarOptSqrIntegrand<StFourierProcessQuadVar, StFourierProductQuadVar>
                                                          (*thisProc,
                                                           thisProd,
                                                           maturity,
                                                           0.0,
                                                           mean));
            }

            return functions;
        }
        catch (exception& e){
            throw ModelException(e, method);
        }
    }

    const DateTime& getStartDate() const {return inst->samples[0].date;}

    /** Post process method for VolSwapFP */
    virtual void postResults(const FourierEngine* model,
                             const Integrator1D*  integrator,
                             const FourierProductIntegrator1D::IntegralArray& integrals,
                             CControl*            control,
                             CResults*            results) {
        static const string method = "VarOptCVFP::postResults";

        try {
#ifdef RG
#else
            double price = integrals[0]/Maths::PI;
            double varCapSqr = integrals[1]*2.0/Maths::PI;
            double varSwapSqr = integrals[2]*2.0/Maths::PI;

            const FourierProcess& process = model->getProcess();
            double varHM;
            if(inst->isFwdStarting()){
                const FwdStFourierProductQuadVar& thisProd = *this;
                const FwdStFourierProcessQuadVar* thisProc = dynamic_cast<const FwdStFourierProcessQuadVar*>(&process);
                if(!thisProc) {
                    throw ModelException(method, "Process does not support StFourierProcessQuadVar or FwdStFourierProcessQuadVar interface" );
                }
                varHM = thisProc->expectation(thisProd, maturity);
            } else {
                const StFourierProductQuadVar& thisProd = *this;
                const StFourierProcessQuadVar* thisProc = dynamic_cast<const StFourierProcessQuadVar*>(&process);
                if(!thisProc) {
                    throw ModelException(method, "Process does not support StFourierProcessQuadVar or FwdStFourierProcessQuadVar interface" );
                }
                varHM = thisProc->expectation(thisProd, maturity);
            }
            //future realized variance or quadratic variation (NORMALIZED BY TIME)
            double futVar = varHM * years;
            //price = price +.5*(varHM * years - strike);
            price = price +.5*(futVar - strike); 
            //double varSwap = varHM * years;
            double varSwap = futVar; 
            varCapSqr = varCapSqr + 0.5*strike*strike - strike*varSwap;
#endif
            // See if naked option price is positive
            if(price <= 0.0) {
                if(price > -0.01) {
                    price = 0.0;
                } else {
                    throw ModelException("(Scaleless) Price of variance  cap is negative " +
                                         Format::toString(price) +
                                         ". Check integrator parameters e.g. tolerance, nbSteps.");
                }
            }

            double varCap = price;

            double histVar = volPast * volPast * pastWeight * totalYears;
            double strikeForSwap = inst->strikeVol * inst->strikeVol * totalYears - histVar;
            double modelSwap = varSwap - strikeForSwap;

            price /= totalYears;
            modelSwap /= totalYears;

            double pv = inst->instSettle->pv(inst->valueDate,
                                             maturity,
                                             inst->discount.get(),
                                             inst->asset.get());

            price *= pv*inst->notional*100.0;
            modelSwap *= pv*inst->notional*100.0;

            if (!inst->dontScaleByStrike) {
                price /= (2.0*inst->strikeVol);
                modelSwap /= (2.0*inst->strikeVol);
            }

            results->storePrice(price, discount->getCcy());

            //store intermediate results used to compute cv_coeff
            results->storeScalarGreek(varSwap,
                                    Results::DEBUG_PACKET,
                                    OutputNameSP(new OutputName(VAR_SWAP)));
            results->storeScalarGreek(varSwapSqr,
                                    Results::DEBUG_PACKET,
                                    OutputNameSP(new OutputName(VAR_SWAP_SQR)));
            results->storeScalarGreek(varCap,
                                    Results::DEBUG_PACKET,
                                    OutputNameSP(new OutputName(VAR_CAP)));
            results->storeScalarGreek(varCapSqr,
                                    Results::DEBUG_PACKET,
                                    OutputNameSP(new OutputName(VAR_CAP_SQR)));
            results->storeScalarGreek(modelSwap,
                                    Results::DEBUG_PACKET,
                                    OutputNameSP(new OutputName(MODEL_SWAP)));
            results->storeScalarGreek(strike,
                                    Results::DEBUG_PACKET,
                                    OutputNameSP(new OutputName(STRIKE)));

            // take care of additional outputs
            if (control && control->isPricing()) {
                OutputRequest* request = control->requestsOutput(OutputRequest::VOL_IN_PAST);
                if (request) {
                    results->storeRequestResult(request, volPast);
                }
                request = control->requestsOutput(OutputRequest::PAST_WEIGHT);
                if (request) {
                    results->storeRequestResult(request, pastWeight);
                }
                request = control->requestsOutput(OutputRequest::VOL_IN_FUTURE);
                if (request) {
                    double volFuture = 0.0;
                    if(inst->numTotalReturns - inst->numPastReturns) {
                        double adj = years*(double)inst->observationsPerYear/(double)(inst->numTotalReturns - inst->numPastReturns);
                        if (inst->subtractMeanVol) {
                            adj *= (double)inst->numTotalReturns/(double)(inst->numTotalReturns -1);
                        }
                        volFuture = sqrt(varHM*adj);
                    }
                    results->storeRequestResult(request, volFuture);
                }
                request = control->requestsOutput(OutputRequest::TOTAL_VOL);
                if (request) {
                    //contributionPast is pastVar normalized by ExpNPast / PPY = volPast * volPast * pastWeight
                    double contributionPast = volPast*volPast*pastWeight;

                    double contributionFuture = 0.0;
                    if (inst->subtractMeanVol) {
                        contributionFuture = futVar * (double)inst->observationsPerYear/(double)(inst->numTotalReturns - 1); 
                    } else {
                        contributionFuture = futVar * (double)inst->observationsPerYear/(double)(inst->numTotalReturns); 
                    }

                    double totalVolSquared = contributionPast + contributionFuture;

                    results->storeRequestResult(request, sqrt(totalVolSquared));
                }
                request = control->requestsOutput(OutputRequest::DISCOUNT_FACTOR);
                if (request) {
                    results->storeRequestResult(request, pv);
                }
                request = control->requestsOutput(OutputRequest::EXPECTED_N);
                if (request) {
                    results->storeRequestResult(request, inst->numTotalReturns);
                }
                // IMPLIED_N. The number of returns 'implied' by samples
                request = control->requestsOutput(OutputRequest::IMPLIED_N);
                if (request) {
                    results->storeRequestResult(request, inst->samples.size()-1);
                }

				// payment_dates
				request = control->requestsOutput(OutputRequest::PAYMENT_DATES);
				if (request) {
					DateTime maturity   = inst->samples.back().date;
					DateTime paymentDate = inst->instSettle->settles(maturity, inst->asset.get());
					DateTimeArray date(1, paymentDate);
					OutputRequestUtil::recordPaymentDates(control,results,&date);
				}

                // DELAY_PRICE
                InstrumentUtil::delayPriceHelper(control,
                                                 results,
                                                 price,
                                                 inst->valueDate,
                                                 inst->discount.get(),
                                                 inst->asset.get(),
                                                 inst->premiumSettle.get());

                // FWD_AT_MAT
                try{
                    InstrumentUtil::recordFwdAtMat(control,
                                                   results,
                                                   maturity,
                                                   inst->valueDate,
                                                   inst->asset.get());
                } catch(exception& e) {
                // Don't rethrow
                ModelException ee = ModelException::addTextToException(e, method + ": Computation of forward price failed");

                IObjectSP fwdObj(new Untweakable(ee));

                results->storeRequestResult(control->requestsOutput(OutputRequest::FWD_AT_MAT),
                                            fwdObj,
                                            OutputNameSP(new OutputName(inst->asset.get()->getName())));
                }
            }
        }

        catch (exception& e){
            throw ModelException(e, method);
        }
    }
};

const string VarOptCVFP::VAR_SWAP       = "varSwap";
const string VarOptCVFP::VAR_SWAP_SQR   = "varSwapSqr";
const string VarOptCVFP::VAR_CAP        = "varCap";
const string VarOptCVFP::VAR_CAP_SQR    = "varCapSqr";
const string VarOptCVFP::MODEL_SWAP     = "modelSwap";
const string VarOptCVFP::STRIKE         = "strike";


////////////////////////////////////////////////////////////////////////////
// FOURIER ENGINE IMPLEMENTATION FOR VARIANE SWAP (NOT CAP)

/** ClosedFormIntegrateLN product */
class VarSwapProduct_FE: virtual public FourierProduct {
public:
    /** Constructor */
    VarSwapProduct_FE(const VarianceSwap* inst);

    /** Price method: computes price and output requests.
        Can be called even for lastDate <= valueDate */
    virtual void price(const FourierEngine* model,
                       Control*             control,
                       Results*             results);
private:
    const VarianceSwap* inst;          //!< Instrument
};


////////////////////////////////////////////////////////////////////////////


/** Implementation of FourierEngine::IntoProduct interface */
FourierProduct* VarianceSwap::createProduct(const FourierEngine* model) const {
    static const string routine("VarianceSwap::createProduct");
    try {
        if(payoffType == "FORWARD") {
            return new VarSwapProduct_FE(this);
        }

        DateTime maturity = samples[samples.size()-1].date;

        if(useCV) {
            return new VarOptCVFP(this, maturity);
        }
	    else {
            return new VarOptFP(this, false, false, maturity);
        }
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

/** Implementation of VarCapCVModel::IntoProduct interface */
VarCapCVModel::Product* VarianceSwap::createProduct(const VarCapCVModel* model) const {
    static const string routine("VarianceSwap::createProduct");

    VarianceSwapSP varCap(copy(this));

    if (CString::equalsIgnoreCase(model->methodology, VarCapCVModel::CONTROL_VARIATE)) {
        varCap->useCV = true;
    }

    return new VarCapCVModel::Product(varCap, varCap->marketAsset.getSP());
}

////////////////////////////////////////////////////////////////////////////


VarSwapProduct_FE::VarSwapProduct_FE(const VarianceSwap* inst):
FourierProduct(inst->asset.get(), inst->valueDate, inst->discount.get(), inst->instSettle.get()),
inst(inst) {}


void VarSwapProduct_FE::price(const FourierEngine* model,
                              Control*             control,
                              Results*             results) {
    static const string routine = "VarSwapProduct_FE::price";
    try {
        // Keep references for easier access
        const DateTime& valueDate = inst->valueDate;

        const CashFlowArray& samples = inst->samples;
        DateTimeArray obsDates;
        DoubleArray   obsSamples;
        for(int iSample = 0; iSample < samples.size(); iSample++) {
            obsDates.push_back(samples[iSample].date);
            obsSamples.push_back(samples[iSample].amount);
        }

        const DateTime& firstDate = obsDates.front();
        const DateTime& lastDate = obsDates.back();
        double strike = inst->strikeVol;
        int expectedN = inst->numTotalReturns;
        int observationsPerYear = inst->observationsPerYear;
        CAssetConstSP asset = inst->asset.getSP();

        // PART 1: past contribution including partly past / partly future contribution
        double pastFloatingVar = 0.0;
        int iStep = 1;
        if(firstDate <= valueDate) {
            // Compute log-returns for past samples
            DoubleArraySP logReturns = VarSwapUtilities::computeHistoricLogReturns(
                asset.get(), obsDates, obsSamples, valueDate, true, !inst->noDivAdj, false);

            // Compute past realized variance
            for(iStep = 0; iStep < obsDates.size() - 1; iStep++) {
                const DateTime& thisStartDate = obsDates[iStep];
                //const DateTime& thisEndDate   = obsDates[iStep + 1];
                if(thisStartDate < valueDate) {
                    // Started in the past
                    double thisContribution = Maths::square((*logReturns)[iStep + 1]);
                    pastFloatingVar += thisContribution;
                } else {
                    // Strictly future
                    break;
                }
            }
        }

        // PART 2: future contribution. Note may be fwd starting
        double futureFloatingVar = 0.0;
        if(valueDate < lastDate) {
            const FourierProcess& process = model->getProcess();
            const StFourierProcessQuadVar* quadVar = dynamic_cast<const StFourierProcessQuadVar*>(&process);
            if(!quadVar) {
                throw ModelException(routine, "Process must implement StFourierProcessQuadVar interface to price Variance Swaps");
            }

            DateTimeArray maturities(1, lastDate);
            if(valueDate < firstDate) {
                // Forward starting case: put the timeline in reverse order
                maturities.push_back(firstDate);
            }

            DoubleArray futureFloatingVarComponents(2);
            for(int iMat = 0; iMat < maturities.size(); iMat++) {
                StFourierProductQuadVar product;
                double annualizedVar = quadVar->expectation(product, maturities[iMat]);
                double time = process.getTimeMetric().yearFrac(valueDate, maturities[iMat]);
                futureFloatingVarComponents[iMat] = annualizedVar * time;

                // Do dividend adjustments
                // VarSwap in stochastic vol is continuous volatility only
                if(inst->noDivAdj) {
                    // No dividend adjustment = use full return = add the effect of divs to vanillas
                    double divVar = VarSwapUtilities::futureDivsVariance(
                        asset.get(), valueDate, valueDate, maturities[iMat]);
                    futureFloatingVarComponents[iMat] += divVar;
                }
            }

            futureFloatingVar = futureFloatingVarComponents[0] -
                futureFloatingVarComponents[1];
        }

        // Compute value of legs
        double floatingLeg =
            ((double)observationsPerYear / (double)expectedN) *
            (pastFloatingVar + futureFloatingVar);

        double fixedLeg = Maths::square(strike);

        // Compute price for non-settled instruments
        double fwdValuedPrice = 100.0 * inst->notional * (floatingLeg - fixedLeg);
        if(!inst->dontScaleByStrike) {
            fwdValuedPrice /=  2.0 * strike;
        }

        double price = 0.0;
        double pv = 1.0;
        DateTime settlement  = instSettle->settles(lastDate, asset.get());
        if((valueDate <= settlement)) {
            // PV and scale
            pv = instSettle->pv(valueDate, lastDate, discount, asset.get());
            price = pv * fwdValuedPrice;
        }

        // Store price
        results->storePrice(price, discount->getCcy());


        // Output Requests now
        if (control && control->isPricing()) {
            // Compute past and future ExpectedN
            VolRequestTime volReq;
            IVolProcessedSP vol(asset->getProcessedVol(&volReq));
            HolidayWrapper assetHols = vol->GetTimeMetric()->getHolidays();
            int pastExpectedN;
            int futureExpectedN;
            VarSwapUtilities::pastAndfutureExpN(obsDates, assetHols, valueDate, expectedN, pastExpectedN, futureExpectedN);

            OutputRequest* request;

            // Compute all vol requests first
            // TOTAL_VOL
            double totalVar = pastFloatingVar + futureFloatingVar;
            double totalvol = sqrt((double)observationsPerYear / expectedN * totalVar);
            request = control->requestsOutput(OutputRequest::TOTAL_VOL);
            if(request) {
                    results->storeRequestResult(request, totalvol);
            }

            // Compute future vol
            double volFuture = 0.0;
            if(futureExpectedN) {
                // There exists future
                volFuture = sqrt((double)(observationsPerYear) / futureExpectedN * futureFloatingVar);
            }

            // VOL_IN_FUTURE
            request = control->requestsOutput(OutputRequest::VOL_IN_FUTURE);
            if(request) {
                results->storeRequestResult(request, volFuture);
            }

            // IND_VOL
            request = control->requestsOutput(OutputRequest::IND_VOL);
            if (request) {
                results->storeRequestResult(request, volFuture);
            }

            // VOL_IN_PAST
            request = control->requestsOutput(OutputRequest::VOL_IN_PAST);
            if (request) {
                double volPast = 0.0;
                if(pastExpectedN) {
                    volPast = sqrt((double)(observationsPerYear) / pastExpectedN * pastFloatingVar);
                }
                results->storeRequestResult(request, volPast);
            }

            // PAST_WEIGHT
            request = control->requestsOutput(OutputRequest::PAST_WEIGHT);
            if (request) {
                double pastWeight = (double)pastExpectedN / (double)expectedN;
                results->storeRequestResult(request, pastWeight);
            }

            //STRIKE_VOL
            request = control->requestsOutput(OutputRequest::STRIKE_VOL);
            if (request) {
               results->storeRequestResult(request, strike);
            }

            //EXPECTED_N
            request = control->requestsOutput(OutputRequest::EXPECTED_N);
            if (request) {
                results->storeRequestResult(request, expectedN);
            }
            // IMPLIED_N. The number of returns 'implied' by samples
            request = control->requestsOutput(OutputRequest::IMPLIED_N);
            if (request) {
                results->storeRequestResult(request, samples.size()-1);
            }

            //DISCOUNT_FACTOR
            request = control->requestsOutput(OutputRequest::DISCOUNT_FACTOR);
            if (request && valueDate <= settlement) {
                results->storeRequestResult(request, pv);
            }

            //PAYMENT_DATES
            request = control->requestsOutput(OutputRequest::PAYMENT_DATES);
            if (request) {
                DateTimeArray date(1, settlement);
                OutputRequestUtil::recordPaymentDates(control, results, &date);
            }

            //KNOWN_CASH_FLOWS
            request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
            if (request) {
               if (lastDate <= valueDate) {
                    CashFlow cf(settlement, fwdValuedPrice);
                   CashFlowArray cfl(1, cf);
                   OutputRequestUtil::recordKnownCashflows(control,
                       results,
                       discount->getCcy(),
                       &cfl);
                }
            }

            //DELAY_PRICE
            InstrumentUtil::delayPriceHelper(control,
                results,
                price,
                valueDate,
                discount,
                asset.get(),
                inst->premiumSettle.get());

            //FWD_AT_MAT
            InstrumentUtil::recordFwdAtMat(control,
                results,
                lastDate,
                valueDate,
                asset.get());
        }
    } catch(exception& e) {
      throw ModelException(e, routine);
    }
}


////////////////////////////////////////////////////////////////////////////


// State var compliant Monte Carlo product
class VarianceSwapMC: public MCProductClient {
private:
    const VarianceSwap* inst;
    double              strikeVol;
    double              notional;
    double              scaleFactor;
    bool                dontScaleByStrike;
    int                 payoffType;
    double              cap;
    double              volPast;

    // state vars
    MCSqrtAnnualQuadVarSP             sqrtAnnualQuadVarGen;      // generator for sqrtAnnualQuadVar
    MCSqrtAnnualQuadVar::IStateVarSP  sqrtAnnualQuadVarSV;       // sqrtAnnualQuadVar state variable
    SVGenDiscFactorSP                    dfGen;                     // generator for discount factors
    SVDiscFactorSP         dfSV;                      // df state variable

    /** Override default method on IMCProduct. This method is called every time
        the path generator is changed (which is, at the moment, when the
        past path generator is created, and then when the future path
        generator is created  */
    virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen){
        static const string routine = "VarianceSwapMC::pathGenUpdated";

        try {
            sqrtAnnualQuadVarSV = sqrtAnnualQuadVarGen->getSqrtAnnualQuadVarSV(newPathGen);
            dfSV = dfGen->getSVDiscFactor(newPathGen);
        } catch (exception& e) {
            throw ModelException(e, routine);
        }
    };

    class VSwapPrices;
    typedef refCountPtr<VSwapPrices> VSwapPricesSP;
    class VSwapPrices: public IMCPricesSimple {
    public:
        enum PricesIndices{
            SWAP_PRICE = 0,
            FUTURE_VOL,
            NB_PRICES
        };

        /** adds supplied price to this set of IMCPrices */
        void add(double price,
                 int    index){
            simplePrices[index]->add(price);
        }

        /** Returns the averaged result, together with the standard error */
        virtual void getResult(double* result,
                               double* resultStdErr) const{
            simplePrices[SWAP_PRICE]->getResult(result, resultStdErr);
        }

        double getTotalVol() const{
            double annualVariance, dummy;
            simplePrices[FUTURE_VOL]->getResult(&annualVariance, &dummy);
            return sqrt(annualVariance);
        }

        /** Returns true if the path, identified by pathIdx, should be
            repriced. If false is returned then there will be no add()
            method called for this path and the IMCPrices object must
            take any appropriate action */
        virtual bool repriceForGreek(int pathIdx){
            return true;
        }

        /** Support product-level caching across tweaks */
        virtual int storagePerPath(IMCProduct* product) const{
            return 0;
        }

        virtual void configureCache(const IntArray& changedAssets){}

        /** Returns a deep copy of this object */
        IMCPrices* clone() const{
            VSwapPricesSP copy(new VSwapPrices());
            copy->simplePrices.resize(simplePrices.size());
            for (unsigned int iPrice = 0; iPrice < copy->simplePrices.size(); ++iPrice){
                IMCPricesSP temp(this->simplePrices[iPrice]->clone());
                copy->simplePrices[iPrice] = DYNAMIC_POINTER_CAST<IMCPricesSimple>(temp);
            }
            return copy.get();
        }

        VSwapPrices(int NbIter,
                    int NbSubSamples):
        simplePrices(NB_PRICES){
            for (unsigned int iPrice = 0; iPrice < simplePrices.size(); ++iPrice){
                simplePrices[iPrice] = IMCPricesSimpleSP(new MCPricesSimple(NbIter, NbSubSamples));
            }
        }

    private:
        VSwapPrices(){}

        /** adds supplied price to this set of IMCPrices */
        virtual void add(double price){
            throw ModelException("IMCPrices::add(double)",
                                 "internal error");
        }

        /** Returns the last price stored. Undefined behaviour if no
            prices yet stored */
        virtual double lastPrice() const{
            throw ModelException("IMCPrices::lastPrice",
                                 "internal error");
        }

        /** On pricing run returns MAX(x, 0.0). It should be used for applying
            the 'final point of optionality' within the product. This allows
            QuickGreeks type of IMCPrices not to apply the max when doing first
            order greeks (this may sound strange but see the docs for why) */
        virtual double maxWithZero(double x) const{
            throw ModelException("IMCPrices::maxWithZero",
                                 "internal error");
        }

        /** Reset this object so that it can be used for the same operation
            again. Normally, a new IMCPrices object is created for each
            pricing run. However, for quick x gamma, it is important to use
            the same one for each pair of assets */
        virtual void reset(){
            throw ModelException("IMCPrices::reset",
                                 "internal error");
        }

        /** Ease cloning */
        virtual IMCPrices* emptyConstructor() const{
            throw ModelException("IMCPrices::emptyConstructor",
                                 "internal error");
        }

        vector<IMCPricesSimpleSP> simplePrices;
    };

public:
    virtual IMCPrices* createOrigPrices(int  nbIter,
                                     int  nbSubSamples,
                                     int  mode) {
        return new VSwapPrices(nbIter, nbSubSamples);
    }

    /** invoked after final simulated path is run. Default does nothing.
        Allows derived classes to store debug information for example */
    virtual void recordExtraOutput(Control*      control,
                                   Results*      results,
                                   const IMCPrices& prices) const{
        const VSwapPrices& myprices = static_cast<const VSwapPrices&>(prices);
        if (control && control->isPricing()) {
            OutputRequest* request =
                control->requestsOutput(OutputRequest::TOTAL_VOL);
            if (request) {
                double volFuture = myprices.getTotalVol();
                results->storeRequestResult(request, volFuture);
            }
            request = control->requestsOutput(OutputRequest::VOL_IN_PAST);
            if (request) {
                results->storeRequestResult(request, volPast);
            }
        }
    }

    // types of spot discretization schemes
    struct PayoffType{
        enum {
            FORWARD = 0,
            FORWARD_CAP,
            NB_ENUMS
        };
    };
    typedef Enum2StringListHelper<PayoffType> PayoffTypeHelper;

    /** Appends 'true' (ie non derived) state variable generators
        required to the supplied collector.*/
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const{
        svCollector->append(sqrtAnnualQuadVarGen.get());
        svCollector->append(dfGen.get());
    }

    // equivalent to InstIntoMCProduct
    VarianceSwapMC(const VarianceSwap*       inst,
                   const SimSeriesConstSP&   simSeries,             // simulation dates
                   const IPastValuesConstSP& mcSqrtAnnualQuadVarPastValues);  // past values

    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&                prices);
};

template<> string nameForType<VarianceSwapMC::PayoffTypeHelper>(VarianceSwapMC::PayoffTypeHelper*){
    return "VarianceSwapMC::PayoffType";
}

template<> string VarianceSwapMC::PayoffTypeHelper::names[VarianceSwapMC::PayoffTypeHelper::EnumList::NB_ENUMS] = {
    "FORWARD",
    "FORWARD_CAP"
};

// equivalent to InstIntoMCProduct
VarianceSwapMC::VarianceSwapMC(const VarianceSwap*       inst,
                               const SimSeriesConstSP&   simSeries,             // simulation dates
                               const IPastValuesConstSP& mcSqrtAnnualQuadVarPastValues):  // past values
MCProductClient(inst->asset.get(),
                inst->valueDate,
                inst->discount.get(),
                IRefLevelSP(IRefLevel::Util::makeZero(inst->valueDate)),
                simSeries,
                IPastValuesSP(IPastValues::Util::makeTrivial(inst->valueDate, 0.0)),
                mcSqrtAnnualQuadVarPastValues,
                inst->instSettle.get(),
                simSeries->getLastDate()),
inst(inst),
strikeVol(inst->strikeVol),
notional(inst->notional),
scaleFactor(100.0),
dontScaleByStrike(inst->dontScaleByStrike),
volPast(0.0),
sqrtAnnualQuadVarGen(new MCSqrtAnnualQuadVar(simSeries)),
dfGen(new SVGenDiscFactor(inst->valueDate,inst->discount.getSP(),
                       inst->instSettle, simSeries->getLastDate())) {
    static const string method("VarianceSwapMC::VarianceSwapMC");
    try{
        if (!dontScaleByStrike) {
            if (Maths::isZero(strikeVol)) {
                throw ModelException(method,
                                     "cannot scale by zero strike! ('dontScaleByStrike' is set to false)");
            }
        }
        payoffType = PayoffTypeHelper::getIndex(inst->payoffType);
        if (payoffType == PayoffTypeHelper::EnumList::FORWARD_CAP){
            cap = inst->cap;
        }
        else if (payoffType == PayoffTypeHelper::EnumList::FORWARD){
            cap = 1.0;
        }
        else{
            throw ModelException(method, "internal error");
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

void VarianceSwapMC::payoff(const IPathGenerator*  pathGen,
                            IMCPrices&                prices) {
    static const string method("VarianceSwapMC::payoff");
    try{
        VSwapPrices& myprices = static_cast<VSwapPrices&>(prices);
        int iAsset = 0;     // only 1 asset
        const SVPath& path = sqrtAnnualQuadVarSV->path(iAsset);
        if (doingPast()){
            volPast = path.begin() == path.end() ? 0.0 : path[path.end() - 1];  // corresponding to today's date
        }
        int iStep = path.end() - 1;     // corresponding to maturity date
        double annualVar = Maths::square(path[iStep]);
        myprices.add(annualVar, VSwapPrices::FUTURE_VOL);
        double payout = notional * scaleFactor * (annualVar - cap * strikeVol * strikeVol);
        if (!dontScaleByStrike) {
            payout /= (2.0 * strikeVol);
        }
        if (payoffType == PayoffTypeHelper::EnumList::FORWARD_CAP){
            payout = Maths::max(payout, 0.0);
        }
        else if (payoffType == PayoffTypeHelper::EnumList::FORWARD){
            // nothing to do
        }
        else{
            throw ModelException(method, "internal error");
        }
#if 1
        payout *= dfSV->firstDF();
#else
        // line below to be reviewed. The issue being that
        // (currently) we can't access values in the future when
        // doing the past. NB apply df before vanillaReprice->store()
        if (!doingPast()){
            payout *= dfSV->path()[0];
        }
#endif
        myprices.add(payout, VSwapPrices::SWAP_PRICE);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* VarianceSwap::createProduct(const MonteCarlo* model) const {
    static const string method("VarianceSwap::createProduct(const MonteCarlo*)");
    // simple simSeries
    // contains first and last sampling dates
    SimSeriesSP simSeries(new SimSeries(1));
    // create sim time line that contains first and last sample dates
    // together with today's date if the latter is not redundant
    DateTimeArray simDates;
    if (valueDate > samples.front().date
        && valueDate < samples.back().date){
        simDates.resize(3);
        simDates[0] = samples.front().date;
        simDates[1] = valueDate;
        simDates[2] = samples.back().date;
    }
    else{
        simDates.resize(2);
        simDates[0] = samples.front().date;
        simDates[1] = samples.back().date;
    }
    simSeries->addDates(simDates);
    // for every sample date + today calculate realized vol
    int nbSamples = samples.size();
    DateTimeArray allDates;
    allDates.reserve(nbSamples + 1);
    DoubleArray allValues;
    allValues.reserve(nbSamples + 1);
    // insert dates smaller than today's
    int iSample = 0;
    for (; iSample < nbSamples && samples[iSample].date < valueDate; ++iSample){
        allDates.push_back(samples[iSample].date);
        int notUsed;
        allValues.push_back(historicalVol(notUsed,
                                          allDates.back(),
                                          samples[iSample].amount));
    }
    // insert today's date if non-redundant
    // and if not before first sampling date
    if (iSample >= nbSamples
        || (iSample !=0 && samples[iSample].date > valueDate)){
        allDates.push_back(valueDate);
        int notUsed;
        allValues.push_back(historicalVol(notUsed));
    }
    // insert dates greater than today's
    for (; iSample < nbSamples; ++iSample){
        allDates.push_back(samples[iSample].date);
        int notUsed;
        double spotPrice = samples[iSample].date == valueDate ?
                               asset->getSpot() :
                               samples[iSample].amount;
        allValues.push_back(historicalVol(notUsed,
                                          allDates.back(),
                                          spotPrice));
    }
    DoubleMatrix allValuesAsMatrix(allValues);
    IPastValuesSP pastValues(
        IPastValues::Util::makeSimple(allDates,
                                      allValuesAsMatrix,
                                      IPastValues::NonNegativePastValueCheck()));
    // if state vars requested
    if (!model->stateVarUsed()){
        throw ModelException(method,
                             "useStateVar needs to be set to true");
    }
    // otherwise, use old methodology
    return new VarianceSwapMC(this, simSeries, pastValues);
}

class  VarianceSwapHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        clazz->setDescription("Variance Swap");
        REGISTER(VarianceSwap, clazz);
        SUPERCLASS(VolVarShell);
        IMPLEMENTS(ISupportVegaMatrixLite);
        IMPLEMENTS(NumericalIntegrationLN::IIntoProduct);
        IMPLEMENTS(ImpliedIntegration::IIntoProduct);
        IMPLEMENTS(ClosedFormIntegrateLN::IIntoProduct);
        IMPLEMENTS(FourierEngine::IIntoProduct);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultVarianceSwap);
        FIELD(useCV,"whether to do CV adjustment");
        FIELD_MAKE_TRANSIENT(useCV);
        FIELD(marketAsset,"asset to store market info");
        FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(marketAsset);
    }

    static IObject* defaultVarianceSwap(){
        return new VarianceSwap();
    }
};

CClassConstSP const VarianceSwap::TYPE = CClass::registerClassLoadMethod(
    "VarianceSwap", typeid(VarianceSwap), VarianceSwapHelper::load);


class VarianceSwapNumerical: virtual public NumericalIntegrationLN::IProduct,
                             virtual public ImpliedIntegration::IProduct,
                             virtual public ClosedFormIntegrateLN::IProduct {

private:
    const VarianceSwap*  vswap; // a reference

    DateTime maturity;
    DateTime settles;
    double   fwd;

public:
    VarianceSwapNumerical(const VarianceSwap* vswap): vswap(vswap) {
        maturity = vswap->samples[vswap->samples.size()-1].date;
        settles  = vswap->instSettle->settles(maturity, vswap->asset.get());
        fwd = vswap->asset->fwdValue(maturity);
    }

    void price(ClosedFormIntegrateLN*   model,
               Control*                 control,
               CResults*                results) {
        static const string method = "ClosedFormIntegrateLN::price";

        // Obtain processed basis from vol
        VarSwapBasis::VarSwapBasisProcSP basisProc;
        try {
            DateTimeArray dates;
            const DateTime& firstSample = vswap->samples.front().date;
            if(firstSample > vswap->valueDate) {
                dates.push_back(firstSample);
            }
            dates.push_back(vswap->samples.back().date);

            VarSwapBasis::VarSwapBasisRequestSP basisRequest(new VarSwapBasis::VarSwapBasisRequest(dates, vswap->discount.getSP()));
            IVolProcessedSP procVol;
            try {
                IVolProcessedSP tmp(vswap->asset->getProcessedVol(basisRequest.get()));
                procVol = tmp;
            } catch(exception&) {
                // Ignore all errors and proceed without basis
                // Could be that VolSpline threw an error like: don't know what VarSwapBasisRequest is
                basisProc = VarSwapBasis::VarSwapBasisProcSP(   );
            }

            if(procVol.get()) {
                VarSwapBasis::VarSwapBasisProcError* error = dynamic_cast<VarSwapBasis::VarSwapBasisProcError*>(procVol.get());
                if(error) {
                    throw ModelException(const_cast<ModelException&>(error->getException()), method);
                }
                basisProc = VarSwapBasis::VarSwapBasisProcSP::dynamicCast(procVol);
            }

            // Get option recorder
            VanillaContractsRecorderSP recorder = VanillaContractsRecorder::createVanillaOptionRecorder(control);

            // check whether future trading time is positive
            ATMVolRequest volReq;
            CVolProcessedBSSP vol(vswap->asset->getProcessedVol(&volReq));
            double t = vol->calcTradingTime(vswap->isFwdStarting() ? vswap->samples.front().date
                : vswap->valueDate, vswap->samples.back().date);
            if (Maths::isZero(t)) {
                vswap->price(0.0, model, control, results);
                return;
            }

            /** now compute varFuture
                more involved than in price(ImpliedIntegration), since ClosedFormIntegrateLN model is specified
                hence, we dont just call utility function futureVar or futureFwdVar */
            double lowStrike = 0.0;
            double highStrike = 0.0;
            int nbSteps = 0;

            // compute var from valueDate to swap->samples.back().date
            double tenor = vol->calcTradingTime(vswap->valueDate, vswap->samples.back().date);
            //DateTime maturity = vswap->samples.back().date);
            DateTime maturity = vswap->samples.back().date;
            double fwd = vswap->asset->fwdValue(maturity);
            model->limits(vswap->asset.get(),vswap->valueDate,maturity,lowStrike,highStrike,nbSteps);

            bool isDoingVegaMatrix = VarianceSwapNumerical::isDoingVegaMatrix(control);

            // varFuture is variance and thus sigma^2 * t
            recorder->setScalingFactor(1.0);
            double varFuture = vswap->varFromCalls(lowStrike,
                                                   highStrike,
                                                   nbSteps,
                                                   maturity,
                                                   fwd,
                                                   model->integrationMethodGet(),
                                                   model->negativeFwdVarAllowed(),
                                                   model->calcAbsPrecision(tenor),
                                                   model->calcRelPrecision(tenor),
                                                   isDoingVegaMatrix,
                                                   basisProc,
                                                   recorder);
            if (vswap->isFwdStarting()) {
                // compute var from valueDate to  vswap->samples.front().date
                tenor = vol->calcTradingTime(vswap->valueDate, vswap->samples.front().date);
                maturity = vswap->samples.front().date;
                fwd = vswap->asset->fwdValue(maturity);
                model->limits(vswap->asset.get(),vswap->valueDate,maturity,lowStrike,highStrike,nbSteps);

                recorder->setScalingFactor(-1.0);
                double varFutureShort = vswap->varFromCalls(lowStrike,
                                                            highStrike,
                                                            nbSteps,
                                                            maturity,
                                                            fwd,
                                                            model->integrationMethodGet(),
                                                            model->negativeFwdVarAllowed(),
                                                            model->calcAbsPrecision(tenor),
                                                            model->calcRelPrecision(tenor),
                                                            isDoingVegaMatrix,
                                                            basisProc,
                                                            recorder);
                // varFuture = variance(T+t) - variance(T)
                varFuture -= varFutureShort;
              }

            if (Maths::isNegative(varFuture)) {
                throw ModelException(method, "future variance (" +
                    Format::toString(varFuture) + ") is negative");
            }
            vswap->price(varFuture, model, control, results, basisProc, recorder);
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    void price(NumericalIntegrationLN* model,
        Control*                control,
        CResults*               results) {
        static const string method = "VarianceSwapNumerical::price";
        try {
            // check whether future trading time is positive
            ATMVolRequest volReq;
            CVolProcessedBSSP vol(vswap->asset->getProcessedVol(&volReq));
            double t = vol->calcTradingTime(vswap->isFwdStarting() ? vswap->samples.front().date
                : vswap->valueDate, vswap->samples.back().date);
            if (Maths::isZero(t)) {
                vswap->price(0.0, model, control, results);
                return;
            }

            // varFuture is variance and thus sigma^2 * t
            double varFuture = 0.0;

            if (vswap->isFwdStarting()) {
                throw ModelException(method, "Cannot price fwd starting var swaps with Numerical Integration model.  Please "
                    "use Implied Integration model.");
            } else {
                varFuture = model->integrate(this);
            }

            if (Maths::isNegative(varFuture)) {
                throw ModelException(method, "future variance (" +
                    Format::toString(varFuture) + ") is negative");
            }

            vswap->price(varFuture, model, control, results);
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    PDFCalculator* pdfCalculator() const {
        static const string method = "VarianceSwapNumerical::pdfCalculator";
        try {
            LinearStrikeTSVolRequest volRequest(0.0,
                                                vswap->valueDate,
                                                maturity,
                                                false);  // note: we do not use the forward starting vol methodology
                                                         // as as the log contract is a European payoff

            PDFRequestLNStrike pdfRequest(&volRequest);

            return vswap->asset->pdfCalculator(&pdfRequest);
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    double payoff(double spot) const {
        return (2.0*log(fwd/spot));
    }

    double centre(const DateTime& date) const {
        return vswap->asset->fwdValue(date);
    }

    double variance(double strike, const DateTime& date) const {
        static const string method = "VarianceSwapNumerical::variance";
        try {
            LinearStrikeTSVolRequest volRequest(strike,
                                                vswap->valueDate,
                                                maturity,
                                                false);  // note: we do not use the forward starting vol methodology
                                                         // as as the log contract is a European payoff

            CVolProcessedBSSP volBS(vswap->asset->getProcessedVol(&volRequest));

            return volBS->CalcVar(vswap->valueDate, date);
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    DateTime time() const {
        return maturity;
    }

    // needed for ImpliedIntegration price below.
    static bool isDoingVegaMatrix(Control* control) {
        //maybe doing old vega matrix....
        SensitivitySP sens = control->getCurrentSensitivity();
        Sensitivity* sensPtr = sens.get();
        bool doingVM = VegaMatrix::TYPE->isInstance(sensPtr) || 
                       VegaProxyMatrix::TYPE->isInstance(sensPtr) ||
                       VegaMatrixLite::TYPE->isInstance(sensPtr);
        return doingVM;
    }

    // ImpliedIntegration methods
    void price(ImpliedIntegration* model,
        Control*            control,
        CResults*           results) {
        static const string method = "ImpliedIntegration::price";
        try {
            // check whether future trading time is positive
            ATMVolRequest volReq;
            CVolProcessedBSSP vol(vswap->asset->getProcessedVol(&volReq));
            double t = vol->calcTradingTime(vswap->isFwdStarting() ? vswap->samples.front().date
                : vswap->valueDate, vswap->samples.back().date);
            if (Maths::isZero(t)) {
                vswap->price(0.0, model, control, results);
                return;
            }

            // varFuture is variance and thus sigma^2 * t
            double varFuture = 0.0;

            /** depending on whether isDoingVegaMatrix is false or true, we use either use
                ImpliedIntegration (false) or ClosedFormIntegrateLN (true) as model; in case
                we use the latter, then the input parameters are hard coded ... */
            if(vswap->isFwdStarting()) {
                varFuture = VarianceSwapUtil::futureFwdVar(vswap->asset.get(),
                                                           vswap->discount.get(),
                                                           vswap->valueDate,
                                                           vswap->samples.front().date,
                                                           vswap->samples.back().date,
                                                           model,
                                                           control);
            } else {
                varFuture = VarianceSwapUtil::futureVar(vswap->asset.get(),
                                                        vswap->discount.get(),
                                                        vswap->valueDate,
                                                        vswap->samples.back().date,
                                                        model,
                                                        control);
            }

            if (Maths::isNegative(varFuture)) {
                throw ModelException(method, "future variance (" +
                    Format::toString(varFuture) + ") is negative");
            }

            vswap->price(varFuture, model, control, results);
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    class Payoff: public Function1DDouble {
    public:
        virtual double operator()(double spot) const {
            return (2.0*log(fwd/spot));
        };

        Payoff(double fwd): fwd(fwd) {}
    private:
        double fwd;
    };

    class Integral: public Function1DDouble {
    public:
        virtual double operator()(double spot) const {
            return (2.0*spot*(log(fwd) - (log(spot) - 1.0)));
        };

        Integral(double fwd): fwd(fwd) {}
    private:
        double fwd;
    };

    /** given a spot level, return the payoff */
    Function1DDouble* payoff() const {
        return new Payoff(fwd);
    }

    /** given a spot level, return integral of the payoff */
    Function1DDouble* indefiniteIntegral() const {
        return new Integral(fwd);
    }

    DateTime today() const {
        return vswap->valueDate;
    }

    LinearImpliedSampler* sampler(
        const DateTime&  date,
        const PDFParams& params) const {
        static const string method = "VarianceSwapNumerical::sampler";
        try {
            LinearStrikeTSVolRequestSP volRequest(
                new LinearStrikeTSVolRequest(0.0,
                                             vswap->valueDate,
                                             maturity,
                                             false)); // note: we do not use the forward starting vol methodology
                                                      // as as the log contract is a European payoff

            PDFRequestLNStrikeSP pdfRequest(new PDFRequestLNStrike(volRequest.get()));

            DateTimeArray datesTo(1, maturity);

            CAssetConstSP asset(vswap->asset.get());

            return new LinearImpliedSampler(asset,
                                            volRequest,
                                            vswap->valueDate,
                                            vswap->valueDate,
                                            datesTo,
                                            params,
                                            pdfRequest);
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }
};

/** Implementation of NumericalIntegrationLN::IntoProduct interface */
NumericalIntegrationLN::IProduct* VarianceSwap::createProduct(NumericalIntegrationLN* model) const {
    return new VarianceSwapNumerical(this);
}

/** Implementation of ImpliedIntegration::IntoProduct interface */
ImpliedIntegration::IProduct* VarianceSwap::createProduct(ImpliedIntegration* model) const {
    return new VarianceSwapNumerical(this);
}

const string ClosedFormIntegrateLN::DISCRETE_INTEGRATION = string("DiscreteIntegration");

const string ClosedFormIntegrateLN::DEFAULT_INTEGRATION_METHOD = string("Default");

const string ClosedFormIntegrateLN::DEFAULT_VEGAMATRIX_STRIKES = string("Default");

const string ClosedFormIntegrateLN::VOLSURFACE_STRIKES = string("VolSurfaceStrikes");
const string ClosedFormIntegrateLN::EQUIDISTANT_STRIKES = string("EquidistantStrikes");
const string ClosedFormIntegrateLN::EQUIDISTANT_VEGA = string("EquidistantVega");

const double ClosedFormIntegrateLN::MIN_TOLERANCE = 1.0e-14;
const double ClosedFormIntegrateLN::MAX_TOLERANCE = 0.0001; // 1 bp in terms of vol

const double ClosedFormIntegrateLN::DEFAULT_ABSOLUTE_VOL_PRECISION = 1.0e-5;
const double ClosedFormIntegrateLN::DEFAULT_RELATIVE_VOL_PRECISION = 1.0e-5;

const string ClosedFormIntegrateLN::TERM_STRUCTURE_DEFAULT = string("Default");
const string ClosedFormIntegrateLN::TERM_STRUCTURE_NONE = string("None");
const string ClosedFormIntegrateLN::TERM_STRUCTURE_ALL = string("All");
const string ClosedFormIntegrateLN::TERM_STRUCTURE_DEFAULT_TENOR = string("1M");

//divMethodology used to compute prices with dividends
const string ClosedFormIntegrateLN::DIV_DEFAULT = string("Default");
const string ClosedFormIntegrateLN::DIV_CONTINUOUS_APPROX = string("ContinuousApprox");
const string ClosedFormIntegrateLN::DIV_BRUTE_FORCE = string("BruteForce");

string ClosedFormIntegrateLN::getDefaultIntegrationMethod() {
    // Returns default method
    return IntFuncInf::TYPE->getName();
}
string ClosedFormIntegrateLN::getDefaultStrikesForVegaMatrix(){
    // Returns default strikes
    return VOLSURFACE_STRIKES;
}
bool ClosedFormIntegrateLN::getUseBasis() const {
    return useBasis;
}

//returns divMethodology used to compute prices with dividends
string ClosedFormIntegrateLN::getDivMethodology() const{    
    return divMethodology;
}

/** Returns a series of dates and weights for the term structure of portfolios */
//using log-forward coefficients (Manos initial implementation)
void ClosedFormIntegrateLN::getPriceTSPortfolios(const CAsset*         asset,
                                                 const DateTime&       valueDate,
                                                 const DateTimeArray&  obsDates,
                                                 DateTimeArray&        datesTS,
                                                 DoubleArray&          weightsTS) const {
    static const string routine = "ClosedFormIntegrateLN::getPriceTSPortfolios";
    try {
        datesTS.resize(0);
        weightsTS.resize(0);

        if(CString::equalsIgnoreCase(termStructMethod, TERM_STRUCTURE_NONE) ||
           obsDates.back() <= valueDate) {
            return;
        }

        DateTimeArray dates;
        if (CString::equalsIgnoreCase(termStructMethod, TERM_STRUCTURE_ALL)) {
            // Compute strictly future dates and fwds at those dates
            dates = valueDate.getFutureDates(obsDates);
            dates.insert(dates.begin(), valueDate);
            if(!dates.size()) {
                return;
            }
        } else {
            // Work with a maturity period
            string tenor = termStructMethod;
            const char *inp = tenor.c_str();
            if(isdigit(*inp)) {
                tenor = "-" + termStructMethod;
            }

            // Start from end and add dates per tenor
            MaturityPeriod frequency(tenor);
            DateTime thisDate = obsDates.back();
            while(thisDate > valueDate) {
                if(thisDate > valueDate) {
                    dates.insert(dates.begin(), thisDate);
                }
                thisDate = frequency.toDate(thisDate);
            }
            dates.insert(dates.begin(), valueDate);
        }

        // First strictly future date and fwd will be dropped eventually
        DoubleArray weights(dates.size());
        DoubleArray fwds(dates.size());
        asset->fwdValue(dates, fwds);

        // Fill in the weights
        for(int iFutStep = 1; iFutStep < dates.size(); iFutStep++) {
            double thisWeight = log(fwds[iFutStep] / fwds[iFutStep - 1]);
            weights[iFutStep] = thisWeight;
        }

        // Drop first entries and return
        dates.erase(dates.begin());
        weights.erase(weights.begin());

        datesTS = dates;
        weightsTS = weights;

    } catch (exception& e) {
        throw ModelException(e, routine);
    }
}

/** Returns a series of dates and weights for the term structure of portfolios */
//using log-pv coeffieints (Gad implementation to account for discrete dividends
void ClosedFormIntegrateLN::getPriceTSPortfolios(const YieldCurve*     discount,
                                                 const DateTime&       valueDate,
                                                 const DateTimeArray&  obsDates,
                                                 DateTimeArray&        datesTS,
                                                 DoubleArray&          weightsTS) const{
    static const string routine = "ClosedFormIntegrateLN::getPriceTSPortfolios";
    try {
        datesTS.resize(0);
        weightsTS.resize(0);

        if(CString::equalsIgnoreCase(termStructMethod, TERM_STRUCTURE_NONE) ||
           obsDates.back() <= valueDate) {
            return;
        }

        DateTimeArray dates;
        if (CString::equalsIgnoreCase(termStructMethod, TERM_STRUCTURE_ALL)) {
            // Compute strictly future dates and fwds at those dates
            dates = valueDate.getFutureDates(obsDates);
            dates.insert(dates.begin(), valueDate);
            if(!dates.size()) {
                return;
            }
        } else {
            // Work with a maturity period
            string tenor = termStructMethod;
            const char *inp = tenor.c_str();
            if(isdigit(*inp)) {
                tenor = "-" + termStructMethod;
            }
            
            // Start from end and add dates per tenor
            MaturityPeriod frequency(tenor);
            DateTime thisDate = obsDates.back();
            while(thisDate > valueDate) {
                if(thisDate > valueDate) {
                    dates.insert(dates.begin(), thisDate);
                }
                thisDate = frequency.toDate(thisDate);
            }
            dates.insert(dates.begin(), valueDate);
        }

        // First strictly future date and fwd will be dropped eventually
        DoubleArray weights(dates.size());
        
        //computes pv Array
        DoubleArray pv(dates.size());
        for(int j=0; j<dates.size(); j++){
            pv[j] = discount->pv(dates[j]);
        }
        
        // Fill in the weights
        for(int iFutStep = 1; iFutStep < dates.size(); iFutStep++) {
            double thisWeight = - log(pv[iFutStep] / pv[iFutStep - 1]);
            weights[iFutStep] = thisWeight;
        }

        // Drop first entries and return
        dates.erase(dates.begin());
        weights.erase(weights.begin());

        datesTS = dates;
        weightsTS = weights;

    } catch (exception& e) {
        throw ModelException(e, routine);
    } 
}

void ClosedFormIntegrateLN::getIntegrator(const Control*  control,
                                          const CAsset*   asset,
                                          const DateTime& valueDate,
                                          const DateTime& maturity,
                                          VanillaContractsRecorderSP recorder,
                                          Integrator1DSP& integrator,
                                          Range&          integrationDomain) const {
    static const string routine = "ClosedFormIntegrateLN::getIntegrator";
    try {
        // Validate Instrument Integration Domain
        const Boundary& lowerInstBound = integrationDomain.getLower();
        const Boundary& upperInstBound = integrationDomain.getUpper();
        double lowerInstStrike;
        if(lowerInstBound.isInfinite()) {
            throw ModelException(routine, "Lower bound for vanilla integration domain cannot be minus infinite");
        } else {
            lowerInstStrike = lowerInstBound.getValue();
            if(Maths::isNegative(lowerInstStrike)) {
                throw ModelException(routine, "Lower bound for vanilla integration domain " +
                    Format::toString(lowerInstStrike) + " cannot be negative");
            }
        }

        // Compute Model Integration Domain
        VolRequestTime volReq;
        IVolProcessedSP vol(asset->getProcessedVol(&volReq));
        double tenor = vol->calcTradingTime(valueDate, maturity);
        double volScale = stdDevInput * sqrt(tenor);
        double lowerModelStrike = exp(-stdDevNbDown * volScale);
        double upperModelStrike = exp(stdDevNbUp * volScale);

        // Come up with an integration domain based on Instrument and Model domains
        // 1) If disjoint then use model
        // 2) If not disjoint then use intersection
        double lowerStrike, upperStrike;
        if(upperInstBound.isInfinite()) {
            // Instrument: semi-infinite, Model: finite
            if(lowerInstStrike > upperModelStrike) {
                // Disjoint so use singleton
                lowerStrike = 1.0;
                upperStrike = 1.0;
            } else {
                // Use intersection
                lowerStrike = Maths::max(lowerModelStrike, lowerInstStrike);
                upperStrike = upperModelStrike;
            }
        } else {
            // Instrument: finite, Model: finite
            double upperInstStrike = upperInstBound.getValue();
            if(lowerInstStrike > upperModelStrike || upperInstStrike < lowerModelStrike) {
                // Disjoint so use singleton
                lowerStrike = 1.0;
                upperStrike = 1.0;
            } else {
                // Use intersection
                lowerStrike = Maths::max(lowerModelStrike, lowerInstStrike);
                upperStrike = Maths::min(upperModelStrike, upperInstStrike);
            }
        }

        // Figure out if we are doing VegaMatrix
        SensitivitySP sens = control->getCurrentSensitivity();
        Sensitivity* sensPtr = sens.get();
        bool doingVegaMatrix = VegaMatrix::TYPE->isInstance(sensPtr);
        bool doingVegaMatrixLite = VegaMatrixLite::TYPE->isInstance(sensPtr);

        // Contrust integrator
        if(doingVegaMatrix || doingVegaMatrixLite || CString::equalsIgnoreCase(integrationMethod, ClosedFormIntegrateLN::DISCRETE_INTEGRATION)) {
            // Figure out number of steps
            double stdDevNbDownFinal = fabs(log(lowerStrike) / volScale);
            double stdDevNbUpFinal = fabs(log(upperStrike) / volScale);
            int nbSteps = (int) ((stdDevNbUpFinal + stdDevNbDownFinal) * stdDevInput * tenor * nbStrikesPerLogStrike);
            // hard coded ... we need a min nb of strikes for short time horizons (sqrt(tenor))
            if (nbSteps < 300) {
                nbSteps = 300;
            }

            // Contruct discrete integrator
            if(doingVegaMatrixLite) {
                // Integrator with Recorder
                integrator = Integrator1DSP(new Trapez1DSimpleRecorder(nbSteps, recorder));
            } else {
                // Integrator without Recorder
                integrator = Integrator1DSP(new Trapez1DSimple(nbSteps));
            }
            // Modify integration range to lie within model limits
            integrationDomain = Range(ClosedBoundary(lowerStrike), ClosedBoundary(upperStrike));
        } else {
            // Create an integrator by name and pass on precision in terms of integrator
            // Leave integration range unaffected
            double absVolPrecisionAdjusted = calcAbsPrecision(tenor);
            double relVolPrecisionAdjusted = calcRelPrecision(tenor);
            integrator = Integrator1D::createIntegrator(integrationMethod,
                                                        absVolPrecisionAdjusted,
                                                        relVolPrecisionAdjusted);
            if(integrationDomain.isOpen()) {
                // Product domain is fine when working on semi-infinite intervals
                // We assume that IntFuncInf does a good job
            } else {
                // Product domain may be difficult to integrate when working on finite intervals
                // IntFuncGeneral etc. don't work well on definite intervals if the function
                // is undersampled
                // Use model bounds
                integrationDomain = Range(ClosedBoundary(lowerStrike), ClosedBoundary(upperStrike));
            }
        }
    } catch (exception& e) {
        throw ModelException(e, routine);
    }
}


/** Implementation of ClosedFormIntegrateLN::IntoProduct interface */
ClosedFormIntegrateLN::IProduct* VarianceSwap::createProduct(ClosedFormIntegrateLN* model) const {
    return new VarianceSwapNumerical(this);
}


bool VolVarSwapLoad()
{
    return (VarianceSwap::TYPE != 0);
}

double VarianceSwap::futureVar(const CAsset*              asset,
                               const YieldCurve*          discount,
                               const DateTime&            valueDate,
                               const DateTime&            volDate,
                               const IModel*              model,
                               Control*                   control) {

    static const string method("VarianceSwap::futureVar");
    try {
        // build a stub variance swap
        VarianceSwap vswap;

        vswap.valueDate    = valueDate;
        vswap.fwdStarting  = false;
        vswap.oneContract  = true;
        vswap.notional     = 1.0;
        vswap.initialSpot  = asset->getSpot();
        vswap.instSettle   = InstrumentSettlementSP(new CashSettlePeriod(0)); // will not affect futureVol
        vswap.asset        = CAssetWrapper(copy(asset));
        vswap.discount     = YieldCurveWrapper(copy(discount));
        vswap.samples.resize(1);
        vswap.samples[0].date = volDate;

        double var = 0.0;

        // check whether we are doing vega matrix ...
        SensitivitySP sens = control->getCurrentSensitivity();
        Sensitivity* sensPtr = sens.get();
        bool vegaMatrix = VegaMatrix::TYPE->isInstance(sensPtr);

        //check which model we have -- ImpliedIntegration or ClosedFormIntegrateLN
        ImpliedIntegration* ii = const_cast<ImpliedIntegration*>(dynamic_cast<const ImpliedIntegration*>(model));
        if(!ii) {
            const ClosedFormIntegrateLN* cfi = dynamic_cast<const ClosedFormIntegrateLN*>(model);
            if(!cfi) { // sthg went wrong .. .
                throw ModelException("Only ImpliedIntegration model or ClosedFormIntegrateLN model allowed.");
            } // now ClosedFormIntegrateLN
            double lowStrike = 0.0;
            double highStrike = 0.0;
            int nbSteps = 0;
            cfi->limits(vswap.asset.get(), vswap.valueDate, volDate, lowStrike, highStrike, nbSteps);
            double fwd = asset->fwdValue(volDate);
            // need a tenor
            ATMVolRequest volReq;
            CVolProcessedBSSP vol(asset->getProcessedVol(&volReq));
            double tenor = vol->calcTradingTime(valueDate, volDate);
            var = vswap.varFromCalls(lowStrike,highStrike,nbSteps,volDate,fwd,
                                     cfi->integrationMethodGet(),
                                     cfi->negativeFwdVarAllowed(),
                                     cfi->calcAbsPrecision(tenor),
                                     cfi->calcRelPrecision(tenor),
                                     vegaMatrix);
        } else { // ImpliedIntegration
            if (vegaMatrix) {
                /* if we are doing vega matrix, then switch from ImpliedIntegration
                    to hard coded version of ClosedFormIntegrateLN, using input from vsn
                    use the existing hard coded values for vega matrix in implied integration */
                ClosedFormIntegrateLN cfModel(0.3,4.0,6.0,100,30,ClosedFormIntegrateLN::DISCRETE_INTEGRATION,
                    ClosedFormIntegrateLN::EQUIDISTANT_STRIKES,sqrt(DBL_EPSILON),sqrt(DBL_EPSILON),false); // these three params are not relevant
                double lowStrike = 0.0;
                double highStrike = 0.0;
                int nbSteps = 0;
                cfModel.limits(vswap.asset.get(), vswap.valueDate, volDate, lowStrike, highStrike,nbSteps);
                // choose 300 as nbSteps in order to be consistent with previous implementation
                var = vswap.varFromCalls(lowStrike,highStrike,300,volDate,asset->fwdValue(volDate),
                                         cfModel.integrationMethodGet(),
                                         true,        // ImpliedIntegration will never fail (true and false have the same impact)
                                         sqrt(DBL_EPSILON), // not relevant
                                         sqrt(DBL_EPSILON),
                                         vegaMatrix); // not relevant
            } else {
                refCountPtr<ImpliedIntegration::IProduct> vsn(vswap.createProduct(ii));
                // vsn = refCountPtr<ImpliedIntegration::IProduct>(); // build a varianceSwapNumerical
                var = ii->integrate(vsn.get()); // get that var, if we are not doing vega matrix
            }
        }
        return var;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


//** VarianceSwapUtil methods */

// Computes integral(startDate,maturity) {sigma^2 dt}.  Note: ignores effect of divs. */
double VarianceSwapUtil::futureFwdVar(
    const CAsset*       asset,
    const YieldCurve*   discount,
    const DateTime&     valueDate,
    const DateTime&     startDate,
    const DateTime&     maturity,
    const IModel*       model,
    Control*            control) {
    static const string method("VarianceSwapUtil::futureFwdVar");
    try {
        // note: we dont need the flush method any more since we now have a samplerMap in ImpliedIntegration
        double varLong  = VarianceSwapUtil::futureVar(asset, discount, valueDate, maturity, model, control);
        double varShort = VarianceSwapUtil::futureVar(asset, discount, valueDate, startDate, model, control);

        double varFwd = varLong - varShort;
        return varFwd;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Utility methods for variance swaps.
Compute integral(valueDate,maturity) {sigma^2 dt}.
Note: ignore effect of divs. */
double VarianceSwapUtil::futureVar(const CAsset*        asset,
                                   const YieldCurve*    discount,
                                   const DateTime&      valueDate,
                                   const DateTime&      volDate,
                                   const IModel*        model,
                                   Control*             control) {
    static const string method("VarianceSwapUtil::futureVar");
    try {
        return VarianceSwap::futureVar(asset, discount, valueDate, volDate, model, control);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Lightweight utility method to price forward starting variance swaps, assuming notional and scaleFactor of 1.0
    called by VarianceIndexForward */
    double VarianceSwapUtil::priceFwdStartingVarSwap(const CAsset*        asset,
        const YieldCurve*    discount,
        double               notional,
        double               scaleFactor,
        bool                 dontScaleByStrike,
        bool                 noDivAdj,
        const DateTime&      valueDate,
        const DateTime&      startDate,
        const DateTime&      maturity,
        double               strikeVol,
        int                  observationsPerYear,
        int                  observationsInSwap,
        const IModel*        model,
        Control*             control)
{
    static const string method("VarianceSwapUtil::priceFwdStartingVarSwap");

    ATMVolRequest volReq;
    CVolProcessedBSSP vol(asset->getProcessedVol(&volReq));
    double t = vol->calcTradingTime(startDate, maturity);

    // need to normalize futureFwdVar return value by t, to compute volSquaredFwd.
    double volSquaredFwd = futureFwdVar(asset, discount, valueDate, startDate, maturity, model, control) / t;

    // add effect for divs.
    if (noDivAdj) {
        volSquaredFwd += VarSwapUtilities::futureDiscreteDivsAdjustment(asset, valueDate, startDate, maturity,
                                                                        observationsPerYear, observationsInSwap);

    }

    return priceVarSwapSimple(volSquaredFwd, strikeVol, notional, scaleFactor, dontScaleByStrike);

}

// returns undiscounted price of a variance swap contract.
double VarianceSwapUtil::priceVarSwapSimple(double totalVar,
                                            double strikeVol,
                                            double notional,
                                            double scaleFactor,
                                            bool dontScaleByStrike)
{
    const string& method = "VarianceSwapUtil::priceVarSwapSimple";
    double result = notional*scaleFactor*(totalVar-strikeVol*strikeVol);

    if (!dontScaleByStrike) {
        if (Maths::isZero(strikeVol)) {
            throw ModelException(method, "Cannot scale by 0.0 strike! ('dontScaleByStrike' is set to false)");
        }
        result /= (2.0*strikeVol);
    }

    return result;
}

DoubleArraySP VarianceSwapUtil::getSensitiveStrikes(OutputNameConstSP       outputName,
                                                    const CAsset*           asset,
                                                    const DateTime&         valueDate,
                                                    const DateTime&         maturityDate,
                                                    const IModel*            model) {
    return VolVarShell::getSensitiveStrikesHelper(outputName,
                                                  valueDate,
                                                  maturityDate,
                                                  asset,
                                                  model);
}


void VarianceSwapUtil::validateBasis(const IModel*        model,
                                     const MarketData*    market,
                                     const CAssetWrapper& asset) {
    static const string routine = "VarianceSwapUtil::validateBasis";

    try {
        // Only allow it for simple assets
        const ClosedFormIntegrateLN* cfModel = dynamic_cast<const ClosedFormIntegrateLN*>(model);
        if(cfModel && cfModel->getUseBasis()) {
            if( XCB::TYPE->isInstance(asset.get()) ||
                Fund::TYPE->isInstance(asset.get()) ||
                ProtAsset::TYPE->isInstance(asset.get()) ||
                ProtEquity::TYPE->isInstance(asset.get()) ) {

                // Require that these assets switch off the basis at the input level
                throw ModelException(routine,
                    "VarSwaps on XCBs, Funds, Ccy Protected and Struck assets cannot use VarSwapBasis. "
                    "Set the useBasis flag to false");
            }
        }
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


/*****************************************/
/** the helper for ClosedFormIntegrateLN */
/*****************************************/

class ClosedFormIntegrateLNHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ClosedFormIntegrateLN, clazz);
        SUPERCLASS(CModelLN);
        EMPTY_SHELL_METHOD(defaultClosedFormIntegrateLN);
        FIELD(stdDevInput, "Std Dev for Strikes");
        FIELD_MAKE_OPTIONAL(stdDevInput);
        FIELD(stdDevNbUp, "Nb of stdDevs on the upside");
        FIELD_MAKE_OPTIONAL(stdDevNbUp);
        FIELD(stdDevNbDown, "Nb of stdDevs on the downside");
        FIELD_MAKE_OPTIONAL(stdDevNbDown);
        FIELD(nbStrikesPerLogStrike, "Nb of strikes per log strike");
        FIELD_MAKE_OPTIONAL(nbStrikesPerLogStrike);
        FIELD(nbStrikesVega, "nb of strikes for vega matrix");
        FIELD_MAKE_OPTIONAL(nbStrikesVega);
        FIELD(integrationMethod, "what integration method");
        FIELD_MAKE_OPTIONAL(integrationMethod);
        FIELD(strikesForVegaMatrix, "what strikes for vega matrix");
        FIELD_MAKE_OPTIONAL(strikesForVegaMatrix);
        FIELD(absVolPrecision, "abs precision of integrator");
        FIELD_MAKE_OPTIONAL(absVolPrecision);
        FIELD(relVolPrecision, "rel precision of integrator");
        FIELD_MAKE_OPTIONAL(relVolPrecision);
        FIELD(useBasis, "whether use cutoff method");
        FIELD_MAKE_OPTIONAL(useBasis);
        FIELD(divMethodology, "divMethodology used to compute prices with dividends");
        FIELD_MAKE_OPTIONAL(divMethodology);
        FIELD(termStructMethod, "Methodology for term structure of portfolios");
        FIELD_MAKE_OPTIONAL(termStructMethod);
        FIELD(names, "names");
        FIELD_MAKE_TRANSIENT(names);
        FIELD(volSurfaceStrikes, "volSurfaceStrikes");
        FIELD_MAKE_TRANSIENT(volSurfaceStrikes);
    }
    // for ClosedFormIntegrateLN::IIntoProduct
    static void loadIntoProduct(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER_INTERFACE(ClosedFormIntegrateLN::IIntoProduct, clazz);
        EXTENDS(Model::IModelIntoProduct);
    }

    static IObject* defaultClosedFormIntegrateLN(){
        return new ClosedFormIntegrateLN();
    }
};

/***********************************************/
/** now some details for ClosedFormIntegrateLN */
/***********************************************/

/** Validation */
void ClosedFormIntegrateLN::validatePop2Object() {
    static const string method = "ClosedFormIntegrateLN::validatePop2Object";
    try {
        // Map default methodology to desired methodology
        if(CString::equalsIgnoreCase(integrationMethod, ClosedFormIntegrateLN::DEFAULT_INTEGRATION_METHOD)) {
            integrationMethod = getDefaultIntegrationMethod();
        }
        if(CString::equalsIgnoreCase(strikesForVegaMatrix, ClosedFormIntegrateLN::DEFAULT_VEGAMATRIX_STRIKES)) {
            strikesForVegaMatrix = getDefaultStrikesForVegaMatrix();
        }
        if (Maths::isNegative(stdDevInput)){
            throw ModelException(method,
                "Input for standard deviation must not be negative, but is " +
                    Format::toString(stdDevInput) + ".");
        }
        if (Maths::isNegative(stdDevNbUp)){
            throw ModelException(method,
                "Input for number of up standard deviations must not be negative, but is " +
                    Format::toString(stdDevNbUp) + ".");
        }
        if (Maths::isNegative(stdDevNbDown)){
            throw ModelException(method,
                "Input for number of down standard deviations must not be negative, but is " +
                    Format::toString(stdDevNbDown) + ".");
        }
        if (Maths::isNegative(absVolPrecision)){
            throw ModelException(method,
                "Input for abs precision must be positive, but is " +
                    Format::toString(absVolPrecision) + ".");
        }
        if (absVolPrecision>MAX_TOLERANCE){ // at least 1 bp in terms of vol
            throw ModelException(method,
                "Input for abs precision must be less than 0.0001, but is " +
                    Format::toString(absVolPrecision) + ".");
        }
        if (Maths::isNegative(relVolPrecision)){
            throw ModelException(method,
                "Input for rel precision must be positive, but is " +
                    Format::toString(relVolPrecision) + ".");
        }
        if (relVolPrecision>MAX_TOLERANCE){
            throw ModelException(method,
                "Input for relVolPrecision must be less than 0.0001, but is " +
                    Format::toString(relVolPrecision) + ".");
        }
        if (   !CString::equalsIgnoreCase(strikesForVegaMatrix, EQUIDISTANT_VEGA)
            && !CString::equalsIgnoreCase(strikesForVegaMatrix, EQUIDISTANT_STRIKES)
            && !CString::equalsIgnoreCase(strikesForVegaMatrix, VOLSURFACE_STRIKES) ) {
            throw ModelException(
                "Wrong input for strikesForVegaMatrix, input is " +
                strikesForVegaMatrix +
                ", but has to be either " +
                DEFAULT_VEGAMATRIX_STRIKES + " , " +
                EQUIDISTANT_VEGA + " , " +
                EQUIDISTANT_STRIKES + " or " +
                VOLSURFACE_STRIKES);
        }
        
        //Validate divMethodology used to compute prices with dividends
        if (   !CString::equalsIgnoreCase(divMethodology, DIV_DEFAULT)
            && !CString::equalsIgnoreCase(divMethodology, DIV_CONTINUOUS_APPROX)
            && !CString::equalsIgnoreCase(divMethodology, DIV_BRUTE_FORCE)
            && !CString::equalsIgnoreCase(divMethodology, string("2"))
            && !CString::equalsIgnoreCase(divMethodology, string("3"))
            && !CString::equalsIgnoreCase(divMethodology, string("4"))
            && !CString::equalsIgnoreCase(divMethodology, string("5"))
            && !CString::equalsIgnoreCase(divMethodology, string("6"))
            && !CString::equalsIgnoreCase(divMethodology, string("7")) ){
            throw ModelException("divMethodology chosen is wrong");
        }

        // Validate term structure methodology
        if(CString::equalsIgnoreCase(termStructMethod, TERM_STRUCTURE_DEFAULT)) {
            termStructMethod = TERM_STRUCTURE_DEFAULT_TENOR;
        }
        if(!CString::equalsIgnoreCase(termStructMethod, TERM_STRUCTURE_NONE) &&
           !CString::equalsIgnoreCase(termStructMethod, TERM_STRUCTURE_ALL)){
            // Try to see if we can construct a Maturity out of this
            try{
                MaturityPeriod frequency(termStructMethod);
            } catch(exception&) {
               throw ModelException(method, "Unrecognized termStructMethod " +
                termStructMethod + ". It has to be one of " +
                TERM_STRUCTURE_NONE + " or " +
                TERM_STRUCTURE_ALL + " or 1D, 1W, 1M etc");
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Constructor
the integrand is call(k)/ k^2 and thus the interesting part is more close to zero,
hence, nb of stddevs is not symmetric */

ClosedFormIntegrateLN::ClosedFormIntegrateLN():CModelLN(TYPE),
stdDevInput(0.3), stdDevNbUp (4.0), stdDevNbDown(8.0), nbStrikesPerLogStrike(100),
nbStrikesVega(30), integrationMethod(DEFAULT_INTEGRATION_METHOD),
strikesForVegaMatrix(DEFAULT_VEGAMATRIX_STRIKES),
absVolPrecision(DEFAULT_ABSOLUTE_VOL_PRECISION),
relVolPrecision(DEFAULT_RELATIVE_VOL_PRECISION),    // default precision for imsl integrator when adjusted
useBasis(true), termStructMethod(TERM_STRUCTURE_DEFAULT),
divMethodology(DIV_DEFAULT){} 

ClosedFormIntegrateLN::ClosedFormIntegrateLN(const string& volType):CModelLN(TYPE, volType),
stdDevInput(0.3), stdDevNbUp (4.0), stdDevNbDown(8.0), nbStrikesPerLogStrike(100),
nbStrikesVega(30), integrationMethod(DEFAULT_INTEGRATION_METHOD),
strikesForVegaMatrix(DEFAULT_VEGAMATRIX_STRIKES),
absVolPrecision(DEFAULT_ABSOLUTE_VOL_PRECISION),
relVolPrecision(DEFAULT_RELATIVE_VOL_PRECISION),    // default precision for imsl integrator when adjusted
useBasis(true), termStructMethod(TERM_STRUCTURE_DEFAULT),
divMethodology(DIV_DEFAULT){
    validatePop2Object();
}

ClosedFormIntegrateLN::ClosedFormIntegrateLN(double        stdDevInput,
                                             double        stdDevNbUp,
                                             double        stdDevNbDown,
                                             int           nbStrikesPerLogStrike,
                                             double        nbStrikesVega,
                                             const string& integrationMethod,
                                             const string& strikesForVegaMatrix,
                                             double        absVolPrecision,
                                             double        relVolPrecision,
                                             bool          useBasis) : CModelLN(TYPE),
stdDevInput(stdDevInput), stdDevNbUp(stdDevNbUp), stdDevNbDown(stdDevNbDown),
nbStrikesPerLogStrike(nbStrikesPerLogStrike), nbStrikesVega(nbStrikesVega), integrationMethod(integrationMethod),
strikesForVegaMatrix(strikesForVegaMatrix), absVolPrecision(absVolPrecision), relVolPrecision(relVolPrecision),useBasis(useBasis),
termStructMethod(TERM_STRUCTURE_DEFAULT),
divMethodology(DIV_DEFAULT){
    validatePop2Object();
}


/** calculate single price and store result in results */
void ClosedFormIntegrateLN::Price(CInstrument*  instrument,
           CControl*     control,
           CResults*     results){
    static const string method = "ClosedFormIntegrateLN::Price";
    if (!IIntoProduct::TYPE->isInstance(instrument)){
        throw ModelException(method, "Instrument of type "+
                             instrument->getClass()->getName() +
                             " does not support ClosedFormIntegrateLN::IntoProduct");
    }
    try {
        if (instrument->priceDeadInstrument(control, results)) {
            return; // done for a dead instrument
        }

        // cast to ClosedFormIntegrateLN::IIntoProduct and create the product
        IIntoProduct& intoProd = dynamic_cast<IIntoProduct&>(*instrument);
        auto_ptr<IProduct> product(intoProd.createProduct(this));

        product->price(this, control, results);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** computes limits for integration, ie lowStrike and highStrike */
const void ClosedFormIntegrateLN::limits(const Asset*       asset,
                                         const DateTime&    today,
                                         const DateTime&    maturity,
                                         double&            lowStrike,
                                         double&            highStrike,
                                         int&               nbSteps) const {
    static const string method = "ClosedFormIntegrateLN::limits";
    try {
        // compute the forward
        double fwd = asset->fwdValue(maturity);
        ATMVolRequest volReq;
        CVolProcessedBSSP vol(asset->getProcessedVol(&volReq));
        double tenor = vol->calcTradingTime(today, maturity);

        // create strikes using fixed std devs and not the market distribution
        // we do this because the market distribution changes between tweaks and we want to keep
        // a fixed domain for integration for stability
        lowStrike = fwd * exp( - stdDevNbDown * stdDevInput * sqrt(tenor));
        highStrike = fwd * exp( stdDevNbUp * stdDevInput * sqrt(tenor));
        // the nb of steps scales with time since the domain increases with time
        nbSteps = (int) ((stdDevNbUp+stdDevNbDown) * stdDevInput * tenor * nbStrikesPerLogStrike);
        // hard coded ... we need a min nb of strikes for short time horizons (sqrt(tenor))
        if (nbSteps<300) {
            nbSteps = 300;
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

double ClosedFormIntegrateLN::calcAbsPrecision(double& tenor) const {
    // transfrom precision in terms of vol to precision in terms of integrator
    return Maths::max(absVolPrecision * absVolPrecision * tenor / 2.0, MIN_TOLERANCE);
}

double ClosedFormIntegrateLN::calcRelPrecision(double& tenor) const {
    // transfrom precision in terms of vol to precision in terms of integrator
    return Maths::max(relVolPrecision * relVolPrecision * tenor / 2.0, MIN_TOLERANCE);
}

DoubleArraySP ClosedFormIntegrateLN::sensitiveStrikes(OutputNameConstSP outputName,
                                                      const Asset*      asset,
                                                      const DateTime&   valueDate,
                                                      const DateTime&   maturity) {
    static const string method = "ClosedFormIntegrateLN::sensitiveStrikes";
    try {
        DoubleArraySP sensStrikes(new DoubleArray(0));

        if (CString::equalsIgnoreCase(strikesForVegaMatrix, VOLSURFACE_STRIKES)) {
            // the following test has to be separated from the previous line otherwise
            // varswap/vaxtoot.xml was crashing due to the poor compiler on NT.opt
            if (names.get()) {
                // loop over names
                int nbNames = names->size();
                for (int iAsset = 0; iAsset < nbNames; iAsset++) {
                    string name = (*names)[iAsset]->idGet(0);
                    if (CString::equalsIgnoreCase(name, outputName->idGet(0))) {
                        int nbStrikes = volSurfaceStrikes[iAsset].size();
                        sensStrikes->resize(nbStrikes);
                        (*sensStrikes) = volSurfaceStrikes[iAsset];
                    }
                }
            }
        } else {
            // some preparation in order to adjust for quanto adjustment
            SensitiveStrikeDescriptor sensStrikeDesc;
            sensStrikeDesc.forwardOnly = false;
            LinearStrikeVolRequestSP volRequest(new LinearStrikeVolRequest(0.0,
                                                                           valueDate,
                                                                           maturity,
                                                                           false));
            if (CString::equalsIgnoreCase(strikesForVegaMatrix, EQUIDISTANT_STRIKES) ) {
                // check strike input here
                if(nbStrikesVega<1){
                    throw ModelException(method,
                        "Input for nb of strikes for vega matrix must be positive, but is " +
                        Format::toString(nbStrikesVega) + ".");
                }
                double lowStrike= 0.0;
                double highStrike = 0.0;
                int nbSteps;
                // get the limits
                limits(asset, valueDate, maturity, lowStrike, highStrike, nbSteps);
                // calc strikes
                for (int i = 0; i <= nbStrikesVega; i++) {
                    double tempStrike = lowStrike + i * (highStrike-lowStrike) / nbStrikesVega;
                    volRequest->setStrike(tempStrike);
                    asset->getSensitiveStrikes(volRequest.get(), outputName, sensStrikeDesc, sensStrikes);
                }
            } else if (CString::equalsIgnoreCase(strikesForVegaMatrix, EQUIDISTANT_VEGA) ) {
                // check strike input here
                if(nbStrikesVega<2){
                    throw ModelException(method,
                        "Input for nb of strikes for vega matrix must be at least 2, but is " +
                            Format::toString(nbStrikesVega) + ".");
                }
                // back out ATM fwd vol for time horizon = maturity
                ATMVolRequest volReq;
                CVolProcessedBSSP vol(asset->getProcessedVol(&volReq));
                double varInput = vol->CalcVar(valueDate, maturity);
                double fwd = asset->fwdValue(maturity);

                for (int i = 0; i < (int)nbStrikesVega; i++) {
                    double percentile = 0.0;
                    if (i==0) {
                        percentile = 0.00001;
                    } else if (i==(int)nbStrikesVega) {
                        percentile = 0.99999;
                    } else {
                        percentile = (double) i / (double) (nbStrikesVega);
                    }
                    double tempStrike = fwd * exp( (-1.0) * varInput / 2.0 + sqrt(varInput) * N1Inverse(percentile) );
                    volRequest->setStrike(tempStrike);
                    asset->getSensitiveStrikes(volRequest.get(), outputName, sensStrikeDesc, sensStrikes);
                }
            } else {
                throw ModelException("Unrecognized vega matrix methodology " + strikesForVegaMatrix);
            }
        }
    return sensStrikes;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

void ClosedFormIntegrateLN::getMarket(const MarketData*  market,
                                      IInstrumentCollectionSP instruments) {
    static const string method = "ClosedFormIntegrateLN::getMarket";
    try {
        if (CString::equalsIgnoreCase(strikesForVegaMatrix,
                                      VOLSURFACE_STRIKES)){
            double arbitraryShiftSize = 0.0001;
            SensControlPerNameSP shift(new VegaMatrix(arbitraryShiftSize));
            OutputNameArrayConstSP sensNames(shift->names(instruments.get()));
            names = sensNames;
            int nbNames = names->size();
            volSurfaceStrikes.resize(nbNames);
            MarketObjectSP getVolSurface;
            for(int iAsset = 0; iAsset < nbNames; iAsset++) {
                try {
                    getVolSurface = market->GetData((*names)[iAsset]->idGet(0), VolSurface::TYPE);
                } catch (exception& e) {
                    throw ModelException(e, method, "Failed to retrieve strikes from vol surface for VegaMatrix for "+
                        (*names)[iAsset]->idGet(0) + ".");
                }
                VolSurfaceSP myVolSurface = VolSurfaceSP::dynamicCast(getVolSurface);
                volSurfaceStrikes[iAsset] = myVolSurface->getStrikes();
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


/** Override default createMDF in order to set the right MDF */
MarketDataFetcherSP ClosedFormIntegrateLN::createMDF() const{

    // Create a MarketDataFetcher that allows VarSwapBasis to go through
    MarketDataFetcherSP fetcher(new MarketDataFetcherLNSpline(getVolType()));
    if (useBasis){
        MDFUtil::setVarSwapBasis(*fetcher, true);
    }
    return fetcher;
}

CClassConstSP const ClosedFormIntegrateLN::TYPE = CClass::registerClassLoadMethod(
    "ClosedFormIntegrateLN", typeid(ClosedFormIntegrateLN), ClosedFormIntegrateLNHelper::load);

CClassConstSP const ClosedFormIntegrateLN::IIntoProduct::TYPE =
CClass::registerInterfaceLoadMethod("ClosedFormIntegrateLN::IIntoProduct",
                                    typeid(ClosedFormIntegrateLN::IIntoProduct),
                                    ClosedFormIntegrateLNHelper::loadIntoProduct);


//////////////////////////////////////////////////////////////////////////
//                                                                      //
//          Credit support for VarianceSwap                             //
//                                                                      //
//  model vol as an asset with no drift and with vol corresponding to   //
//  vol of vol                                                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

class VarianceSwapCreditSupport : virtual public Generic1FactorCreditSupport
{
public:

    //////////////////////////////////////////////////////////////////////////
    //                                                                      //
    // Local helper class to use an "asset" to mimic variance               //
    // the fwd is 1.0, the vol is the vol of variance due to vol of vol and //
    // mean reversion. the mean reversion aspect of dynamics itself is not  //
    // simulated                                                            //
    //                                                                      //
    //////////////////////////////////////////////////////////////////////////

     class VarAssetForCredit;

    // helper class for VolProcessed
    class VarAssetVolProcForCredit : public CVolProcessedBS
    {
    public:
        static CClassConstSP const TYPE;

        VarAssetVolProcForCredit(const VarAssetForCredit *varAsset=0)
            : CVolProcessedBS(TYPE), varAsset(varAsset){};

        string getName() const { return varAsset->getName(); }

        double calcTradingTime(const DateTime &date1, const DateTime &date2) const
        {   return varAsset->vol->getTimeMetric().yearFrac(date1, date2); }

        TimeMetricConstSP GetTimeMetric() const {
            return TimeMetricConstSP::attachToRef(&varAsset->vol->getTimeMetric());
        }

        double CalcVar(const DateTime &date1, const DateTime &date2) const
        {
            double volOfVol = varAsset->vol->getVolOfVol();
            double meanRev = varAsset->vol->getMeanReversion();
            const TimeMetric &metric = varAsset->vol->getTimeMetric();
            double maxT = metric.yearFrac(varAsset->baseDate, varAsset->maxDate);
            double tau1 = metric.yearFrac(varAsset->baseDate, date1);
            double vol1 = calcVolOfVar(volOfVol, meanRev, maxT, Maths::min(maxT, tau1));
            double tau2 = metric.yearFrac(varAsset->baseDate, date2);
            double vol2 = calcVolOfVar(volOfVol, meanRev, maxT, Maths::min(maxT, tau2));
            return vol2 * vol2 * tau2 - vol1 * vol1 * tau1;
        }

        double CalcVol(const DateTime& date1, const DateTime& date2) const
        {
            double dt = calcTradingTime(date1, date2);
            if( Maths::isZero(dt) ) dt = 1.e-6;
            return sqrt(CalcVar(date1, date2)/dt);
        }

        void populateCompositeVol(CompositeVol* compositeVol) const
        { throw ModelException("VarAssetVolProcForCredit::populateCompositeVol", "Not supported"); }

        static void load(CClassSP& clazz)
        {
            clazz->setPrivate(); // make invisible to EAS/spreadsheet
            REGISTER(VarianceSwapCreditSupport::VarAssetVolProcForCredit, clazz);
            SUPERCLASS(CVolProcessedBS);
            EMPTY_SHELL_METHOD(defaultVarAssetVolProcForCredit);
        }
	    static IObject* defaultVarAssetVolProcForCredit()
        { return new VarAssetVolProcForCredit(); }

        // simple function to get vol of variance with stochastic vol, include sqrt(dtau) factor
        static inline double f(double x){ return (x<1e-6?1:((1-exp(-x))/x)); }

        static double calcVolOfVar(double volOfVol, double meanRev, double T, double tau)
        {
            double dtau = T - tau;
            double x = meanRev * tau;
            double f1 = f(x);
            double f2 = f(2.0 * x);
            double f3 = f(meanRev * dtau);
            double a = (x<1e-4)?(1.0/3.0):((1 - 2.0*f1 + f2)/x/x);
            double b = f3 * f1 * f1;
            double c = f3 * f3 * f2;

            double factor = (tau * tau * a + tau * dtau * b + dtau * dtau * c);
            return volOfVol * sqrt(factor) / T;
        }

    private:
        const VarAssetForCredit *varAsset; // $unregistered
    };

    // helper class for asset
    class VarAssetForCredit : public CAsset
    {
    public:
        static CClassConstSP const TYPE;

        VarAssetForCredit() : CAsset(TYPE){};

        VarAssetForCredit(smartConstPtr<VolSV> vol, const DateTime baseDate, const DateTime maxDate)
            : CAsset(TYPE), vol(vol), baseDate(baseDate), maxDate(maxDate) {};

        string getName() const { return vol->getName(); }

        double getSpot() const { return 1.0; }

        // the IMarketObservable interface for retrieving a single sample
        double pastValue(const DateTime&            sampleDate,
                        const ObservationType*      obsType,
                        const ObservationSource*    source,
                        const FixingType*           fixType,
                        const IObservationOverride* overrides,
                        const SamplingConvention*   sampleRule) const{
            throw ModelException("VarAssetForCredit::pastValue",
                        "Cannot implement this method for VarAssetForCredit");
        }

        // IMarketObservable - retrieve a single observation date
        // Returns false if obs is to be omitted
        bool observationDate(const DateTime&           sampleDate,
                             const ObservationSource*  source,
                             const SamplingConvention* sampleRule,
                             DateTime*                 obsDate) const {
            throw ModelException("VarAssetForCredit::observationDate",
                        "Cannot implement this method for VarAssetForCredit");
        }

        // the IMarketObservable interface for retrieving past samples events
        double addPastSampleEvent(const DateTime&             sampleDate,
                                const ObservationType*      obsType,
                                const ObservationSource*    source,
                                const FixingType*           fixType,
                                const IObservationOverride* overrides,
                                const SamplingConvention*   sampleRule,
                                PastSamplesCollector*        collector) const {
            throw ModelException("VarAssetForCredit::addPastSampleEvent", 
                        "Cannot implement this method for VarAssetForCredit");
        }

        // the IMarketObservable interface for
        // is the given date a holiday for the relevant source
        bool isHoliday(const DateTime& sampleDate,
                       const ObservationSource*   source) const{
            throw ModelException("VarAssetForCredit::isHoliday",
                        "Cannot implement this method for VarAssetForCredit");
        }

        double fwdValue(const DateTime& date) const { return 1.0; }

        CVolProcessed* getProcessedVol(const CVolRequest* volRequest) const
        { return (CVolProcessed *)(new VarAssetVolProcForCredit(this)); }

        void fwdValue(const DateTimeArray &dateList, const FwdValueAlgorithm &algo,
                    CDoubleArray& result) const
        { throw ModelException("VarAssetForCredit::fwdValue", "FwdValueAlgorithm interface not supported"); }

        string getYCName() const
        { throw ModelException("VarAssetForCredit::getYCName", "Not supported"); }

        DateTime settleDate(const DateTime& tradeDate) const
        { throw ModelException("VarAssetForCredit::settleDate", "Not supported"); }

        PDFCalculator* pdfCalculator(const PDFRequest* request) const
        { throw ModelException("VarAssetForCredit::pdfCalculator", "Not supported"); }

        static void load(CClassSP& clazz)
        {
            clazz->setPrivate(); // make invisible to EAS/spreadsheet
            REGISTER(VarianceSwapCreditSupport::VarAssetForCredit, clazz);
            SUPERCLASS(CAsset);
            EMPTY_SHELL_METHOD(defaultVarAssetForCredit);
        }
	    static IObject* defaultVarAssetForCredit()
        { return new VarAssetForCredit(); }

    public:
        smartConstPtr<VolSV> vol; // $unregistered
        DateTime baseDate; // $unregistered
        DateTime maxDate; // $unregistered
    };

public:

    Generic1Factor* getInst() const {
        return (Generic1Factor *)(instr.get());
    }

    // return the variance asset, overwrite the default Generic1FactorCreditSupport implimentation
    virtual CreditUndSP getUnderlier() const
    {
        return CreditUndSP( new AssetUnd( varAsset.get() ) );
    }

    /** preprocess instrument for a given set of path dates */
    virtual void preProcess(const DateTimeArray& dates,
							const DoubleArray& atmFwd,
							const DoubleArray& atmVar);

    /** calculate values for a given path */
    virtual void calcPathValues(
        DoubleArray& results,
        const DateTimeArray& dates,
        const double* spots,
        double spotRef);

    /** return model for this instrument */
    virtual IModelSP getModel();

    /** return instrument's last exposure date */
    virtual DateTime getInstLastExposureDate() const;

    VarianceSwapCreditSupport(CInstrument* inst, CMarketDataSP market);

private:

    IModelSP model;
    VarianceSwap* instrOrig;
    smartPtr<VarianceSwap> instr;

    smartPtr<VarAssetForCredit> varAsset;

    // preprocessed values to speed up calculation
    // MTM(t) = pv(t, T) * (fwd - strike), fwd0 = fwd * pv(0, T), same for strike
    DateTime    maturity;
    double      fwd0, strike0;
    DoubleArray pvs;
};

VarianceSwapCreditSupport::VarianceSwapCreditSupport(
        CInstrument* inst, CMarketDataSP market):
Generic1FactorCreditSupport()
{
    static const string method = "VarianceSwapCreditSupport::VarianceSwapCreditSupport";

    // model created in the constructor
    model = IModelSP(new ClosedFormIntegrateLN("VolPreferred"));
    model->validatePop2Object();

    // call model->getInstrumentAndModelMarket to initiate market data selection
    model->getInstrumentAndModelMarket(market.get(), inst);

    // keep original
    instrOrig = dynamic_cast<VarianceSwap*>(inst);
    instrOrig->Validate();
    // copy instrument
    copyInstrument(instr, instrOrig);

    // we don't handle fwd starting case
    if( instr->samples.front().date.getDate() > instr->getValueDate().getDate() )
        throw ModelException(method, "Forward start sample not supported");
    if (instr->payoffType != "FORWARD")
        throw ModelException(method, "Only support payoffType FORWARD." );

    // get the stochastic vol info
    ATMVolRequest volReq;
    CVolProcessedBSSP vol(Generic1FactorCreditSupport::getAsset()->getProcessedVol(&volReq));
    try {
        VolSVSP volSV = VolSVSP::dynamicCast(market->GetData(vol->getName(), VolSV::TYPE));
        volSV->getMarket(model.get(), market.get(), vol->getName());
        DateTime baseDate = instr->getValueDate();
        DateTime maxDate = instr->samples.back().date;
        varAsset = smartPtr<VarAssetForCredit>(new VarAssetForCredit(volSV, baseDate, maxDate));
    }
    catch ( exception & )
    {
        throw ModelException(method, "Require VolSJ instance for "
            + vol->getName() + " in the market data");
    }
}

/** return model for this instrument */
IModelSP VarianceSwapCreditSupport::getModel()
{
    return model;
}

/** return instrument's last exposure date */
DateTime VarianceSwapCreditSupport::getInstLastExposureDate() const
{
    return instr->endDate(0);
}


/** preprocess instrument for a given set of path dates */
void VarianceSwapCreditSupport::preProcess(const DateTimeArray& dates,
   									       const DoubleArray& atmFwd,
									       const DoubleArray& atmVar)
{
    static const string method = "VarianceSwapCreditSupport::preProcess";

    maturity = instr->samples.back().date;

    // no need to do any processing if past maturity
    // we ignore the diff between maturity and settlement
    if( dates[0].isGreater(maturity) ) return;

    static const double VOL_TWEAK = 0.005;
    int i, numDate = dates.size();
    pvs.resize(numDate);

    // get current MTM and future vol
    OutputRequestSP request(new OutputRequest(OutputRequest::VOL_IN_FUTURE));
    OutputRequestArrayConstSP requests(new OutputRequestArray(1, request));
    const string& requestPacketName = request->getPacketName();
    const string& requestName = request->getRequestName();
    OutputNameSP requestOutName(new OutputName(requestName));
    smartPtr<Control> control(new Control(SensitivityArrayConstSP(   ), requests,0,""));
    CResults result;

    model->Price(instr.get(), control.get(), &result);
    double price = result.retrievePrice();
    double volFuture = result.retrieveScalarGreek(requestPacketName, requestOutName);

    const TimeMetric &metric = varAsset->vol->getTimeMetric();
    double t = metric.yearFrac(varAsset->baseDate, varAsset->maxDate);

    // tweak future vol to compute delta and effective strike
    // this way, we don't need to explicitly deal with volPast and scaling factors
    double volUp = volFuture + VOL_TWEAK;
    instr->price(volUp * volUp * t, model.get(), control.get(), &result); // price needs varFuture as input
    double priceUp = result.retrievePrice();
    double volDn = volFuture - VOL_TWEAK;
    instr->price(volDn * volDn * t, model.get(), control.get(), &result); // price needs varFuture as input
    double priceDn = result.retrievePrice();

    double slope = (priceUp - priceDn)/(volUp * volUp - volDn * volDn);
    if (Maths::isNegative(slope) == Maths::isPositive(instr->notional)){
        throw ModelException(method,
                             "Internal error. VarianceSwap price must increase "
                             "with future vol (if remove sign of notional)");
    }

    fwd0 = slope * volFuture * volFuture;
    strike0 =  fwd0 - price;

    // get pv to valueDate to scale fwd and strike for exposures
    DateTime today = instr->getValueDate();
    for(i=0; i<numDate; i++)
    {
        if( dates[i].isGreater(maturity) ) break;
        pvs[i] = instr->discount->pv(today, dates[i]);
    }

}



/** calculate values for a given path */
void VarianceSwapCreditSupport::calcPathValues(DoubleArray& results,
                                               const DateTimeArray& dates,
                                               const double* spots,
                                               double spotRef)
{
    // this product does not care for forward starting or otherwise
    int i = 0;
    for (; i<dates.size(); i++)
    {
        if( dates[i].isGreater(maturity) ) break;

        results[i] = (fwd0 * spots[i] - strike0)/pvs[i];
    }
    for (; i<dates.size(); i++) results[i] = 0.0;

}

CreditSupportSP VarianceSwap::createCreditSupport(CMarketDataSP market)
{
    return CreditSupportSP(new VarianceSwapCreditSupport(this, market));
}

CClassConstSP const VarianceSwapCreditSupport::VarAssetVolProcForCredit::TYPE = CClass::registerClassLoadMethod(
    "VarianceSwapCreditSupport::VarAssetVolProcForCredit", typeid(VarianceSwapCreditSupport::VarAssetVolProcForCredit), load);

CClassConstSP const VarianceSwapCreditSupport::VarAssetForCredit::TYPE = CClass::registerClassLoadMethod(
    "VarianceSwapCreditSupport::VarAssetForCredit", typeid(VarianceSwapCreditSupport::VarAssetForCredit), load);


// *****************************************************************/
// ************ VANILLA VARSWAP CLASS DECLARATIONS *****************/
// *****************************************************************/
// dependent on VarianceSwap, VolVarShell etc.

class VanillaVarSwap;
typedef smartPtr<VanillaVarSwap> VanillaVarSwapSP;
typedef smartConstPtr<VanillaVarSwap> VanillaVarSwapConstSP;
typedef array<VanillaVarSwapSP, VanillaVarSwap> VanillaVarSwapArray;
typedef smartPtr<VanillaVarSwapArray> VanillaVarSwapArraySP;

DEFINE_TEMPLATE_TYPE(VanillaVarSwapArray);

// Flow variance swap
class VanillaVarSwap: public CInstrument,
                      virtual public VanVSModel::IIntoProduct,
                      virtual public Theta::Shift,
                      virtual public LastSensDate,
                      virtual public ISensitiveStrikes,
                      virtual public PastSamplesEvent::IEventHandler,
                      virtual public ISupportVegaMatrixLite{

public:
    static CClassConstSP const TYPE;
    friend class VanillaVSHelper;
    friend class VanillaVSProd;
    friend class VSWConvert;
    friend class VarSwapAddin;

    /** Support VEGA_MATRIX_LITE */
    virtual void avoidVegaMatrixLite(const IModel* model);
    
    void validatePop2Object();

    /** Called once before the initial pricing */
    virtual void Validate();

    /** Returns the value date (aka today) the instrument is currently
        pricing for */
    virtual DateTime getValueDate() const ;

    /** retrieve market data */
    void GetMarket(const IModel*       model,
                   const CMarketDataSP market);

    /** Implementation of VanVSModel::IntoProduct interface */
    virtual VanVSModel::IProduct* createProduct(const VanVSModel* model) const;

    virtual bool priceDeadInstrument(CControl* control, CResults* results) const;

    // Builds a regular VarianceSwap (or cap) from this
    VarianceSwapSP convert(bool toCap) const;

    int numPastReturns() const;

    virtual bool sensShift(Theta* shift);

    virtual DateTime endDate(const Sensitivity* sensControl) const;

    /** indicates whether VEGA_MATRIX is sensible for this instrument */
    virtual bool avoidVegaMatrix(const IModel* model);

    /** returns all strikes on the vol surface to which
        this instrument is sensitive */
    virtual DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                              const IModel*      model);

    /** Returns the name of the instrument's discount currency. */
    virtual string discountYieldCurveName() const;

    static double volFromVarAdjustment(double futureVarAdjustment, VanillaVarSwapConstSP inst) {
        // Convert to Fair Vol adjustment
        double strike = inst->strike;
        double factor = (double)(inst->observationsPerYear) / (double)(inst->expectedN);
        double volAdjustment = sqrt(Maths::square(strike) + factor * futureVarAdjustment) - strike;
        return volAdjustment;
    }

    // implementation of PastSamplesEvent::IEventHandler interface
    void getEvents(const PastSamplesEvent* samples, IModel* model, 
                   const DateTime& eventDate, EventResults* events) const {
        static const string method("GenericNFBase::getEvents");

        try {
            PastSamplesCollectorSP collector =
                        PastSamplesCollectorSP(new PastSamplesCollector(eventDate));            

            FixingTypeSP fixType(new AssetFixType(asset->getTrueName()));

            CashFlowArraySP smpls = this->samples;
            
            // Fill in historic samples
            for(int i=0; i<smpls->size() &&
                        (*smpls)[i].date <= eventDate; i++) {
                asset->addPastSampleEvent((*smpls)[i].date,
                                          ((*observationBuilder->obsTypes())[i]).get(),
                                          assetHistorySourceObject.get(),
                                          fixType.get(),
                                          0, // overrides use eventually
                                          isdaDateAdjust.get(),
                                          collector.get());
            }

            // now put all the events in the EventsResults object
            collector->finaliseResults(events);

        } catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** Stochastic Interest rates adjustment for Variance Swap */
    static double stochRatesVolAdjustment(
        double rateVolatility,
        double rateMeanReversion,
        double equityRateCorrelation,
        double equityVol,
        TimeMetricConstSP metric,
        VanillaVarSwapConstSP inst) {

        static const string method = "VanillaVarSwap::stochRatesVSwapVolAdjustment";
        try {
            double timeToMat = metric->yearFrac(inst->firstDate, inst->lastDate);
            double rateMeanReversionSq = Maths::square(rateMeanReversion);
            double phi = exp(- rateMeanReversion * timeToMat);

            // Convexity adjustment under T-Fwd measure
            double rateConvxAdj = Maths::square(rateVolatility) / rateMeanReversionSq *
                (timeToMat - 2.0 * (1.0 - phi) / rateMeanReversion + (1.0 - Maths::square(phi)) / (2.0 * rateMeanReversion));
            // Correlation adjustment due to moving from Q to T-Fwd measure
            double corrAdj = equityRateCorrelation * equityVol * rateVolatility / rateMeanReversionSq *
                (rateMeanReversion * timeToMat - 1 + phi);

            // Implied Total Variance adjustment (i.e. not divided by T)
            double futureVarAdjustment = -2.0 * (rateConvxAdj + corrAdj);
            double volAdjustment = volFromVarAdjustment(futureVarAdjustment, inst);
            return volAdjustment;
        } catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** Jumps adjustment for Variance Swap */
    static double jumpsVolAdjustment(
        double jumpsPerYear,
        double logSpotJumpMean,
        double logSpotJumpStdDev,
        TimeMetricConstSP metric,
        VanillaVarSwapConstSP inst) {

        static const string method = "VanillaVarSwap::jumpsVolAdjustment";
        try {
            // Assuming compound Poisson jumps with crash rate lambda and jump size Y. Then
            // VSwap = LogContract + 2 * lambda * E[1 + Y + 0.5*Y^2 - exp(Y)]
            // When Y ~ Normal then we get a closed form

            double timeToMat = metric->yearFrac(inst->firstDate, inst->lastDate);

            // Implied Total Variance adjustment (i.e. not divided by T)
            double jumpVar = Maths::square(logSpotJumpStdDev);
            double futureVarAdjustment = 2.0 * jumpsPerYear * timeToMat * (
                1.0 + logSpotJumpMean + 0.5 * (Maths::square(logSpotJumpMean) + jumpVar) - exp(logSpotJumpMean + jumpVar / 2.0));
            double volAdjustment = volFromVarAdjustment(futureVarAdjustment, inst);
            return volAdjustment;
        } catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    static VanillaVarSwapSP makeParSwap(const CAssetWrapper& asset,
                                        const YieldCurveWrapper& discount,
                                        const DateTime& firstDate,
                                        const DateTime& lastDate,
                                        bool dividendAdjusted,
                                        InstrumentSettlementSP settlement,
                                        const HolidayWrapper& assetHols,
                                        CMarketDataSP market,
                                        IModelSP model,
                                        bool isOption,
                                        double cap,
                                        ScenarioSP scenario) {
        static const string method = "VanillaVarSwap::makeParSwap";
        try {
            // Fill in expectedN  and strike
            HolidayWrapper holsTmp = assetHols;
            holsTmp.getData(model.get(), market.get());
            int expectedN = holsTmp->businessDaysDiff(firstDate, lastDate);

            // Price copy of instrument to get VOL_IN_FUTURE
            OutputRequestArraySP outputs(new OutputRequestArray(0));
            outputs->push_back(OutputRequestSP(new OutputRequest(OutputRequest::VOL_IN_FUTURE)));
            CControlSP ctrl(new Control(
                SensitivityArrayConstSP(   ),
                outputs,
                false,
                ""));

            CInstrumentSP inst(new VanillaVarSwap(asset, discount, firstDate, lastDate, dividendAdjusted,
                settlement, assetHols, expectedN, 0.2, isOption, cap));
            inst->validatePop2Object();
            CResultsSP results(model->go(inst, scenario, ctrl, market));

            // Override strike of original instrument and call validate
            IObjectConstSP obj = results-> retrieveRequestResult(OutputRequest::VOL_IN_FUTURE);
            if(Untweakable::TYPE->isInstance(obj)) {
                const Untweakable* tmp = dynamic_cast<const Untweakable*>(obj.get());
                throw ModelException(method, tmp->getMessage());
            }
            CDoubleConstSP vol(CDoubleConstSP::dynamicCast(obj));
            double parStrike = vol->doubleValue();

            VanillaVarSwapSP parInst(new VanillaVarSwap(asset, discount, firstDate, lastDate, dividendAdjusted,
                settlement, assetHols, expectedN, parStrike, isOption, cap));

            return parInst;
        } catch(exception& e) {
            throw ModelException(e, method);
        }
    }

private:
    // constructor
    VanillaVarSwap(): CInstrument(TYPE), expectedN(0),
                                         subtractMeanVol(false),
		                                 divAdjOnExDate(false),
                                         scaleByStrike(true),
                                         capOnly(false),
                                         cap(2.5),
                                         sampleInterval(1),
                                         firstSample(0.0),
                                         validateExpectedN(true)
                                         {}

    /** Reduced constructor for par swap-only instrument on vanilla asset */
    VanillaVarSwap(const CAssetWrapper& asset,
                   const YieldCurveWrapper& discount,
                   const DateTime& firstDate,
                   const DateTime& lastDate,
                   bool dividendAdjusted,
                   InstrumentSettlementSP settlement,
                   const HolidayWrapper& assetHols,
                   int expectedN,
                   double strike,
                   bool isOption,
                   double cap):
    CInstrument(TYPE), asset(asset), ccyTreatment(CAsset::CCY_TREATMENT_NONE), discount(discount),
    fwdStarting(false), strike(strike), expectedN(expectedN), observationsPerYear(252), subtractMeanVol(false),
    dividendAdjusted(dividendAdjusted), divAdjOnExDate(false),
    scaleByStrike(true), isCapped(isOption), capOnly(isOption), cap(cap),
    instSettle(settlement), premiumSettle(settlement),
    firstDate(firstDate), lastDate(lastDate), sampleTimeRule("EOD"), sampleInterval(1), firstSample(0.0), assetHols(assetHols), validateExpectedN(true) {
        static const string method = "VanillaVarSwap::VanillaVarSwap";

        // Create empty asset history
        CClassConstSP clazz = CClass::forName(typeid(AssetHistory));
        IObjectSP assetHistTmp(clazz->newInstance());
        assetHistory = AssetHistoryWrapper(AssetHistorySP::dynamicCast(assetHistTmp));

        validatePop2Object();
        try {

        } catch(exception& e) {
            throw ModelException(e, method);
        }
    }

    // not implemented
    VanillaVarSwap(const VanillaVarSwap& rhs);
    VanillaVarSwap& operator=(const VanillaVarSwap& rhs);

    // grab relevant samples from cache
    //void grabSamples();

    // Returns the first sample needed in the AssetHistory for pricing
    //DateTime firstSampleNeeded();

    // market data
    DateTime                valueDate;
    CAssetWrapper           asset;
    CAssetWrapper           capAsset;       // Needed incase market data differs between VarSwap model and cap model
    string                  ccyTreatment;
    YieldCurveWrapper       discount;
    AssetHistoryWrapper     assetHistory;           // superceded by assetHistorySource

    // product data
    bool                    fwdStarting;            // input ignored - set in theta roll depending on valueDate <> firstDate
    //bool                    oneContract;            // always true
    double                  strike;
    int                     expectedN;              // Expected number of returns over life of contract. Set on initial instrument creation.
    int                     observationsPerYear;
    bool                    subtractMeanVol;
    bool                    dividendAdjusted;
    bool                    divAdjOnExDate;
    bool                    scaleByStrike;
    bool                    isCapped;               // (isCapped, capOnly) = (false, false): price vswap only; (false, true): validated against
    bool                    capOnly;                // (isCapped, capOnly) = (true, false): price vswap and cap; (true, true): price cap only
    double                  cap;                    // Cap expressed as multiplier by strike vol - this is different to the cap on the
                                                    // VolVarShell which is a multiplier by strike variance
    InstrumentSettlementSP  instSettle;
    InstrumentSettlementSP  premiumSettle;

    // Sample definition data - start, end, frequency and timing rule
    // No holidays are used in determination of samples used.
    // Instead all samples are used - if there is a holiday it will contribute zero variance
    DateTime                firstDate;          // time superceded by startObsType
    DateTime                lastDate;            // time superceded by endObsType
    string                  sampleTimeRule;     // superceded by obsType
    int                     sampleInterval;     // better to use a MaturityPeriod?
    double                  firstSample;        // at time of trade


    // New fields for AssetHistory. Used in construction of an ObservationBuilderDaily
    string                  startObsType;
    string                  obsType;
    string                  endObsType;
    string                  assetHistorySource;
    HolidayWrapper          assetHols;          // Instrument-Level holidays - sampling now excludes these.
                                                // Also used for numFutureReturns and numPastReturns
    SamplingConventionSP    isdaDateAdjust;               // How/whether to roll/omit if sample date is a holiday

    // internal data
    CashFlowArraySP         samples;
    IObservationBuilderSP   observationBuilder;
    ObservationSourceSP     assetHistorySourceObject; //!< object for asset history source
    bool                    validateExpectedN;
};

void VanillaVarSwap::avoidVegaMatrixLite(const IModel* model) {
    static const string method = "VanillaVarSwap::avoidVegaMatrixLite";
    
    if (avoidVegaMatrix(model)) {
        // Check basic VEGA_MATRIX
        throw ModelException(method, "Instrument does not support VEGA_MATRIX and hence not VEGA_MATRIX_LITE");
    } else if (!SimpleEquity::TYPE->isInstance(asset.get())){
        // Allow LITE only for SimpleEquity    
        throw ModelException(method, "Only SimpleEquity underlyings supported");
    } else if(isCapped) {
        if(capOnly) {
            // Disallow Caps
            throw ModelException (method, "CapOnly payoffs not supported for VEGA_MATRIX_LITE");
        } else {
            // Disallow Capped Swaps
            throw ModelException (method, "Capped Swap payoffs not supported for VEGA_MATRIX_LITE");
        }
    } else {
        // Make sure we have ClosedFormIntegrateLN
        const VanVSModel* umbrella = dynamic_cast<const VanVSModel*>(model);
        if (!umbrella) {
            throw ModelException("VEGA_MATRIX_LITE is only supported for model VanVSModel.");
        }
        if (!ClosedFormIntegrateLN::TYPE->isInstance(umbrella->varSwapModel)) {
            throw ModelException("VEGA_MATRIX_LITE is only supported for ClosedFormIntegrateLN.");
        }
    }
}



void VanillaVarSwap::validatePop2Object() {
    static const string method = "VanillaVarSwap::validatePop2Object";
    try {

        if (firstDate >= lastDate) {
            throw ModelException(method, "firstDate (" + firstDate.toString() +
                                 ") is equal or after the lastDate (" + lastDate.toString()
                                 + ")");
        }

        if (!isCapped && capOnly) {
            // Can't be uncapped and asking to price cap only
            throw ModelException(method, "Price for the cap only is requested but the instrument is not capped");
        }

        if (isCapped && !Maths::isPositive(cap)) {
            throw ModelException(method, "Variance swap is capped but has a cap <= 0");
        }

        if (subtractMeanVol) {
            // Subtracting mean vol not supported on vanilla instrument
            throw ModelException(method, "Subtract Mean Vol flag must be false on the Vanilla Variance Swap");
        }

        if (sampleInterval != 1) {
            throw ModelException(method, "A sample interval of " + Format::toString(sampleInterval)
                                          + " is not presently supported");
        }

        capAsset = CAssetWrapper(asset.getName());

        // Validate that firstDate and lastDate are not weekends
        if (firstDate.isWeekend()) {
            throw ModelException(method, "Start date (" + firstDate.toString() +
                                         ") can't be on a weekend");
        }
        if (lastDate.isWeekend()) {
            throw ModelException(method, "End date (" + lastDate.toString() +
                                         ") can't be on a weekend");
        }

        // If at least one of startObsType, obsType, endObsType, assetHistorySource are passed
        // then validate that all are passed
        if (!(startObsType.empty() &&
            endObsType.empty() &&
            obsType.empty() &&
            assetHistorySource.empty())) {
                
            if (startObsType.empty() ||
                endObsType.empty() ||
                obsType.empty() ||
                assetHistorySource.empty()) {

                    throw ModelException(method, "None or all of startObsType, obsType, endObsType, "
                                                  "assetHistorySource must be passed");

            }
        }

        if (!isdaDateAdjust) {
            isdaDateAdjust = SamplingConventionSP(new UnadjustedConvention()); 
        }
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

void VanillaVarSwap::Validate() {
static const string method = "VanillaVarSwap::Validate";
    try {
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

/** Returns the value date (aka today) the instrument is currently
pricing for */
DateTime VanillaVarSwap::getValueDate() const {
    return valueDate;
}

/** retrieve market data */
// We assume that the model for the varswap is passed in the composite model
// (irrespective of isCapped/capOnly flags)
void VanillaVarSwap::GetMarket(const IModel*       model,
                               const CMarketDataSP market) {
static const string method = "VanillaVarSwap::GetMarket";
    try {

        const VanVSModel* compositeModel = dynamic_cast<const VanVSModel*>(model);

        if (!compositeModel) {
            throw ModelException(method, "Passed model must of type VanVSModel");
        }

        const IModel* varSwapModel = compositeModel->varSwapModel.get();
        if (!varSwapModel) {
            throw ModelException(method, "No variance swap model has been supplied");
        }

        // If capped, populate capModel and capAsset -
        // the assumption is that only the asset market data may differ between models
        if (isCapped) {
            if (compositeModel->capModel.get()) {
                VarCapCVModel*     tmpModel = dynamic_cast<VarCapCVModel*>(compositeModel->capModel.get());
                if(tmpModel){
                    tmpModel->volModel->getMarket(market.get(),
                                                  IInstrumentCollection::singleton(this));
                    CAsset::getAssetMarketData(tmpModel->volModel.get(), market.get(), ccyTreatment,
                                               discount, capAsset);
                } else {
                    compositeModel->capModel->getMarket(market.get(),
                                                        IInstrumentCollection::singleton(this));
                    CAsset::getAssetMarketData(compositeModel->capModel.get(), market.get(), ccyTreatment,
                                               discount, capAsset);
                }
            } else {
                throw ModelException(method, "Model to price cap must be provided if variance swap is capped");
            }
        }

        market->GetReferenceDate(valueDate);

        discount.getData(varSwapModel, market);
        instSettle->getMarket(varSwapModel, market.get());

        // Populate asset
        CAsset::getAssetMarketData(varSwapModel, market.get(), ccyTreatment,
                                   discount, asset);

        if (premiumSettle.get()) {
            premiumSettle->getMarket(varSwapModel, market.get());
        }

       
        VarianceSwapUtil::validateBasis(varSwapModel, market.get(), asset);

        if (assetHols.getName() == "") {

            // Default to hols on asset (these are settlement hols!)
            // Holiday is immutable so ok to const_cast
            HolidaySP settleHols(const_cast<Holiday *>(AssetUtil::getHoliday(asset.get()).get()));

            assetHols.setObject(settleHols);

        } else {
            assetHols.getData(model, market);
        }

        // Build a transient ObservationBuilderDaily from the following:
        // startDate = firstDate, endDate = lastDate
        // excludeWeekends = assetHols.useWeekends
        // excludeDates = assetHols.dates
        // startObsType, endObsType, obsType

        // For AssetHistory backwards-compatibility. Default new fields as appropriate.
        if (startObsType.empty() &&
            endObsType.empty() &&
            obsType.empty() &&
            assetHistorySource.empty()) {
            // New fields not populated
            if (sampleTimeRule == DateTime::START_OF_DAY) {
                obsType = "Open";
            } else {
                obsType = "Close";                
            }

            if (lastDate.getTime() == DateTime::START_OF_DAY_TIME) {
                endObsType = "Open";
            } else if (lastDate.getTime() == DateTime::END_OF_DAY_TIME - 1) {
                endObsType = "OSPOption";
            } else {
                endObsType = "Close";
            }

            if (firstDate.getTime() == DateTime::START_OF_DAY_TIME) {
                startObsType = "Open";
            } else if (firstDate.getTime() == DateTime::END_OF_DAY_TIME - 1) {
                startObsType = "OSPOption";
            } else {
                startObsType = "Close";
            }

            // Only attempt to populate the assetHistory object from the market cache if it is actually needed for past samples
            // Otherwise default an empty AssetHistory object. 
            if (valueDate >= firstDate)
            {
               assetHistory.getData(varSwapModel, market.get());
            }
            else {
                // with a matching name exists in the market cache
                MarketObjectSP emptyAssetHist = MarketObjectSP(dynamic_cast<MarketObject *>(
                                                        new AssetHistory("dummy", 0)));
                assetHistory.setObject(emptyAssetHist);
            }

            assetHistorySource = assetHistory->getSource();
        }

        // populate an assetHistorySourceObject
        assetHistorySourceObject 
            = ObservationSourceSP(new ObservationSource(assetHistorySource));

        bool excludeWeekends;
        DateTimeArray excludeDates = *(assetHols->toALIB(excludeWeekends).get()); // bad(ish)

        observationBuilder = IObservationBuilderSP(
                                new ObservationBuilderDaily(firstDate,
                                                            lastDate,
                                                            excludeWeekends,
                                                            excludeDates,
                                                            startObsType,
                                                            obsType,
                                                            endObsType));

        DateTimeArraySP obsDates = observationBuilder->dateList();
        DoubleArraySP obsSamples = DoubleArraySP(new DoubleArray(obsDates->size()));

        // Default expectedN if not passed
        if (!expectedN) {
            expectedN = obsSamples->size()-1;
        }

        // Validate expectedN against the number of samples
        if (validateExpectedN) {
            if (obsSamples->size()-1 != expectedN) {
                throw ModelException(method, 
                    "Number of returns based on instrument holidays (" + Format::toString(obsSamples->size()-1) +
                    ") is not equal to expectedN (" + Format::toString(expectedN) + ")");
            }
        }

		// need to grab the samples we want from the asset history at
		// this point as for any later time shifts we need to put in
		// what tomorrow's sample would be.
        // note we need the first sample in the asset history cache (irrespective of
        // whether the firstSample override is used).
        VarSwapUtilities::populateSamples(asset.get(), assetHistorySourceObject.get(),
                                          valueDate, *obsDates, *observationBuilder->obsTypes(), 
                                          *obsSamples, isdaDateAdjust);

        samples = CashFlowArraySP(new CashFlowArray(0));
        int i=0;
        for (i=0; i<obsSamples->size(); i++) {
            samples->push_back(CashFlow((*obsDates)[i], (*obsSamples)[i]));
        }

        // Override first sample if specified by instrument
        if(!Maths::isZero(firstSample)) {
            samples->front() = CashFlow((*samples)[0].date, firstSample);
        }

        //grabSamples();
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

// Returns the first sample needed in the AssetHistory for pricing
//DateTime VanillaVarSwap::firstSampleNeeded() {
//    DateTime startDate;
//    if (!firstSample) {
//        startDate = firstDate;
//    } else {
//        startDate = DateTime(firstDate.getDate() + sampleInterval, DateTime::timeConvert(sampleTimeRule));
//    }
//    return startDate;
//}

bool VanillaVarSwap::priceDeadInstrument(CControl* control, CResults* results) const {
    static const string method = "VanillaVarSwap::priceDeadInstrument";
    try {
        //
        //               valueDate >= lastDate

        // Preferably allow the real VarianceSwap to take care of this
        return false;

    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}


// Date when to stop tweaking
DateTime VanillaVarSwap::endDate(const Sensitivity* sensControl) const {
    static const string method = "VanillaVarSwap::endDate";
    try {
        const DateTime& instEnd  = instSettle->settles(lastDate, asset.get());
        const DateTime& assetEnd = asset->settleDate(lastDate);
        const DateTime& end = instEnd.isGreater(assetEnd) ? instEnd : assetEnd;
        return end;
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

// Indicates whether VEGA_MATRIX is sensible for this instrument
bool VanillaVarSwap::avoidVegaMatrix(const IModel* model) {
    static const string method = "VanillaVarSwap::avoidVegaMatrix";
    try {
        const VanVSModel* umbrella = dynamic_cast<const VanVSModel*>(model);

        if (!umbrella) {
            throw ModelException("Model must be of type VanVSModel for VanillaVarSwap.");
        }

        if (ImpliedIntegration::TYPE->isInstance(umbrella->varSwapModel)) {
            return false;
        } else if (ClosedFormIntegrateLN::TYPE->isInstance(umbrella->varSwapModel)) {
            return false;
        } else {
            return true;
        }
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

// Returns all strikes on the vol surface to which this instrument is sensitive
DoubleArraySP VanillaVarSwap::getSensitiveStrikes(OutputNameConstSP outputName,
                                                  const IModel*      model) {
    static const string method = "VanillaVarSwap::getSensitiveStrikes";
    try {
        if (avoidVegaMatrix(model)) {
            throw ModelException(method, "VEGA_MATRIX is not valid for this instrument.");
        }

        const VanVSModel* umbrella = dynamic_cast<const VanVSModel*>(model);

        if (!umbrella) {
            throw ModelException("Model must be of type VanVSModel for VanillaVarSwap.");
        }

        return VarianceSwapUtil::getSensitiveStrikes(outputName,
                                                     asset.get(),
                                                     valueDate,
                                                     lastDate,
                                                     umbrella->varSwapModel.get());

    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

VarianceSwapSP VanillaVarSwap::convert(bool toCap) const {
    static const string method = "VanillaVarSwap::convert";
    try {
        string payoffType = "FORWARD";
        if (toCap) {
            payoffType += "_CAP";
        }

        // Translate the cap from a multiplier by strike vol to a multiplier by strike variance
        double varMultiplier = cap * cap;

        VarianceSwapSP vswap(new VarianceSwap(valueDate,           firstDate,
                                              fwdStarting,         firstSample,        ccyTreatment,
                                              instSettle,          premiumSettle,
                                              asset,               discount,
                                              (*samples.get()),    strike,
                                              observationsPerYear, subtractMeanVol,
                                              payoffType,          !scaleByStrike,
                                              varMultiplier,       !dividendAdjusted,  divAdjOnExDate,
                                              true,                numPastReturns(),
                                              expectedN));

        if (toCap) {
            // Different asset is used incase of model-contigent market data
            vswap->asset = capAsset;
        }

        // Validate the real var swap
        vswap->validatePop2Object();
        vswap->Validate();

        if (toCap) {
            // for control variate
            vswap->marketAsset = asset;
        }

        return vswap;
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

// Calculates the number of returns up to the next sample date
// (where next sample date >= value date)
// Function will need modifying to support frequency != 1
// This is used for PAST_WEIGHT dependent outputs (not FV)
int VanillaVarSwap::numPastReturns() const {
    static const string method = "VanillaVarSwap::numPastReturns";
    try {

        int numPastReturns = 0;
        if (valueDate >= lastDate) {
            numPastReturns = expectedN;
        }
        else if (valueDate > firstDate) {

            // Find next sample date
            int i = 0;
            while (i < samples->size() && (*samples)[i].date < valueDate) {
                i++;
            }

            // expected number of returns from next sample date
            int numFutureReturns = assetHols->businessDaysDiff((*samples)[i].date, lastDate);

            // number of returns to next sample is defined as:
            // (total num returns expected on instrument creation) - (expected num returns from next sample)
            numPastReturns = expectedN - numFutureReturns;
            numPastReturns = Maths::max(numPastReturns, 1);
        }

        return numPastReturns;

    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

bool VanillaVarSwap::sensShift(Theta* shift) {
    static const string method("VanillaVarSwap::sensShift (theta)");
    try {
        DateTime newDate = shift->rollDate(valueDate);

        if (!samples) {
            throw ModelException(method, "samples not grabbed");
        }

        // fill vol sample point if needed
        double spot    = 0.0;
        bool   useSpot = !shift->useAssetFwds();

        if (useSpot) {
            spot = asset->getSpot();
        }

        bool pastRollDate = false;
        int  i = 0;

        while (!pastRollDate && i < samples->size() ) {
            if (((*samples)[i].date.isGreater(valueDate) &&
                 !(*samples)[i].date.isGreater(newDate)) ||
                ((*samples)[i].date.equals(valueDate)    &&
                 Maths::isZero((*samples)[i].amount))) {
                if (useSpot) {
                    (*samples)[i].amount = spot;
                }
                else {
                    (*samples)[i].amount = asset->fwdValue((*samples)[i].date);
                }
            }
            else if ((*samples)[i].date.isGreater(newDate)) {
                pastRollDate = true;
            }
            ++i;
        }

        /* If fwd start date falls between value date and theta date
           have to set first sample (populated in initialSpot on VarianceSwap construction).*/

        // fwdStarting flag is not checked anymore in validatePop2Object, so make sure that it has the correct value
        fwdStarting = (samples->front().date > valueDate);
        if ( fwdStarting                           &&
             firstDate.isGreaterOrEqual(valueDate) &&
             newDate.isGreaterOrEqual(firstDate)   )
        {
            // returns the fwd price if shift is a theta fs
            firstSample = samples->front().amount;

            // not fwd starting anymore
            fwdStarting = false;
        }

        valueDate = newDate;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
    return true; // our components have theta type sensitivity
}

// Superceded by VarSwapUtilities::populateSamples()
//void VanillaVarSwap::grabSamples() {
//    static const string method = "VanillaVarSwap::grabSamples";
//    try {
//        // Whether first sample is specified at the instrument level
//        bool overrideFirstSample = !Maths::isZero(firstSample);
//
//        // Generate the sample dates from the sample schedule
//        samples = CashFlowArraySP(new CashFlowArray(0));
//
//        if(!overrideFirstSample) {
//            // start date
//            samples->push_back(CashFlow(firstDate, 0));
//        }
//
//        // intermediate dates
//        DateTime date(firstDate.getDate() + sampleInterval, DateTime::timeConvert(sampleTimeRule));
//
//		const Holiday* hols = assetHols.operator->();
//
//        while (date.getDate() < lastDate.getDate()) {
//            if (hols->isBusinessDay(date)) {
//                CashFlow dateOnly(date, 0);
//                samples->push_back(dateOnly);
//            }
//            date = date.rollDate(sampleInterval);
//        }
//
//        // end date
//        samples->push_back(CashFlow(lastDate, 0));
////        int useMyCode = 1;
//
//        // Populate sample values
////        SimpleEquity* simpEq = dynamic_cast<SimpleEquity*>(asset.get());
////        if (!useMyCode || !simpEq) {
//            assetHistory.getSP()->getSamples(samples, valueDate);
////        } else {
//            // RIGHT NOW FOR TRUE ASSET HISTORY
//            // ONLY IMPLEMENTED FOR SIMPLEEQUITY AT THE MOMENT
///*            ObservationTypeSP obsType(new ObservationExact());
//            ObservationSourceSP source = IMarketObservable::getDefaultObsSource();
//            FixingTypeSP fixType(new AssetFixingType(simpEq->getTrueName()));
//            SamplingConventionSP sampleRule(new UnadjustedConvention()); 
//            for (int i =0; i < samples->size() &&  
//                        (*samples)[i].date <= valueDate; i++) {
//                (*samples)[i].amount = simpEq->pastValue((*samples)[i].date,
//                                                         obsType.get(),
//                                                         source.get(),
//                                                         fixType.get(),
//                                                         0,
//                                                         sampleRule.get());
//            }
//        }*/
//
//        // Populate first sample if specified by instrument
//        if(overrideFirstSample) {
//            samples->insert(samples->begin(), CashFlow(firstDate, firstSample));
//        }
//    }
//    catch (exception& e){
//        throw ModelException(e, method);
//    }
//}


/** Returns the name of the instrument's discount currency */
string VanillaVarSwap::discountYieldCurveName() const {
    return discount.getName();
}


class VanillaVSProd: virtual public VanVSModel::IProduct {
private:

    const VanillaVarSwap* inst; // a reference

public:

    VanillaVSProd(const VanillaVarSwap* instrument): inst(instrument){}

    // This is the method responsible for pricing the VarSwap and the cap and aggregating
    // the results
    void price(VanVSModel*       model,
               Control*          control,
               CResults*         results) const;
};

// *****************************************************************/
// ************ VANILLA VARSWAP CLASS DEFINITIONS ******************/
// *****************************************************************/

const string VanVSModel::METHOD_NONE     = "none";

VanVSModel::IProduct* VanillaVarSwap::createProduct(const VanVSModel* model) const {
    return new VanillaVSProd(this);
}


class  VanillaVSHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        clazz->setDescription("Vanilla Variance Swap");
        REGISTER(VanillaVarSwap, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(VanVSModel::IIntoProduct);
        IMPLEMENTS(LastSensDate);
        IMPLEMENTS(ISensitiveStrikes);
        IMPLEMENTS(Theta::Shift);
        IMPLEMENTS(PastSamplesEvent::IEventHandler);
        IMPLEMENTS(ISupportVegaMatrixLite);
        EMPTY_SHELL_METHOD(defaultVanillaVarSwap);
        FIELD(valueDate, "Value date");
        FIELD_MAKE_OPTIONAL(valueDate);
        FIELD(asset, "Variance Swap underlying");
        FIELD(capAsset, "Variance Swap cap underlying");
        FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(capAsset);
        FIELD(ccyTreatment, "Currency Treatment");
        FIELD(discount, "Discount curve");
        FIELD(assetHistory, "Underlying sample history. Superceded by assetHistorySource");    // superceded by assetHistorySource
        FIELD_MAKE_OPTIONAL(assetHistory);                  
        FIELD(fwdStarting, "Ignored");
        FIELD_MAKE_OPTIONAL(fwdStarting);
        FIELD(strike, "Contract strike in vol unit");
        FIELD(expectedN, "Expected number of observations over life of contract");
        FIELD_MAKE_OPTIONAL(expectedN);
        FIELD(observationsPerYear, "Observations per year");
        FIELD(subtractMeanVol, "whether mean vol is subtracted");
        FIELD_MAKE_OPTIONAL(subtractMeanVol);               //default to false
        FIELD(dividendAdjusted, "whether to adjust for dividends");
		FIELD(divAdjOnExDate, "Whether to add to ex-date sample or subtract from previous sample");
		FIELD_MAKE_OPTIONAL(divAdjOnExDate);
        FIELD(scaleByStrike, "whether to divide by 2K");
        FIELD_MAKE_OPTIONAL(scaleByStrike);                 //default to true
        FIELD(isCapped, "whether the variance swap is capped");
        FIELD(capOnly, "whether to price the cap only");
        FIELD_MAKE_OPTIONAL(capOnly);
        FIELD(cap, "Cap expressed as a multiplier by strike vol");
        FIELD_MAKE_OPTIONAL(cap);
        FIELD(instSettle, "Instrument settlement at maturity");
        FIELD(premiumSettle, "Premium settlement");
        FIELD_MAKE_OPTIONAL(premiumSettle);
        FIELD(firstDate, "First sample date");       // time superceded by startObsType
        FIELD(lastDate, "Last sample date");         // time superceded by startObsType
        FIELD(sampleTimeRule, "Sample timing rule. Superceded by obsType"); // superceded by obsType
        FIELD_MAKE_OPTIONAL(sampleTimeRule);
        FIELD(sampleInterval, "Sampling interval in calendar days"); // must be 1
        FIELD_MAKE_OPTIONAL(sampleInterval);
        FIELD(firstSample, "Initial sample");
        FIELD_MAKE_OPTIONAL(firstSample);
        FIELD(startObsType, "Observation type for start date");
        FIELD_MAKE_OPTIONAL(startObsType);
        FIELD(obsType, "Observation type for intermediate dates");
        FIELD_MAKE_OPTIONAL(obsType);
        FIELD(endObsType, "Observation type for end date");
        FIELD_MAKE_OPTIONAL(endObsType);
        FIELD(assetHistorySource, "Asset History Source");
        FIELD_MAKE_OPTIONAL(assetHistorySource);
        FIELD_NO_DESC(observationBuilder);
        FIELD_MAKE_TRANSIENT(observationBuilder);
        FIELD_NO_DESC(samples);
        FIELD_MAKE_TRANSIENT(samples);
        FIELD(assetHols, "Instrument-level holidays");          // used in ObservationBuilderDaily
        FIELD_MAKE_OPTIONAL(assetHols);
        FIELD(isdaDateAdjust, "ISDA Convention for sampling on holidays");
        FIELD_MAKE_OPTIONAL(isdaDateAdjust);
        FIELD(assetHistorySourceObject, "Asset history source object");
        FIELD_MAKE_TRANSIENT(assetHistorySourceObject);
        FIELD(validateExpectedN, "Validate expectedN against sample schedule");
        FIELD_MAKE_OPTIONAL(validateExpectedN);
    }

    static IObject* defaultVanillaVarSwap(){
        return new VanillaVarSwap();
    }
};

VanVSModel::VanVSModel(IModelSP varSwapModel):CModel(TYPE), varSwapModel(varSwapModel), capModel(0) {
    static const string method = "VanVSModel::VanVSModel";
    try {
        validatePop2Object();
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

void VanVSModel::validatePop2Object() {
    static const string method = "VanVSModel::validatePop2Object";
    try {
        if (!(NumericalIntegrationLN::TYPE->isInstance(varSwapModel) ||
                ImpliedIntegration::TYPE->isInstance(varSwapModel)  ||
                ClosedFormIntegrateLN::TYPE->isInstance(varSwapModel) ||
                FourierEngine::TYPE->isInstance(varSwapModel) ||
                MonteCarlo::TYPE->isInstance(varSwapModel))) {
            throw ModelException(method, "varSwapModel must be of type NumericalIntegrationLN, ImpliedIntegration, "
                                              "ClosedFormIntegrateLN, FourierEngine or MonteCarlo");
        }
        if(capModel.get()&& !methodology.empty()){
            if(CString::equalsIgnoreCase(methodology, VarCapCVModel::CONTROL_VARIATE) ||
                CString::equalsIgnoreCase(methodology, VarCapCVModel::MULTIPLICATIVE_SCALING) ||
                CString::equalsIgnoreCase(methodology, VarCapCVModel::ADJUST_MEAN_VOL)){
				useCV = true;
			} else if(CString::equalsIgnoreCase(methodology, VarCapCVModel::DEFAULT)){
				methodology = VarCapCVModel::ADJUST_MEAN_VOL;
				useCV = true;
			} else if(CString::equalsIgnoreCase(methodology, VanVSModel::METHOD_NONE)){
				// do nothing
			} else {
				throw ModelException(method, "methodology must be none, default, control_variate, "
												"multiplicative_scaling or adjust_mean_vol.");
			}
        }
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
};

// inherited from CModel
void VanVSModel::Price(CInstrument*  instrument,
                       CControl*     control,
                       CResults*     results)
{
    static const string method = "VanVSModel::Price";
    IProduct* product = 0;
    try {
        if (!IIntoProduct::TYPE->isInstance(instrument)){
            throw ModelException("Instrument of type "+
                                 instrument->getClass()->getName() +
                                 " does not support VanVSModel::IntoProduct");
        }
        if (instrument->priceDeadInstrument(control, results)) {
            return; // done for a dead instrument
        }
        // cast to VanVSModel::IIntoProduct
        IIntoProduct& intoProd = dynamic_cast<IIntoProduct&>(*instrument);

        // create and price the product
        product = intoProd.createProduct(this);
        product->price(this, control, results);
        delete product;
    }
    catch (exception& e) {
        delete product;
        throw ModelException(e, method);
    }
};


void VanVSModel::setDomesticYCName (string discountYieldCurveName) const {
    // Call parent
    CModel::setDomesticYCName(discountYieldCurveName);

    // Call method on sub-models
    if(varSwapModel.get()) {
        varSwapModel->setDomesticYCName(discountYieldCurveName);
    }
    if(capModel.get()) {
        capModel->setDomesticYCName(discountYieldCurveName);
    }
}


//needed?
MarketObjectSP VanVSModel::GetMarket(const MarketData*    market,
                                     const string&        name,
                                     const CClassConstSP& type) const {
    static const string method = "VanVSModel::GetMarket";
    try {
        return varSwapModel->GetMarket(market, name, type);
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

void VanVSModel::getMarket(const MarketData*  market,
                           IInstrumentCollectionSP instruments) {
    static const string method = "VanVSModel::getMarket";
    try {
        varSwapModel->getMarket(market, instruments);
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

IModel::WantsRiskMapping VanVSModel::wantsRiskMapping() const {
    WantsRiskMapping s = varSwapModel->wantsRiskMapping(),
                     c = !capModel ? riskMappingIrrelevant :
                                     capModel->wantsRiskMapping();

    return (s == riskMappingDisallowed || c == riskMappingDisallowed) ?
                riskMappingDisallowed :
           (s == riskMappingAllowed || c == riskMappingAllowed) ?
                riskMappingAllowed :
           riskMappingIrrelevant;
}

// This is the method responsible for pricing the VarSwap and the cap and aggregating
// the results
void VanillaVSProd::price(VanVSModel*       model,
                          Control*          control,
                          CResults*         results) const {
    static const string method = "VanVSModel::price";
    try {

        if (inst->isCapped && !(model->capModel)) {
            throw ModelException(method, "Variance swap is capped but model for "
                                         "pricing cap has not been provided.");
        }

        VarianceSwapSP vswap;

        if (!inst->capOnly) {
            // Price the var swap

            vswap = inst->convert(false);

            model->varSwapModel->Price(vswap.get(), control, results);

            if (!inst->isCapped) {
                // Var swap only
                results->storePrice(results->retrievePrice(), results->getCcyName());
            }

            if (control && control->isPricing()) {

			    // Store price of var swap as an output request
			    OutputRequest* request = control->requestsOutput(OutputRequest::VARSWAP_FV);
			    if (request) {
			        results->storeRequestResult(request, results->retrievePrice());
                }
            }
        }

        if (inst->isCapped) {
            // Price the cap

            vswap = inst->convert(true);

            // Control
            CControlSP ctrl(copy(control));
            CControl* ctrlToUse = 0;
            if (inst->capOnly) {
                // Use original control
                ctrlToUse = control;
            } else {
                // Copy original control & rest it
                ctrlToUse = ctrl.get();
                if(ctrlToUse) {
                    ctrlToUse->reset();  // want this as if nothing has happened
                }
            }

            // Results
            CResultsSP capResults(new Results);
            CResults*  resultsToUse = capResults.get();
            if (inst->capOnly) {
                // Price cap only, override resultsToUse with main results object
                resultsToUse = results;
            }

			// Price using appropriate results and control
			if (model->useCV) {
				model->capCVModel = VarCapCVModelSP(new VarCapCVModel(model->capModel,
															model->varSwapModel,
															model->methodology));
				model->capCVModel->Price(vswap.get(), ctrlToUse, resultsToUse);
			}else{
				model->capModel->Price(vswap.get(), ctrlToUse, resultsToUse);
			}

            // Store price of var cap as an output request
            if (control && control->isPricing()) {
                OutputRequest* request = control->requestsOutput(OutputRequest::VARCAP_FV);
                if (request) {
                    results->storeRequestResult(request, resultsToUse->retrievePrice());
                }
            }

            if (!inst->capOnly) {
                // Aggregate the cap and vswap results -
                // cap is defined as long call on variance so need to subtract results
#if 0
                results->add(capResults.get(),
                             CControlSP(control),
                             -1.0,      // scaleFactor
                             true);     // sameInstrument
#endif

                // accept the fact that the cap doesn't add anything worthwhile to
                // results other than price, to avoid killing ourselves over
                // aggregation of results
                double value = results->retrievePrice() - capResults->retrievePrice();
                results->storePrice(value, results->getCcyName());

                // Aggregate cashflows
                try {
                    if (control && control->isPricing()) {
                        OutputRequest* request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
                        if(request) {
                            results->add(OutputRequest::KNOWN_CASHFLOWS, capResults.get(), -1.0);
                        }
                    }
                } catch(exception&) {
                    // Ignore errors and forget about aggregration
                }
            }
        }
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
};

// Type registration
CClassConstSP const VanVSModel::TYPE =
    CClass::registerClassLoadMethod("VanVSModel", typeid(VanVSModel),
    VanVSModel::load);

CClassConstSP const VanillaVarSwap::TYPE = CClass::registerClassLoadMethod(
    "VanillaVarSwap", typeid(VanillaVarSwap), VanillaVSHelper::load);

CClassConstSP const VanVSModel::IIntoProduct::TYPE =
CClass::registerInterfaceLoadMethod("VanVSModel::IIntoProduct",
                                    typeid(VanVSModel::IIntoProduct), 0);

// IMS compatible VanVSModel shell for ClosedFormIntegrateLN (varswap)
// and FourierEngine (cap)
class VanVSModelCFIN: public VanVSModel {
public:
    static CClassConstSP const TYPE;

private:
    VanVSModelCFIN():VanVSModel(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(VanVSModelCFIN, clazz);
        SUPERCLASS(VanVSModel);
        EMPTY_SHELL_METHOD(defaultVanVSModelCFIN);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultVanVSModelCFIN(){
        return new VanVSModelCFIN();
    }
};

CClassConstSP const VanVSModelCFIN::TYPE =
CClass::registerClassLoadMethod(
    "VanVSModelCFIN", typeid(VanVSModelCFIN), load);

// IMS compatible VanVSModel shell for ImpliedIntegration (varswap)
// and FourierEngine (cap)
class VanVSModelIMIG: public VanVSModel {
public:
    static CClassConstSP const TYPE;

private:
    VanVSModelIMIG():VanVSModel(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(VanVSModelIMIG, clazz);
        SUPERCLASS(VanVSModel);
        EMPTY_SHELL_METHOD(defaultVanVSModelIMIG);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultVanVSModelIMIG(){
        return new VanVSModelIMIG();
    }
};

CClassConstSP const VanVSModelIMIG::TYPE =
CClass::registerClassLoadMethod(
    "VanVSModelIMIG", typeid(VanVSModelIMIG), load);

// IMS compatible VanVSModel shell for ClosedFormIntegrateLN (varswap)
// and FourierEngineDefault (cap)
class BSSwapSVCap: public VanVSModel {
public:
    static CClassConstSP const TYPE;

private:
    BSSwapSVCap():VanVSModel(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(BSSwapSVCap, clazz);
        SUPERCLASS(VanVSModel);
        EMPTY_SHELL_METHOD(defaultBSSwapSVCap);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultBSSwapSVCap(){
        return new BSSwapSVCap();
    }
};

CClassConstSP const BSSwapSVCap::TYPE =
CClass::registerClassLoadMethod(
    "BSSwapSVCap", typeid(BSSwapSVCap), load);

// IMS compatible VanVSModel shell for FourierEngineDefault (varswap)
// and FourierEngineDefault (cap)
class SVSwapSVCap: public VanVSModel {
public:
    static CClassConstSP const TYPE;

private:
    SVSwapSVCap():VanVSModel(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(SVSwapSVCap, clazz);
        SUPERCLASS(VanVSModel);
        EMPTY_SHELL_METHOD(defaultSVSwapSVCap);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultSVSwapSVCap(){
        return new SVSwapSVCap();
    }
};

CClassConstSP const SVSwapSVCap::TYPE =
CClass::registerClassLoadMethod(
    "SVSwapSVCap", typeid(SVSwapSVCap), load);

// IMS compatible VanVSModel shell for ClosedFormIntegrateLN (varswap)
// and FourierEngineVSCurve (cap)
class BSSwapVSCurveCap: public VanVSModel {
public:
    static CClassConstSP const TYPE;

private:
    BSSwapVSCurveCap():VanVSModel(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(BSSwapVSCurveCap, clazz);
        SUPERCLASS(VanVSModel);
        EMPTY_SHELL_METHOD(defaultBSSwapVSCurveCap);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultBSSwapVSCurveCap(){
        return new BSSwapVSCurveCap();
    }
};

CClassConstSP const BSSwapVSCurveCap::TYPE =
CClass::registerClassLoadMethod(
    "BSSwapVSCurveCap", typeid(BSSwapVSCurveCap), load);

//-----------------------------------------------------------------------
// Converts a VSW (VanillaVarSwap) into a VAX (VarianceSwap) for use by the credit offline
class VSWConvert: public CObject{
public:
    static CClassConstSP const TYPE;

    VanillaVarSwapSP  vsw;          // VanillaVarSwap
    CMarketDataSP     market;       // AssetHistory needed for sample schedule generation etc.
    IModelSP          model;        // Model needed for GetMarket()

    /** set an object in a data dictionary */
    static IObjectSP convert(VSWConvert* params){
        return params->createVAX();
    }

    IObjectSP createVAX() const {
        static const string method = "VSWConvert::createVAX";
        try {
            vsw->GetMarket(model.get(), market);

            vsw->Validate();

            if (vsw->isCapped && !vsw->capOnly) {
                throw ModelException(method, "Passed inst cannot be a capped variance swap (i.e. have isCapped=true and capOnly=false)");
            }

            return IObjectSP(vsw->convert(vsw->isCapped));
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    VSWConvert(): CObject(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(VSWConvert, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCreate);
        FIELD(vsw,   "VanillaVarSwap instrument");
        FIELD(market,   "Market data");
        FIELD(model, "Model");
        Addin::registerClassObjectMethod(
            "VSW_CONVERT",
            Addin::UTILITIES,
            "Convert a VanillaVarSwap (VSW) into a VarianceSwap (VAX)",
            TYPE,
            false,
            Addin::returnHandle,
            (Addin::ObjMethod*)convert);
    }

    static IObject* defaultCreate(){
        return new VSWConvert();
    }
};

CClassConstSP const VSWConvert::TYPE =
CClass::registerClassLoadMethod("VSWConvert", typeid(VSWConvert), load);

/* VarSwapCallsPutsAddin */
class VarSwapCallsPutsAddin: public CObject{
public:
    static CClassConstSP const TYPE;

    MarketDataSP                market;
    CAssetWrapper               asset;
    DateTime                    valueDate;
    DateTime                    maturity;
    DoubleArray                 relativeStrikes;
    string                      volType; // optional, default = "VolPreferred"
    bool                        exact; // optional, default = false (true = IntFuncInf on interval relativeStrikes)

private:
	ObjectArraySP varSwapPutsCalls() {
        static const string method = "VarSwapCallsPutsAddin::varSwapPutsCalls";
        try {
            MarketDataFetcherSP mdf(new MarketDataFetcherLN(volType));
            NonPricingModel model(mdf);
            asset.getData(&model, market);

            int nbStrikes = relativeStrikes.size();

            refCountPtr<Range> interval;
            interval = refCountPtr<Range>(new Range(ClosedBoundary(relativeStrikes[0]), ClosedBoundary(relativeStrikes[nbStrikes-1])));
            double fwd = asset->fwdValue(maturity);

            // LinearStrikeVolRequestSP volRequest(new LinearStrikeVolRequest(0.0,valueDate,maturity,false));
            // VarianceSwap::VarSwapIntegrand myIntegrand(*interval, fwd, valueDate, maturity, asset, volRequest);
            VanillaContractsRecorderSP recorder(new VanillaContractsRecorder(1.0, false));
            VarianceSwap::VarSwapIntegrandWeightConstSP optionWeights(new 
                VarianceSwap::VarSwapIntegrandWeight(*interval));
            StaticReplicationIntegrandSP integrand(new StaticReplicationIntegrand(
                optionWeights, recorder, asset.getSP(), valueDate, maturity, false));
            const Function1DDouble& myIntegrand = *integrand;

            ObjectArraySP result(new ObjectArray(2));
            // compute operator(relativeStrikes)
            DoubleArraySP operatorValues(new DoubleArray(nbStrikes));
            for (int i=0; i<nbStrikes; i++) {
                (*operatorValues)[i] = myIntegrand(relativeStrikes[i]);
            }
            // compute integral for range given by relativeStrikes
            Integrator1DSP myIntegrator = Integrator1D::createIntegrator("IntFuncGeneral", DBL_EPSILON, DBL_EPSILON);
            // wrap a doulbe into an object (alternative would be to create doublearray of length one)
            IObjectSP integral(CDouble::create(myIntegrator->integrate(myIntegrand)));
            // combine results and return them
            (*result)[0] = operatorValues;
            (*result)[1] = integral;
            return result;

        } catch (exception& e){
            throw ModelException(e, method);
        }
    }

    static IObjectSP callVarSwapPutsCalls(VarSwapCallsPutsAddin* params) {
        return params->varSwapPutsCalls();
    }

    /** for reflection */
    VarSwapCallsPutsAddin():
    CObject(TYPE),
    volType("VolPreferred"),
    exact(false) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(VarSwapCallsPutsAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultVarSwapCallsPutsAddin);
        FIELD(market, "");
        FIELD(asset, "");
        FIELD(valueDate, "");
        FIELD(maturity, "");
        FIELD(relativeStrikes, "");
        FIELD(volType, "");
        FIELD_MAKE_OPTIONAL(volType);
        FIELD(exact, "");
        FIELD_MAKE_OPTIONAL(exact);

        Addin::registerClassObjectMethod("VOLVAR_SWAP_CALLSPUTS",
                                          Addin::MARKET,
                                          "returns value of operator at relative strike",
                                          TYPE,
                                          false,
                                          Addin::returnHandle,
                                          (Addin::ObjMethod*)callVarSwapPutsCalls);
    }

    static IObject* defaultVarSwapCallsPutsAddin(){
        return new VarSwapCallsPutsAddin();
    }
};

CClassConstSP const VarSwapCallsPutsAddin::TYPE =
CClass::registerClassLoadMethod("VarSwapCallsPutsAddin", typeid(VarSwapCallsPutsAddin), load);

/* GenerateSamplesAddin */
class GenerateSamplesAddin: public CObject{
public:
    static CClassConstSP const TYPE;

    DateTime        startDate;
    DateTime        endDate;
    HolidaySP       holidays;
    string          frequency;
    string          badDayConvention;
    string          timeRule;

private:
	IObjectSP generateSamples() {
        static const string method = "GenerateSamplesAddin::generateSamples";
        try {
            StubSP stub(StubFactory::make("None"));
            BadDayConventionSP badDay(BadDayConventionFactory::make(badDayConvention));
            DayCountConventionSP dayCount(DayCountConventionFactory::make("Act/360"));

            MaturityPeriod periodRef(frequency);
            int frequencyCount;
            string frequencyInterval;
            periodRef.decompose(frequencyCount, frequencyInterval);

            int time = DateTime::timeConvert(timeRule);
            DateTime firstSample(startDate.getDate(), time);
            DateTime lastSample(endDate.getDate(), time);

            CashFlowArray samples = SwapTool::cashflows(firstSample,            // startDate
                                                        lastSample,             // maturity
                                                        stub.get(),             // stub
                                                        true,                   // stubAtEnd
                                                        badDay.get(),           // accrualBadDayConv
                                                        badDay.get(),           // payBadDayConv
                                                        holidays.get(),         // holidays
                                                        true,                   // keepStartDate (false => startDate is removed if amount is zero)
                                                        false,                  // addPrincipal (true => principal (1.0) is added to amount at endDate)
                                                        0.0,                    // rate (interest rate)
                                                        frequencyCount,         // count
                                                        frequencyInterval,      // period
                                                        dayCount.get());        // dayCountConv (for dicounting of amounts)
            return IObjectSP(samples.clone());

        } catch (exception& e){
            throw ModelException(e, method);
        }
    }

    static IObjectSP callGenerateSamples(GenerateSamplesAddin* params) {
        return params->generateSamples();
    }

    /** for reflection */
    GenerateSamplesAddin():
    CObject(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(GenerateSamplesAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultGenerateSamplesAddin);
        FIELD(startDate, "start date");
        FIELD(endDate, "end date");
        FIELD(holidays, "holidays");
        FIELD(frequency, "sample frequency (as string)");
        FIELD(badDayConvention, "bad day convention (as string)");
        FIELD(timeRule, "time rule for sample schedule");

        Addin::registerClassObjectMethod("VAR_SWAP_GENERATE_SAMPLES",
                                          Addin::MARKET,
                                          "returns sample schedule for variance swap",
                                          TYPE,
                                          false,
                                          Addin::returnHandle,
                                          (Addin::ObjMethod*)callGenerateSamples);
    }

    static IObject* defaultGenerateSamplesAddin(){
        return new GenerateSamplesAddin();
    }
};

CClassConstSP const GenerateSamplesAddin::TYPE =
CClass::registerClassLoadMethod("GenerateSamplesAddin", typeid(GenerateSamplesAddin), load);


/////////////////////////////////////////////////////////////////////////////
// Control Variates  /////////////
/////////////////////////////////////////////////////////////////////////////

/**  VarianceSwapLeastSquareFit objFunc class */
class VarianceSwapLeastSquareFit: public Calibrator::ObjFuncLeastSquare{
public:
    static CClassConstSP const TYPE;

    void getMarket(MarketData* market){}

    IObjectSP getAdjustableGroup(){
        return IObjectSP::attachToRef(this);
    }

    int getNbFuncs() const{
        return 1;
    }

    void calcValue(CDoubleArray& funcvals) const {
        static const string method = "VarianceSwapLeastSquareFit::calcValue";
        try{
            try{
                if(!converted){
                    convert();
                }
                volModel->Price(vswap.get(), control.get(), results.get());
            }
            catch(exception& e){
                throw ModelException(method,"Failed to compute model price of future variance");
            }
            funcvals[0] = results->retrievePrice() - marketSwap;
        }
        catch(exception& e){
            throw ModelException(e, method);
        }
    }

    // price model price of the original var cap
    double priceVarCap(CControl* extControl,
                       CResults* extResults) const {
        static const string method = "VarianceSwapLeastSquareFit::priceVarCap";
        try{
            if(converted){
                vswap->strikeVol = strikeVol;
                vswap->dontScaleByStrike = dontScaleByStrike;
                converted = false;
            }
            volModel->Price(vswap.get(), extControl, extResults);
        }
        catch(exception& e){
            throw ModelException(method,"Failed to compute model price of future variance");
        }
        return extResults->retrievePrice();
    }

    // for reflection
    VarianceSwapLeastSquareFit():
    Calibrator::ObjFuncLeastSquare(TYPE){}

    VarianceSwapLeastSquareFit(IModelSP            volModel,
                               IModelSP            fwdModel,
                               VarianceSwapSP      vswap,
                               CAssetSP            marketAsset):
    Calibrator::ObjFuncLeastSquare(TYPE),
    volModel(volModel),
    fwdModel(fwdModel),
    vswap(vswap),
    marketAsset(marketAsset){
        validatePop2Object();
    }

private:
    virtual void validatePop2Object(){
        static const string method = "VarianceSwapLeastSquareFit::validatePop2Object";
        try {
            //temp control and results
            control = CControlSP(new CControl());
            results = CResultsSP(new CResults());

            if(!FourierEngine::TYPE->isInstance(volModel)){
                throw ModelException(method, "volModel must be Fourier Engine");
            }

            if(!CString::equalsIgnoreCase(vswap->payoffType, "FORWARD_CAP")){
                throw ModelException(method, "the instrument must be a variance cap");
            }

            convert();

            // calc market price of var cap
            string      payoffType = vswap->payoffType;
            CAssetSP    tmpAsset(vswap->asset.getSP());

            vswap->payoffType = "FORWARD";
            (vswap->asset).setObject(CAssetSP(marketAsset));

            try {
                fwdModel->Price(vswap.get(), control.get(), results.get());
            }
            catch (exception& e){
                string message = "Failed to compute market price of future variance";
                throw ModelException::addTextToException(e, message);
            }
            marketSwap = results->retrievePrice();

            vswap->payoffType = payoffType;
            (vswap->asset).setObject(tmpAsset);

        } catch (exception& e){
            throw ModelException(e, method);
        }
    }

    // convert into var swap
    void convert() const {
        strikeVol = vswap->strikeVol;
        dontScaleByStrike = vswap->dontScaleByStrike;
        vswap->strikeVol = 0.0;
        vswap->dontScaleByStrike = true;
        converted = true;
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(VarianceSwapLeastSquareFit, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultVarianceSwapLeastSquareFit);

		FIELD(volModel, "");
        FIELD(fwdModel, "");
        FIELD(vswap, "");
        FIELD(marketAsset, "");

        // transients
        FIELD(control, "");
        FIELD_MAKE_TRANSIENT(control);
        FIELD(results, "");
        FIELD_MAKE_TRANSIENT(results);
        FIELD(strikeVol, "");
        FIELD_MAKE_TRANSIENT(strikeVol);
        FIELD(dontScaleByStrike, "");
        FIELD_MAKE_TRANSIENT(dontScaleByStrike);
        FIELD(marketSwap, "");
        FIELD_MAKE_TRANSIENT(marketSwap);
        FIELD(converted, "");
        FIELD_MAKE_TRANSIENT(converted);
    }

    static IObject* defaultVarianceSwapLeastSquareFit(){
        return new VarianceSwapLeastSquareFit();
    }

    //registered fields
    IModelSP            volModel;       // FourierEngine
    IModelSP            fwdModel;       // ClosedFormIntegrateLN
    VarianceSwapSP      vswap;          // Variance cap
    CAssetSP            marketAsset;    // Asset with VolSurface

    // transcient fields
    CControlSP          control;
    CResultsSP          results;
    mutable double      strikeVol;
    mutable bool        dontScaleByStrike;
    mutable double      marketSwap;
    mutable bool        converted;      // true -> strikeVol = 0, dontScaleByStrike = true
};
typedef smartPtr<VarianceSwapLeastSquareFit> VarianceSwapLeastSquareFitSP;

CClassConstSP const VarianceSwapLeastSquareFit::TYPE =
CClass::registerClassLoadMethod("VarianceSwapLeastSquareFit", typeid(VarianceSwapLeastSquareFit), load);


void VarCapCVModel::validatePop2Object(){
    if (CString::equalsIgnoreCase(methodology, DEFAULT)){
        methodology = ADJUST_MEAN_VOL;
    }
}

// registration, invoked when class is 'loaded'
void VarCapCVModel::load(CClassSP& clazz){
    clazz->setPublic();
    REGISTER(VarCapCVModel, clazz);
    SUPERCLASS(CModel);
    EMPTY_SHELL_METHOD(defaultVarCapCVModel);
    FIELD(volModel,"Model to price Variance Cap");
    FIELD(fwdModel,"Model to price Variance Swap");
    FIELD(methodology,"Methodology used to adjust variance cap to be consistent with variance swap");
    FIELD_MAKE_OPTIONAL(methodology);
}

// for VarCapCVModel::IIntoProduct
void VarCapCVModel::loadIntoProduct(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER_INTERFACE(VarCapCVModel::IIntoProduct, clazz);
    EXTENDS(Model::IModelIntoProduct);
}

IObject* VarCapCVModel::defaultVarCapCVModel(){
    return new VarCapCVModel();
}

// constructor
VarCapCVModel::VarCapCVModel():CModel(TYPE), methodology(DEFAULT) {}

VarCapCVModel::VarCapCVModel(IModelSP  volModel,
                             IModelSP  fwdModel,
                             string    methodology):
CModel(TYPE),
volModel(volModel),
fwdModel(fwdModel),
methodology(methodology) {}

void VarCapCVModel::Price(CInstrument* instrument,
                       CControl*    control,
                       CResults*    results) {
    static const string method = "VarCapCVModel::Price";

    try{
        // Get Product from instrument
        if (!IIntoProduct::TYPE->isInstance(instrument)){
            throw ModelException(method, "Instrument of type "+
                                 instrument->getClass()->getName() +
                                 " does not support VarCapCVModel::IntoProduct");
        }

        IIntoProduct& intoProd = dynamic_cast<IIntoProduct&>(*instrument);
        refCountPtr<Product> prod(intoProd.createProduct(this));

        double price = 0.0;
        CResultsSP tmpRes(new Results());
        VarianceSwap* vswap = (prod->varCap).get();

        // valueDate >= matDate is taken care of here
        if(vswap->priceDeadInstrument(control, results)){
            return; // dead instrument priced
        }

        // adjustment ...
        if (CString::equalsIgnoreCase(methodology, CONTROL_VARIATE)){
            // Price components now
            try{
                volModel->Price(vswap, control, results);   // Essentially FourierEngine
            }
            catch(exception& e){
                string message = "Failed to compute model price of variance swap";
                throw ModelException::addTextToException(e, message);
            }

            CAssetSP   tmpAsset(vswap->asset.getSP());
            string type = vswap->payoffType;
            (vswap->asset).setObject(prod->marketAsset);
            vswap->payoffType = "FORWARD";
            try {
                fwdModel->Price(vswap, control, tmpRes.get());
            }
            catch (exception& e){
                string message = "Failed to compute market price of variance swap";
                throw ModelException::addTextToException(e, message);
            }
            (vswap->asset).setObject(tmpAsset);
            vswap->payoffType = type;
            double marketSwap = tmpRes->retrievePrice();

            // Extract information from results
            OutputNameSP    nameVarSwap(new OutputName(VarOptCVFP::VAR_SWAP));
            OutputNameSP    nameVarSwapSqr(new OutputName(VarOptCVFP::VAR_SWAP_SQR));
            OutputNameSP    nameVarCap(new OutputName(VarOptCVFP::VAR_CAP));
            OutputNameSP    nameVarCapSqr(new OutputName(VarOptCVFP::VAR_CAP_SQR));
            OutputNameSP    nameModelSwap(new OutputName(VarOptCVFP::MODEL_SWAP));
            OutputNameSP    nameStrike(new OutputName(VarOptCVFP::STRIKE));

            price                   = results->retrievePrice();
            double varSwap          = results->retrieveScalarGreek(Results::DEBUG_PACKET,nameVarSwap);
            double varSwapSqr       = results->retrieveScalarGreek(Results::DEBUG_PACKET,nameVarSwapSqr);
            double varCap           = results->retrieveScalarGreek(Results::DEBUG_PACKET,nameVarCap);
            double varCapSqr        = results->retrieveScalarGreek(Results::DEBUG_PACKET,nameVarCapSqr);
            double modelSwap        = results->retrieveScalarGreek(Results::DEBUG_PACKET,nameModelSwap);
            double strike           = results->retrieveScalarGreek(Results::DEBUG_PACKET,nameStrike);

            // Combine them to get the price
            double cv_coeff = (varCapSqr-varCap*(varSwap-strike))/(varSwapSqr-varSwap*varSwap);

            // Price and PV
            price += cv_coeff*(marketSwap - modelSwap);

            // See if option price is positive
            if(price <= 0.0) {
                if(price > -0.01) {
                    price = 0.0;
                } else {
                    throw ModelException(method, "Control Variate Adjusted price for Variance Cap is negative.");
                }
            }

            // take care of additional outputs
            if (control && control->isPricing()) {
                OutputRequest* request = control->requestsOutput(OutputRequest::CV_COEFF);
                if (request) {
                    results->storeRequestResult(request, cv_coeff);
                }
             }
        } else if (CString::equalsIgnoreCase(methodology, MULTIPLICATIVE_SCALING)){

            if(Maths::isZero(vswap->strikeVol)){
                //effectively variance swap
                CAssetSP   tmpAsset(vswap->asset.getSP());
                string type = vswap->payoffType;
                (vswap->asset).setObject(prod->marketAsset);
                vswap->payoffType = "FORWARD";
                try {
                    fwdModel->Price(vswap, control, tmpRes.get());
                }
                catch (exception& e){
                    string message = "Failed to compute market price of variance swap";
                    throw ModelException::addTextToException(e, message);
                }
                (vswap->asset).setObject(tmpAsset);
                vswap->payoffType = type;
                price = tmpRes->retrievePrice();

                //simply to have requests reported
                try{
                    volModel->Price(vswap, control, results);
                }
                catch(exception& e){
                    string message = "Failed to compute model price of future variance";
                    throw ModelException::addTextToException(e, message);
                }
            } else {
                int     numReturns;
                double  volPast    = vswap->historicalVol(numReturns);
                double  pastWeight = (double)numReturns/(double)(vswap->numTotalReturns);

                // model price and market price of future variance
                double strikeVol = vswap->strikeVol;
                double cap = vswap->cap;
                bool dontScaleByStrike = vswap->dontScaleByStrike;
                vswap->strikeVol = volPast * sqrt(pastWeight);
                vswap->cap = 1.0;
                vswap->dontScaleByStrike = true;
                try{
                    volModel->Price(vswap, control, results);
                }
                catch(exception& e){
                    string message = "Failed to compute model price of future variance";
                    throw ModelException::addTextToException(e, message);
                }
                double modelSwap = results->retrievePrice();
                vswap->cap = cap;

                CAssetSP   tmpAsset(vswap->asset.getSP());
                string type = vswap->payoffType;
                (vswap->asset).setObject(prod->marketAsset);
                vswap->payoffType = "FORWARD";
                try {
                    fwdModel->Price(vswap, control, tmpRes.get());
                }
                catch (exception& e){
                    string message = "Failed to compute market price of future variance";
                    throw ModelException::addTextToException(e, message);
                }
                double marketSwap = tmpRes->retrievePrice();
                (vswap->asset).setObject(tmpAsset);
                vswap->payoffType = type;
                vswap->strikeVol = strikeVol;
                vswap->dontScaleByStrike = dontScaleByStrike;

                // adjust cap
                double ratio = marketSwap / modelSwap;
                vswap->cap = cap / ratio + volPast*volPast*pastWeight/strikeVol/strikeVol*(1.0-1.0/ratio);
                try{
                    volModel->Price(vswap, control, results);
                }
                catch(exception& e){
                    string message = "Failed to price variance cap";
                    throw ModelException::addTextToException(e, message);
                }
                price = results->retrievePrice()*ratio;
                vswap->cap = cap;
            }
        } else if (CString::equalsIgnoreCase(methodology, ADJUST_MEAN_VOL)) {
			// currency protection is not supported
			if(CString::equalsIgnoreCase(vswap->ccyTreatment, CAsset::CCY_TREATMENT_PROTECTED)||
				CString::equalsIgnoreCase(vswap->ccyTreatment, CAsset::CCY_TREATMENT_STRUCK)) {
                throw ModelException(method, "Currency treatment P or S is not supported by adjust_mean_vol method.");
            }

            // calibrator
            Calibrator          calibrator(LevenbergMarquardtSP(new LevenbergMarquardt()));

            // objFunc
            VarianceSwapLeastSquareFit      objFunc(volModel,
                                                    fwdModel,
                                                    VarianceSwapSP(copy(vswap)),
                                                    prod->marketAsset);

            // ids
            string  name;
            EquityBase*     equityAsset = dynamic_cast<EquityBase*>(prod->marketAsset.get());
            if(equityAsset){
                name = equityAsset->getVolName();
            } else {
                name = prod->marketAsset->getName();
            }
            string volType;
            if(FourierEngine::TYPE->isInstance(volModel.get())){
                FourierEngine*   tmpModel = dynamic_cast<FourierEngine*>(volModel.get());
                if(FourierProcessSV::TYPE->isInstance(tmpModel->getProcess())){
                    volType = "VolSV";
                } else if (FourierProcessSVJ::TYPE->isInstance(tmpModel->getProcess())){
                    volType = "VolSVJ";
                } else if (FourierProcessSVCJ::TYPE->isInstance(tmpModel->getProcess())){
                    volType = "VolSVCJ";
                } else {
                    throw ModelException(method,
                        "the fourier process must be FourierProcessSV,FourierProcessSVJ or FourierProcessSVCJ");
                }
            } else {
                throw ModelException(method, "the model must be FourierEngine");
            }

            Calibrator::InstanceIDArray ids;
            ids.reserve(1);
            ids.push_back(Calibrator::InstanceIDSP(
                new Calibrator::InstanceIDDb(volType,
                                             name,
                                             "meanVol",
                                             false,
                                             0.2,
                                             0)));

		    // run the calibrator
		    calibrator.run(objFunc, ids);

            price = objFunc.priceVarCap(control, results);

        } else {
            throw ModelException(method, "Methodology must be control_variate, multiplicative_scaling or adjust_mean_vol");
        }

        // Store results
        results->storePrice(price, prod->varCap->discount->getCcy());

        if (CString::equalsIgnoreCase(methodology, CONTROL_VARIATE)||CString::equalsIgnoreCase(methodology,MULTIPLICATIVE_SCALING)){
            //report "market" VOL_IN_FUTURE and TOTAL_VOL
            if (control && control->isPricing()) {
                 OutputRequest*  request = control->requestsOutput(OutputRequest::VOL_IN_FUTURE);
                 if (request) {
                     const CDouble* requestResult = dynamic_cast<const CDouble*>(tmpRes->retrieveRequestResult("VOL_IN_FUTURE").get());
                     if(requestResult) {
                         results->storeRequestResult(request, requestResult->doubleValue());
                     }
                 }
                 request = control->requestsOutput(OutputRequest::TOTAL_VOL);
                 if (request) {
                     const CDouble* requestResult = dynamic_cast<const CDouble*>(tmpRes->retrieveRequestResult("TOTAL_VOL").get());
                     if(requestResult) {
                         results->storeRequestResult(request, requestResult->doubleValue());
                     }
                 }
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

IModel::WantsRiskMapping VarCapCVModel::wantsRiskMapping() const {
    return (volModel->wantsRiskMapping() == riskMappingDisallowed ||
            fwdModel->wantsRiskMapping() == riskMappingDisallowed) ?
                riskMappingDisallowed :
           (volModel->wantsRiskMapping() == riskMappingAllowed ||
            fwdModel->wantsRiskMapping() == riskMappingAllowed) ?
                riskMappingAllowed :
           riskMappingIrrelevant;
}

// convenience wrapper for variance caps
class VarianceCap: public VarianceSwap {

public:
    static CClassConstSP const TYPE;

private:
    VarianceCap() :VarianceSwap(TYPE) {
        payoffType = "FORWARD_CAP";
    }

    static void load(CClassSP& clazz) {
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VarianceCap, clazz);
        SUPERCLASS(VarianceSwap);
        EMPTY_SHELL_METHOD(defaultVarianceCap);
    }
    static IObject* defaultVarianceCap() {
        return new VarianceCap();
    }
};

CClassConstSP const VarianceCap::TYPE = CClass::registerClassLoadMethod(
    "VarianceCap", typeid(VarianceCap), VarianceCap::load);

//////////////////////////////////////////////////////////////////////////////////////////////


/** Default VanillaVarSwap Maker */
class DefaultParVanillaVarSwapMaker: public CObject,
                                     virtual public ClientRunnable{
public:
    static CClassConstSP const TYPE;

    // EdrAction version of addin
    IObjectSP run() {
        return IObjectSP(makeParSwaps());
    }

    /** Creates par swap instruments for input maturities */
    VanillaVarSwapArraySP makeParSwaps() {
        int nbInsts = maturities.size();
        VanillaVarSwapArraySP insts(new VanillaVarSwapArray(nbInsts));
        for(int iInst = 0; iInst < nbInsts; iInst++) {
            (*insts)[iInst] = VanillaVarSwap::makeParSwap(
                asset,
                discount,
                firstDate,
                maturities[iInst],
                dividendAdjusted,
                settlement,
                assetHols,
                market,
                model,
                isOption,
                cap,
                scenario);
        }

        return insts;
    }

    /** Full constructor */
    DefaultParVanillaVarSwapMaker(const CAssetWrapper&     asset,
                                  const YieldCurveWrapper& discount,
                                  const DateTimeArray&     maturities,
                                  bool                     dividendAdjusted,
                                  const HolidayWrapper&    assetHols,
                                  CMarketDataSP            market,
                                  bool                     isOption,
                                  double                   cap):
    CObject(TYPE), asset(asset), discount(discount), maturities(maturities), dividendAdjusted(dividendAdjusted),
    assetHols(assetHols), market(market), isOption(isOption), cap(cap) {
        validatePop2Object();
    }

private:
    CAssetWrapper           asset;
    YieldCurveWrapper       discount;
    DateTime                firstDate;
    DateTimeArray           maturities;
    bool                    dividendAdjusted;
    InstrumentSettlementSP  settlement;
    HolidayWrapper          assetHols;
    CMarketDataSP           market;
    IModelSP                model;
    bool                    isOption;
    double                  cap;
    ScenarioSP              scenario;

    /** for reflection */
    DefaultParVanillaVarSwapMaker(): CObject(TYPE), isOption(false), cap(0.0) {
        // Assume a swap
        IModelSP swapModel(new ClosedFormIntegrateLN("VolPreferred"));
        model = VanVSModelSP(new VanVSModel(swapModel));
    }


    void validatePop2Object() {
        // VarSwap starting tonight
        if(firstDate.empty()) {
            DateTime valueDate = market->GetReferenceDate();
            firstDate = DateTime(valueDate.getDate(), maturities[0].getTime());
        }

        if( !settlement) {
            settlement = InstrumentSettlementSP(new CashSettlePeriod(0));
        }
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic();
        clazz->setDescription("Create default par variance swaps");
        REGISTER(DefaultParVanillaVarSwapMaker, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultDefaultParVanillaVarSwapMaker);
        IMPLEMENTS(ClientRunnable);
        FIELD(asset, "Asset");
        FIELD(discount, "Discount Curve");
        FIELD(firstDate, "First sample");
        FIELD_MAKE_OPTIONAL(firstDate);
        FIELD(maturities, "Maturities");
        FIELD(dividendAdjusted, "Dividend Adjusted");
        FIELD(settlement, "Settlement");
        FIELD_MAKE_OPTIONAL(settlement);
        FIELD(assetHols, "Asset holidays");
        FIELD(market, "Market data");
        // Optional
        FIELD(model, "Model");
        FIELD_MAKE_OPTIONAL(model);
        FIELD(isOption, "Variance option");
        FIELD_MAKE_OPTIONAL(isOption);
        FIELD(cap, "Cap Value");
        FIELD_MAKE_OPTIONAL(cap);
        FIELD(scenario, "Scenarios");
        FIELD_MAKE_OPTIONAL(scenario);


        Addin::registerObjectMethod("DEFAULT_PAR_VAN_VARSWAP",
                                    Addin::RISK,
                                    "Creates a default par vanilla var swap",
                                    true,
                                    Addin::returnHandle,
                                    &DefaultParVanillaVarSwapMaker::makeParSwaps);
    }

    static IObject* defaultDefaultParVanillaVarSwapMaker(){
        return new DefaultParVanillaVarSwapMaker();
    }
};

CClassConstSP const DefaultParVanillaVarSwapMaker::TYPE =
CClass::registerClassLoadMethod("DefaultParVanillaVarSwapMaker", typeid(DefaultParVanillaVarSwapMaker), load);


//////////////////////////////////////////////////////////////////////////////////////////////


/** VarSwap Addin */
class VarSwapAddin: public CObject{
public:
    static CClassConstSP const TYPE;

    /** Computes par var swap vol for input maturities */
    DoubleArraySP getStrikes() {
        int nbInsts = insts->size();
        DoubleArraySP vols(new DoubleArray(nbInsts));
        for(int iInst = 0; iInst < nbInsts; iInst++) {
            (*vols)[iInst] = (*insts)[iInst]->strike;
        }
        return vols;
    }

    VarSwapAddin(VanillaVarSwapArraySP insts): CObject(TYPE), insts(insts) {}

private:
    VanillaVarSwapArraySP insts;

    /** for reflection */
    VarSwapAddin(): CObject(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(VarSwapAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultVarSwapAddin);
        FIELD(insts, "Variance Swaps");

        Addin::registerObjectMethod("VARSWAP_STRIKE",
                                    Addin::RISK,
                                    "Computes the par vol of the vanilla variance swap",
                                    true,
                                    Addin::expandSimple,
                                    &VarSwapAddin::getStrikes);
    }

    static IObject* defaultVarSwapAddin(){
        return new VarSwapAddin();
    }
};

CClassConstSP const VarSwapAddin::TYPE =
CClass::registerClassLoadMethod("VarSwapAddin", typeid(VarSwapAddin), load);


//////////////////////////////////////////////////////////////////////////////////////////////


/* Stochastic Rates Adjustment */
class VSwapStochRatesAdjustment: public CObject{
public:
    static CClassConstSP const TYPE;

private:
    VanillaVarSwapArraySP insts;
    TimeMetricSP          timeMetric;
    double                rateVolatility;
    double                rateMeanReversion;
    double                equityRateCorrelation;
    DoubleArraySP         equityVol;

    DoubleArraySP stochRatesAdjustment() {
        static const string method = "VSwapStochRatesAdjustment::stochRatesAdjustment";
        try {
            // Compute adjustment per maturity
            int nbDates = insts->size();
            DoubleArraySP output(new DoubleArray(nbDates));
            for(int iDate = 0; iDate < nbDates; iDate++) {
                double volAdjustment = VanillaVarSwap::stochRatesVolAdjustment(
                    rateVolatility,
                    rateMeanReversion,
                    equityRateCorrelation,
                    (*equityVol)[iDate],
                    timeMetric,
                    (*insts)[iDate]);
                (*output)[iDate] = volAdjustment;
            }

            return output;
        } catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** for reflection */
    VSwapStochRatesAdjustment(): CObject(TYPE) {}

    void validatePop2Object() {
        static const string method = "VSwapStochRatesAdjustment::stochRatesAdjustment";
        try {
            if(Maths::isNegative(rateVolatility)) {
                throw ModelException(
                    method, "rate volatility must be non-negative. Value is " + Format::toString(rateVolatility));
            }
            if(Maths::isNegative(rateMeanReversion)) {
                throw ModelException(
                    method, "rate mean reversion must be non-negative. Value is " + Format::toString(rateMeanReversion));
            }
            if(equityRateCorrelation < -1 || equityRateCorrelation > 1) {
                throw ModelException(
                    method, "equity rate correlation must be in [-1, 1]. Value is " + Format::toString(equityRateCorrelation));
            }
            for(int iDate = 0; iDate < equityVol->size(); iDate++) {
                if(Maths::isNegative((*equityVol)[iDate])) {
                    throw ModelException(
                        method, "equity volatility must be non-negative. Value for maturity " +
                        Format::toString(iDate + 1) + " is " + Format::toString((*equityVol)[iDate]));
                }
            }
        } catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(VSwapStochRatesAdjustment, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultVSwapStochRatesAdjustment);

        FIELD(insts, "Variance Swap instruments");
        FIELD(timeMetric, "Time metric");
        FIELD(rateVolatility, "Interest rate volatility");
        FIELD(rateMeanReversion, "Rate mean reversion");
        FIELD(equityRateCorrelation, "Equity - rate correlation");
        FIELD(equityVol, "Equity volatility");

        Addin::registerObjectMethod("VAR_SWAP_STOCHASTIC_RATES_ADJUSTMENT",
                                    Addin::RISK,
                                    "VSwap Adjustment for stochastic rates",
                                    true,
                                    Addin::expandSimple,
                                    &VSwapStochRatesAdjustment::stochRatesAdjustment);
    }

    static IObject* defaultVSwapStochRatesAdjustment(){
        return new VSwapStochRatesAdjustment();
    }
};

CClassConstSP const VSwapStochRatesAdjustment::TYPE =
CClass::registerClassLoadMethod("VSwapStochRatesAdjustment", typeid(VSwapStochRatesAdjustment), load);


//////////////////////////////////////////////////////////////////////////////////////////////


/* Stochastic Rates Adjustment */
class VSwapJumpsAdjustment: public CObject{
public:
    static CClassConstSP const TYPE;

private:
    VanillaVarSwapArraySP insts;
    TimeMetricSP          timeMetric;
    double                jumpPerYear;
    double                logSpotJumpMean;
    double                logSpotJumpStdDev;

    DoubleArraySP jumpsAdjustment() {
        static const string method = "VSwapJumpsAdjustment::jumpsAdjustment";
        try {
            // Compute adjustment per maturity
            int nbDates = insts->size();
            DoubleArraySP output(new DoubleArray(nbDates));
            for(int iDate = 0; iDate < nbDates; iDate++) {
                double volAdjustment = VanillaVarSwap::jumpsVolAdjustment(
                    jumpPerYear,
                    logSpotJumpMean,
                    logSpotJumpStdDev,
                    timeMetric,
                    (*insts)[iDate]);
                (*output)[iDate] = volAdjustment;
            }

            return output;
        } catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** for reflection */
    VSwapJumpsAdjustment(): CObject(TYPE) {}

    void validatePop2Object() {
        static const string method = "VSwapJumpsAdjustment::stochRatesAdjustment";
        try {
            if(Maths::isNegative(jumpPerYear)) {
                throw ModelException(
                    method, "Number of jumps per year must be non-negative. Value is " + Format::toString(jumpPerYear));
            }
            if(Maths::isNegative(logSpotJumpStdDev)) {
                throw ModelException(
                    method, "Std dev of jumps must be non-negative. Value is " + Format::toString(logSpotJumpStdDev));
            }
        } catch (exception& e) {
            throw ModelException(e, method);
        }
    }


    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(VSwapJumpsAdjustment, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultVSwapJumpsAdjustment);

        FIELD(insts, "Variance Swap instruments");
        FIELD(timeMetric, "Time metric");
        FIELD(jumpPerYear, "Jumps per year");
        FIELD(logSpotJumpMean, "Mean of jumps in log spot");
        FIELD(logSpotJumpStdDev, "StdDev of jumps in log spot");

        Addin::registerObjectMethod("VAR_SWAP_JUMPS_ADJUSTMENT",
                                    Addin::RISK,
                                    "VSwap Adjustment for jumps",
                                    true,
                                    Addin::expandSimple,
                                    &VSwapJumpsAdjustment::jumpsAdjustment);
    }

    static IObject* defaultVSwapJumpsAdjustment(){
        return new VSwapJumpsAdjustment();
    }
};

CClassConstSP const VSwapJumpsAdjustment::TYPE =
CClass::registerClassLoadMethod("VSwapJumpsAdjustment", typeid(VSwapJumpsAdjustment), load);


DRLIB_END_NAMESPACE
