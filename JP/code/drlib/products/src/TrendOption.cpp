//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : TrendOption.cpp
//
//   Description : Option on the number of times the return on an asset changes direction
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : November 25, 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/TrendOption.hpp"
#include "edginc/Vanilla.hpp"
#include "edginc/AtMaturity.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"



DRLIB_BEGIN_NAMESPACE

const string TrendOption::SPREAD_NONE = "NONE";
const string TrendOption::SPREAD_UP = "UP";
const string TrendOption::SPREAD_DOWN = "DOWN";
const string TrendOption::SPREAD_SYMMETRIC = "SYMMETRIC";
const double TrendOption::MINIMUM_ALPHA = 0.001;

double binModel(DoubleArraySP odds, DoubleArraySP payoff);
double binModel(DoubleArraySP odds, DoubleArraySP payoff)
{
    int i;
    int j;
    double pu;
    double pd;
    DoubleArraySP value(new DoubleArray());
    
    value->resize(payoff->size());

    for (j=0; j<payoff->size(); j++) {
        (*value)[j] = (*payoff)[j];
    }

    for (i=odds->size()-1; i>=0; i--) {
        pu = (*odds)[i];
        pd = 1 - pu;
        
        for (j=0; j<=i; j++) {
            (*value)[j] = (*value)[j+1] * pu + (*value)[j] * pd;
        }
    }

    return (*value)[0];
}



////////////////////
void TrendOption::Validate()
{
    static const string method = "TrendOption::Validate";
    int i;


    Generic1Factor::validate();
    
    // Monitor dates must be increasing
    DateTime::ensureIncreasing(*(monitorDates.get()),
                               "monitorDates",
                               true);

    // historic monitor prices must be set.
    i = 0;
    while (i < monitorDates->size() && (*monitorDates)[i] < valueDate) {
        if (i >= histMonSamples->size() || (*histMonSamples)[i] <= 0.) {
            throw ModelException(method, "at least one historic monitoring sample has not been set");
        }
        i++;
    }
    

    // check spread types
    if (!(CString::equalsIgnoreCase(spreadType, TrendOption::SPREAD_NONE) ||
          CString::equalsIgnoreCase(spreadType, TrendOption::SPREAD_UP) ||
          CString::equalsIgnoreCase(spreadType, TrendOption::SPREAD_DOWN) ||
          CString::equalsIgnoreCase(spreadType, TrendOption::SPREAD_SYMMETRIC))) {

        throw ModelException(method, "Unsupported spread type (" + spreadType +
             "). Valid spread types are NONE, UP, DOWN, and SYMMETRIC.");
    }

    if (fwdStarting) {
        if (startDate > (*monitorDates)[0]) {
            throw ModelException(method, "Start date can't be after 1st monitoring date");
        }

        if (valueDate >= startDate && initialSpot <= 0.) {
            throw ModelException(method, "Initial spot must be set if valuing on or after start date");
        }
    }
    
    return;
}

bool TrendOption::avoidVegaMatrix(const IModel* model) 
{
    // return true if all of the options are historic.
    if ((*monitorDates)[monitorDates->size()-1] <= valueDate) {
        return true;
    }


    return false;
}


// Returns the strike for vega matrix. Throws if there's more than 1.
DoubleArraySP TrendOption::getSensitiveStrikes(OutputNameConstSP outputName,
                                               const IModel*       model)
{

    DoubleArraySP sensStrikes = DoubleArraySP(new DoubleArray(0));
    DoubleArraySP sensStrikesHigh = DoubleArraySP(new DoubleArray(0));
    DoubleArraySP sensStrikesLow = DoubleArraySP(new DoubleArray(0));

    double lowStrike;
    double highStrike;

    if (avoidVegaMatrix(model)) {
        throw ModelException("TrendOption::getSensitiveStrikes", 
                             "VEGA_MATRIX is not valid for this instrument");
    }

    // get start date for vol interpolation
    DateTime imntStartDate = fwdStarting ? startDate:valueDate;

    // create a vol request object to be passed on

    SensitiveStrikeDescriptor sensStrikeDesc;
    sensStrikeDesc.forwardOnly = false;

    int i;
    int j;
    for (i=1; i<monitorDates->size(); i++) {
        if ((*monitorDates)[i] > valueDate) {
            
            if (valueDate >= (*monitorDates)[i-1]) {
                // the already started option
                getDigStrikes((*histMonSamples)[i-1], &lowStrike, &highStrike);

                LinearStrikeVolRequestSP volRequest(new LinearStrikeVolRequest(highStrike, 
                                                               (*monitorDates)[i-1], 
                                                               (*monitorDates)[i],
                                                               false)); // fwdStarting
                
                asset->getSensitiveStrikes(volRequest.get(), outputName, 
                                           sensStrikeDesc, sensStrikesHigh);
                for (j=0; j<sensStrikesHigh->size(); j++) {
                    sensStrikes->push_back((*sensStrikesHigh)[j]);
                }
                volRequest->setStrike(lowStrike);
                asset->getSensitiveStrikes(volRequest.get(), outputName, 
                                           sensStrikeDesc, sensStrikesLow);
                for (j=0; j<sensStrikesLow->size(); j++) {
                    sensStrikes->push_back((*sensStrikesLow)[j]);
                }

            } else {
                getDigStrikes(1, &lowStrike, &highStrike);

                LinearStrikeVolRequestSP volRequest(new LinearStrikeVolRequest(highStrike, 
                                                               (*monitorDates)[i-1], 
                                                               (*monitorDates)[i],
                                                               true)); // fwdStarting
                
                asset->getSensitiveStrikes(volRequest.get(), outputName, 
                                           sensStrikeDesc, sensStrikesHigh);
                for (j=0; j<sensStrikesHigh->size(); j++) {
                    sensStrikes->push_back((*sensStrikesHigh)[j]);
                }
                volRequest->setStrike(lowStrike);
                asset->getSensitiveStrikes(volRequest.get(), outputName, 
                                           sensStrikeDesc, sensStrikesLow);
                for (j=0; j<sensStrikesLow->size(); j++) {
                    sensStrikes->push_back((*sensStrikesLow)[j]);
                }
            }
        }
    }

    return sensStrikes;
}

void TrendOption::validatePop2Object()
{
    static const string method("TrendOption::validatePop2Object");
    // one time initialization 
    if (histMonSamples->size() < monitorDates->size()) {
        histMonSamples->resize(monitorDates->size());
    }
}


// get the high and low strikes for the digital approximation. Also returns the strikeSpread for
// scaling the callSpread result. Not that for forward starting the strikes are percents but the
// strikeSpread is scaled by the spot at start so that it's always an absolute number.
void TrendOption::getDigStrikes(double strike, double *lowStrike, double *highStrike) const
{
    static const string method("TrendOption::getDigStrikes");

    if (CString::equalsIgnoreCase(spreadType, TrendOption::SPREAD_NONE)) {
        *lowStrike = strike*(1. - TrendOption::MINIMUM_ALPHA);
        *highStrike = strike*(1. + TrendOption::MINIMUM_ALPHA);
    } else if (CString::equalsIgnoreCase(spreadType, TrendOption::SPREAD_UP)) {
        *lowStrike = strike;
        *highStrike = strike*(1. + strikeSpread);
    } else if (CString::equalsIgnoreCase(spreadType, TrendOption::SPREAD_DOWN)) {
        *lowStrike = strike*(1. - strikeSpread);
        *highStrike = strike;
    } else if (CString::equalsIgnoreCase(spreadType, TrendOption::SPREAD_SYMMETRIC)) {
        *lowStrike = strike*(1. - strikeSpread/2.);
        *highStrike = strike*(1. + strikeSpread/2.);
    } else {
        throw ModelException(method, "Unsupported spread type (" + spreadType +
             "). Valid spread types are NONE, UP, DOWN, and SYMMETRIC.");
    }

    return;
}

/** when to stop tweaking */
DateTime TrendOption::endDate(const Sensitivity* ) const {
    return (*monitorDates)[monitorDates->size()-1];
}

bool TrendOption::priceDeadInstrument(CControl* , CResults * results) const
{// to do :
    return false;
}


//////////////////////////////////////////////////////////
/** product class for payoff */
//////////////////////////////////////////////////////////
class TrendOptionClosedForm: public CClosedFormLN::IProduct{
private:
    const TrendOption*  inst; // a reference

public:
    TrendOptionClosedForm(const TrendOption* inst): inst(inst){}

    void price(CClosedFormLN*   model,
               Control*         control, 
               CResults*        results) const;
protected:

};



// the main entry point for model
void TrendOptionClosedForm::price(CClosedFormLN*   ,
                                   Control*        control, 
                                   CResults*       results) const
{
    static const string method = "TrendOptionClosedForm::price";
    try {
        int i;
        int numPast = 0; // the number of stock obsevations in the past
        DoubleArraySP oddsStockUp(new DoubleArray());
        InstrumentSettlementSP matSettle = InstrumentSettlementSP(new AtMaturity());
        
        
        oddsStockUp->resize(inst->monitorDates->size()-1);

        // -- past observations. Set the odds to one for the direction the stock went.
        
        // first count the number of past observaltions
        for (i=0; i<inst->monitorDates->size(); i++) {
            if ((*inst->monitorDates)[i] <= inst->valueDate) {
                numPast++;
            }
        }

        // set the odds for as many as we can
        for (i=0; i<numPast-1; i++) {
            if ((*inst->histMonSamples)[i+1] > (*inst->histMonSamples)[i]) {
                (*oddsStockUp)[i] = 1.;
            } else {
                (*oddsStockUp)[i] = 0.;
            }
        }
        
        // -- future observations
        double highStrike;
        double lowStrike;
        bool fwdStarting;
        double strike;
        double digPrice;
        DateTime startDate;
        DateTime matDate;

        int numPastIntervals;
        if (numPast == 0) {
            numPastIntervals = 0;
        } else {
            numPastIntervals = numPast-1;
        }

        for (i=numPastIntervals; i<inst->monitorDates->size()-1; i++) {

            if (numPast == 0) { // first interval yet to start
                fwdStarting = true;
                strike = 1.;
            } else if (i == numPastIntervals) { // in the middle of an interval
                fwdStarting = false;
                strike = (*inst->histMonSamples)[i];
            } else { // future intervals
                fwdStarting = true;
                strike = 1.;
            }
            startDate = (*inst->monitorDates)[i];
            matDate = (*inst->monitorDates)[i+1];


            inst->getDigStrikes(strike, &lowStrike, &highStrike);

            digPrice = CVanilla::priceSpread(
                inst->valueDate, // valueDate
                startDate,  // startDate
                matDate, // matDate
                true, // isCall
                fwdStarting, // fwdStarting,
                true, // oneContract,
                1, // notional,
                1, // initialSpot,
                lowStrike, // lowStrike,
                highStrike, // highStrike,
                matSettle.get(),
                inst->asset.get(),
                inst->discount.get()) / (highStrike - lowStrike);
            
            if (fwdStarting == true) {
                (*oddsStockUp)[i] = digPrice/inst->asset->fwdValue((*inst->monitorDates)[i]);
            } else {
                (*oddsStockUp)[i] = digPrice;
            }
        }
         
        // Compute the odds of trend changes
        DoubleArraySP oddsTrendChange(new DoubleArray());
        oddsTrendChange->resize(inst->monitorDates->size()-2);

        for (i=0; i<oddsTrendChange->size(); i++) {
            (*oddsTrendChange)[i] = 1. - (*oddsStockUp)[i]*(*oddsStockUp)[i+1] 
                - (1.-(*oddsStockUp)[i])*(1.-(*oddsStockUp)[i+1]);
        }

        // set the payoff
        DoubleArraySP payoff(new DoubleArray());
        payoff->resize(inst->monitorDates->size()-1);
        for (i=0; i<payoff->size(); i++) {
            if (inst->isCall == true) {
                (*payoff)[i] = inst->notional*Maths::max(inst->trendMult*i - inst->strike, 0.);
            } else {
                (*payoff)[i] = inst->notional*Maths::max(inst->strike - inst->trendMult*i, 0.);
            }
        }
        
        double returnPrice = binModel(oddsTrendChange, payoff);
        
        returnPrice *= inst->instSettle->pv(
            (*inst->monitorDates)[inst->monitorDates->size()-1],
            inst->discount.get(), 
            inst->asset.get());

        results->storePrice(returnPrice, inst->discount->getCcy());

        // add any requested outputs such as FWD_AT_MAT, DELAY_PRICE
        inst->addRequests(control,
                          results,
                          returnPrice,
                          inst->endDate(NULL));
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Implementation of ClosedForm::IntoProduct interface */
CClosedFormLN::IProduct* TrendOption::createProduct(
    CClosedFormLN* ) const{

    return new TrendOptionClosedForm(this);
}

/** Rolls the value date for theta */
bool TrendOption::sensShift(Theta* shift)
{    
    int i;
    DateTime newDate = shift->rollDate(valueDate);

    if (valueDate == newDate) {
        // only set stuff if it was not set before
        for (i=0; i<monitorDates->size(); i++) {
            if (newDate == (*monitorDates)[i] && Maths::isZero((*histMonSamples)[i])) {
                (*histMonSamples)[i] = asset->getThetaSpotOnDate(shift, (*monitorDates)[i]);
            }
        }

    } else {
              
        for (i=0; i<monitorDates->size(); i++) {
            if (newDate.isGreaterOrEqual((*monitorDates)[i]) && (*monitorDates)[i].isGreater(valueDate)) {
                (*histMonSamples)[i] = asset->getThetaSpotOnDate(shift, (*monitorDates)[i]);                    
            }
        }
    }

    // roll today 
    valueDate = newDate;
    
    return true;
};

// for reflection
TrendOption::TrendOption(): Generic1Factor(TYPE){
    trendMult = 1.;
}
    

class TrendOptionHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(TrendOption, clazz);
        SUPERCLASS(Generic1Factor);
        IMPLEMENTS(CClosedFormLN::IIntoProduct);
        IMPLEMENTS(LastSensDate);
        EMPTY_SHELL_METHOD(defaultTrendOption);
        FIELD(isCall, "is the option a call as opposed to a put");
        FIELD(strike, "the strike of the option");
        FIELD(trendMult, "multiplies the number of trend switches in the payoff. redundant");
        FIELD_MAKE_OPTIONAL(trendMult);
        FIELD(monitorDates, "one date for eack sample");
        FIELD(histMonSamples, "level of u/l on hist sample dates");
        FIELD(spreadType, "increasing, decreasing, or none");
        FIELD(strikeSpread, "spread for barrier call spread approx");
    }

    static IObject* defaultTrendOption(){
        return new TrendOption();
    }
};

CClassConstSP const TrendOption::TYPE = CClass::registerClassLoadMethod(
    "TrendOption", typeid(TrendOption), TrendOptionHelper::load);
   
bool TrendOptionLoad()
{
    return (TrendOption::TYPE != 0);
}




DRLIB_END_NAMESPACE

