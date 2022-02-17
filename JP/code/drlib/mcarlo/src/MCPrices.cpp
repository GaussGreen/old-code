#include "edginc/config.hpp"
#include "edginc/Sensitivity.hpp"
#include "edginc/TwoSidedDeriv.hpp"
#include "edginc/Maths.hpp"
#include "edginc/MCPrices.hpp"
#include "edginc/MCProduct.hpp"

DRLIB_BEGIN_NAMESPACE

// IMCPrices
/** Returns the suggested mode for handling specified sensitivity.
    This is either NONE, FIRST_ORDER, or SECOND_ORDER only */
IQuickGreeks::QuickGreekType IQuickGreeks::defaultMode(
    const Sensitivity* sens){
    IQuickGreeks::QuickGreekType mode;
    // first of all can we do quick greeks for this sensitivity?
    // (exclude all discrete type of shifts eg theta)
    if (sens->discreteShift()){
        mode = NONE;
    } else {
        // the next question is whether this sensitivity is calculated by
        // a two sided shift
        if (dynamic_cast<const ITwoSidedDeriv*>(sens) != 0){
            mode = SECOND_ORDER;
        } else {
            mode = FIRST_ORDER;
        }
    }
    return mode;
}

/** Returns true if 'quick greeks' is applicable for this
    sensitivity */
bool IQuickGreeks::doQuickGreeks(const Sensitivity* sens){
    // exclude all discrete type of shifts eg theta
    return (!sens->discreteShift());
}

IMCPrices::~IMCPrices() {}

MCPricesSimple::MCPricesSimple(
    int             NbIter,
    int             NbSubSamples):
        NbIter(NbIter), 
        NbSubSamples(NbSubSamples), 
        numPerSubSample(NbIter/NbSubSamples)
{
    reset();
}

MCPricesSimple::~MCPricesSimple(){}

/** deep copy 
    Allocate via emptyConstructor so the correct memory footprint 
    is created.  */
IMCPrices* MCPricesSimple::clone() const{
    MCPricesSimple* copy = dynamic_cast<MCPricesSimple*>(emptyConstructor());
    copy->NbIter = NbIter;
    copy->NbSubSamples = NbSubSamples;
    copy->numPerSubSample = numPerSubSample;
    copy->subSampPos = subSampPos;
    copy->sampleSum = sampleSum;
    copy->sumSoFar = sumSoFar;
    copy->sumSqrSoFar = sumSqrSoFar;
    copy->lastPriceStored = lastPriceStored;
    return copy;
}

/** Keeps running total of sum and square of sum (wrt sample) */
void MCPricesSimple::add(double   price) {
    lastPriceStored = price;
    sampleSum += price;
    subSampPos++;
    if (subSampPos == numPerSubSample){
        subSampPos = 0;
        double volatile sampleMean = sampleSum/numPerSubSample;
        sumSoFar += sampleMean; // update running total (of sum)
        sumSqrSoFar += sampleMean * sampleMean; // update total of sum^2
        sampleSum = 0.0;
    }
}
/** Returns the last price passed to the add() method. Undefined
    behaviour if no prices yet stored */
double MCPricesSimple::lastPrice() const{
    return lastPriceStored;
}


// bool MCPricesSimple::repriceForGreek(int pathIdx) {
//     return true;
// }

void MCPricesSimple::getResult(double* result,
                                        double* resultStdErr) const{
    if (subSampPos != 0){
        throw ModelException("MCPricesSimple::getResult", "Internal" 
                             " error - not all prices have been supplied");
    }
    double volatile sum = sumSoFar/(double)NbSubSamples;
    double volatile sumSqr = sumSqrSoFar/(double)NbSubSamples;
    *result = sum;
    if (NbSubSamples > 1) {
        double volatile tmp = sum * sum; /* force read/write to memory - more
                                            consistent numbers */
        double stdErr = (sumSqr - tmp)/((double)NbSubSamples-1.0);
        *resultStdErr = stdErr < 0.0? 0.0: sqrt(stdErr);
    } else {
        *resultStdErr = 0.0;
    }
}

/** Clears out SumSubPrices and resets iSubSample */
void MCPricesSimple::reset(){
    subSampPos = 0;
    sampleSum = 0.0;
    sumSoFar = 0.0;
    sumSqrSoFar = 0.0;
}

/**  returns MAX(x, 0.0). */
double MCPricesSimple::maxWithZero(double x) const{
    return Maths::max(x,0.0);
}

IMCPrices* MCPricesSimple::emptyConstructor() const {
    return new MCPricesSimple(1, 1);
}

int MCPricesSimple::storagePerPath(IMCProduct* product) const {
    return 0;
}

void MCPricesSimple::configureCache(const IntArray& changedAssets) {
    ; // empty
}


/** Constructor: First two params as per MCPricesSimple. The reprice
    parameter is used to control which paths are skipped. Note that
    a reference is taken to the reprice object. The caller must ensure
    that the right data is passed to it and that it is updated 
    correctly for each greek */
MCPricesGeneral::MCPricesGeneral(
    int               nbIter,
    int               nbSubSamples,
    const IRepriceSP& reprice):
        MCPricesSimple(nbIter, nbSubSamples), mode(NONE), 
        path(0), reprice(reprice){}

IMCPrices* MCPricesGeneral::emptyConstructor() const 
{
    return new MCPricesGeneral(1, 1, IRepriceSP());
}
        
/** Returns a deep copy of this object */
IMCPrices* MCPricesGeneral::clone() const
{
    // parent first s
    MCPricesGeneral& copy = dynamic_cast<MCPricesGeneral&>(*MCPricesSimple::clone());
    // then our fields
    IReprice* repriceCopy = !reprice? 0: 
        &dynamic_cast<IReprice&>(*reprice->clone());
    IRepriceSP repriceSP(repriceCopy);
    copy.reprice = repriceSP;
    copy.mode = mode;
    copy.path = path;
    return &copy;
}

/** Sets mode. Does not support CROSS_GAMMA. tolerance is 
    threshold used to decide when to skip paths for second order
    greeks */
void MCPricesGeneral::setMode(IQuickGreeks::QuickGreekType mode){
    switch (mode){
    case NONE:
    case FIRST_ORDER:
    case SECOND_ORDER:
    case CROSS_GAMMA:
        this->mode = mode;
        break;
    default:
        throw ModelException("MCPricesGeneral::setMode",
                             "Unrecogonised mode");
    }
}

/** Returns the reprice object */
IRepriceSP MCPricesGeneral::getReprice(){
    return reprice;
}

/** Uses internal Reprice object to determine whether path can
    be skipped or not. If the path is to be skipped add() is called
    with the relevant value */
bool MCPricesGeneral::repriceForGreek(int pathIdx){
    if (!reprice){
        // shouldn't get here unless mode == NONE, but just in case ...
        return true;
    }
    switch (mode){
    case NONE:
        return true;
    case FIRST_ORDER:
    {
        if (reprice->firstDerivZero(pathIdx)){
            add(reprice->originalPrice(pathIdx));
            return false;
        }
        return true;
    }
    case SECOND_ORDER:
    {
        if (reprice->firstNumericalDerivZero(pathIdx)){
            add(reprice->originalPrice(pathIdx));
            return false;
        }
        return true;
    }
    case CROSS_GAMMA:
    {
        bool crossNumericalDerivZero = 
            reprice->crossNumericalDerivZero(pathIdx);
        if (crossNumericalDerivZero){
            MCPricesSimple::add(0.0); // we report change in price
            path++; // must increment path as we don't call our add method
        }
        return (!crossNumericalDerivZero);
    }
    default:
        return true; // shouldn't get here
    }
}
    
/** adds supplied price to this set of IMCPrices. Running totals of
    sum and sum^2 (wrt sample) are kept. Note that this must be the
    real price (ie after any MAX etc) */
void MCPricesGeneral::add(double price){
    if (mode == CROSS_GAMMA){
        // for cross gamma we report change in price
        double origPrice = reprice->originalPrice(path);
        path++;
        double newPrice = price - origPrice;
        MCPricesSimple::add(newPrice);
        // and correct lastPriceStored
        lastPriceStored = price;
    } else {
        MCPricesSimple::add(price);
    }
}

/** Resets internal variables ready for another pricing run */
void MCPricesGeneral::reset(){
    MCPricesSimple::reset();
    path = 0;
}

/**  returns MAX(x, 0.0) unless doing first order quick greeks
     in which case x is returned */
double MCPricesGeneral::maxWithZero(double x) const{
    if (mode == FIRST_ORDER){
        return x;
    }
    return Maths::max(x,0.0);
}


MCPricesGeneral::~MCPricesGeneral(){}

DRLIB_END_NAMESPACE


