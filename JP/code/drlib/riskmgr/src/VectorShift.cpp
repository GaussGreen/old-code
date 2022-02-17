//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VectorShift.cpp
//
//   Description : Tweaking for things with a name and a set of expiries ie
//                 rho pointwise, vega pointwise etc
//
//   Author      : Mark A Robson
//
//   Date        : 5 March 2001
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#define QLIB_VECTORSHIFT_CPP
#include "edginc/VectorShift.hpp"
#include "edginc/ExpiryWindow.hpp"
#include "edginc/ExpiryResult.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/Results.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/SensMgr.hpp"


DRLIB_BEGIN_NAMESPACE

/** Slightly more specialised ScalarShift where the shifting
    information is qualified by an expiry (eg vega pointwise).
    Note still abstract */

VectorShift::~VectorShift(){}

void VectorShift::setExpiryToTweak(ExpiryConstSP expiry){
    /* Store the expiry in the cachedExpiries array and set the expiriesAreSet
     * flag so that we do not overwrite them */
    cachedExpiries.reset(new ExpiryArray(1, ExpirySP::constCast(expiry)));
    expiriesAreSet = true;
}

void VectorShift::setExpiriesToTweak(ExpiryArrayConstSP expiries){
    /* Store the expiry in the cachedExpiries array and set the expiriesAreSet
     * flag so that we do not overwrite them */
    cachedExpiries.reset(new ExpiryArray(*expiries));
    expiriesAreSet = true;
}

/** implements a one sided vector derivative for each instance of the
    market data which is sensitive to this SensControl */
void VectorShift::calculate(TweakGroup*  tweakGroup,
                            CResults*    results){
    calculate(tweakGroup, results, results);
}

/** As above but uses alternative results set for base price */
void VectorShift::calculate(TweakGroup*  tweakGroup,
                            CResults*    resultsForTweak,
                            CResults*    resultsForBasePrice){
    try{
        // get list of names to calculate result for. Remove blanks
        OutputNameArrayConstSP names(this->names(tweakGroup, resultsForTweak));

        if (names->empty()) {
            resultsForTweak->storeNotApplicable(this);            
        }
    
        // see if the instrument has a last sens date method
        LastSensDate* lsd = 
            dynamic_cast<LastSensDate*>(tweakGroup->getInstrument());
        DateTime      endDate;
        DateTime      valueDate;
        if (lsd) {
            valueDate = tweakGroup->getInstrument()->getValueDate();
        }
        for (int idx = 0; idx < names->size(); idx++){
            // store the name of what we want to shift
            setMarketDataName((*names)[idx]);
            /* skip over where result has been calculated already */
            if (!resultsForTweak->exists(this)){
                double origShiftSize = getShiftSize(); // save this
                try {
                    if (lsd) {
                        // end date can be on a per name basis
                        endDate = lsd->endDate(this);
                    }
                    // now get expiries for which we should tweak
                    cachedExpiries = getExpiries(tweakGroup);
                    // and then what shifts to make for each expiry
                    DoubleArrayConstSP tweakShiftSizes = 
                        calculateTweakSizes(tweakGroup, (*names)[idx], 
                                            *cachedExpiries);
                    // create room for storing the results
                    ExpiryResultArraySP tweaks(new ExpiryResultArray(
                        cachedExpiries->size()));
                    // then loop over the expiries
                    bool expired = false;
                    for (int jdx = 0; jdx < cachedExpiries->size(); jdx++){
                        // store the expiry/shiftSize which we want to tweak
                        setExpiry((*cachedExpiries)[jdx]);
                        setShiftSize((*tweakShiftSizes)[jdx]);
                        // calculate sens (if not expired)
                        double firstDeriv = expired? 0.0:
                            calcOneSidedFirstDeriv(tweakGroup, 
                                                   resultsForBasePrice);
                        // store result in array
                        (*tweaks)[jdx] = ExpiryResult(expiry, firstDeriv);

                        // do we need to tweak anymore ?
                        if (lsd) {
                            expired = expiry->toDate(valueDate).
                                isGreater(endDate);
                        }
                    }
                    setShiftSize(origShiftSize); // restore value
                    // and store it
                    resultsForTweak->storeGreek(tweaks, this);
                }
                catch (exception& e) {
                    setShiftSize(origShiftSize); // restore value
                    resultsForTweak->storeGreek(IObjectSP(
                                                    new Untweakable(e)), this);
                }
            }
        }
    } catch (exception& e){
        resultsForTweak->storeGreek(IObjectSP(new Untweakable(e)), 
                                    getSensOutputName(),
                                    OutputNameSP(new OutputName("")));
    }
}

/** Calculate shift sizes for given name and set of expiries. Default 
    implementation returns DoubleArray(expiries.size(), getShiftSize()) */
DoubleArrayConstSP VectorShift::calculateTweakSizes(
    IObject*           tweakGroup, 
    OutputNameConstSP  name,
    const ExpiryArray& expiries){
    return DoubleArrayConstSP(new DoubleArray(expiries.size(), getShiftSize()));
}

/** Returns the expiry which is currently being tweaked */
ExpiryConstSP VectorShift::getExpiry() const {
    if ( !expiry ) {
       throw ModelException("expiry is Null",  "VectorShift::getExpiry");
    }
    return expiry;
}

ExpiryWindowConstSP VectorShift::getExpiryWindow() const {
    return ExpiryWindow::around(getExpiries(), getExpiry());
}

/** sets the expiry which is currently being tweaked to a certain value*/
void VectorShift::setExpiry(ExpiryConstSP newExpiry) {
    expiry = newExpiry;
}

bool VectorShift::getExpiriesAreSet() const {
    return expiriesAreSet;
}

/** Returns the expiries which are to be tweaked */
ExpiryArrayConstSP VectorShift::getExpiries(const IObject* tweakGroup){
    if (expiriesAreSet){
        // If expiriesAreSet, the cachedExpiries contains the expiries to be tweaked
        return cachedExpiries;
    }
    IObjectConstSP expiriesObj = SensMgrConst(tweakGroup).qualifier(this);
    ExpiryArrayConstSP expiries = ExpiryArrayConstSP::dynamicCast(expiriesObj);
    return expiries;
}


/** Returns the expiries which are to be tweaked. Only valid once
    the tweaking process has started. ie value is derived during
    the tweaking */
ExpiryArrayConstSP VectorShift::getExpiries()const{
    if (!cachedExpiries ) {
       throw ModelException("expiries are Null",  "VectorShift::getExpiries");
    }
    return cachedExpiries;
}

VectorShift::VectorShift(const CClassConstSP& clazz,
                         const string&        outputName,
                         const double&        shiftSize): 
    ScalarShift(clazz, outputName, shiftSize),
    expiriesAreSet(false) {}

/** for reflection */
VectorShift::VectorShift(const CClassConstSP& clazz,
                         const string&        sensName):
    ScalarShift(clazz, sensName),
    expiriesAreSet(false) {}

class VectorShiftHelper{
public:
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(VectorShift, clazz);
        SUPERCLASS(ScalarShift);
        FIELD(expiry, "expiry");
        FIELD_MAKE_TRANSIENT(expiry);
        FIELD(cachedExpiries, "cachedExpiries");
        FIELD_MAKE_TRANSIENT(cachedExpiries);
        FIELD(expiriesAreSet, "expiriesAreSet");
        FIELD_MAKE_TRANSIENT(expiriesAreSet);

    }
};

CClassConstSP const VectorShift::TYPE = CClass::registerClassLoadMethod(
    "VectorShift", typeid(VectorShift), VectorShiftHelper::load);

DRLIB_END_NAMESPACE

