//---------------------------------------------------------------------------
//
//   Group       : Convertibles Research
//
//   Filename    : ResetSchedule.cpp
//
//   Description : A reset schedule class
//
//   Author      : André Segger
//
//   Date        : 07 January 2003
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/ResetSchedule.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Addin.hpp"

DRLIB_BEGIN_NAMESPACE

ResetSchedule::ResetSchedule(const DateTime&      valueDate,
                             const DateTimeArray& dates,
                             const DoubleArray&   minResetStrikes,
                             const DoubleArray&   maxResetStrikes) : CObject(TYPE),
                                                                   resetDates(new DateTimeArray(dates)),
                                                                   minResetStrike(new DoubleArray(minResetStrikes)),
                                                                   maxResetStrike(new DoubleArray(maxResetStrikes)),
                                                                   valueDate(valueDate) {
    resetType = "UP";
    resetTypes = CStringArraySP(new CStringArray(resetDates->size()));
    for (int i=0 ; i<resetTypes->size() ; ++i) {
        (*resetTypes)[i] = "UP";
    }
    validatePop2Object();
    fwdStarting = false;    
    initLevel = CashFlowArray(0);
}

ResetSchedule::~ResetSchedule() {}

int ResetSchedule::length() const {
    return resetDates->getLength();
}


void ResetSchedule::validatePop2Object() {
    static const string method = "ResetSchedule::validatePop2Object";
    try {
        int i;

        if (!resetDates) {
            throw ModelException(method, "Reset dates must not be Null");
        }
        if (!resetLevel) {
            throw ModelException(method, "Reset levels must not be Null");
        }
        if (!minResetStrike) {
            throw ModelException(method, "Min reset strikes must not be Null");
        }

        if ( resetDates->size() != minResetStrike->size() ) {
           throw ModelException(method,
                                "Number of reset dates (" +
                                Format::toString(resetDates->size()) +
                                ") must be the same as number of min reset strikes (" +
                                Format::toString(minResetStrike->size()) + ").");
        }

        if ( !!maxResetStrike && resetDates->size() < maxResetStrike->size() ) {
           throw ModelException(method,
                                "Number of reset dates (" +
                                Format::toString(resetDates->size()) +
                                ") must be the same as number of max reset strikes (" +
                                Format::toString(maxResetStrike->size()) + ").");
        }

        if ( !(resetDates->size() == resetLevel->size())) {
           throw ModelException(method,
                                "Number of reset dates (" +
                                Format::toString(resetDates->size()) +
                                ") must be the same as number of reset levels (" +
                                Format::toString(resetLevel->size()) + ").");
        }

        if ( !!resetDates ) {
           for (i = 1; i < length(); i++) {
               if (!(*resetDates)[i].isGreater((*resetDates)[i-1])) {
                   throw ModelException(method,
                                        "Reset dates not in increasing order (" +
                                        (*resetDates)[i].toString() + " <= " +
                                        (*resetDates)[i-1].toString() + ")");
               }
           }
        }

        if ( !!resetTypes && resetTypes->size() > 0) {
            if (!(resetTypes->size() == resetDates->size())) {
               throw ModelException(method,
                                "Number of reset dates (" +
                                Format::toString(resetDates->size()) +
                                ") must be the same as number of reset types (" +
                                Format::toString(resetTypes->size()) + ").");
            }
        }

        if (!!receiveIntrinsic) {
            if (!(receiveIntrinsic->size() == resetDates->size())) {
               throw ModelException(method,
                                "Number of reset dates (" +
                                Format::toString(resetDates->size()) +
                                ") must be the same as number of receive intrinsic flags (" +
                                Format::toString(receiveIntrinsic->size()) + ").");
            }
        } else {
            receiveIntrinsic = BoolArraySP(new BoolArray(resetDates->size()));
            for (i=0 ; i<receiveIntrinsic->size() ; ++i) {
                (*receiveIntrinsic)[i] = false;
            }
        }

        // create a max reset strike if necessary
        if ( !maxResetStrike ) {
            maxResetStrike = DoubleArraySP(new DoubleArray(resetDates->size()));
        } else {
            // check that all max reset strikes greater or equal to min reset strikes
           for (i = 0; i < maxResetStrike->size(); ++i) {
               if ((*minResetStrike)[i] > (*maxResetStrike)[i]) {
                   throw ModelException(method, "Min strikes must not be greater than max strikes!");
               }
               if (!Maths::isPositive((*minResetStrike)[i])) {
                   throw ModelException(method, "Min strikes must be strictly positive!");
               }
           }
        }

        // create a parity array if necessary
        if ( !parity) {
            parity = DoubleArraySP(new DoubleArray(resetDates->size()));
            for (i=0 ; i<parity->size() ; ++i) {
                (*parity)[i] = 1.0;
            }

        } else {
            if (resetDates->size() > parity->size()) {
                // extend the parity array if necessary
                int numStrikes = parity->size();
                parity->resize(resetDates->size());
                for (i=numStrikes ; i<resetDates->size() ; ++i) {
                    (*parity)[i] = 1.0;
                }
            } else {
                if (resetDates->size() < parity->size()) {
                    throw ModelException(method, "Must not have more parity values than reset dates!");
                }
            }
        }


        // extend the max reset strike array if necessary
        if ( resetDates->size() > maxResetStrike->size() ) {
            int numStrikes = maxResetStrike->size();
            maxResetStrike->resize(resetDates->size());
            for (i=numStrikes ; i<maxResetStrike->size() ; ++i) {
                (*maxResetStrike)[i] = 0.0;
            }
        }

        // for fwdStarting
        if (fwdStarting){
            if (initLevel.size() != 1)
	            throw ModelException(method, "The size of initLevel array must be one row.  Not available Average-In, yet.");
            if (resetDates->size() > 0) {
                if (initLevel[0].date > (*resetDates)[0]) {
                    throw ModelException(method, "startDate (last date in initLevel) must be earlier than the first reset dates!");
                }
            } else {
                throw ModelException(method, "Reset dates must be provided if the instrument is forward starting");
            }
        }
    }
    catch (exception& e) {
        throw ModelException(&e, method);
    }
}

string ResetSchedule::getResetType(const DateTime& date)
{
    string type = "";
    if (!!resetTypes && resetTypes->size() > 0) {
        int i = 0;
        while ( i < resetDates->size() &&
                date >= (*resetDates)[i] ) {
            ++i;
        }

        if (i >= resetDates->size()-1) {
            type = (*resetTypes)[resetDates->size()-1];
        } else {
            type = (*resetTypes)[i];
        }
    } else {
        type = resetType;
    }
    return type;
}

bool ResetSchedule::getReceiveIntrinsic(const DateTime& date) {
    static const string method = "ResetSchedule::getReceiveIntrinsic";
    try {

        for (int i = 0; i<resetDates->size(); i++) {
            if ((*resetDates)[i] == date) {
                return (*receiveIntrinsic)[i];
            }
        }

        throw ModelException(method, "Date passed is not a reset date");

    }
    catch (exception &e) {
        throw ModelException(e, method);
    }

}

bool ResetSchedule::isFlat(const DateTime& date)
{
    bool flat = true;
    if (!!resetTypes && resetTypes->size() > 0) {
        string type = "";
        int i = 0;
        while ( i < resetDates->size() &&
                date >= (*resetDates)[i] ) {
            ++i;
        }

        if (i >= resetDates->size()-1) {
            --i;
        }

        flat = true;
        for (int j=1 ; j<=i ;++j) {
            if ((*resetTypes)[j] != (*resetTypes)[j-1]) {
                flat = false;
                break;
            }
        }
    }

    return flat;
}


ResetSchedule::ResetSchedule(): CObject(TYPE) {
    fwdStarting = false;
}

/** return date list (deep copy) */
DateTimeArray ResetSchedule::getDates() const{
    return *resetDates;
}

/** returns reference to date array */
const DateTimeArray &ResetSchedule::getDateArray() const
{
    return *resetDates;
}

/** return value list (deep copy) */
DoubleArray ResetSchedule::getMinResetStrikes() const{
    return *minResetStrike;
}

/** return max reset strike list (deep copy) - bit garbage ... do I need this? */
DoubleArray ResetSchedule::getMaxResetStrikes() const{
    if (!maxResetStrike) {
        DoubleArray dblArray(0);
        return dblArray;
    } else {
        return *maxResetStrike;
    }
}

/** return the parity at which the bond gets reset (deep copy) */
DoubleArray ResetSchedule::getParity() const {
    if (!parity) {
        DoubleArray dblArray(0);
        return dblArray;
    } else {
        return *parity;
    }
}



/** get reset information for a particular date */
bool ResetSchedule::hasReset(const DateTime& date) const
{
    bool found = false;
    for(int i=0;i<resetDates->size();++i) {
       if ( (*resetDates)[i] == date) {
          found = true;
          break;
       }
    }

    return found;
}

// Get conversion ratio following reset based on spot and reset parameters
// faceValue, maxStrike, minStrike, spot quantities passed must be in consistent ccy (payoff or underlying)
// Returns Ratio = FaceValue / (spot * resetPremium) subject to minStrike<=spot*resetPremium<=maxStrike
 double ResetSchedule::getResetRatio(double faceValue,    
                                           double maxStrike,    // Max conversion price following reset 
                                           double minStrike,    // Min conversion price following reset 
                                           double resetPremium,   // Multiple of spot to give new conversion price after reset (otherwise known as parity)
                                           double spot) {
    static const string method = "ResetSchedule::getResetRatio";
    try {
        double targetConversionPrice = spot*resetPremium;
        targetConversionPrice = max(minStrike, targetConversionPrice);
        targetConversionPrice = min(maxStrike, targetConversionPrice);
        double ratio = faceValue / targetConversionPrice ;
        return ratio;
    }
    catch (exception &e) {
        throw ModelException(e, method);
    }
}

/** get current conversion ratio based on the initial conv ratio and a history of fixings */
double ResetSchedule::getCurrentConversionRatio(const double&   initialConvRatio,
                                                const DateTime& valueDate,
                                                const double&   faceValue)
{
    static const string method = "ResetSchedule::getCurrentConversionRatio";
    try {
        double  currentConvRatio = initialConvRatio;
        int     i                = 0;
        double targetCR;
        while ( i < resetDates->size() &&
            valueDate >= (*resetDates)[i] ) {

            // Validate that the past reset level is populated
            if ((valueDate > (*resetDates)[i]) && Maths::isZero((*resetLevel)[i])) {
                throw ModelException(method,
                        "Past reset level for " + (*resetDates)[i].toString() +
                        " must be populated");
            }

            // assumes all strikes were using the same fx rate
            // parity[] is the "multiplier by spot to give new conversion price at reset"
            // minResetStrike/maxResetStrike are the min/max allowable conversion price following reset

            targetCR = ResetSchedule::getResetRatio(faceValue, 
                                     (*maxResetStrike)[i],
                                     (*minResetStrike)[i],
                                     (*parity)[i],
                                     (*resetLevel)[i]);

            if ((!!resetTypes && resetTypes->size() > 0)) {
                if ((*resetTypes)[i] == "UP") {
                    if (targetCR > currentConvRatio) {
                        currentConvRatio = targetCR;
                    }
                } else {
                    currentConvRatio = targetCR;
                }
            } else {
                if (resetType == "UP") {
                    if (targetCR > currentConvRatio) {
                        currentConvRatio = targetCR;
                    }
                } else {
                    currentConvRatio = targetCR;
                }
            }
            ++i;
        }

        return currentConvRatio;
    }
    catch (exception &e) {
        throw ModelException(e, method);
    }
}

double ResetSchedule::getMaxResetStrike(const DateTime& date) const
{
    double maxStrike = DBL_MAX;
    bool   found     = false;
    int    idx       = 0;
    int    i;

    for( i=0 ; i<resetDates->size() ; ++i ) {
       if ( (*resetDates)[i] == date) {
          found = true;
          idx   = i;
          break;
       }
    }

    if ( found ) {
        if (!!resetTypes && resetTypes->size() > 0) {
            maxStrike = (*maxResetStrike)[i];
        } else {
            if (resetType == "UP") {
                for( i=0 ; i<=idx ; ++i ) {
                    if ( !!resetLevel && resetLevel->size() > 0 ) {
                        if ( (*resetDates)[i] <= valueDate  &&
                             (*resetLevel)[i] < maxStrike   &&
                             Maths::isPositive((*resetLevel)[i]) ) {
                            maxStrike = (*resetLevel)[i];
                        }
                    }
                }
                if ( Maths::isPositive((*maxResetStrike)[idx])) {
                    maxStrike =  Maths::min(maxStrike,(*maxResetStrike)[idx]);
                }

            } else {
                maxStrike = (*maxResetStrike)[i];
            }
        }
    }

    return maxStrike;
}

double ResetSchedule::getMinResetStrike(const DateTime& date) const
{
    double minStrike = 0.0;
    bool   found     = false;
    int    idx       = 0;
    int    i;

    for( i=0 ; i<resetDates->size() ; ++i ) {
       if ( (*resetDates)[i] == date) {
          found = true;
          idx   = i;
          break;
       }
    }

    if ( found ) {
        minStrike = (*minResetStrike)[i];
    }

    return minStrike;
}

double ResetSchedule::getParity(const DateTime& date) const {
    double parityLevel = 0.0;
    bool   found     = false;
    int    idx       = 0;
    int    i;

    for( i=0 ; i<resetDates->size() ; ++i ) {
       if ( (*resetDates)[i] == date) {
          found = true;
          idx   = i;
          break;
       }
    }
    if ( found ) {
        parityLevel = (*parity)[i];
    }
    return parityLevel;
}


void ResetSchedule::rollDate(const DateTime& oldValueDate, const DateTime& newValueDate,
                             const double newSpot, const bool isCcyStruck, const double newFX)
{
    int i;
    valueDate = newValueDate;
    DateTime startDate = fwdStarting? initLevel[0].date : valueDate;    
    if (fwdStarting && startDate > oldValueDate && startDate <= valueDate){
        initLevel[0].amount = newSpot;
    }
    else{
        for( i=0 ; i<resetDates->size() ; ++i ) {
           if ( (*resetDates)[i] > oldValueDate && (*resetDates)[i] <= valueDate) {
               (*resetLevel)[i] = newSpot;
               if ( isCcyStruck && !!resetFX) {
                   (*resetFX)[i]    = newFX;
               }
           } else if ( (*resetDates)[i] >= oldValueDate && (*resetDates)[i] <= valueDate) {
               // roll-to-now: set the reset level, if it hasn't been set yet
               if (!Maths::isPositive((*resetLevel)[i])) {
                   (*resetLevel)[i] = newSpot;
                   if ( isCcyStruck && !!resetFX) {
                       (*resetFX)[i]    = newFX;
                   }
               };
           }
        }
    }
}

void ResetSchedule::setValueDate(const DateTime& valueDate)
{
    this->valueDate = valueDate;
}

void ResetSchedule::preProcessSchedule(const bool isCcyStruck, const double& fxSpot, DoubleArraySP fwdFXs,
                                       const double& initialConvPrice)
{
    static const string method = "ResetSchedule::preProcessSchedule";
    try {
        // currency translate all values into bond currency
        int i;
        if ( isCcyStruck ) {
            // need to do some validation here, since I don't have the asset info in validatePop2Obj
            if (!resetFX) {
                throw ModelException(method, "Reset fx rates must not be Null");
            }
            if (resetFX->size() != resetDates->size()) {
                throw ModelException(method,
                                    "number of dates (" +
                                    Format::toString(resetDates->size()) +
                                    ") must be the same as number of past FX fixings (" +
                                    Format::toString(resetFX->size()) + ").");
            }

            // try/catch clause to avoid Solaris compiler bug
            try {
                // fx-translate the schedule
                for (i=0 ; i<resetDates->size() ; ++i) {
                    if ( (*resetDates)[i] <= valueDate ) {

                        // Validate that the past resetFX is populated
                        if (Maths::isZero((*resetFX)[i])) {
                            throw ModelException(method,
                                "Past resetFX value for " + (*resetDates)[i].toString()
                                + " must be populated");
                        }
                        (*minResetStrike)[i] *= (*resetFX)[i];
                        (*maxResetStrike)[i] *= (*resetFX)[i];
                        (*resetLevel)[i]     *= (*resetFX)[i];
                    }
                    else {
                        (*minResetStrike)[i] *= (*fwdFXs)[i];
                        (*maxResetStrike)[i] *= (*fwdFXs)[i];
                    }
                }
            }
            catch (exception &e) {
                throw ModelException(e, method);
            }
        }
        // if the maxStrikes are 0.0, default to the initialConvPrice
        for (i=0 ; i<resetDates->size() ; ++i) {
            if (Maths::isZero((*maxResetStrike)[i])) {
                (*maxResetStrike)[i] = initialConvPrice;
            }
        }
    }
    catch (exception &e) {
        throw ModelException(e, method);
    }
}

bool ResetSchedule::isStarted(DateTime valDate){
    if (fwdStarting)
        return initLevel[initLevel.size()-1].date < valDate;
    else
        return false;
}

DateTime ResetSchedule::getStartDate(DateTime valDate){
    static const string method = "ResetSchedule::getStartDate";
    try{
        DateTime startDate;
        startDate = initLevel[initLevel.size()-1].date;
        // if it's already started, return valDate
        if(isStarted(valDate))
            startDate = valDate;
        return startDate;
    }
    catch (exception &e) {
        throw ModelException(e, method);
    }
}

double ResetSchedule::getStartLevel(DateTime valDate, double altLvl){
    static const string method = "ResetSchedule::getStartLevel";
    try{
        double startLevel = altLvl; // return altLvl if it's not started.
        if (isStarted(valDate)){
            // check the past sampling are given.  
            for (int i=0; i<initLevel.size(); i++){
                if (initLevel[i].amount <= 0.0 && valDate > initLevel[i].date){
                    throw ModelException(method,
                                    "Past sampled value of initLevel [" + initLevel[i].date.toString()
                                    + "] must be more than 0.");
                }
            }
            // assuming one size (not average in)
            startLevel = initLevel[initLevel.size()-1].amount;
        }
        return startLevel;
    }
    catch (exception &e) {
        throw ModelException(e, method);
    }
}

bool ResetSchedule::isFwdStart()
{
        return fwdStarting;
}

void ResetSchedule::scaleLevels(const double scaleFactor)
{
    static const string method = "ResetSchedule::scaleLevel";
    try {
        int i;
        for (i=0 ; i<resetDates->size() ; ++i) {
            (*minResetStrike)[i] *= scaleFactor;
            (*maxResetStrike)[i] *= scaleFactor;
        }
    }
    catch (exception &e) {
        throw ModelException(e, method);
    }
}

class ResetScheduleHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ResetSchedule, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultResetSchedule);
        FIELD(resetDates,        "reset dates");
        FIELD(minResetStrike,    "minimum conversion price following reset");
        FIELD(maxResetStrike,    "maximum conversion price following reset");
        FIELD_MAKE_OPTIONAL(maxResetStrike);
        // parity is a misleading name. When < 1.0 will give a parity after reset > 1.0.
        FIELD(parity,            "multiplier by spot to give new conversion price after reset");
        FIELD_MAKE_OPTIONAL(parity);
        FIELD(resetLevel,        "past reset levels");
        FIELD(resetFX,           "past reset fx rates");
        FIELD_MAKE_OPTIONAL(resetFX);
        FIELD(resetTypes,        "reset types");
        FIELD_MAKE_OPTIONAL(resetTypes);
        FIELD(resetType,  "UP = reset conversion ratio up only, UP_DOWN = always reset conversion ratio");
        FIELD_MAKE_OPTIONAL(resetType);
        FIELD(valueDate,  "Value Date");
        FIELD_MAKE_OPTIONAL(valueDate);
        FIELD(fwdStarting,  "TRUE for forward starting, the initial conversion will determined in future.");
        FIELD_MAKE_OPTIONAL(fwdStarting);
        FIELD(initLevel,  "fwd Start Dates and Spot Level (It's allow only one date for time being!!)");
        FIELD_MAKE_OPTIONAL(initLevel);
        FIELD(receiveIntrinsic, "Receive max(parity-faceValue,0) on reset date");
        FIELD_MAKE_OPTIONAL(receiveIntrinsic);
        Addin::registerConstructor("RESET_SCHEDULE",
                                   Addin::UTILITIES,
                                   "Creates a Reset Schedule",
                                   ResetSchedule::TYPE);
    }

    static IObject* defaultResetSchedule(){
        return new ResetSchedule();
    }
};

CClassConstSP const ResetSchedule::TYPE = CClass::registerClassLoadMethod(
    "ResetSchedule", typeid(ResetSchedule), ResetScheduleHelper::load);

DRLIB_END_NAMESPACE

