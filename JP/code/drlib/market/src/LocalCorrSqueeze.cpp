//----------------------------------------------------------------------------
//
//   Group       : QR Equities London
//
//   Filename    : LocalCorrSqueeze.cpp
//
//   Description : Holds market data for local correlation skew model
//
//   Author      : Eva Strasser
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#ifndef LOC_CORR_SQUEEZE_CPP
#define LOC_CORR_SQUEEZE_CPP

#include "edginc/LocalCorrSqueeze.hpp"
#include "edginc/Maths.hpp"
#include "edginc/MarketDataFetcher.hpp"
#include "edginc/MDFUtil.hpp"
#include "edginc/Nrfns.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/MarketData.hpp"
#include "edginc/TimeMetric.hpp"

DRLIB_BEGIN_NAMESPACE

/** Interface for market data for local correlation skew model */

CClassConstSP const LocalCorrSqueeze::TYPE = CClass::registerClassLoadMethod(
    "LocalCorrSqueeze", typeid(LocalCorrSqueeze), LocalCorrSqueeze::load);

const string LocalCorrSqueeze::LOCAL_CORR_SQUEEZE_STEP = "STEP";
const string LocalCorrSqueeze::LOCAL_CORR_SQUEEZE_LINEAR = "LINEAR";

LocalCorrSqueeze::~LocalCorrSqueeze() {} 

void LocalCorrSqueeze::validatePop2Object() {
    static const string method("LocalCorrSqueeze::validatePop2Object");
    try {
        /** populate interpolationTypeInt */
        if (CString::equalsIgnoreCase(interpolationType, LOCAL_CORR_SQUEEZE_STEP)) {
            interpolationTypeInt = 1; // STEP
        } else if (CString::equalsIgnoreCase(interpolationType, LOCAL_CORR_SQUEEZE_LINEAR)) {
            interpolationTypeInt = 2; // LINEAR
        } else {
            throw ModelException(method, "Squeeze interpolation type has to be either "
                    + LocalCorrSqueeze::LOCAL_CORR_SQUEEZE_STEP + " or " 
                    + LocalCorrSqueeze::LOCAL_CORR_SQUEEZE_LINEAR + " but is " 
                    + interpolationType + ".");
        }        

        /** STEP -- size of strike vector + 1 should equal nb of columns of squeeze matrix */
        if (interpolationTypeInt == 1) {
            if (!(strikes->size() + 1 == squeezes->numCols()))  {
                throw ModelException(method, "For piecewise constant squeezes, "
                    "size of strikes +1 has to equal nb of columns of squeeze matrix. \n"
                    "This is not the case for category " + name + ".");
            }
        }
        
        /** LINEAR -- size of strike vector should equal nb of columns of squeeze matrix */
        if (interpolationTypeInt == 2) {
            if (!(strikes->size() == squeezes->numCols())) {
                throw ModelException(method, "For piecewise linear squeezes, "
                "size of strikes has to equal nb of columns of squeeze matrix. \n"
                "This is not the case for category " + name + ".");
            }
        }

        /** expiries and nb rows of matrix have to coincide */
        if (!(expiries->size() == squeezes->numRows())) {
            throw ModelException(method, "Size of expiries has to equal nb of rows of squeeze matrix. \n" 
                "This is not the case for category " + name + ".");
        }        
        
        /** check that strikes are (strictly) increasing */        
        for (int iStrike=1; iStrike<strikes->size(); iStrike++) {
            if ( !Maths::isPositive( (*strikes)[iStrike]-(*strikes)[iStrike-1] )) {
                throw ModelException(method, "Strikes must be strictly increasing, but this "
                    "is not the case for region " + name + ".");
            }
        }

        /** check that strikes are in reasonable range */
        if (Maths::isPositive(strikes->back() - 5.0)) {
            throw ModelException(method, "Maximum strike should not be higher than +5.0, but amounts to "
                + Format::toString(strikes->back()) 
                + " in region " 
                + name + ".");
        }
        if (Maths::isNegative(strikes->front() + 5.0)) {
            throw ModelException(method, "Minimum strike should not be lower than -5.0, but amounts to "
                + Format::toString(strikes->front()) 
                + " in region " 
                + name + ".");
        }
        
        /** check that squeezes are between -1 and +1 */
        minSqueeze =  1.0;
        maxSqueeze =  -1.0;        
        for (int iCol=0; iCol<squeezes->numCols(); iCol++) {            
            for (int iRow=0; iRow<squeezes->numRows(); iRow++) {
                if (!Maths::isNegative( abs((*squeezes)[iCol][iRow]) - 1.0 )) {
                    throw ModelException(method, 
                        "Squeezes must be between -1 and +1, but squeeze in row " 
                        + Format::toString(iRow) 
                        + " and column " 
                        + Format::toString(iCol) 
                        + " equals " 
                        + Format::toString((*squeezes)[iCol][iRow])
                        + " for region " + name + ".");                
                }            
                if ( Maths::isNegative( (*squeezes)[iCol][iRow] - minSqueeze) ) {
                    minSqueeze = (*squeezes)[iCol][iRow];
                }
                if ( Maths::isPositive( (*squeezes)[iCol][iRow] - maxSqueeze) ) {
                    maxSqueeze = (*squeezes)[iCol][iRow];
                }
            }
        }      
        
        /** set up array of strike indices used for sensitivities */
        if (interpolationTypeInt == 1) { // STEP
            strikeIndices = DoubleArraySP(new DoubleArray(squeezes->numCols())); // copy strikes
            for (int iStrike = 0; iStrike < squeezes->numCols(); iStrike++) {
                (*strikeIndices)[iStrike] = double(iStrike);
            }            
        } else { // interpolationTypeInt == 2 LINEAR
            strikeIndices = DoubleArraySP(new DoubleArray(*strikes));
        }        
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

string LocalCorrSqueeze::getName() const {
    return name;
}

double LocalCorrSqueeze::getMinSqueeze() const {
    return minSqueeze;
}

double LocalCorrSqueeze::getMaxSqueeze() const {
    return maxSqueeze;
}

bool LocalCorrSqueeze::isZeroObject() const {
    return isZeroObj;
}

double LocalCorrSqueeze::computeSqueeze(double marketFactor,
                                        double tradingTime) const {
    static const string method("LocalCorrSqueeze::getSqueeze");    
    try {         
        if (isZeroObj) {
            return 0.0;
        } 
        /** find index for interpolator */
        static unsigned long timeIdx;
        locate(&(*expiriesInTradingTime)[0]-1, expiriesInTradingTime->size(), tradingTime, &timeIdx);            
        if (expiries->size() == timeIdx) {
            timeIdx -= 1;
        }
        if (interpolationTypeInt == 1) {
            static unsigned long position;        
            locate(&(*strikes)[0]-1, strikes->size(), marketFactor, &position);
            return (*squeezes)[position][timeIdx];        
        } else { // if (interpolationTypeInt == 2) 
            return linearInterpSqueeze[timeIdx]->value(marketFactor); 
        } 
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

/** override default implementation and set up linear interpolator for squeezes (LINEAR) */
void LocalCorrSqueeze::getMarket(const IModel* model, const MarketData* market) {
    static const string method("LocalCorrSqueeze::getMarket");
    try {
        /** convert array of expiries into array of (dummy) trading time */
        /** get ref date from market */
    	DateTime refDate;
	    market->GetReferenceDate(refDate);
        /** check that expiries are increasing */
        if (!Expiry::isIncreasing(refDate, expiries.get())) {
            throw ModelException(method, "Expiries for LocalCorrSqueeze object " + name + " are not increasing!");
        }
        /** set up very basic dummy time metric */
        string name = "Holiday";
	    const HolidaySP holiday(new Holiday(name, DateTimeArray(0), false));
	    TimeMetricSP metric(new TimeMetric(0.05, holiday.get()));
        /** carry out conversion */
	    expiriesInTradingTime = DoubleArraySP(new DoubleArray(expiries->size()));
	    for (int bm = 0; bm < expiries->size(); bm++) {
		    DateTime date = (*expiries)[bm]->toDate(refDate); 		    
		    (*expiriesInTradingTime)[bm] = metric->yearFrac(refDate,date);
	    }
       /** set up array of interpolators */
	   if (interpolationTypeInt == 2) {
           linearInterpSqueeze = InterpolantArray(expiries->size());
           for (int iExpiry=0; iExpiry<expiries->size(); iExpiry++) {
               DoubleArray squeezeArray(squeezes->numCols());
                for (int i=0; i<squeezes->numCols(); i++) {
                    squeezeArray[i] = (*squeezes)[i][iExpiry];
                }                
                linearInterpSqueeze[iExpiry] = LinearInterpolantSP(new 
		            LinearInterpolant(*strikes.get(), squeezeArray));
           }           
        } else {
            linearInterpSqueeze = InterpolantArray(0);
        }        
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

                      
LocalCorrSqueeze::LocalCorrSqueeze(): 
MarketObject(TYPE), 
isZeroObj(false),
name(""),
minSqueeze(0.0),
maxSqueeze(0.0) {}

IObject* LocalCorrSqueeze::defaultLocalCorrSqueeze() {
    return new LocalCorrSqueeze();
}

/** Support for MATRIX sensitivity */
string LocalCorrSqueeze::sensName(const LocalCorrExpiryAndStrikewise*) const {
    return name; 
}

ExpiryAndStrikeArrayConstSP LocalCorrSqueeze::sensQualifiers(const LocalCorrExpiryAndStrikewise*) const {
    return ExpiryAndStrike::series(expiries, *strikeIndices);
}

TweakOutcome LocalCorrSqueeze::sensShift(const PropertyTweak<LocalCorrExpiryAndStrikewise>& tweak) {
    static const string method = "LocalCorrSqueeze::sensShift";
    try {
        /** do nthg if it is a dummy object */
        if (isZeroObj) {
            return TweakOutcome(tweak.coefficient, false); 
        }
        /** find the idx for the expiry in question */
        int iExpiry = 0;
        for (iExpiry=0; iExpiry<expiries->size(); iExpiry++) {
            if ( (*expiries)[iExpiry]->equals(tweak.qualifier->expiry.get()) ) {
                break;
            }
        }
        /** find the idx for the strike in question */
        int iStrike = 0;
        for (iStrike = 0; iStrike<strikeIndices->size(); iStrike++) {
            if ( Maths::isZero((*strikeIndices)[iStrike]-tweak.qualifier->strike) ) {
                break;
            }
        }
        /** shift the corresponding squeeze and put a collar on top */
        if (!Maths::isZero(tweak.coefficient)) {
            (*squeezes)[iStrike][iExpiry] += tweak.coefficient;
            (*squeezes)[iStrike][iExpiry] = min(max((*squeezes)[iStrike][iExpiry], -1.0),1.0); // COLLAR
        }        

        if (interpolationTypeInt == 2) {           
            DoubleArray squeezeArray(squeezes->numCols());
            for (int i=0; i<squeezes->numCols(); i++) {
                squeezeArray[i] = (*squeezes)[i][iExpiry];
            }                
            linearInterpSqueeze[iExpiry] = LinearInterpolantSP(new 
		        LinearInterpolant(*strikes.get(), squeezeArray));                   
        }
        return TweakOutcome(tweak.coefficient, false); 
    } catch (exception &e) {
        throw ModelException(e, method, "LocalCorrSqueezeExpiryAndStrikewise tweak failed for " + getName());
    }
}


/** Support for POINTWISE sensitivity */
string LocalCorrSqueeze::sensName(const LocalCorrExpiry*) const {
    return name; 
}

ExpiryWindowArrayConstSP LocalCorrSqueeze::sensQualifiers(const LocalCorrExpiry*) const{
    return ExpiryWindow::series(expiries);
}

TweakOutcome LocalCorrSqueeze::sensShift(const PropertyTweak<LocalCorrExpiry>& tweak) {
    static const string method = "LocalCorrSqueeze::sensShift";
    try {
        ASSERT(!!tweak.tag);
        /** do nthg if it is a dummy object */
        if (isZeroObj) {
            return TweakOutcome(tweak.coefficient, false); 
        }
        /** find the idx for the expiry in question */
        int iExpiry = 0; 
        for (iExpiry=0; iExpiry<expiries->size(); iExpiry++) {
            if ( (*expiries)[iExpiry]->equals(tweak.qualifier->expiry.get()) ) {
                break;
            }
        }
        /** shift the corresponding squeezes and put a collar on top */        
        if (!Maths::isZero(tweak.coefficient)) {
            if (tweak.tag->op == LocalCorrExpiry::PARALLEL_SHIFT)  {
                for (int iStrike=0; iStrike<squeezes->numCols(); iStrike++) {
                    (*squeezes)[iStrike][iExpiry] += tweak.coefficient;
                    (*squeezes)[iStrike][iExpiry] = min(max((*squeezes)[iStrike][iExpiry], -1.0),1.0); // COLLAR
                }
            } else {
                if (interpolationTypeInt == 1) { // STEP
                    bool passZero = false;                    
                    for (int iStrike = 0; iStrike<squeezes->numCols(); iStrike++) {                        
                        if (Maths::isNegative((*strikes)[iStrike])) {
                            (*squeezes)[iStrike][iExpiry] -= tweak.coefficient * (*strikes)[iStrike];
                        } else if (!Maths::isNegative((*strikes)[iStrike]) && !passZero) {
                            // do nthg for two in a row
                            passZero = true;
                            iStrike++;
                        } else { 
                            (*squeezes)[iStrike][iExpiry] -= tweak.coefficient * (*strikes)[iStrike-1];
                        }                        
                    } 
                } else { // (interpolationType == 2) LINEAR
                    for (int iStrike = 0; iStrike<squeezes->numCols(); iStrike++) {
                        (*squeezes)[iStrike][iExpiry] -= tweak.coefficient * (*strikes)[iStrike];
                    }
                }
            }
        }
        // set up new interpolant
        if (interpolationTypeInt == 2) {
            DoubleArray squeezeArray(squeezes->numCols());
            for (int i=0; i<squeezes->numCols(); i++) {
                squeezeArray[i] = (*squeezes)[i][iExpiry];
            }                
            linearInterpSqueeze[iExpiry] = LinearInterpolantSP(new 
		        LinearInterpolant(*strikes.get(), squeezeArray));
        }        
        return TweakOutcome(tweak.coefficient, false); 
    } catch (exception &e) {
        throw ModelException(e, method, "LocalCorrSqueezeExpiry tweak failed for " + getName());
    }
}


/** Support for PARALLEL sensitivity */
string LocalCorrSqueeze::sensName(const LocalCorrVoid*) const {
    return name; /** this is the name of the object */
}

TweakOutcome LocalCorrSqueeze::sensShift(const PropertyTweak<LocalCorrVoid>& tweak) {
    static const string method = "LocalCorrSqueeze::sensShift";
    try {
        ASSERT(!!tweak.tag);
        if (isZeroObj) {
            return TweakOutcome(tweak.coefficient, false); 
        }
        if (!Maths::isZero(tweak.coefficient)) {
            if (tweak.tag->op == LocalCorrVoid::PARALLEL_SHIFT) {
                for (int iSqueeze = 0; iSqueeze<squeezes->numCols(); iSqueeze++) {
                    for (int iExpiry = 0; iExpiry<squeezes->numRows(); iExpiry++) {
                        (*squeezes)[iSqueeze][iExpiry] += tweak.coefficient;
                    }
                }
            } else {
               if (interpolationTypeInt == 1) { // STEP
                    bool passZero = false;                    
                    for (int iStrike = 0; iStrike<squeezes->numCols(); iStrike++) {
                        for (int iExpiry = 0; iExpiry<squeezes->numRows(); iExpiry++) {
                            if (Maths::isNegative((*strikes)[iStrike])) {
                                (*squeezes)[iStrike][iExpiry] -= tweak.coefficient * (*strikes)[iStrike];
                            } else if (!Maths::isNegative((*strikes)[iStrike]) && !passZero) {
                                // do nthg for two in a row
                                passZero = true;
                                iStrike++;
                            } else { 
                                (*squeezes)[iStrike][iExpiry] -= tweak.coefficient * (*strikes)[iStrike-1];
                            }
                        }
                    } 
                } else { // (interpolationType == 2) LINEAR
                    for (int iStrike = 0; iStrike<squeezes->numCols(); iStrike++) {
                        for (int iExpiry = 0; iExpiry<squeezes->numRows(); iExpiry++) {
                            (*squeezes)[iStrike][iExpiry] -= tweak.coefficient * (*strikes)[iStrike];
                        }
                    }
                }
            }
        }
        
        // set up new interpolant
        if (interpolationTypeInt == 2) {           
           for (int iExpiry=0; iExpiry<expiries->size(); iExpiry++) {
               DoubleArray squeezeArray(squeezes->numCols());
                for (int i=0; i<squeezes->numCols(); i++) {
                    squeezeArray[i] = (*squeezes)[i][iExpiry];
                }                
                linearInterpSqueeze[iExpiry] = LinearInterpolantSP(new 
		            LinearInterpolant(*strikes.get(), squeezeArray));
           }           
        }
        return TweakOutcome(tweak.coefficient, false); 
    } catch (exception &e) {
        throw ModelException(e, method, "LocalCorrSqueezeVoid tweak failed for " + getName());
    }
}

LocalCorrSqueeze::LocalCorrSqueeze(bool isZeroObj): 
MarketObject(TYPE), 
isZeroObj(isZeroObj), 
name(""),
minSqueeze(0.0), 
maxSqueeze(0.0) {}

void LocalCorrSqueeze::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(LocalCorrSqueeze, clazz);
    SUPERCLASS(MarketObject);    
    IMPLEMENTS(ITweakableWithRespectTo<LocalCorrExpiryAndStrikewise>);
    IMPLEMENTS(ITweakableWithRespectTo<LocalCorrExpiry>);
    IMPLEMENTS(ITweakableWithRespectTo<LocalCorrVoid>);
    EMPTY_SHELL_METHOD(defaultLocalCorrSqueeze);
    FIELD(name, "name of region");    
    FIELD(expiries, "vector of expiries");    
    FIELD(strikes, "vector of  strikes");    
    FIELD(squeezes, "vector of vector of squeezes");    
    FIELD(interpolationType, "string identifyer for interpolation type");    
    FIELD(minSqueeze, "min squeeze used");
    FIELD_MAKE_TRANSIENT(minSqueeze);
    FIELD(maxSqueeze, "max squeeze used");
    FIELD_MAKE_TRANSIENT(maxSqueeze);
    FIELD(isZeroObj, "is zero obj");
    FIELD_MAKE_TRANSIENT(isZeroObj);
    FIELD(expiriesInTradingTime, "");
    FIELD_MAKE_TRANSIENT(expiriesInTradingTime);
    FIELD(strikeIndices, "");
    FIELD_MAKE_TRANSIENT(strikeIndices);
    FIELD(interpolationTypeInt, "");
    FIELD_MAKE_TRANSIENT(interpolationTypeInt);
    FIELD(linearInterpSqueeze, "");
    FIELD_MAKE_TRANSIENT(linearInterpSqueeze);
    MarketDataFetcher::setDefaultRetrievalMode(LocalCorrSqueeze::TYPE, false, NULL);
}

DEFINE_TEMPLATE_TYPE(LocalCorrSqueezeArray);

DRLIB_END_NAMESPACE
#endif // LOC_CORR_SQUEEZE_CPP

