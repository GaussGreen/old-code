//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CorrelationTerm.cpp
//
//   Description : Holds implied correlation parameters for term structure
//
//   Author      : Oliver Brockhaus
//
//   Date        : 05 Dec 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_CORRELATIONTERM_CPP
#include "edginc/CorrelationTerm.hpp"
#include "edginc/Maths.hpp"
#include "edginc/MaturityTimePeriod.hpp"
#include "edginc/Addin.hpp"
#include "edginc/CorrelationCategory.hpp"
#include "edginc/Format.hpp"
#include "edginc/NonPricingModel.hpp"
#include "edginc/MarketDataFetcher.hpp"
#include "edginc/MDFUtil.hpp"
#include "edginc/MarketData.hpp"

DRLIB_BEGIN_NAMESPACE

// place holder for corr term data
CorrTermData::CorrTermData(ExpiryArraySP    benchmarkTermExpiries,
                           ExpiryArraySP    shortTermExpiries,
                           ExpiryArraySP    longTermExpiries,
                           CDoubleMatrixSP  correlationMatrix,
                           CDoubleMatrixSP  shortTermSqueeze,
                           CDoubleMatrixSP  longTermSqueeze) :
benchmarkTermExpiries(benchmarkTermExpiries), 
shortTermExpiries(shortTermExpiries), 
longTermExpiries(longTermExpiries), 
correlationMatrix(correlationMatrix), 
shortTermSqueeze(shortTermSqueeze),
longTermSqueeze(longTermSqueeze) {}


CorrelationTerm::~CorrelationTerm(){}

const string CorrelationTerm::SHORT_TERM = "1M";
const string CorrelationTerm::LONG_TERM = "10Y";
const int CorrelationTerm::TERM_TIME = DateTime::END_OF_DAY_TIME;
const double CorrelationTerm::EIGEN_VALUE_FLOOR = 0.0000000001;
const double CorrelationTerm::MAX_SQ_ERROR = 0.0025;

/** Validation */
void CorrelationTerm::validatePop2Object(){
    static const string method("CorrelationTerm::validatePop2Object");
    try {
        if (name.empty()){
            name = category1 < category2? (category1+ "_" + category2): (category2 + "_" + category1);
        } else if (name == category1 || name == category2){
                throw ModelException(method, "Correlation Category's name must be different to"
                                        " either category's name");
        }    
        if ( Maths::isPositive(fabs(corrShortTermSqueeze) - 1.0) ){
            throw ModelException(method, "corrShortTermSqueeze "
                + Format::toString(corrShortTermSqueeze)
                + " for "
                + name
                + " must be between -1.0 and 1.0");
        }
        if ( Maths::isPositive(fabs(corrLongTermSqueeze) - 1.0) ){
            throw ModelException(method, "corrLongTermSqueeze "
                + Format::toString(corrLongTermSqueeze)
                + " for " 
                + name 
                + " must be between -1.0 and 1.0");
        }
        if (corrShortTermExpiry->equals(corrLongTermExpiry.get())) {
            throw ModelException(method, "expiries for "+name+
                                    " must be different");
        }
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

void CorrelationTerm::getMarket(const IModel* model, const MarketData* market) {
    static const string method("CorrelationTerm::getMarket");
    try {
        DateTime valueDate;
        market->GetReferenceDate(valueDate);
        
        ExpiryArraySP checkExpiries(new ExpiryArray(1, corrShortTermExpiry));
        checkExpiries->push_back(corrLongTermExpiry);
        if (!Expiry::isIncreasing(valueDate, checkExpiries.get())) {
            throw ModelException(method, "short term expiry "
                + corrShortTermExpiry->toDate(valueDate).toString() 
                + " must be strictly less than long term expiry "
                + corrLongTermExpiry->toDate(valueDate).toString() 
                + " for correlation term object " + this->getName());
        }

    } catch (exception& e) {
        throw ModelException(e, method);
    }
}


/** Initialises this piece of market data - 
    records the pair of names idenitfying the correlation */
void CorrelationTerm::initialise(MarketData* market){
    string marketCorrelationTermName;
    bool   hasCorrelationTerm = market->hasCorrelationTermData(category1, category2);

    if ( hasCorrelationTerm ) {
         marketCorrelationTermName = market->getCorrelationTermName(category1, category2);
    }

    if ( !hasCorrelationTerm || marketCorrelationTermName == name ) {
        market->setCorrelationTermName(category1, category2, name);
    } else {
        throw ModelException("CorrelationTerm::initialise", 
                 "There is already a correlationTerm in the cache for assets " 
                 + category1 + " and " + category2 + " with name " + marketCorrelationTermName + ".\n"
                 " Cannot add a correlationTerm for the same categorys with a different name (" + name + ").");
    }
}

double CorrelationTerm::getCorrShortTermSqueeze() const {
    return corrShortTermSqueeze;
}

ExpirySP CorrelationTerm::getShortTermExpiry() const {
    return corrShortTermExpiry;
}

double CorrelationTerm::getCorrLongTermSqueeze() const {
    return corrLongTermSqueeze;
}

ExpirySP CorrelationTerm::getLongTermExpiry() const {
    return corrLongTermExpiry;
}

// Sensitivity methods

/** Returns the name of the protected equity - used to determine
    whether to tweak the object */
string CorrelationTerm::sensName(ShortTermSqueezeTweak* shift) const {
    return getName();
}

string CorrelationTerm::sensName(LongTermSqueezeTweak* shift) const {
    return getName();
}

/** Shifts the object using given shift */    
bool CorrelationTerm::sensShift(ShortTermSqueezeTweak* shift, 
                                bool useShiftSizeSign) {
    /** dont tweak zero objects */
    if (isZeroObject) {
        return false;
    }
    // must set the intitial value
    shift->setInitialValue(corrShortTermSqueeze);
    
    double shiftSize = shift->getShiftSize();
 
    // then shift (zero squeeze is positive squeeze)
    if (!useShiftSizeSign || Maths::isNegative(corrShortTermSqueeze)) {
        corrShortTermSqueeze += shiftSize; 
    } else {
        corrShortTermSqueeze -= shiftSize;
    }
    corrShortTermSqueeze = min(max(corrShortTermSqueeze, -1.0),1.0); // COLLAR
    return false; // none of our components has a SHORT_TERM_SQUEEZE_TWEAK type sensitivity
}

bool CorrelationTerm::sensShift(ShortTermSqueezeTweak* shift) {
    return sensShift(shift, shift->overrideShiftSizeSign());
}

bool CorrelationTerm::sensShift(LongTermSqueezeTweak* shift,
                                bool useShiftSizeSign) {

    /** dont tweak zero objects */
    if (isZeroObject) {
        return false;
    }
    // must set the intitial value
    shift->setInitialValue(corrLongTermSqueeze);
    
    double shiftSize = shift->getShiftSize();
 
    // then shift (zero squeeze is positive squeeze)
    if (!useShiftSizeSign || Maths::isNegative(corrLongTermSqueeze)) {
        corrLongTermSqueeze += shiftSize; 
    } else {
        corrLongTermSqueeze -= shiftSize;
    }
    corrLongTermSqueeze = min(max(corrLongTermSqueeze, -1.0),1.0); // COLLAR
    return false; // none of our components has a LONG_TERM_SQUEEZE_TWEAK type sensitivity
}

bool CorrelationTerm::sensShift(LongTermSqueezeTweak* shift) {
    return sensShift(shift, shift->overrideShiftSizeSign());
}


string CorrelationTerm::getName() const {
    return name;
}

CorrelationTerm::CorrelationTerm() : MarketObject(TYPE),
    corrShortTermExpiry(ExpirySP(new MaturityTimePeriod(SHORT_TERM, TERM_TIME))),
    corrLongTermExpiry(ExpirySP(new MaturityTimePeriod(LONG_TERM, TERM_TIME))),
    isZeroObject(false) {}

/** constructor -- in order to create dummy correlationTerm object */
CorrelationTerm::CorrelationTerm(bool isZeroObject) : MarketObject(TYPE), 
    name(""),  category1(""), category2(""), 
    corrShortTermSqueeze(0.0), 
    corrShortTermExpiry(ExpirySP(new MaturityTimePeriod(SHORT_TERM, TERM_TIME))),
    corrLongTermSqueeze(0.0), 
    corrLongTermExpiry(ExpirySP(new MaturityTimePeriod(LONG_TERM, TERM_TIME))),
    isZeroObject(isZeroObject) {}

class CorrelationTermHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CorrelationTerm, clazz);
        SUPERCLASS(MarketObject);
        IMPLEMENTS(ShortTermSqueezeTweak::IShift);
        IMPLEMENTS(LongTermSqueezeTweak::IShift);
        EMPTY_SHELL_METHOD(defaultCorrelationTerm);
        FIELD(name, "name for correlation term");
        FIELD_MAKE_OPTIONAL(name);
        FIELD(category1,      "Name of Category 1");
        FIELD(category2,      "Name of Category 2");
        FIELD(corrShortTermSqueeze, "correlation short squeeze ");
        FIELD(corrLongTermSqueeze, "correlation long squeeze ");
        FIELD(corrShortTermExpiry, "correlation short date");
        FIELD_MAKE_OPTIONAL(corrShortTermExpiry);
        FIELD(corrLongTermExpiry, "correlation long date");
        FIELD_MAKE_OPTIONAL(corrLongTermExpiry);
        FIELD(isZeroObject, "whether or not its a zero object");
        FIELD_MAKE_TRANSIENT(isZeroObject);
        // by default don't get CorrelationTerm
        MarketDataFetcher::setDefaultRetrievalMode(CorrelationTerm::TYPE,
                                                   false, NULL);
    }

    static IObject* defaultCorrelationTerm(){
        return new CorrelationTerm();
    }
};

CClassConstSP const CorrelationTerm::TYPE = CClass::registerClassLoadMethod(
    "CorrelationTerm", typeid(CorrelationTerm), CorrelationTermHelper::load);

//support for arrays
DEFINE_TEMPLATE_TYPE(CorrelationTermArray);

///// support for wrapper class /////

// definition of TYPE for MarketWrapper template class
DEFINE_TEMPLATE_TYPE(CorrelationTermWrapper);

/** Casts array element to an IObject */
IObjectConstSP arrayObjectCast<CorrelationTermWrapper>::toIObject(
    const CorrelationTermWrapper& value){
    IObjectConstSP objValue(IObjectConstSP::attachToRef(&value));
    return objValue;
}

/** Casts array element to an IObject */
IObjectSP arrayObjectCast<CorrelationTermWrapper>::toIObject(CorrelationTermWrapper& value){
    return IObjectSP::attachToRef(&value);
}

/** Turns the IObjectSP into a CorrelationTerm */
CorrelationTermWrapper arrayObjectCast<CorrelationTermWrapper>::fromIObject(IObjectSP& value){
    CorrelationTermWrapper *ptr = dynamic_cast<CorrelationTermWrapper *>(value.get());
    if (!ptr){
        throw ModelException("arrayObjectCast::fromIObject", "Object is not"
                             " a CorrelationTermWrapper but a "+
                             value->getClass()->getName());
    }
    return *ptr;
}

DEFINE_TEMPLATE_TYPE(CorrelationTermWrapperArray);

bool CorrelationTermLoad() {
    return (CorrelationTerm::TYPE != 0);
}

/** Pull out the forward correlation matrix between two dates from the market data */
DoubleMatrixArraySP CorrelationTerm::CorrelationTermMatrix (
    const DateTime&             valueDate, 
    const DateTimeArray&        toDates,
    const DoubleMatrix&         fwdVarAtDates,    
    const TimeMetricArray&      timeMetricArray,
    CorrTermDataSP              data,
    const TCalcType             calcType,
    double                      eigenValueFloor,
    double                      maxSqError,
    const DoubleArray&          barriers,
    const BoolArray&            skipFwdCorrelation) {

    static const string method("CorrelationTerm::CorrelationTermMatrix");

    try{                 
        /** dimension checking I ... */
        int pos, iAsset, jAsset, iDate;
        int numAssets = data->correlationMatrix->numCols();
        int numDates = toDates.size();
        double sqError = 0.0;

        int numCombinations = (numAssets * numAssets - numAssets) / 2;
        /** temporary check, needs to be fixed in the long run -> in DMGT already */
        if(skipFwdCorrelation.size() != numCombinations) {
            throw ModelException(method, "Market Data mismatch, GaussTerm used, but Gauss provided");
        }

        if( (data->correlationMatrix->numRows() != numAssets)
            || (data->shortTermSqueeze->numRows() != numAssets)
            || (data->shortTermSqueeze->numCols() != numAssets)
            || (data->longTermSqueeze->numRows() != numAssets)
            || (data->longTermSqueeze->numCols() != numAssets) ) {
            throw ModelException(method,
                "Num rows correlation = " + 
                Format::toString(data->correlationMatrix->numRows()) +
                ", num rows corrShortTermSqueeze = " + 
                Format::toString(data->shortTermSqueeze->numRows()) +
                ",\n num cols corrShortTermSqueeze = " + 
                Format::toString(data->shortTermSqueeze->numCols()) +
                ", num rows corrLongTermSqueeze = " + 
                Format::toString(data->longTermSqueeze->numRows()) +
                ",\n num cols corrLongTermSqueeze = " + 
                Format::toString(data->longTermSqueeze->numCols()) +
                " should all be equal to" +
                " num cols correlation = " + Format::toString(numAssets) );
        }

        /** dimension checking II ... fwdVarAtDates */
        if( calcType == forward ) { // need vars only for forward correlation calculation
            if( fwdVarAtDates.numCols() != numAssets ) {
                throw ModelException(method,"Num variances = " 
                                            + Format::toString(fwdVarAtDates.numCols()) 
                                            + " different from num assets = "
                                            + Format::toString(numAssets) );
            }
            if ( fwdVarAtDates.numRows() != numDates ) {
                throw ModelException(method,"Num variances = " 
                                          + Format::toString(fwdVarAtDates.numRows()) 
                                          + " different from num dates = "
                                          + Format::toString(numDates) );
            }
        }

        /** SPOT CORRELATIONS (var not needed here) */
        DoubleMatrixArraySP spotResult = DoubleMatrixArraySP(new DoubleMatrixArray(numDates));
        DoubleMatrixArraySP fwdResult = DoubleMatrixArraySP(new DoubleMatrixArray(numDates));
        DoubleArray spotCorrelations(numDates);
        for( iDate=0; iDate<numDates; iDate++ ) {
            (*spotResult)[iDate] = CDoubleMatrixSP(new DoubleMatrix(numAssets,numAssets));
            (*fwdResult)[iDate] = CDoubleMatrixSP(new DoubleMatrix(numAssets,numAssets));
        }
        DoubleArrayArray timeMetricFraction((numAssets*numAssets-numAssets)/2);
        pos = 0;
        for( iAsset=0; iAsset<numAssets; iAsset++ ) {
            for( iDate=0; iDate<numDates; iDate++ ) {
                (*(*spotResult)[iDate])[iAsset][iAsset] = 1;
                (*(*fwdResult)[iDate])[iAsset][iAsset] = 1;
            }
            for( jAsset=iAsset+1; jAsset<numAssets; jAsset++, pos++ ) {
                spotCorrelations = spotCorrelation(
                    (*data->correlationMatrix)[iAsset][jAsset],
                    (*data->benchmarkTermExpiries)[pos],
                    (*data->shortTermSqueeze)[iAsset][jAsset],
                    (*data->shortTermExpiries)[pos],
                    (*data->longTermSqueeze)[iAsset][jAsset],
                    (*data->longTermExpiries)[pos],
                    valueDate,
                    toDates,
                    timeMetricArray[iAsset],
                    timeMetricArray[jAsset]);
                for( iDate=0; iDate<numDates; iDate++ ) { // fill matrix
                    (*(*spotResult)[iDate])[iAsset][jAsset] = spotCorrelations[iDate];
                    (*(*spotResult)[iDate])[jAsset][iAsset] = spotCorrelations[iDate];
                }
            }
        }        
        /** SPOT CORRELATIONS -- check whether or not positive definite */
        for( iDate = 0; iDate < numDates; iDate++ ) {
            sqError = 0.0;
            DoubleMatrix corrTerm = (*spotResult)[iDate]->symmToCorrel(&sqError, eigenValueFloor);
            if( Maths::isPositive(sqError - maxSqError) ) {
                throw ModelException(method, 
                    "Average distance to positive definite modification of correlation matrix is " 
                        + Format::toString(sqrt(sqError)) + " on " 
                        + toDates[iDate].toString() 
                        + " and should not be more than "
                        + Format::toString(sqrt(maxSqError)) + ".");
            }
            (*spotResult)[iDate] = CDoubleMatrixSP(copy(&corrTerm));
        }

        /** FORWARD CORRELATIONS -- attention: vars are _fwd_ vars */
        if( calcType == forward ) {
            /** toDates[0]: fwd corr = spot corr, just check admissibility */ 
            DoubleMatrix corrTerm = (*spotResult)[0]->symmToCorrel(&sqError, eigenValueFloor);
            (*fwdResult)[0] = CDoubleMatrixSP(copy(&corrTerm));
        
            DoubleArray thisTotalVar(numAssets);
            DoubleArray thisTotalSqrtVar(numAssets);
            DoubleArray prevTotalSqrtVar(numAssets);
            BoolArray   skipCalcFwdCorr(numAssets);

             /** initialization */
            for (iAsset=0; iAsset<numAssets; iAsset++) {
                thisTotalVar[iAsset] = fwdVarAtDates[iAsset][0];
            }

            for( iDate=1; iDate<numDates; iDate++ ) {                
                /** vols */                
                for( iAsset=0; iAsset<numAssets; iAsset++ ) {                
                    // floor at zero, but should not be necessary
                    prevTotalSqrtVar[iAsset] = sqrt(max(thisTotalVar[iAsset],0.0));
                    thisTotalVar[iAsset] += fwdVarAtDates[iAsset][iDate]; 
                    thisTotalSqrtVar[iAsset] = sqrt(max(thisTotalVar[iAsset],0.0));                    
                    int nbBusDaysDiff = timeMetricArray[iAsset]->getHolidays()->
                        businessDaysDiff(toDates[iDate-1],toDates[iDate]);
                    if (nbBusDaysDiff > 0) {
                        skipCalcFwdCorr[iAsset] = false; 
                    } else {
                        skipCalcFwdCorr[iAsset] = true; 
                    }
                }
                /** forward correlations */
                pos = 0;
                for( iAsset=0; iAsset<numAssets; iAsset++ ) {
                    for( jAsset=iAsset+1; jAsset<numAssets; jAsset++, pos++ ) {
                        double divisor = sqrt(max(fwdVarAtDates[iAsset][iDate],0.0))
                            * sqrt(max(fwdVarAtDates[jAsset][iDate],0.0));
                        if ( Maths::isZero(divisor) ||
                                skipCalcFwdCorr[iAsset] || 
                                skipCalcFwdCorr[jAsset] ) {
                            (*(*fwdResult)[iDate])[iAsset][jAsset] = 
                                (*(*fwdResult)[iDate-1])[iAsset][jAsset];
                        } else if (skipFwdCorrelation[pos] == true) {
                            (*(*fwdResult)[iDate])[iAsset][jAsset] = 
                                (*(*spotResult)[iDate])[iAsset][jAsset];
                        } else {
                            (*(*fwdResult)[iDate])[iAsset][jAsset] = 
                                ( thisTotalSqrtVar[iAsset] * thisTotalSqrtVar[jAsset] * (*(*spotResult)[iDate])[iAsset][jAsset] 
                                    - prevTotalSqrtVar[iAsset] * prevTotalSqrtVar[jAsset] * (*(*spotResult)[iDate-1])[iAsset][jAsset] )
                                / divisor;
                        }
                        (*(*fwdResult)[iDate])[jAsset][iAsset]
                            = (*(*fwdResult)[iDate])[iAsset][jAsset];
                    }
                }

                // MODIFY CORRS
                sqError = 0.0;
                corrTerm = (*fwdResult)[iDate]->symmToCorrel(&sqError, eigenValueFloor);
                if (barriers.size()<2)
                    throw ModelException(method, "The array of barriers is empty.");

                if( Maths::isPositive(sqError - barriers[0]) &&
                    Maths::isPositive(barriers[1] - sqError)) {
                    corrTerm = *(*fwdResult)[iDate-1];
                } else if( Maths::isPositive(sqError - barriers[1]) ) {
                    throw ModelException(method, 
                        "Average distance to positive definite modification of correlation matrix is " 
                                                 + Format::toString(sqrt(sqError)) + " on " 
                                                 + toDates[iDate].toString() 
                                                 + " and should not be more than "
                                                 + Format::toString(sqrt(barriers[1])) + ".");
                }
                (*fwdResult)[iDate] = CDoubleMatrixSP(copy(&corrTerm)); // upper triangular matrix
            }
            return fwdResult;          
        } else {
            return spotResult;
        }
    } catch (exception& e){
        throw ModelException(e, "CorrelationTerm::CorrelationTermMatrix");
    }
}

/** Help function for CorrelationTermMatrix: parametrisation of term structure */
DoubleArray CorrelationTerm::spotCorrelation(
    double                      refCorrelation,
    ExpirySP                    corrBenchmarkExpiry,
    double                      corrShortTermSqueeze,
    ExpirySP                    corrShortTermExpiry,
    double                      corrLongTermSqueeze,
    ExpirySP                    corrLongTermExpiry,
    const DateTime&             valueDate,
    const DateTimeArray&        toDates,
    const TimeMetricConstSP&    timeMetricOne,
    const TimeMetricConstSP&    timeMetricTwo) {

    static const string method("CorrelationTerm::spotCorrelation");

    // todo: in case of error, somehow return asset resp category info

    try{
        int iDate = 0;
        int numDates = toDates.size();
        DoubleArray spotCorrelations(numDates);

        double shortCorr = refCorrelation 
            + corrShortTermSqueeze * ( Maths::isNegative(corrShortTermSqueeze) ?
            Maths::max(0.0,refCorrelation) : (1.0 - refCorrelation) );
        double longCorr = refCorrelation 
            + corrLongTermSqueeze * ( Maths::isNegative(corrLongTermSqueeze) ?
            Maths::max(0.0,refCorrelation) : (1.0 - refCorrelation) );

        /** four cases: 1) t[0] = t[1] = t[2]   2) t[0] = t[1] < t[2]
                        3) t[0] < t[1] = t[2]   4) t[0] < t[1] < t[2] */

        DoubleArray t(3);
        t[0] = sqrt(timeMetricOne->yearFrac(valueDate, corrShortTermExpiry->toDate(valueDate))
            * timeMetricTwo->yearFrac(valueDate, corrShortTermExpiry->toDate(valueDate)));
        t[1] = sqrt(timeMetricOne->yearFrac(valueDate, corrBenchmarkExpiry->toDate(valueDate))
            * timeMetricTwo->yearFrac(valueDate, corrBenchmarkExpiry->toDate(valueDate)));
        t[2] = sqrt(timeMetricOne->yearFrac(valueDate, corrLongTermExpiry->toDate(valueDate))
            * timeMetricTwo->yearFrac(valueDate, corrLongTermExpiry->toDate(valueDate)));
        
        if ( Maths::isZero(t[1]-t[0]) ) {
            throw ModelException(method, "Zero time fraction between short-term and reference.");
        }
        if ( Maths::isZero(t[2]-t[1]) ) {
            throw ModelException(method, "Zero time fraction between reference and long-term.");
        }

        /*if ( Maths::isZero(t[1]-t[0]) && Maths::isZero(t[2]-t[1]) ) {
            if ( !Maths::isZero(longCorr-refCorrelation)) {
                throw ModelException(method, "Zero time fraction between reference and long-term but non-zero squeeze factor.");
            } else if ( !Maths::isZero(refCorrelation-shortCorr) ) {
                throw ModelException(method, "Zero time fraction between short-term and reference but non-zero squeeze factor.");
            }
            // just assign refCorrelation to each date
            for( iDate=0; iDate<numDates; iDate++ ) {
                spotCorrelations[iDate] = refCorrelation; 
            }
        } else if ( Maths::isZero(t[1]-t[0]) && Maths::isPositive(t[2]-t[1]) ) {
            if (!Maths::isZero(refCorrelation-shortCorr)) {
                throw ModelException(method, "Zero time fraction between short-term and reference but non-zero squeeze factor.");
            }
            double  nu = log(2.0) / (t[2]-t[0]); // mu = \infty, mu/(mu-nu) = 1
            CDoubleMatrixSP rhoTransform = CDoubleMatrixSP(new DoubleMatrix(2,2));
            for (iDate=0; iDate<2; iDate++) {
                (*rhoTransform)[iDate][0] = exp(-nu*t[iDate]);
                (*rhoTransform)[iDate][1] = 1.0 - exp(-nu*t[iDate]);
            }
            
            DoubleArray rhoIn(2);
            rhoIn[0] = shortCorr; // equals as well refCorrelation
            rhoIn[1] = longCorr; 
            DoubleArray rhoOut = rhoTransform->solve(rhoIn); // rhoOut[0] = \theta_0, rhoOut[1] = \rho_\infty

            for( iDate=0; iDate<numDates; iDate++ ) {
                double time = sqrt(timeMetricOne->yearFrac(valueDate,toDates[iDate]) 
                    * timeMetricTwo->yearFrac(valueDate,toDates[iDate])); 
                spotCorrelations[iDate] = 
                    Maths::isNegative(time) ? rhoOut[0] : // limit case for time = 0 equals \hat \theta_0
                        rhoOut[1] + exp(-nu*time) * (rhoOut[0] - rhoOut[1]); 
                spotCorrelations[iDate] = Maths::max(Maths::min(spotCorrelations[iDate],1.0),-1.0);
            }
        } else if ( Maths::isPositive(t[1]-t[0]) && Maths::isZero(t[2]-t[1]) ) {
            if (!Maths::isZero(longCorr-refCorrelation)) {
                throw ModelException(method, "Zero time fraction between reference and long-term but non-zero squeeze factor.");
            }
            // NO SOLUTION
        }*/
        
        double  mu = log(2.0) / (t[1]-t[0]), nu = log(2.0) / (t[2]-t[0]);
        CDoubleMatrixSP rhoTransform = CDoubleMatrixSP(new DoubleMatrix(3,3));
        for( iDate=0; iDate<3; iDate++ ) {
            (*rhoTransform)[iDate][0] = exp(-mu*t[iDate]);
            (*rhoTransform)[iDate][1] = mu / (mu-nu) * (exp(-nu*t[iDate]) - (*rhoTransform)[iDate][0]);
            (*rhoTransform)[iDate][2] = 1.0 - (*rhoTransform)[iDate][0] - (*rhoTransform)[iDate][1];
        }

        DoubleArray rhoIn(3);
        rhoIn[0] = shortCorr;
        rhoIn[1] = refCorrelation;
        rhoIn[2] = longCorr;
        DoubleArray rhoOut = rhoTransform->solve(rhoIn);

        for( iDate=0; iDate<numDates; iDate++ ) {            
            double time = sqrt(timeMetricOne->yearFrac(valueDate,toDates[iDate])
                *timeMetricTwo->yearFrac(valueDate,toDates[iDate]));
            spotCorrelations[iDate] = 
                Maths::isNegative(time) ? rhoOut[0] : // limit case for time = 0 equals \hat \rho_0
                    rhoOut[2] + exp(-mu*time) * (rhoOut[0] - rhoOut[2]) 
                    + mu/(mu-nu) * (exp(-nu*time)-exp(-mu*time)) * (rhoOut[1] - rhoOut[2]);            
            spotCorrelations[iDate] = Maths::max(Maths::min(spotCorrelations[iDate],1.0),-1.0);
        }
        return spotCorrelations;

    } catch (exception& e){
        throw ModelException(e, "CorrelationTerm::spotCorrelation");
    }
}

/** Computes short and long term squeeze matrices of dimension numAssets 
    and populates expiry arrays from array of correlation term objects 
    => static method */
CorrTermDataSP CorrelationTerm::getCorrelationTermSqueezesAndExpiries (
                                     const DateTime&                refDate,                  
                                     int                            numAssets,
                                     const CorrelationCommonArray&  corrObjects,      
                                     const CorrelationTermArray&    corrTermArray) {
    static const string method("Correlation::getCorrelationTermSqueezesAndExpiries");
    try {
        /** allocate some space */
        CDoubleMatrixSP correlationMatrix(new DoubleMatrix(numAssets, numAssets));
        CDoubleMatrixSP shortTermSqueeze(new DoubleMatrix(numAssets, numAssets));
        CDoubleMatrixSP longTermSqueeze(new DoubleMatrix(numAssets, numAssets));
        int len = (numAssets * numAssets - numAssets)/2;
        
        ExpiryArraySP benchmarkTermExpiries(new ExpiryArray(len));
        ExpiryArraySP shortTermExpiries(new ExpiryArray(len));
        ExpiryArraySP longTermExpiries(new ExpiryArray(len));
        if (corrTermArray.size() > 0) {               
            int pos = 0;
            int corrPos = 0;
            for (int i = 0; i < numAssets; i++) {
                (*correlationMatrix)[i][i] = 1.0;
                for (int j = i+1; j < numAssets; j++, pos++, ++corrPos) {

                    if (!corrObjects[corrPos].get()) {
                        /** fill everything with zeros ... */
                        (*correlationMatrix)[i][j] = 0.0;
                        (*correlationMatrix)[j][i] = (*correlationMatrix)[i][j];
                    
                        (*shortTermSqueeze)[i][j] = 0.0;
                        (*shortTermSqueeze)[j][i] = (*shortTermSqueeze)[i][j];

                        (*longTermSqueeze)[i][j] = 0.0;
                        (*longTermSqueeze)[j][i] = (*longTermSqueeze)[i][j];

                        (*shortTermExpiries)[pos] = 
                            ExpirySP(new MaturityTimePeriod(SHORT_TERM, TERM_TIME));
                        (*benchmarkTermExpiries)[pos] = 
                            ExpirySP(new MaturityTimePeriod(Correlation::BENCHMARK_EXPIRY, TERM_TIME));
                        (*longTermExpiries)[pos] = 
                            ExpirySP(new MaturityTimePeriod(LONG_TERM, TERM_TIME));
                    } else {
                        // skip self-correlations
                        // ES: isnt that redundant, given the way corrObjects is constructed?
                        while (corrObjects[corrPos]->getAsset1Name() == corrObjects[corrPos]->getAsset2Name())
                            ++ corrPos;

                        if (! Correlation::TYPE->isInstance(corrObjects[corrPos].get()))
                        {
                            // TO DO: handle the case the correlation is post calibrated
                            throw ModelException(method, "correlation is calibrated");
                        }
                        CorrelationSP p = CorrelationSP::dynamicCast(corrObjects[corrPos]);

                        (*benchmarkTermExpiries)[pos] = p->getCorrExpiry(); 
                        (*shortTermExpiries)[pos] = corrTermArray[pos]->getShortTermExpiry();
                        (*longTermExpiries)[pos] = corrTermArray[pos]->getLongTermExpiry();

                        // another validation ... 
                        ExpiryArraySP checkExpiries(new ExpiryArray(1, (*shortTermExpiries)[pos]));
                        checkExpiries->push_back((*benchmarkTermExpiries)[pos]);
                        checkExpiries->push_back((*longTermExpiries)[pos]);

                        if (!Expiry::isIncreasing(refDate, checkExpiries.get())) {
                            throw ModelException(method, "benchmark term expiry "
                                + (*benchmarkTermExpiries)[pos]->toDate(refDate).toString() 
                                + " must be strictly greater than short term expiry "
                                + (*shortTermExpiries)[pos]->toDate(refDate).toString()
                                + " and strictly less than long term expiry " 
                                + (*longTermExpiries)[pos]->toDate(refDate).toString()
                                + " for " + corrObjects[pos]->getName());
                        }

                        if ((*shortTermExpiries)[pos]->equals((*benchmarkTermExpiries)[pos].get()) ||
                            (*longTermExpiries)[pos]->equals((*benchmarkTermExpiries)[pos].get()) ) {
                            throw ModelException(method, "Correlation benchmark expiry " +
                                (*benchmarkTermExpiries)[pos]->toString() + 
                                " of correlation " +
                                corrObjects[pos]->getName() +
                                " must not coincide with either " +
                                " of the two expiries of the respective correlation category " +
                                corrTermArray[pos]->getName() );
                        }

                        (*correlationMatrix)[i][j] = p->getCorrelation();
                        (*correlationMatrix)[j][i] = (*correlationMatrix)[i][j];
                    
                        (*shortTermSqueeze)[i][j] = corrTermArray[pos]->getCorrShortTermSqueeze();
                        (*shortTermSqueeze)[j][i] = (*shortTermSqueeze)[i][j];

                        (*longTermSqueeze)[i][j] = corrTermArray[pos]->getCorrLongTermSqueeze();
                        (*longTermSqueeze)[j][i] = (*longTermSqueeze)[i][j];
                    }
                }
            } 
        } else {
            shortTermExpiries = ExpiryArraySP(new ExpiryArray(len, 
                ExpirySP(new MaturityTimePeriod(SHORT_TERM, TERM_TIME))));
            benchmarkTermExpiries = ExpiryArraySP(new ExpiryArray(len, 
                ExpirySP(new MaturityTimePeriod(Correlation::BENCHMARK_EXPIRY, TERM_TIME))));
            longTermExpiries = ExpiryArraySP(new ExpiryArray(len, 
                ExpirySP(new MaturityTimePeriod(LONG_TERM, TERM_TIME))));
        }
        return CorrTermDataSP(new CorrTermData(benchmarkTermExpiries,
                                               shortTermExpiries,
                                               longTermExpiries,
                                               correlationMatrix,
                                               shortTermSqueeze,
                                               longTermSqueeze));
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

/* ComputeCorrTermStructureAddin */
class ComputeCorrTermStructureAddin: public CObject{
public:
    static CClassConstSP const TYPE;

   	MarketDataSP        market;
    StringArray         assets;
    DateTimeArray       toDates;    
    DateTimeSP          corrSwapMaturity;
    TimeMetricSP        timeMetric; // could potentially make optional and retrieve from market, if available

private:
	IObjectSP computeCorrTerm(){
        static const string method = "ComputeCorrTermStructureAddin::computeCorrTerm";
        try {
            int numAssets = assets.size();
                       
            StringArraySP resultCorrNames(new StringArray((numAssets * numAssets - numAssets)/2));
            DoubleArrayArraySP resultCorrTerm(new DoubleArrayArray((numAssets * numAssets - numAssets)/2));

            // create correlation matrix
            DoubleMatrix refCorrelation(numAssets, numAssets);
            CorrelationCommonArray corrObjects((numAssets * numAssets - numAssets)/2);
            CorrelationTermArray corrTermArray((numAssets * numAssets - numAssets)/2);
            ExpiryArray corrBenchmarkExpiries((numAssets * numAssets - numAssets)/2);

            int pos = 0;
            int i = 0;
            TimeMetricArray timeMetricArray(numAssets);

            string volType = "VolPreferred";
            MarketDataFetcherSP mdf(new MarketDataFetcher()); 
            if (corrSwapMaturity.get()) {                
                mdf->setCorrSwapExpiry(corrSwapMaturity);
            }
            MDFUtil::setCorrTermMode(*mdf, true);        
            NonPricingModel model(mdf);

            for (i = 0; i < numAssets; i++) {
                refCorrelation[i][i] = 1.0;
                
                MarketObjectSP mo1 = 
                    market->GetData(assets[i], CorrelationCategory::TYPE);
                CorrelationCategorySP category1 = CorrelationCategorySP::dynamicCast(mo1);

                for (int j = i + 1; j < numAssets; j++, pos++) {
                    string corrName = market->getCorrelationName(assets[i],assets[j]);
                    (*resultCorrNames)[pos] = corrName; // save name

                    MarketObjectSP corrObj(market->GetData(corrName, CorrelationBase::TYPE));
                    CorrelationBaseSP corrBase = CorrelationBaseSP::dynamicCast(corrObj);
                    corrObjects[pos] = CorrelationSP::dynamicCast(corrBase);
                    (*corrObjects[pos]).getMarket(&model,market.get());
                    if (! Correlation::TYPE->isInstance(corrObjects[pos].get()))
                    {
                        // TO DO: handle the case the correlation is post calibrated
                        throw ModelException(method, "correlation is calibrated");
                    }
                    CorrelationSP p = CorrelationSP::dynamicCast(corrObjects[pos]);
                    double corr = p->getCorrelation();
                    refCorrelation[i][j] = corr;
                    refCorrelation[j][i] = corr;
                    corrBenchmarkExpiries[pos] = p->getCorrExpiry();

                    MarketObjectSP mo2 = 
                        market->GetData(assets[j], CorrelationCategory::TYPE);
                    CorrelationCategorySP category2 = CorrelationCategorySP::dynamicCast(mo2);

                    const string& corrTermName = 
                        market->getCorrelationTermName(category1->getCategoryName(), category2->getCategoryName());
                    MarketObjectSP moTerm = market->GetData(corrTermName, CorrelationTerm::TYPE);
                    corrTermArray[pos] = CorrelationTermSP::dynamicCast(moTerm);
                }
                timeMetricArray[i] = timeMetric;
            }

            DoubleArrayArray timeMetricFraction((numAssets*numAssets-numAssets)/2);

            DateTime valueDate;
            market->GetReferenceDate(valueDate);
            
            CorrTermDataSP data = 
                CorrelationTerm::getCorrelationTermSqueezesAndExpiries(valueDate,
                                                                       numAssets,
                                                                       corrObjects, 
                                                                       corrTermArray);            

            // ultimately want to call spotCorrelation
            pos = 0; 
            for (i = 0; i < numAssets; i++) {
                for (int j = i + 1; j < numAssets; j++, pos++) {
                    DoubleArray spotCorrelations(toDates.size());
                    spotCorrelations = 
                        CorrelationTerm::spotCorrelation((*data->correlationMatrix)[i][j],
                                                         (*data->benchmarkTermExpiries)[pos],
                                                         (*data->shortTermSqueeze)[i][j],
                                                         (*data->shortTermExpiries)[pos],
                                                         (*data->longTermSqueeze)[i][j],
                                                         (*data->longTermExpiries)[pos],
                                                         valueDate,
                                                         toDates,
                                                         timeMetricArray[i],
                                                         timeMetricArray[j]);
                    (*resultCorrTerm)[pos] = spotCorrelations;
                }
            }
            ObjectArraySP result(new ObjectArray(2));
            (*result)[0] = resultCorrNames;
            (*result)[1] = resultCorrTerm;
            return result;

        } catch (exception& e){
            throw ModelException(e, method);
        }
    }

    static IObjectSP callComputeCorrTerm(ComputeCorrTermStructureAddin* params) {
        return params->computeCorrTerm();
    }

    /** constructor -> ugly, since hardcoded, needs to be modified */
    ComputeCorrTermStructureAddin() : CObject(TYPE) {
        if (!timeMetric.get()) {
            string name = "Holiday";
            const HolidaySP holiday(new Holiday(name, DateTimeArray(0), false));
            timeMetric = TimeMetricSP(new TimeMetric(0.05, holiday.get()));
        }
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ComputeCorrTermStructureAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultComputeCorrTermStructureAddin);
		FIELD(market, "");
        FIELD(assets, "");
        FIELD(toDates, "");
        FIELD(corrSwapMaturity, "");
        FIELD_MAKE_OPTIONAL(corrSwapMaturity);
        FIELD(timeMetric, "");
        FIELD_MAKE_OPTIONAL(timeMetric); // if in market, then not necessary to pass extra


        Addin::registerClassObjectMethod("COMPUTE_CORR_TERM_STRUCTURE",
                                          Addin::MARKET,
                                          "returns correlation term structure",
                                          TYPE,
                                          false,
                                          Addin::expandMulti,
                                          (Addin::ObjMethod*)callComputeCorrTerm);
    }

    static IObject* defaultComputeCorrTermStructureAddin(){
        return new ComputeCorrTermStructureAddin();
    }
};

CClassConstSP const ComputeCorrTermStructureAddin::TYPE =
CClass::registerClassLoadMethod("ComputeCorrTermStructureAddin", typeid(ComputeCorrTermStructureAddin), load);


DRLIB_END_NAMESPACE
