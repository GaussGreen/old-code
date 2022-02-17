//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PDFDefaultLNStrike.cpp
//
//   Description : Default Implementation of PDFCalculator for log normal vols
//
//   Author      : Andrew J Swain
//
//   Date        : 15 March 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/PDFDefaultLNStrike.hpp"
#include "edginc/Black.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

PDFDefaultLNStrike::~PDFDefaultLNStrike(){}

/** Populates supplied lostrikes and histrikes with strikes * (1.0 - epsilon)
    and strikes * (1.0 + epsilon) where epsilon is modified so that none 
    of the strikes overlap */
void PDFDefaultLNStrike::constructStrikes(
    const CSliceDouble&     strikes,
    double&                 epsilon, // note modified
    CSliceDouble&           lostrikes,
    CSliceDouble&           histrikes) {
    int n = strikes.size(); // number of strikes
    int i;
    for(i = 0; i < n; i++) {
        lostrikes[i] = strikes[i] * (1.0 - epsilon);
        histrikes[i] = strikes[i] * (1.0 + epsilon);
    }

    double thisOverlap, maxOverlap = 0.0;
    int j;
    for(i = 0; i < n-1; i++) {
        thisOverlap = lostrikes[i+1] - histrikes[i];
        if(Maths::isNegative(thisOverlap - maxOverlap)) {
            // record maxOverlap and its position
            maxOverlap = thisOverlap;
            j = i;
        }
    }
    
    // if they overlap adjust epsilon so that maxOverlap is 0.0
    if(Maths::isNegative(maxOverlap)) {
        epsilon = (strikes[j+1]-strikes[j])/(strikes[j+1]+strikes[j]);
        for(i = 0; i < n; i++) {
            lostrikes[i] = strikes[i] * (1.0 - epsilon);
            histrikes[i] = strikes[i] * (1.0 + epsilon);
        }
    }
}


/** validates that the strikes for specified maturity are ok */
void PDFDefaultLNStrike::validateStrikes(const CSliceDouble& strikes, 
                                         const DateTime&     maturity) {
    for (int i = 0; i < strikes.size(); i++) {
        if (Maths::isNegative(strikes[i]) ||
            (i > 0 && !Maths::isPositive(strikes[i] - strikes[i-1])))
            throw ModelException("validateStrikes", "Strikes for maturity " + 
                                 maturity.toString() + " are not increasing" +
                                 " (" + Format::toString(strikes[i-1])+") vs ("+ 
                                 Format::toString(strikes[i]) + ")");
    }
}

/** awaiting spec */
void PDFDefaultLNStrike::probabilities(const DoubleArray& strikes,
                                       const DateTime&    maturity,
                                       DoubleArray&       probs) const {
    // map through wider interface
    DoubleArray       strikesCopy(strikes); // to avoid const problems
    DateTimeArray     maturities(1, maturity);
    CSliceDouble      myStrikes(&strikesCopy[0], strikes.size());
    CSliceDouble      myProbs(&probs[0], probs.size());
    probabilities(myStrikes, maturities, myProbs);
}

/* same as next spreads method but takes in strikes as a DoubleArray.
   Just calls next spreads() method having mapped data types */
void PDFDefaultLNStrike::spreads(const DoubleArray& strikes,
                                 const DateTime&    maturity,
                                 double&            adjEpsilon,
                                 DoubleArray&       probs) const {
    // map through wider interface
    DoubleArray       strikesCopy(strikes); // to avoid const problems
    DateTimeArray     maturities(1, maturity);
    CSliceDouble      myStrikes(&strikesCopy[0], strikes.size());
    CSliceDouble      myProbs(&probs[0], probs.size());
    spreads(myStrikes, maturities, adjEpsilon, myProbs);
}

/** Caclulates spreads using specified strike(1+epsilon) and
    strike(1-epsilon) */
void PDFDefaultLNStrike::spreads(const CLatticeDouble&    strikes,
                                 const DateTimeArray&     maturities,
                                 double&                  adjEpsilon,
                                 CLatticeDouble&          spread) const {
    static const string method = "PDFDefaultLNStrike::spreads";
    try{
        // calculate forwards at maturities
        DoubleArray fwds(maturities.size());
        asset->fwdValue(maturities, fwds);
        // then construct two CLatticeDoubles of lo and hi strikes
        CLatticeDouble loStrikes(strikes.sizes());
        CLatticeDouble hiStrikes(strikes.sizes());
        // then loop through dates setting up strikes
        int iStep;
        for(iStep = 0; iStep < maturities.size(); iStep++) {
            if (strikes[iStep].size() != spread[iStep].size()) {
                throw ModelException(method, "Size mismatch between strikes ("+
                                     Format::toString(strikes[iStep].size()) +
                                     ") and probabilities ("+ 
                                     Format::toString(spread[iStep].size())+
                                     ")" + " at maturity " + 
                                     maturities[iStep].toString());
            }
            validateStrikes(strikes[iStep], maturities[iStep]);
            constructStrikes(strikes[iStep], adjEpsilon,
                             loStrikes[iStep], hiStrikes[iStep]);
        }
        // pass through to core routine
        spreads(loStrikes, hiStrikes, maturities, fwds, spread);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Caclulates spreads using specified lo and hi strikes */
void PDFDefaultLNStrike::spreads(
    const CLatticeDouble&    loStrikes,
    const CLatticeDouble&    hiStrikes,
    const DateTimeArray&     maturities,
    const DoubleArray&       fwds, // at maturities
    CLatticeDouble&          spread) const {
    static const string method = "PDFDefaultLNStrike::spreads";
    try{
        // set up our dates for the callSpread method
        DateTimeArray volDates(2);
        volDates[0] = fwdStarting? lnRequest->getStartDateTime(): valueDate;
        DoubleArray results(1); // reserve memory for results
        VolRequestLNStrikeSP request(copy(lnRequest->getVolRequest()));

        // allow negative fwd variance
        request->allowNegativeFwdVar(true);
        bool set = false;

        // keep copy to restore later on
        bool prevFwdStarting = fwdStarting;
        DateTime prevVolDateZero = volDates[0];
        double prevSpotAtStart = spotAtStart;

        for(int iStep = 0; iStep < maturities.size(); iStep++) {
            if(volFwdStarting && maturities[iStep].
               isGreater(request->getStartDate()) && !set) {
                fwdStarting = true;
                volDates[0] = request->getStartDate();
                spotAtStart = asset->fwdValue(volDates[0]);
                set = true;
            }
            volDates[1] = maturities[iStep];
            bool adjustRequest = !fwdStarting && volFwdStarting;

            for (int i = 0; i < loStrikes[iStep].size(); i++) {
                double loStrike = loStrikes[iStep][i];
                double hiStrike = hiStrikes[iStep][i];

                request->setStrike(adjustRequest? 
                                   loStrike/spotAtStart: loStrike);
                // NB Using optimal form of CalcVar (at least for VolSurface)
                CVolProcessedBSSP procVol(asset->
                                          getProcessedVol(request.get()));
                
                procVol->CalcVar(volDates, CVolProcessedBS::forward, results);

                double lovar = results[0];
                procVol.reset(); /* avoid holding two processed vols at the
                                    same time - performance improvement
                                    depending on processed vol implementation*/

                request->setStrike(adjustRequest?
                                   hiStrike/spotAtStart: hiStrike);
                procVol = CVolProcessedBSSP(asset->
                                            getProcessedVol(request.get()));
                procVol->CalcVar(volDates, CVolProcessedBS::forward, results);
                double hivar = results[0];
                
                if (fwdStarting) {
                    loStrike *= spotAtStart;
                    hiStrike *= spotAtStart;
                }

                double lo = Black::price(true, fwds[iStep], 
                                         loStrike, 1.0, lovar);
                double hi = Black::price(true, fwds[iStep], 
                                         hiStrike, 1.0, hivar);
                
                spread[iStep][i] = lo-hi;
            }
        }
        
        fwdStarting = prevFwdStarting;
        volDates[0] = prevVolDateZero;
        spotAtStart = prevSpotAtStart;

    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** awaiting spec */
void PDFDefaultLNStrike::probabilities(
    const CLatticeDouble&    strikes,
    const DateTimeArray&     maturities,
    CLatticeDouble&          probs) const {
    static const string method = "PDFDefaultLNStrike::probabilities";

    bool ok = true; // used to allow all errors to be reported
    string failedStrikes;

    try{
        double adjEpsilon = epsilon; // work on copy
        // calculate spreads
        spreads(strikes, maturities, adjEpsilon, probs);
        // work out probabilities
        const VolRequestLNStrike* volRequest = lnRequest->getVolRequest();
        const DateTime& startDate = volRequest->getStartDate();

        bool set = false;
        double prevSpotAtStart = spotAtStart;

        for(int iStep = 0; iStep < maturities.size(); iStep++) {
            bool relativeStrikes = fwdStarting || 
                                   !fwdStarting && volFwdStarting && maturities[iStep].isGreater(startDate);
            if(volFwdStarting && maturities[iStep].isGreater(startDate) && !set) {
                spotAtStart = asset->fwdValue(startDate);
                set = true;
            }

            
            for (int i = 0; i < strikes[iStep].size(); i++) {
                // calculate probability
                probs[iStep][i] /= 2*adjEpsilon*strikes[iStep][i] * (
                    relativeStrikes? spotAtStart: 1.0);
                // validate probabilities
                if(i > 0 && probs[iStep][i-1] < probs[iStep][i] -
                   accuracy * adjEpsilon * adjEpsilon) {
                    ok = false;
                    failedStrikes += "Increasing probabilities " + 
                        Format::toString(probs[iStep][i-1])+ " " +
                        Format::toString(probs[iStep][i]) + 
                        " at strikes " + 
                        Format::toString(strikes[iStep][i-1]) + " " +
                        Format::toString(strikes[iStep][i]) + 
                        " at maturity " + 
                        maturities[iStep].toString() + ". \n";
                }

                if(probs[iStep][i] < - accuracy * adjEpsilon * adjEpsilon) {
                    ok = false;
                    failedStrikes += "Probability " +  
                        Format::toString(probs[iStep][i]) + 
                        " < 0% at strike " + 
                        Format::toString(strikes[iStep][i]) +
                        " at maturity " + 
                        maturities[iStep].toString() + ". \n";
                }
                
                if(probs[iStep][i] > 1 + accuracy * adjEpsilon * adjEpsilon) {
                    ok = false;
                    failedStrikes += "Probability " +
                        Format::toString(probs[iStep][i]) + 
                        " > 100% at strike " + 
                        Format::toString(strikes[iStep][i]) +
                        " at maturity " + 
                        maturities[iStep].toString() + ". \n";
                }
            }
        }
        spotAtStart = prevSpotAtStart;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
    if(!ok) {
        throw PDFCalculatorException(method, 
                                     "Asset " + asset->getName() + 
                                     ": \n" + failedStrikes);
    }
}


/** awaiting spec */
void PDFDefaultLNStrike::localDensity(const DoubleArray& strikes,
                                      const DateTime&    maturity,
                                      DoubleArray&       density) const {
    static const string method = "PDFDefaultLNStrike::localDensity";

    bool ok = true; // used to allow all errors to be reported
    string failedStrikes;

    try {
        // need to do two sets of spreads, first at strike-epsilon/2 the
        // other at strike+epsilon/2
        // So start by constructing these strikes
        DoubleArray       strikesCopy(strikes); // to avoid const problems
        CSliceDouble      myStrikes(&strikesCopy[0], strikes.size());
        CSliceDouble      loStrikes(strikes.size());
        CSliceDouble      hiStrikes(strikes.size());
        validateStrikes(myStrikes, maturity);
        double adjEpsilon = epsilon; // work on copy
        constructStrikes(myStrikes, adjEpsilon, loStrikes, hiStrikes);
        DateTimeArray maturities(1, maturity);
        // calculate forwards at maturities
        DoubleArray fwds(maturities.size());
        asset->fwdValue(maturities, fwds);
        // calculate spreads at loStrikes
        CSliceDouble  loSpread(strikes.size());
        spreads(loStrikes, myStrikes, maturities, fwds, loSpread);
        // and then at hi strikes
        CSliceDouble  hiSpread(strikes.size());
        spreads(myStrikes, hiStrikes, maturities, fwds, hiSpread);

        // turn spreads into local densitiy
        for (int i = 0; i < strikes.size(); i++) {
            double k = fwdStarting? strikes[i]*spotAtStart: strikes[i];;
            double spread = loSpread[i]-hiSpread[i];
            density[i] = spread/(adjEpsilon*adjEpsilon*k*k);
            if (density[i] < - accuracy * adjEpsilon * adjEpsilon) {
                ok = false;
                failedStrikes += ", " + Format::toString(strikes[i+1]);
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
    if(!ok) {
        throw PDFCalculatorException(method, 
                                     "Negative density at strikes " + 
                                     failedStrikes + " for asset " + 
                                     asset->getName() +
                                     " on date " + 
                                     maturity.toString());
    }
}

/** awaiting spec */
void PDFDefaultLNStrike::integratedDensity(const DoubleArray& strikes,
                                           const DateTime&    maturity,
                                           DoubleArray&       density) const {
    static const string method = "PDFDefaultLNStrike::integratedDensity";
    if (strikes.size() < 3) {
        throw ModelException(method, "need at least two strikes");
    }

    try{
        // calculate mid strikes
        CSliceDouble midStrikes(strikes.size()-1);
        for (int i = 1; i < strikes.size(); i++){
            midStrikes[i-1] = (strikes[i-1] + strikes[i])/2.0;
        }
        // set up datetime array
        DateTimeArray maturities(1, maturity);
        // work on a copy of epsilon
        double adjEpsilon = epsilon;
        // calculate spreads for these strikes
        CSliceDouble spread(midStrikes.size());
        spreads(midStrikes, maturities, adjEpsilon, spread); 

        double sum = 0.0;
        // do density excluding first point
        string failedStrikes;
        bool ok = true;
        double scaleFactor = fwdStarting? spotAtStart: 1.0;
        for (int j = 1; j < strikes.size(); j++){
            if (j < strikes.size()-1){
                density[j] = 
                    (spread[j-1]/(2.0*adjEpsilon*scaleFactor*midStrikes[j-1]))-
                    (spread[j]/(2.0*adjEpsilon*scaleFactor*midStrikes[j]));
            } else {
                density[j] = 
                    spread[j-1]/(2.0*adjEpsilon*scaleFactor*midStrikes[j-1]);
            }
            if (density[j] < - accuracy * adjEpsilon * adjEpsilon) {
                ok = false;
                failedStrikes += ", " + Format::toString(j+1);
            }
            sum += density[j];
        }

        // and 1st density
        density[0] = 1.0 - sum;
        if(!ok) {
            throw PDFCalculatorException(
                method, 
                "Negative integrated density at strikes " + 
                failedStrikes + " for asset " + asset->getName() +
                " on date " + 
                maturity.toString());
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** calculates transient fields */
void PDFDefaultLNStrike::validatePop2Object(){
    epsilon  = lnRequest->getCallSpreadWidth();
    accuracy = lnRequest->getAccuracy();
    const VolRequestLNStrike* volRequest = lnRequest->getVolRequest();
    const DateTime& startDate = lnRequest->getStartDateTime();
    fwdStarting = startDate.isGreater(valueDate);
    volFwdStarting = volRequest->getStartDate().isGreater(valueDate);
    if (fwdStarting || volFwdStarting) {
        spotAtStart = asset->fwdValue(startDate);
    } else {
        spotAtStart = 0.0;
    }
}

/** constructor - to do alter type of request to PDFRequestLNStrike */
PDFDefaultLNStrike::PDFDefaultLNStrike(
    const DateTime&    valueDate,
    const CAsset*      asset,
    const PDFRequest*  request):
    PDFCalculator(TYPE), valueDate(valueDate), asset(asset){
    if (!PDFRequestLNStrike::TYPE->isInstance(request)) {
        throw ModelException("PDFDefaultLNStrike::PDFDefaultLNStrike",
                             "Default implementation only works for LN "
                             "strike requests - given " +
                             request->getClass()->getName());
    }
    
    const PDFRequestLNStrike& requestLN = 
        dynamic_cast<const PDFRequestLNStrike&>(*request);
    // attach to ref but take copy if object is on the stack
    lnRequest = PDFRequestLNStrikeConstSP(copyIfRef(&requestLN));
    validatePop2Object();
}

PDFDefaultLNStrike::PDFDefaultLNStrike(
    const CClassConstSP&             clazz,
    const DateTime&                  valueDate,
    const CAssetConstSP&             asset,
    const PDFRequestLNStrikeConstSP& lnRequest):
    PDFCalculator(clazz), valueDate(valueDate), asset(asset),
    lnRequest(lnRequest){
    validatePop2Object();
}

PDFDefaultLNStrike::PDFDefaultLNStrike():PDFCalculator(TYPE) {}

PDFDefaultLNStrike::PDFDefaultLNStrike(
    const CClassConstSP& clazz):PDFCalculator(clazz) {}

class PDFDefaultLNStrike::Helper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(PDFDefaultLNStrike, clazz);
        SUPERCLASS(PDFCalculator);
        EMPTY_SHELL_METHOD(createDefault);
        FIELD(valueDate, "valueDate");
        FIELD(asset, "asset");
        FIELD(lnRequest, "request");
        FIELD(epsilon, "epsilon");
        FIELD_MAKE_TRANSIENT(epsilon);
        FIELD(accuracy, "accuracy");
        FIELD_MAKE_TRANSIENT(accuracy);
//        FIELD(startDate, "startDate");
//        FIELD_MAKE_TRANSIENT(startDate);
        FIELD(fwdStarting, "fwdStarting");
        FIELD_MAKE_TRANSIENT(fwdStarting);
        FIELD(volFwdStarting, "volFwdStarting");
        FIELD_MAKE_TRANSIENT(volFwdStarting);
        FIELD(spotAtStart, "spotAtStart");
        FIELD_MAKE_TRANSIENT(spotAtStart);
    }
     
    static IObject* createDefault(){
        return new PDFDefaultLNStrike();
    }    
};

CClassConstSP const PDFDefaultLNStrike::TYPE = CClass::registerClassLoadMethod(
"PDFDefaultLNStrike", typeid(PDFDefaultLNStrike), 
PDFDefaultLNStrike::Helper::load);

DRLIB_END_NAMESPACE
