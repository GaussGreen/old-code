//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolSpline.cpp
//
//   Description : Vol Spline
//
//   Author      : Regis Guichard
//
//   Date        : 04 April 02
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolSpline.hpp"
#include "edginc/Nrfns.hpp"

DRLIB_BEGIN_NAMESPACE

const double VolSpline::exphalf = exp(0.5);
const double VolSpline::default_upperSigmaUnscaled = 1.0;

/** Our VolParam class. This class exists primarily to keep the
    ComputeImpVol and spotVolSurfaceFromStrikes methods out of the
    CVolBaseParamSurface interface */
class VSVolParam: public CVolParam{
public:
    static CClassConstSP const TYPE;
    // constructor
    VSVolParam(): CVolParam(TYPE){}

    virtual void ComputeImpVol(const CVolBase*          vol,
                               const CLatticeDouble&    strikes,
                               const DateTimeArray&     maturities,
                               CLatticeDouble&          impV) const{
        // turn the vol into what we must have
        const VolSpline* myVol = 
            static_cast<const VolSpline *>(vol);
        // then just pass through the parameterised vol
        myVol->ComputeImpVol(strikes, maturities, impV);
    }

    /*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
    virtual VolSurface* spotVolSurfaceFromStrikes(
        const CVolBase*       vol,
        const CDoubleArray&   strikes) const{
        // turn the vol into what we must have
        const VolSpline* myVol = 
            static_cast<const VolSpline *>(vol);
        // then just pass through the parameterised vol
        return myVol->spotVolSurfaceFromStrikes(strikes);
    }

private:
    static void load(CClassSP& clazz){
        REGISTER(VSVolParam, clazz);
        SUPERCLASS(CVolParam);
        EMPTY_SHELL_METHOD(defaultCtor);
    }
    static IObject* defaultCtor(){
        return new VSVolParam();
    }
};

CClassConstSP const VSVolParam::TYPE =
CClass::registerClassLoadMethod("VSVolParam", typeid(VSVolParam), VSVolParam::load);

/** method that builds a CVolParam. */
CVolParam* VolSpline::createVolParam() const{
    return new VSVolParam();
}

string VolSpline::sensName(DeltaSurface* shift) const{
    return getName();
}

bool VolSpline::sensShift(DeltaSurface* shift){
    double shiftSize = shift->getShiftSize();
    if (!Maths::isZero(shiftSize)){
        // only bother if non zero 
        // what to do - spot moves to spot + shift
        // simply move strike ref k -> k (S+dS)/S
        double spot    = shift->getSpot();
        double newSpot = spot * (1.0 + shiftSize);

        strikeRef *= newSpot/spot;

        // and rebuild cache by shifting backbone ourselves
        VolSurfaceSP bbvol(copy(getBackboneSurface()));
        bbvol->sensShift(shift);
        buildCache(bbvol.get());
    }
    // DON'T tweak the backbone !!!
    return false; 
}
    
void VolSpline::ComputeImpVol(const CLatticeDouble&      strikes,
                              const DateTimeArray&       maturities,
                              CLatticeDouble&            impV) const {
    static const string routine("VolSpline::ComputeImpVol");

    if ((maturities.size() != strikes.size()) ||
        (maturities.size() != impV.size())) {
        throw ModelException(routine, "Size mismatch between strikes ("+ 
                             Format::toString(strikes.size()) +
                             "), maturities ("+ 
                             Format::toString(maturities.size())+
                             ") and impV ("+ 
                             Format::toString(impV.size())+ ")");
    }  
    
    int iBackboneDateStart = 0;
    for (int iMat = 0; iMat < maturities.size(); iMat++) {
        int nbStrikes = strikes[iMat].size();
        if (nbStrikes != impV[iMat].size()){
            throw ModelException(routine, "Size mismatch between strikes"
                                 " & maturities for Mat " +
                                 maturities[iMat].toString() +
                                 " (n "+ Format::toString(iMat) + ")");
        }
        int iBackboneDate = iBackboneDateStart;
        while(iBackboneDate < nbBackboneDates && (*backboneDates)[iBackboneDate] < maturities[iMat]) {
            ++iBackboneDate;
        }
        iBackboneDateStart = Maths::max(0, iBackboneDate - 1);
        int iBackboneDateEnd = Maths::min(iBackboneDate, nbBackboneDates - 1);
        CDoubleMatrixSP matrixsp(spotVolMatrixFromStrikes(&strikes[iMat][0],
                                                          nbStrikes,
                                                          iBackboneDateStart,
                                                          iBackboneDateEnd));
        CDoubleMatrix& matrix = *matrixsp;
        if (iBackboneDateEnd == iBackboneDateStart){
            int iStrike = 0;
            for (; iStrike < nbStrikes; ++iStrike){
                impV[iMat][iStrike] = matrix[iStrike][0];
            }
        }
        else{
            double yearFrac = timeMetric->yearFrac(baseDate, maturities[iMat]);
            int iStrike = 0;
            for (; iStrike < nbStrikes; ++iStrike){
                double varStart = Maths::square(matrix[iStrike][0] * (*sqrTimeFracs)[iBackboneDateStart]);
                double varEnd = Maths::square(matrix[iStrike][1] * (*sqrTimeFracs)[iBackboneDateEnd]);
                double yearFracDiff = (*timeFracs)[iBackboneDateEnd] - (*timeFracs)[iBackboneDateStart];
                double slope = (varEnd - varStart) / yearFracDiff;
                double var = varStart + slope * (yearFrac - (*timeFracs)[iBackboneDateStart]);
                if (Maths::isNegative(var)){
                    throw ModelException(routine,
                                         "Negative variance at maturity date " + 
                                         maturities[iMat].toString() + 
                                         " and at strike " +
                                         Format::toString(strikes[iMat][iStrike]));
                }
                else if(Maths::isZero(var)){
                    impV[iMat][iStrike] = 0.0;
                } 
                else{
                    impV[iMat][iStrike] = sqrt(var / yearFrac);
                }
            }
        }
    }
}

double VolSpline::scalessStrike(double strike, 
                                double strikeRefVol,
                                double sqrTimeFrac) const{
    return ((strike - strikeRef) / (strikeRefVol * sqrTimeFrac * strikeRef));
}

double VolSpline::leftExtrapolate(double strike, int iMat) const{
    return ((*backboneVols)[iMat][0] + 
            (*lowerSlopes)[iMat] * 
            (strike - (*backboneStrikes)[0]));
}

double VolSpline::rightExtrapolate(double strike, 
                                   double upperStrike,
                                   double upperVol,
                                   double upperSlope,
                                   double upperCurvature,
                                   double strikeRefVol,
                                   double sqrTimeFrac) const{
    double z = scalessStrike(strike, 
                             strikeRefVol,
                             sqrTimeFrac);
    double rightz = scalessStrike(upperStrike,
                                  strikeRefVol,
                                  sqrTimeFrac);

//    double upperSigma = upperSigmaUnscaled * sqrTimeFrac;
    double upperSigma = upperSigmaUnscaled / sqrTimeFrac;
    double x = (z - rightz) / upperSigma;

    double phi_0 = upperVol;
    double x_dash_invert = upperSigma * strikeRef * strikeRefVol * sqrTimeFrac;
    double phi_dash_0 = x_dash_invert * upperSlope;
    double phi_dash_dash_0 = Maths::square(x_dash_invert) * upperCurvature;

    double a = phi_0 + phi_dash_0 + phi_dash_dash_0;
    double b = - phi_dash_dash_0;
    double c = - exphalf * phi_dash_0;

    return ((a + b * exp(-0.5 * x * x) + c * exp(-0.5 * Maths::square(x + 1.0))));
}

// do some caching to avoid memory reallocation
CDoubleMatrixSP VolSpline::createVolMatrix(int nbStrikes,
                                           int nbDates) const{
    if (!volMatrix || volMatrix->numCols() != nbStrikes
        || volMatrix->numRows() != nbDates){
        volMatrix = CDoubleMatrixSP(new CDoubleMatrix(nbStrikes,
                                                      nbDates));
    }
    return volMatrix;
}


/*** Creates a "Vol Matrix" on provided Strikes and on provided
     backbone date indices */
CDoubleMatrixSP VolSpline::spotVolMatrixFromStrikes(
    const double*         strikes,
    int                   strikes_size,
    int                   iBackboneDateStart,
    int                   iBackboneDateEnd) const{
    static const string routine("VolSpline::spotVolMatrixFromStrikes");

    try{
        if (iBackboneDateStart < 0
            || iBackboneDateStart > iBackboneDateEnd
            || iBackboneDateEnd >= nbBackboneDates){
            throw ModelException(routine,
                                 Format::toString("Must have 0 <= iBackboneDateStart <= iBackboneDateEnd < nbBackboneDates;\n"
                                                  "got %ld, %ld and %ld respectively",
                                                  iBackboneDateStart,
                                                  iBackboneDateEnd,
                                                  nbBackboneDates));
        }
        CDoubleMatrixSP matrixsp(createVolMatrix(strikes_size,
                                                 iBackboneDateEnd - iBackboneDateStart + 1));
        CDoubleMatrix& matrix = *matrixsp;
        int iMat = 0;
        int iBackboneDate = iBackboneDateStart;
        for (; iBackboneDate <= iBackboneDateEnd; iBackboneDate++, iMat++) {
            int iStrike = 0;
            for (; iStrike < strikes_size; ++iStrike) {
                if(nbBackboneStrikes > 1) {                
                    // left hand side
                    if (Maths::isNegative(strikes[iStrike] - (*backboneStrikes)[0])) {  // strikes[iStrike] < (*backboneStrikes)[0])
                        matrix[iStrike][iMat] = leftExtrapolate(strikes[iStrike], iBackboneDate);
                    }
                    // right hand side
                    else if (!Maths::isNegative(strikes[iStrike] - (*backboneStrikes)[nbBackboneStrikes - 1])) {    // strikes[iStrike] >= (*backboneStrikes)[nbBackboneStrikes - 1]
                        matrix[iStrike][iMat] = rightExtrapolate(strikes[iStrike], 
                                                                 (*upperStrikes)[iBackboneDate],
                                                                 (*upperVols)[iBackboneDate],
                                                                 (*upperSlopes)[iBackboneDate],
                                                                 (*upperCurvatures)[iBackboneDate],
                                                                 (*strikeRefVols)[iBackboneDate],
                                                                 (*sqrTimeFracs)[iBackboneDate]);
                    }
                    // middle
                    else {
                        matrix[iStrike][iMat] = splines[iBackboneDate]->valueWithGuess(strikes[iStrike],
                                                                                       splineStrikeIndices[iMat],
                                                                                       0);
                    }
                }
                else{
                    matrix[iStrike][iMat] = (*backboneVols)[iBackboneDate][0];
                }
            }
        }
        return matrixsp;
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

/*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
VolSurface* VolSpline::spotVolSurfaceFromStrikes(
    const CDoubleArray&   strikes) const{
    static const string routine("VolSpline::spotVolSurfaceFromStrikes");
    try{
        CDoubleMatrixSP matrix(spotVolMatrixFromStrikes(&strikes[0],
                                                        strikes.size(),
                                                        0,
                                                        nbBackboneDates-1));
        /** for performance use special constructor */
        return new VolSurface(getBackboneSurface(), strikes, *matrix);
    } catch (exception& e){
        string message = "Failed constructing surface from strike " + Format::toString(strikes[0]);
        throw ModelException(e, routine, message);
    }
}

VolSpline::VolSpline():
    CVolBaseParamSurface(TYPE),
    strikeRef(0.0),
    upperSigmaUnscaled(default_upperSigmaUnscaled),
    nbBackboneDates(0),
    nbBackboneStrikes(0) {}

VolSpline::VolSpline(const VolSurface& volSurface,
                     double strikeRef,
                     double upperSigmaUnscaled):
    CVolBaseParamSurface(TYPE, volSurface, CVolParamSP(new VSVolParam())),
    strikeRef(strikeRef),
    upperSigmaUnscaled(upperSigmaUnscaled) {
    buildCache();
}

void VolSpline::validatePop2Object(){
    static const string method = "VolSpline::validatePop2Object";

    if (Maths::isNegative(strikeRef)){
        throw ModelException(method, "strikeRef must be non-negative. (NB: Will be internally set if zero)");
    }
    if (!Maths::isPositive(upperSigmaUnscaled)){
        throw ModelException(method, "upperSigmaUnscaled must be non-negative");
    }
}

/*** Build the parameterised vol and cache any values **/
void VolSpline::buildCache() {
    /* Some minimum amount of stuff to cache */
    buildCache(getBackboneSurface());
}

/*** Build the parameterised vol and cache any values **/
void VolSpline::buildCache(const VolSurface* backbone) {
    backboneStrikes = DoubleArraySP(copy(&backbone->getStrikes()));
    nbBackboneStrikes = backboneStrikes->size();
    backboneDates = DateTimeArraySP(copy(&backbone->getDates()));
    nbBackboneDates = backboneDates->size();
    backboneVols = CDoubleMatrixSP(copy(&backbone->getVolMatrix()));
    backboneVols->transpose(); // transpose so as to be able get a ref to array of backboneVols at given maturity
    baseDate = backbone->getBaseDate();
    timeMetric = backbone->getTimeMetric();
    /* This is weak. When VolSurface is splined into VolSpline at getMarket time, I have 
       no (easy!) way to access the asset's spot price (and therefore can't set strikeRef
       to spot price). So, I have to do something rather ugly here... */
    if (!Maths::isPositive(strikeRef)){
        strikeRef = 0.0;
        int iStrike = 0;
        for (; iStrike < nbBackboneStrikes; ++iStrike){
            strikeRef += (*backboneStrikes)[iStrike];
        }
        strikeRef /= nbBackboneStrikes;     // average of backbone strikes
    }
    // backbone stuff
    sqrTimeFracs = DoubleArraySP(new DoubleArray(nbBackboneDates));
    timeFracs = DoubleArraySP(new DoubleArray(nbBackboneDates));
    int iMat = 0;
    for (; iMat < nbBackboneDates; iMat++) {
        // get vol at strike ref and time fracs
        (*timeFracs)[iMat] = timeMetric->yearFrac(baseDate, (*backboneDates)[iMat]);
        (*sqrTimeFracs)[iMat] = sqrt((*timeFracs)[iMat]);
    }
    /* If there's only one strike, there isn't much else to cache */
    if (nbBackboneStrikes == 1){
        return;
    }
    /* Stuff needed for left hand side / lower extrapolation */
    lowerSlopes  = DoubleArraySP(new DoubleArray(nbBackboneDates));
    /* Stuff needed for right hand side / upper extrapolation */
    upperStrikes = DoubleArraySP(new DoubleArray(nbBackboneDates));
    upperVols = DoubleArraySP(new DoubleArray(nbBackboneDates));
    upperSlopes = DoubleArraySP(new DoubleArray(nbBackboneDates));
    upperCurvatures = DoubleArraySP(new DoubleArray(nbBackboneDates));
    /* Back bone stuff */
    strikeRefVols = DoubleArraySP(new DoubleArray(nbBackboneDates));
    splines.resize(nbBackboneDates);
    splineStrikeIndices = IntArray(nbBackboneDates, -1);
    /* Loop over maturities */
    for (iMat = 0; iMat < nbBackboneDates; iMat++) {
        // get right hand side strike and vol
        (*upperStrikes)[iMat] = (*backboneStrikes)[nbBackboneStrikes-1];
        (*upperVols)[iMat] = (*backboneVols)[iMat][nbBackboneStrikes-1];
        // get slope from linear fit
        if(nbBackboneStrikes == 2) {
            (*lowerSlopes)[iMat]  = ((*backboneVols)[iMat][0]  - (*backboneVols)[iMat][1]) / 
                ((*backboneStrikes)[0]     - (*backboneStrikes)[1]);
            (*lowerSlopes)[iMat]  = Maths::min((*lowerSlopes)[iMat], 0.0);

            (*upperSlopes)[iMat] = ((*backboneVols)[iMat][nbBackboneStrikes-1]  - (*backboneVols)[iMat][nbBackboneStrikes-2]) / 
                ((*backboneStrikes)[nbBackboneStrikes-1]     - (*backboneStrikes)[nbBackboneStrikes-2]);
        } 
        // get slope from quadratic fit
        else { 
            DoubleArray leftCoeffs(3);
            DoubleArray rightCoeffs(3);
            polcof(&(*backboneStrikes)[0], 
                   &(*backboneVols)[iMat][0], 
                   2, 
                   &leftCoeffs[0]);
            polcof(&(*backboneStrikes)[nbBackboneStrikes-3], 
                   &(*backboneVols)[iMat][nbBackboneStrikes-3], 
                   2, 
                   &rightCoeffs[0]);
            (*lowerSlopes)[iMat]  = 2.0 * leftCoeffs[2]  * (*backboneStrikes)[0] + leftCoeffs[1];
            (*lowerSlopes)[iMat]  = Maths::min((*lowerSlopes)[iMat], 0.0);

            (*upperSlopes)[iMat] = 2.0 * rightCoeffs[2] * (*backboneStrikes)[nbBackboneStrikes-1] + rightCoeffs[1];
        }
        // construct spline
        NRSpline spliner(1, (*lowerSlopes)[iMat],   // constant first deriv on left
                         1, (*upperSlopes)[iMat]);  // constant first deriv on right
        splines[iMat] = NRSpline::InterpolantSP::constCast(
            NRSpline::computeInterp(spliner,
                                    &(*backboneStrikes)[0],
                                    (*backboneVols)[iMat],
                                    nbBackboneStrikes));
        // get vol at strike ref and time fracs
        (*strikeRefVols)[iMat] = splines[iMat]->value(strikeRef);
        // get right hand side second derivative
        (*upperCurvatures)[iMat] = splines[iMat]->value((*upperStrikes)[iMat] - DBL_EPSILON, 2);
    }
}


class VolSplineHelper{        
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VolSpline, clazz);
        SUPERCLASS(CVolBaseParamSurface);
        IMPLEMENTS(IVolatilityBS);
        IMPLEMENTS(IVolatilityDVF);
        IMPLEMENTS(DeltaSurface::IShift);
        EMPTY_SHELL_METHOD(defaultCtor);

        FIELD(strikeRef, "Backbone Strike");
        FIELD_MAKE_OPTIONAL(strikeRef);        
        FIELD(upperSigmaUnscaled, "Transition speed from spline to flat.");        
        FIELD_MAKE_OPTIONAL(upperSigmaUnscaled);        

        FIELD(timeMetric, "");
        FIELD_MAKE_TRANSIENT(timeMetric);
        FIELD(baseDate, "");
        FIELD_MAKE_TRANSIENT(baseDate);
        FIELD(backboneDates, "");
        FIELD_MAKE_TRANSIENT(backboneDates);
        FIELD(nbBackboneDates, "");
        FIELD_MAKE_TRANSIENT(nbBackboneDates);
        FIELD(backboneStrikes, "");
        FIELD_MAKE_TRANSIENT(backboneStrikes);
        FIELD(nbBackboneStrikes, "");
        FIELD_MAKE_TRANSIENT(nbBackboneStrikes);
        FIELD(backboneVols, "");
        FIELD_MAKE_TRANSIENT(backboneVols);
        FIELD(lowerSlopes, "");
        FIELD_MAKE_TRANSIENT(lowerSlopes);
        FIELD(upperStrikes, "");
        FIELD_MAKE_TRANSIENT(upperStrikes);
        FIELD(upperVols, "");
        FIELD_MAKE_TRANSIENT(upperVols);
        FIELD(upperSlopes, "");
        FIELD_MAKE_TRANSIENT(upperSlopes);
        FIELD(upperCurvatures, "");
        FIELD_MAKE_TRANSIENT(upperCurvatures);
        FIELD(strikeRefVols, "");
        FIELD_MAKE_TRANSIENT(strikeRefVols);
        FIELD(timeFracs, "");
        FIELD_MAKE_TRANSIENT(timeFracs);
        FIELD(sqrTimeFracs, "");
        FIELD_MAKE_TRANSIENT(sqrTimeFracs);
        FIELD(splines, "");
        FIELD_MAKE_TRANSIENT(splines);
        FIELD(splineStrikeIndices, "");
        FIELD_MAKE_TRANSIENT(splineStrikeIndices);
        FIELD(volMatrix, "");
        FIELD_MAKE_TRANSIENT(volMatrix);
    }

    static IObject* defaultCtor(){
        return new VolSpline();
    }
};

CClassConstSP const VolSpline::TYPE =
CClass::registerClassLoadMethod("VolSpline", typeid(VolSpline), VolSplineHelper::load);

// external symbol to allow class to be forced to be linked in
bool VolSplineLinkIn(){
    return (VolSpline::TYPE != 0);
}

DRLIB_END_NAMESPACE


