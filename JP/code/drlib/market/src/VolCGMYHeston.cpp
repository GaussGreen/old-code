//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolCGMYHeston.cpp
//
//   Description : 
//
//   Date        : 11 November 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolCGMYHeston.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

VolCGMYHeston::~VolCGMYHeston(){
    // trying to shut up compiler
}
void CGMYe::validatePop2Object(){
    static const string routine = "CGMYe::validatePop2Object";
    try {
        double gammaYNegative = Maths::SpecialFunc::gamma(-yNegative);
        cNegativeTimesGammaYNegative = cNegative * gammaYNegative;
        double gammaYPositive = Maths::SpecialFunc::gamma(-yPositive);
        cPositiveTimesGammaYPositive = cPositive * gammaYPositive;
        mToThePowerYPositive = pow(m, yPositive);
        gToThePowerYNegative = pow(g, yNegative);
        sqGaussianVol = gaussianVol * gaussianVol;
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}

Complex CGMYe::calcLaplaceExponent(double         tau,
                                   const Complex& u) const{   
    static const string routine = "CGMYe::calcLaplaceExponent";    
    // Carr, Geman, Madan, Yor: Stochastic Volatility for Levy processes page 10
    double u_re = u.real();
    if(!Maths::isPositive(m - u_re) || !Maths::isPositive(g + u_re) ) {
        throw ModelException(routine, "The real part of the frequency u (" 
                                      + Format::toString(u_re) 
                                      + ") must be in the strip (-g, m) = ("
                                      + Format::toString(-g) 
                                      + ", " + Format::toString(m) + ")");
    }

    Complex cumulant = cPositiveTimesGammaYPositive * (pow(m - u, yPositive) - mToThePowerYPositive) 
                       + cNegativeTimesGammaYNegative * (pow(g + u, yNegative) - gToThePowerYNegative) 
                       + (u * u) * (0.5 * sqGaussianVol);
    return (cumulant * tau);
}

CGMYe::CGMYe(double  cNegative,
             double  cPositive,
             double  g,
             double  m,
             double  yNegative,
             double  yPositive,
             double  gaussianVol):
CObject(TYPE),
cNegative(cNegative),
cPositive(cPositive),
g(g),
m(m),
yNegative(yNegative),
yPositive(yPositive),
gaussianVol(gaussianVol){
    static const string routine = "CGMYe::CGMYe";
    try {
        validatePop2Object();
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}

CGMYe::CGMYe():
CObject(TYPE),
cNegative(0.0),
cPositive(0.0),
g(0.0),
m(0.0),
yNegative(0.0),
yPositive(0.0),
gaussianVol(0.0),
cNegativeTimesGammaYNegative(0.0),
cPositiveTimesGammaYPositive(0.0),
mToThePowerYPositive(0.0),
gToThePowerYNegative(0.0),
sqGaussianVol(0.0){}

IObject* CGMYe::defaultCtor(){
    return new CGMYe();
}

void CGMYe::load(CClassSP& clazz){
    REGISTER(CGMYe, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultCtor);
    FIELD(cNegative, "cNegative");
    FIELD(cPositive, "cPositive");
    FIELD(g, "g");
    FIELD(m, "m");
    FIELD(yNegative, "yNegative");
    FIELD(yPositive, "yPositive");
    FIELD(gaussianVol, "gaussianVol");
    FIELD(cNegativeTimesGammaYNegative, "");
    FIELD_MAKE_TRANSIENT(cNegativeTimesGammaYNegative);
    FIELD(cPositiveTimesGammaYPositive, "");
    FIELD_MAKE_TRANSIENT(cPositiveTimesGammaYPositive);
    FIELD(mToThePowerYPositive, "");
    FIELD_MAKE_TRANSIENT(mToThePowerYPositive);
    FIELD(gToThePowerYNegative, "");
    FIELD_MAKE_TRANSIENT(gToThePowerYNegative);
    FIELD(sqGaussianVol, "");
    FIELD_MAKE_TRANSIENT(sqGaussianVol);
}

CClassConstSP const CGMYe::TYPE =
    CClass::registerClassLoadMethod("CGMYe", typeid(CGMYe), load);

void VolCGMYHeston::VolCGMYHestonParam::ComputeImpVol(const CVolBase*          vol,
                                                      const CLatticeDouble&    strikes,
                                                      const DateTimeArray&     maturities,
                                                      CLatticeDouble&          impV) const{
    // turn the vol into what we must have
    const VolCGMYHeston* myVol = static_cast<const VolCGMYHeston *>(vol);
    // then just pass through the parameterised vol
    myVol->ComputeImpVol(strikes, maturities, impV);
}

/*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
VolSurface* VolCGMYHeston::VolCGMYHestonParam::spotVolSurfaceFromStrikes(
    const CVolBase*       vol,
    const CDoubleArray&   strikes) const{
    // turn the vol into what we must have
    const VolCGMYHeston* myVol = 
        static_cast<const VolCGMYHeston *>(vol);
    // then just pass through the parameterised vol
    return myVol->spotVolSurfaceFromStrikes(strikes);
}

void VolCGMYHeston::validatePop2Object(){
    try {
        Calibrator::IAdjustable::checkRange(this);
        // build CGMYe
        update();
    }
    catch(exception& e){
        throw ModelException(e, "VolCGMYHeston::validatePop2Object");
    }
}

void VolCGMYHeston::ComputeImpVol(const CLatticeDouble&      strikes,
                           const DateTimeArray&       maturities,
                           CLatticeDouble&            impV) const {
    static const string routine("VolCGMYHeston::ComputeImpVol");
    throw ModelException(routine, "Not supported");
    if ((maturities.size() != strikes.size()) ||
        (maturities.size() != impV.size())) {
        throw ModelException(routine, "Size mismatch between strikes ("+ 
                             Format::toString(strikes.size()) +
                             "), maturities ("+ 
                             Format::toString(maturities.size())+
                             ") and impV ("+ 
                             Format::toString(impV.size())+ ")");
    }
    
    for (int iMat = 0; iMat < maturities.size(); iMat++) {
        if (strikes[iMat].size() != impV[iMat].size()){
            throw ModelException(routine, "Size mismatch between strikes"
                                 " & maturities for Mat " +
                                 maturities[iMat].toString() +
                                 " (n "+ Format::toString(iMat) + ")");
        }
        for (int iStrike = 0; iStrike < strikes[iMat].size(); iStrike ++) {
            impV[iMat][iStrike] = 0.0;
        }
    }
}


/** calculate a simple expected variance */
double VolCGMYHeston::FutureVariance(double mat) const{
    return 0.0;
}

    
/*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
VolSurface* VolCGMYHeston::spotVolSurfaceFromStrikes(
    const CDoubleArray&   strikes) const{
    static const string routine("VolCGMYHeston::spotVolSurfaceFromStrikes");
    try{
        throw ModelException(routine, "Not supported");

        const VolSurface* backbone = getBackboneSurface();
        const DateTimeArray& dates = backbone->getDates();
        CDoubleMatrix matrix(strikes.size(),
                             dates.size());
        for (int iStrike = 0; iStrike < strikes.size(); iStrike++) {
            for (int iMat = 0; iMat < dates.size(); iMat++) {
                matrix[iStrike][iMat] = 0.0;
            }
        }
        /** for performance need constructor that takes in
            cached values (to do) */
        VolSurface* volSurf = 
            new VolSurface(getName(),
                           timeMetric.get(),
                           strikes,
                           matrix,
                           backbone->getExpiries().get(),
                           baseDate);
        return volSurf;
    } catch (exception& e){
        throw ModelException(e, routine);
    }
} 

/** Needed for IAdjustable interface. Returns market data name for vol */
string VolCGMYHeston::getName() const{
    return CVolBaseParamSurface::getName();
}

/** Called after adjustments have been made to fields */
void VolCGMYHeston::fieldsUpdated(const CFieldArray& fields){
    update();
}

/** Called after adjustments have been made to fields */
void VolCGMYHeston::update(){
    static const string routine = "VolCGMYHeston::update";
    try {
        // refresh cgmye
        cgmye = CGMYeSP(new CGMYe(cNegative,
                                  cPositive,
                                  g,
                                  m,
                                  yNegative,
                                  yPositive,
                                  gaussianVol));
        heston = HestonSP(new Heston(initialVol,
                                     meanVol,
                                     meanReversRate,
                                     volVol,
                                     0.0)); // correlation
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}

/** method that builds a CVolParam. */
CVolParam* VolCGMYHeston::createVolParam() const{
    return new VolCGMYHestonParam();
}

IObject* VolCGMYHeston::defaultCtor(){
    return new VolCGMYHeston();
}

/** Invoked when Class is 'loaded' */
void VolCGMYHeston::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(VolCGMYHeston, clazz);
    SUPERCLASS(VolBaseParam);
//    IMPLEMENTS(IVolatilityBS);
//    IMPLEMENTS(IVolatilityDVF);
    IMPLEMENTS(Calibrator::IAdjustable);    
    IMPLEMENTS(IDynamicsParameter);    
    EMPTY_SHELL_METHOD(defaultCtor);

    // Clock
    FIELD(initialVol, "initialVol");
    FIELD(meanVol, "meanVol");
    FIELD(meanReversRate, "meanReversRate");
    FIELD(volVol, "volVol");
    // Levy
    FIELD(cNegative, "C negative");
    FIELD(cPositive, "C positive");
    FIELD(g, "G");
    FIELD(m, "m");
    FIELD(yNegative, "Y negative");
    FIELD(yPositive, "Y positive");
    FIELD(gaussianVol, "gaussianVol");

    // transient        
    FIELD(baseDate, "");
    FIELD_MAKE_TRANSIENT(baseDate);
    FIELD(cgmye, "");
    FIELD_MAKE_TRANSIENT(cgmye);
    FIELD(heston, "");
    FIELD_MAKE_TRANSIENT(heston);

    // add our fields and their ranges to central list
    // clock
    Calibrator::IAdjustable::registerField(
        clazz, "initialVol",
        new Range(ClosedBoundary(0.0), Infinity(Infinity::Plus)));
    Calibrator::IAdjustable::registerField(
        clazz, "meanVol", 
        new Range(ClosedBoundary(0.0), Infinity(Infinity::Plus)));
    Calibrator::IAdjustable::registerField(
        clazz, "meanReversRate",
        new Range(ClosedBoundary(0.0), Infinity(Infinity::Plus)));
    Calibrator::IAdjustable::registerField(
        clazz, "volVol", 
        new Range(ClosedBoundary(0.0), Infinity(Infinity::Plus)));
    Calibrator::IAdjustable::registerField(
        clazz, "meanReversRate",
        new Range(ClosedBoundary(0.0), Infinity(Infinity::Plus)));
    // Levy
    Calibrator::IAdjustable::registerField(
        clazz, "cNegative",
        new Range(ClosedBoundary(0.0), Infinity(Infinity::Plus)));
    Calibrator::IAdjustable::registerField(
        clazz, "cPositive", 
        new Range(ClosedBoundary(0.0), Infinity(Infinity::Plus)));
    Calibrator::IAdjustable::registerField(
        clazz, "g",
        new Range(ClosedBoundary(0.0), Infinity(Infinity::Plus)));
    Calibrator::IAdjustable::registerField(
        clazz, "m", 
        new Range(ClosedBoundary(0.0), Infinity(Infinity::Plus)));
    Calibrator::IAdjustable::registerField(
        clazz, "yPositive",
        new Range(Infinity(Infinity::Minus), OpenBoundary(2.0)));
    Calibrator::IAdjustable::registerField(
        clazz, "yNegative", 
        new Range(Infinity(Infinity::Minus), OpenBoundary(2.0)));
    Calibrator::IAdjustable::registerField(
        clazz, "gaussianVol", 
        new Range(ClosedBoundary(0.0), Infinity(Infinity::Plus)));
}

VolCGMYHeston::VolCGMYHeston(): VolBaseParam(TYPE),
                                initialVol(0.2),
                                volVol(0.4),
                                meanVol(0.2),
                                meanReversRate(1.0),
                                cNegative(0.35),
                                cPositive(0.35),
                                g(0.88),
                                m(20),
                                yNegative(0.4),
                                yPositive(0.4),
                                gaussianVol(0.1) {}

/*** Build the parameterised vol and cache any values **/
void VolCGMYHeston::buildCache() {
    const VolSurface* backbone = getBackboneSurface();
    baseDate = backbone->getBaseDate();
}

//////////// class FourierProcessCGMYHeston ///////////////////

/* Started Log Return */
Complex VolCGMYHeston::scalelessCumulant(const StFourierProcessLogRtn& process,
                                         const StFourierProductLogRtn& product, 
                                         const Complex& z, 
                                         const DateTime& matDate) const{
    static const string method = "VolCGMYHeston::scalelessCumulant";
    try{
        double tau = timeMetric->yearFrac(baseDate,
                                          matDate);

        Complex alpha, beta;
        try{    // compiler problem
            calcCumulantComponents(tau,
                                   z,
                                   *cgmye,
                                   *heston,
                                   alpha,
                                   beta);
        }
        catch(...){
            throw;
        }
        
        return (alpha + beta * Maths::square(initialVol));
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

/* Fwd starting Log Return  */
Complex VolCGMYHeston::scalelessCumulant(const FwdStFourierProcessLogRtn& process,
                                         const FwdStFourierProductLogRtn& product, 
                                         const Complex& z, 
                                         const DateTime& matDate) const{
    static const string method = "VolCGMYHeston::scalelessCumulant";
    try{
        /* Time fractions from today till start of option and then from
           start of option till maturity */
        double t = timeMetric->yearFrac(baseDate,
                                        product.getStartDate());
        double tau = timeMetric->yearFrac(product.getStartDate(),
                                          matDate);

        /* From start date to maturity */
        Complex alpha_tT, beta_tT;
        try{    // compiler problem
            calcCumulantComponents(tau,
                                   z,
                                   *cgmye,
                                   *heston,
                                   alpha_tT,
                                   beta_tT);
        }
        catch(...){
            throw;
        }

        /* Then, from today till start date */
        Complex alpha_t, beta_t;
        try{    // compiler problem
            heston->calcAlphaBeta(t,
                                  0.0,
                                  beta_tT,
                                  0.0,
                                  alpha_t,
                                  beta_t);
        }
        catch(...){
            throw;
        }
        
        return (alpha_t + alpha_tT + beta_t * Maths::square(initialVol));
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

void VolCGMYHeston::calcCumulantComponents(double          tau,    // tau = T - t (in years)
                                           const Complex&  u,
                                           const CGMYe&    cgmye,
                                           const Heston&   heston,
                                           Complex&        alpha,
                                           Complex&        beta) {
    static const string method = "FourierProcessCGMYHeston::calcCumulantComponents";    
    try {    
        double LevyTime = 1.0;
        Complex levyCumulant = cgmye.calcLaplaceExponent(LevyTime, u);
        Complex levyCumulantAtOne = cgmye.calcLaplaceExponent(LevyTime, 1.0);
        Complex z = levyCumulant - u * levyCumulantAtOne;
        
        // Evaluate the Laplace transform of the integrated variance at z
        heston.calcAlphaBeta(tau,
                             0.0,
                             0.0,                                  
                             z,
                             alpha,
                             beta);
    } catch(exception& e) {
        throw ModelException(e, method);
    }
}

CClassConstSP const VolCGMYHeston::TYPE =
CClass::registerClassLoadMethod("VolCGMYHeston", typeid(VolCGMYHeston), load);

CClassConstSP const VolCGMYHeston::VolCGMYHestonParam::TYPE =
CClass::registerClassLoadMethod("VolCGMYHeston::VolCGMYHestonParam", typeid(VolCGMYHestonParam), load);

/* external symbol to allow class to be forced to be linked in */
bool VolCGMYHestonLinkIn(){
    return (VolCGMYHeston::TYPE != 0);
}

DRLIB_END_NAMESPACE
