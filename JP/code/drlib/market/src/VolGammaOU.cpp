//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolGammaOU.cpp
//
//   Description : 
//
//   Date        : 01 July 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolGammaOU.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE
VolGammaOU::Helper::Helper(double meanVol,
                           double meanReversRate,
                           double volVol,                                  
                           double correlation):
VolOUHelper(TYPE, meanReversRate, correlation),
alpha(0.0),     // does not matter if volVol or meanVol == 0.0
nu(0.0){
    if(!Maths::isZero(meanVol) && !Maths::isZero(volVol)){
        alpha = Maths::square(meanVol / volVol);
        nu    = Maths::square(volVol * alpha);
    }
}

#if 0
void VolGammaOU::Helper::mapFromVarsToVols(const double& alpha,
                                           const double& nu,
                                           double&       meanVol,
                                           double&       volVol);
#endif

/** Computes lambda * \int_0^t k(f(s)) where
    k(x) = nu * x / (alpha - x), alpha > x and    
    f(s) = c1 + c2 exp(-lambda(t-s)) */
Complex VolGammaOU::Helper::integral(double         tau,
                                     const Complex& c1,
                                     const Complex& c2){    
    static const string method = "VolGammaOU::Helper::integral";    
    try {
        if(Maths::isZero(nu)) {
            return 0.0;
        }

        Complex alphaStar   = alpha - c1;
        double  effTime     = lambda * tau;
        double  expMEffTime = exp(-effTime);

        if(Maths::isPositive(c1.real() + c2.real() * expMEffTime - alpha)) { 
            throw ModelException(method,
                                 Format::toString("(c1, c2) = "
                                                  "(%f + i %f, %f + i %f) "
                                 "is outside domain of definition",
                                 c1.real(),
                                 c1.imag(),
                                 c2.real(),
                                 c2.imag()));
        }

        return ( nu * (c1 * effTime + alpha * log((alphaStar - c2 * expMEffTime) / (alphaStar - c2))) 
                 / alphaStar );
    } catch(exception& e) {
        throw ModelException(e, method);
    }
}

/** Returns log E_0 exp(u Z(1)) = nu * x / (alpha - x) for alpha > x*/
Complex VolGammaOU::Helper::cumulantBDLP(const Complex& u){
    static const string method = "VolGammaOU::Helper::cumulantBDLP";    
    try {
        if(Maths::isZero(nu)) {
            return 0.0;
        }
    
        if(Maths::isPositive(u.real() - alpha)) { 
            throw ModelException(method,
                                 Format::toString("u = (%f + i %f) "
                                 "is outside domain of definition",
                                 u.real(),
                                 u.imag()));
        }

        return (nu * u / (alpha - u) );
    } catch(exception& e) {
        throw ModelException(e, method);
    }
}

/** Invoked when Class is 'loaded' */
void VolGammaOU::Helper::load(CClassSP& clazz){
    REGISTER(Helper, clazz);
    SUPERCLASS(VolOUHelper);
    EMPTY_SHELL_METHOD(defaultCtor);
    FIELD(alpha, "alpha");
    FIELD(nu, "nu");
}

VolGammaOU::Helper::Helper():
VolOUHelper(TYPE),
alpha(0.0),
nu(0.0){}

IObject* VolGammaOU::Helper::defaultCtor(){
    return new Helper();
}

CClassConstSP const VolGammaOU::Helper::TYPE =
CClass::registerClassLoadMethod("VolGammaOU::Helper", typeid(VolGammaOU::Helper), load);

#if 0
void VolGammaOU::Helper::mapFromVarsToVols(const double& alpha,
                                                 const double& nu,
                                                 double&       meanVol,
                                                 double&       volVol) {

    static const string method = "VolGammaOU::Helper::mapFromVarsToVols";
    
    try {
        if(Maths::isZero(nu)) {
            meanVol= 0.0;
            volVol = 0.0;
            return;
        }

        meanVol = sqrt(nu / alpha);
        volVol  = sqrt(nu / Maths::square(alpha));

    } catch(exception& e) {
        throw ModelException(e, method);
    }

}
#endif

// VolGammaOU
void VolGammaOU::VolGammaOUParam::ComputeImpVol(
                                        const CVolBase*          vol,
                                        const CLatticeDouble&    strikes,
                                        const DateTimeArray&     maturities,
                                        CLatticeDouble&          impV) const{
    // turn the vol into what we must have
    const VolGammaOU* myVol = 
        static_cast<const VolGammaOU *>(vol);
    // then just pass through the parameterised vol
    myVol->ComputeImpVol(strikes, maturities, impV);
}

/*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
VolSurface* VolGammaOU::VolGammaOUParam::spotVolSurfaceFromStrikes(
    const CVolBase*       vol,
    const CDoubleArray&   strikes) const{
    // turn the vol into what we must have
    const VolGammaOU* myVol = 
        static_cast<const VolGammaOU *>(vol);
    // then just pass through the parameterised vol
    return myVol->spotVolSurfaceFromStrikes(strikes);
}

void VolGammaOU::validatePop2Object(){
    try {
        Calibrator::IAdjustable::checkRange(this);
        // build Helper
        update();
    }
    catch(exception& e){
        throw ModelException(e, "VolGammaOU::validatePop2Object");
    }
}

void VolGammaOU::ComputeImpVol(const CLatticeDouble&    strikes,
                               const DateTimeArray&     maturities,
                               CLatticeDouble&          impV) const {
    static const string routine("VolGammaOU::ComputeImpVol");
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

#if 0
/** calculate a simple expected variance */
double VolGammaOU::FutureVariance(double mat) const{
    double mVar = meanVol * meanVol;
    double effTime = meanReversRate * mat;
    double expMEffTime = exp(-effTime);
            
    return ( (1.0 - expMEffTime) * (initialVol * initialVol - mVar) /
             effTime + mVar);
}
#endif

/*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
VolSurface* VolGammaOU::spotVolSurfaceFromStrikes(
    const CDoubleArray&   strikes) const{
    static const string routine("VolGammaOU::spotVolSurfaceFromStrikes");
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
string VolGammaOU::getName() const{
    return CVolBaseParamSurface::getName();
}

/** Called after adjustments have been made to fields (eg calibrator) */
void VolGammaOU::fieldsUpdated(const CFieldArray& fields){
    update();
}

/** Called after adjustments have been made to fields */
void VolGammaOU::update(){
    helper = HelperSP(new Helper(meanVol,
                                 meanReversRate,
                                 volVol,                                  
                                 correlation));
}

/** method that builds a CVolParam. */
CVolParam* VolGammaOU::createVolParam() const{
    return new VolGammaOUParam();
}

IObject* VolGammaOU::defaultCtor(){
    return new VolGammaOU();
}

/** Invoked when Class is 'loaded' */
void VolGammaOU::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(VolGammaOU, clazz);
    SUPERCLASS(VolBaseParam);
//    IMPLEMENTS(IVolatilityBS);
//    IMPLEMENTS(IVolatilityDVF);
    IMPLEMENTS(IDynamicsParameter);
    IMPLEMENTS(Calibrator::IAdjustable);
    IMPLEMENTS(VolAJDSuper::ISuperposable);
    EMPTY_SHELL_METHOD(defaultCtor);
    FIELD(initialVol, "initialVol");
    FIELD(meanVol, "meanVol");
    FIELD(meanReversRate, "meanReversRate");
    FIELD(volVol, "volVol");
    FIELD(correlation, "correlation");

    // transient        
    FIELD(baseDate, "");
    FIELD_MAKE_TRANSIENT(baseDate);
    FIELD(helper, "");
    FIELD_MAKE_TRANSIENT(helper);

    // add our fields and their ranges to central list
    Calibrator::IAdjustable::registerField(
        clazz, "initialVol",
        new Range(ClosedBoundary(0.0), Infinity(Infinity::Plus)));
    Calibrator::IAdjustable::registerField(
        clazz, "correlation", 
        new Range(Infinity(Infinity::Minus), ClosedBoundary(0.0)));
    Calibrator::IAdjustable::registerField(
        clazz, "volVol", 
        new Range(ClosedBoundary(0.0), Infinity(Infinity::Plus)));
    Calibrator::IAdjustable::registerField(
        clazz, "meanVol", 
        new Range(ClosedBoundary(0.0), Infinity(Infinity::Plus)));
    Calibrator::IAdjustable::registerField(
        clazz, "meanReversRate",
        new Range(ClosedBoundary(0.0), Infinity(Infinity::Plus)));
}

VolGammaOU::VolGammaOU(): VolBaseParam(TYPE),
                          initialVol(0.2),
                          correlation(-0.7),
                          volVol(0.4),
                          meanVol(0.2),
                          meanReversRate(1.0){}

VolGammaOU::~VolGammaOU(){
    // trying to shup the compiler
}

/*** Build the parameterised vol and cache any values **/
void VolGammaOU::buildCache() {
    const VolSurface* backbone = getBackboneSurface();
    baseDate = backbone->getBaseDate();
}

// SUPERPOSABLE
/** Calculates the components alpha and betas that appear in 
    the exponent of the time-t joint Bilateral Laplace
    of X_T = (Y_T,..., X^i_T,...) where Y_T = ln(S_T / F(0, T)) 
    is the 'weighted' dimension-less log spot at time T and 
    the X^i_T's are any additional factors.
    The Bilateral Laplace transform is of the form
        exp(alpha(t, T, u) + sum_i betas^i(t, T, u) * X^i_t) */
void VolGammaOU::calcCumulantComponents(const DateTime&     fromDate,
                                        const DateTime&     toDate,
                                        const ComplexArray& u,
                                        double              weight,
                                        Complex&            alpha,
                                        ComplexArray&       betas) const{
    static const string method("VolSV::calcCumulantComponents");
    try{
        betas[0] = u[0];
        double tau = timeMetric->yearFrac(fromDate,
                                          toDate);
        helper->calcJointCumulantComponents(tau,
                                            u[0],     // log spot
                                            u[1],     // instantaneous variance
                                            0.0,      // integrated variance
                                            weight,
                                            alpha,
                                            betas[1]);
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

/** Returns the initial values Y_0,..., X^i_0,...
    Naturally, Y_0 = 0. */
double VolGammaOU::getInitialValue(int i) const{
    if (i == 1){
        return (initialVol * initialVol);
    }
    return 0.0;
}

/** Returns the nber of factors. The first factor is understood to be 
    the log spot. */
int VolGammaOU::getNbFactors() const{
    return 2;
}

void VolGammaOU::getMarket(const IModel*     model, 
                           const MarketData* market, 
                           const string&     name){
    VolBaseParam::getMarket(model, market, name);
}

//////////// class FourierProcessGammaOU ///////////////////
// Characteristic function
void VolGammaOU::calcCumulantComponents(
    double         tau,    // tau = T - t (in years)
    const Complex& z,
    Complex&       comp1,
    Complex&       comp2) const {
    static const string method = "VolGammaOU::calcCumulantComponents";
    try{        
        helper->calcJointCumulantComponents(tau,
                                            z,
                                            0.0,
                                            0.0,
                                            comp1,
                                            comp2);
    } catch(exception& e){
        throw ModelException(e, method);
    }
}

/* Started Log Return */
Complex VolGammaOU::scalelessCumulant(const StFourierProcessLogRtn& process,
                                      const StFourierProductLogRtn& product, 
                                      const Complex& z, 
                                      const DateTime& matDate) const{

    static const string method = "VolGammaOU::scalelessCumulant";
    try{
        double tau = timeMetric->yearFrac(baseDate,
                                          matDate);

        Complex comp1, comp2;
        calcCumulantComponents(tau,
                               z,
                               comp1,
                               comp2);
        
        return (comp1 + comp2 * Maths::square(initialVol));
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

/* Fwd starting Log Return  */
Complex VolGammaOU::scalelessCumulant(
    const FwdStFourierProcessLogRtn& process,
    const FwdStFourierProductLogRtn& product, 
    const Complex& z, 
    const DateTime& matDate) const{
    static const string method = "VolGammaOU::scalelessCumulant";
    try{
        /* Time fractions from today till start of option and then from
           start of option till maturity */
        double t = timeMetric->yearFrac(baseDate,
                                        product.getStartDate());
        double tau = timeMetric->yearFrac(product.getStartDate(),
                                          matDate);

        /* Calculate comp1, comp2 components from start date till maturity */
        Complex comp1_tT, comp2_tT;
        helper->calcJointCumulantComponents(tau,
                                            z,
                                            0.0,
                                            0.0,
                                            comp1_tT,
                                            comp2_tT);
        
        /* Then, from today till start date */
        Complex comp1_t, comp2_t;
        helper->calcJointCumulantComponents(t,
                                            0.0,
                                            comp2_tT,
                                            0.0,
                                            comp1_t,
                                            comp2_t);
        
        return (comp1_t + comp1_tT + comp2_t * Maths::square(initialVol));
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

/* Started integrated variance */
Complex VolGammaOU::cumulant(const StFourierProcessIntVar& process,
                             const StFourierProductIntVar& product, 
                             const Complex z, 
                             const DateTime& matDate) const {    
    static const string method = "VolGammaOU::cumulant";

    try{
        double tau = timeMetric->yearFrac(baseDate,
                                          matDate);

        Complex comp1, comp2;
        helper->calcJointCumulantComponents(tau,
                                            0.0,
                                            0.0,
                                            z,
                                            comp1,
                                            comp2);
        
        return (comp1 + comp2 * Maths::square(initialVol));
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}
    
/* Forward Starting integrated variance */
Complex VolGammaOU::cumulant(const FwdStFourierProcessIntVar& process,
                             const FwdStFourierProductIntVar& product, 
                             const Complex z, 
                             const DateTime& matDate) const {        
    static const string method = "VolGammaOU::cumulant";
    try{
        /* Time fractions from today till start of option and then from
           start of option till maturity */
        double t = timeMetric->yearFrac(baseDate,
                                        product.getStartDate());
        double tau = timeMetric->yearFrac(product.getStartDate(),
                                          matDate);

        /* Calculate alpha, beta component from start date till maturity */
        Complex comp1_tT, comp2_tT;          
        helper->calcJointCumulantComponents(tau,
                                            0.0,
                                            0.0,
                                            z,
                                            comp1_tT,
                                            comp2_tT);

        /* Then, from today till start date */
        Complex comp1_t, comp2_t;
        helper->calcJointCumulantComponents(t,
                                            0.0,
                                            0.0,
                                            comp2_tT,
                                            comp1_t,
                                            comp2_t);

        return (comp1_t + comp1_tT + comp2_t * Maths::square(initialVol));
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}
   
CClassConstSP const VolGammaOU::TYPE =
CClass::registerClassLoadMethod("VolGammaOU", typeid(VolGammaOU), load);

CClassConstSP const VolGammaOU::VolGammaOUParam::TYPE =
CClass::registerClassLoadMethod("VolGammaOU::VolGammaOUParam", typeid(VolGammaOUParam), load);

/* external symbol to allow class to be forced to be linked in */
bool VolGammaOULinkIn(){
    return (VolGammaOU::TYPE != 0);
}

DRLIB_END_NAMESPACE
