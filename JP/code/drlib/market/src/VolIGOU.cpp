//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolIGOU.cpp
//
//   Description : 
//
//   Date        : 07 November 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolIGOU.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

VolIGOU::~VolIGOU(){
    // trying to shut up compiler
}

VolIGOU::Helper::Helper(double meanVol,
                        double meanReversRate,
                        double volVol,                                  
                        double correlation):
VolOUHelper(TYPE, meanReversRate, correlation),
gamma(0.0),     // does not matter if volVol or meanVol == 0.0
delta(0.0){
    if(!Maths::isZero(meanVol) && !Maths::isZero(volVol)) {
        gamma = meanVol / volVol;
        delta = Maths::square(meanVol) * gamma;
    }
}

#if 0
void VolIGOU::Helper::mapFromVarsToVols(const double& gamma,
                                        const double& delta,
                                        double&       meanVol,
                                        double&       volVol);
#endif

/** Computes lambda * \int_0^t k(f(s)) where
    k(x) = delta * x / sqrt(gamma^2 - 2 * x), gamma^2 - 2 * x > 0 and    
    f(s) = c1 + c2 exp(-lambda(t-s)) */
Complex VolIGOU::Helper::integral(double         tau,
                                  const Complex& c1,
                                  const Complex& c2) {    
    static const string method = "VolIGOU::Helper::integral";
    try {
        if(Maths::isZero(delta)) {
            return 0.0;
        }

        double effTime     = lambda * tau;
        double expMEffTime = exp(-effTime);
        double gammaSq     = Maths::square(gamma);
    
        Complex gammaSqMinusTwoCOne = gammaSq - 2.0 * c1;
        Complex sqrtGgammaSqMinusTwoCOne = sqrt(gammaSqMinusTwoCOne);

        // part1
        Complex part1 = sqrt(gammaSqMinusTwoCOne - 2.0 * c2 * expMEffTime) - 
                        sqrt(gammaSqMinusTwoCOne - 2.0 * c2);

        // part2
        Complex ratio = (sqrtGgammaSqMinusTwoCOne + sqrt(gammaSqMinusTwoCOne - 2.0 * c2 * expMEffTime)) / 
                        (sqrtGgammaSqMinusTwoCOne + sqrt(gammaSqMinusTwoCOne - 2.0 * c2));

    
        Complex part2 = 2.0 * c1 * (effTime / 2.0 + log(ratio)) / sqrtGgammaSqMinusTwoCOne;

        return (delta * (part1 + part2));
    } catch(exception& e) {
        throw ModelException(e, method);
    }
}

/** Returns log E_0 exp(u Z(1)) = nu * x / (alpha - x) for alpha > x*/
Complex VolIGOU::Helper::cumulantBDLP(const Complex& u) {
    static const string method = "VolIGOU::Helper::cumulantBDLP";    
    try {
        if(Maths::isZero(delta)) {
            return 0.0;
        }
    
        if(!Maths::isPositive(Maths::square(gamma) - 2.0 * u.real())) { 
            throw ModelException(method,
                                 Format::toString("u = (%f + i %f) "
                                 "is outside domain of definition",
                                 u.real(),
                                 u.imag()));
        }
        return (delta * u / sqrt(Maths::square(gamma) - 2.0 * u) );
    } catch(exception& e) {
        throw ModelException(e, method);
    }
}

Complex VolIGOU::Helper::atan(Complex x){
    Complex i(0.0, 1.0);
    return (- i * atanh(i * x));
}

Complex VolIGOU::Helper::atanh(Complex x){
    return (0.5 * log((1.0 + x) / (1.0 - x)) );
}

/** Invoked when Class is 'loaded' */
void VolIGOU::Helper::load(CClassSP& clazz){
    REGISTER(Helper, clazz);
    SUPERCLASS(VolOUHelper);
    EMPTY_SHELL_METHOD(defaultCtor);
    FIELD(gamma, "gamma");
    FIELD(delta, "delta");
}

VolIGOU::Helper::Helper():
VolOUHelper(TYPE),
gamma(0.0),
delta(0.0){}

IObject* VolIGOU::Helper::defaultCtor(){
    return new Helper();
}

CClassConstSP const VolIGOU::Helper::TYPE =
CClass::registerClassLoadMethod("VolIGOU::Helper", typeid(VolIGOU::Helper), load);

#if 0
void VolIGOU::Helper::mapFromVarsToVols(const double& gamma,
                                  const double& delta,
                                  double&       meanVol,
                                  double&       volVol) {

    static const string method = "VolIGOU::Helper::mapFromVarsToVols";
    
    try {
        if(Maths::isZero(delta)) {
            meanVol= 0.0;
            volVol = 0.0;
            return;
        }

        meanVol = sqrt(delta / gamma);
        volVol  = sqrt(delta / (gamma * Maths::square(gamma)));

    } catch(exception& e) {
        throw ModelException(e, method);
    }

}
#endif

// VolIGOU class
void VolIGOU::VolIGOUParam::ComputeImpVol(const CVolBase*          vol,
                                          const CLatticeDouble&    strikes,
                                          const DateTimeArray&     maturities,
                                          CLatticeDouble&          impV) const{
    // turn the vol into what we must have
    const VolIGOU* myVol = static_cast<const VolIGOU *>(vol);
    // then just pass through the parameterised vol
    myVol->ComputeImpVol(strikes, maturities, impV);
}

/*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
VolSurface* VolIGOU::VolIGOUParam::spotVolSurfaceFromStrikes(
    const CVolBase*       vol,
    const CDoubleArray&   strikes) const{
    // turn the vol into what we must have
    const VolIGOU* myVol = 
        static_cast<const VolIGOU *>(vol);
    // then just pass through the parameterised vol
    return myVol->spotVolSurfaceFromStrikes(strikes);
}

void VolIGOU::validatePop2Object(){
    try {
        Calibrator::IAdjustable::checkRange(this);
        // build Helper
        update();
    }
    catch(exception& e){
        throw ModelException(e, "VolIGOU::validatePop2Object");
    }
}

void VolIGOU::ComputeImpVol(const CLatticeDouble& strikes,
                            const DateTimeArray&  maturities,
                            CLatticeDouble&       impV) const {
    static const string routine("VolIGOU::ComputeImpVol");
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
double VolIGOU::FutureVariance(double mat) const{
    /*the expected change in log(S) under a jump */
    return 0.0;
}

/** Needed for IAdjustable interface. Returns market data name for vol */
string VolIGOU::getName() const{
    return CVolBaseParamSurface::getName();
}

/** Called after adjustments have been made to fields (eg calibrator) */
void VolIGOU::fieldsUpdated(const CFieldArray& fields){
    update();
}

/** Called after adjustments have been made to fields */
void VolIGOU::update(){
    helper = HelperSP(new Helper(meanVol,
                                 meanReversRate,
                                 volVol,                                  
                                 correlation));
}

/** method that builds a CVolParam. */
CVolParam* VolIGOU::createVolParam() const{
    return new VolIGOUParam();
}

IObject* VolIGOU::defaultCtor(){
    return new VolIGOU();
}
    
/*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
VolSurface* VolIGOU::spotVolSurfaceFromStrikes(
    const CDoubleArray&   strikes) const{
    static const string routine("VolIGOU::spotVolSurfaceFromStrikes");
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

/** Invoked when Class is 'loaded' */
void VolIGOU::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(VolIGOU, clazz);
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

VolIGOU::VolIGOU(): VolBaseParam(TYPE),
                    initialVol(0.2),
                    correlation(-0.7),
                    volVol(0.4),
                    meanVol(0.2),
                    meanReversRate(1.0){}

/*** Build the parameterised vol and cache any values **/
/*** Build the parameterised vol and cache any values **/
void VolIGOU::buildCache() {
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
void VolIGOU::calcCumulantComponents(const DateTime&     fromDate,
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
double VolIGOU::getInitialValue(int i) const{
    if (i == 1){
        return (initialVol * initialVol);
    }
    return 0.0;
}

/** Returns the nber of factors. The first factor is understood to be 
    the log spot. */
int VolIGOU::getNbFactors() const{
    return 2;
}

void VolIGOU::getMarket(const IModel*     model, 
                           const MarketData* market, 
                           const string&     name){
    VolBaseParam::getMarket(model, market, name);
}

/* Started Log Return */
Complex VolIGOU::scalelessCumulant(const StFourierProcessLogRtn& process,
                                   const StFourierProductLogRtn& product, 
                                   const Complex& z, 
                                   const DateTime& matDate) const{
    static const string method = "VolIGOU::scalelessCumulant";
    try{
        double tau = timeMetric->yearFrac(baseDate,
                                          matDate);

        Complex comp1, comp2;
        helper->calcJointCumulantComponents(tau,
                                            z,
                                            0.0,
                                            0.0,
                                            comp1,
                                            comp2);

        return (comp1 + comp2 * Maths::square(initialVol));
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

/* Fwd starting Log Return  */
Complex VolIGOU::scalelessCumulant(const FwdStFourierProcessLogRtn& process,
                                   const FwdStFourierProductLogRtn& product, 
                                   const Complex& z, 
                                   const DateTime& matDate) const{
    static const string method = "VolIGOU::scalelessCumulant";
    try{
        /* Time fractions from today till start of option and then from
           start of option till maturity */
        double t = timeMetric->yearFrac(baseDate,
                                        product.getStartDate());
        double tau = timeMetric->yearFrac(product.getStartDate(),
                                          matDate);

        /* Calculate comp1, comp2 from start date till maturity */
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
Complex VolIGOU::cumulant(const StFourierProcessIntVar& process,
                          const StFourierProductIntVar& product, 
                          const Complex z, 
                          const DateTime& matDate) const {    
    static const string method = "VolIGOU::cumulant";
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
Complex VolIGOU::cumulant(const FwdStFourierProcessIntVar& process,
                          const FwdStFourierProductIntVar& product, 
                          const Complex z, 
                          const DateTime& matDate) const {        
    static const string method = "VolIGOU::cumulant";
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
   
CClassConstSP const VolIGOU::TYPE =
CClass::registerClassLoadMethod("VolIGOU", typeid(VolIGOU), load);

CClassConstSP const VolIGOU::VolIGOUParam::TYPE =
CClass::registerClassLoadMethod("VolIGOU::VolIGOUParam", typeid(VolIGOUParam), load);

/* external symbol to allow class to be forced to be linked in */
bool VolIGOULinkIn(){
    return (VolIGOU::TYPE != 0);
}

DRLIB_END_NAMESPACE
