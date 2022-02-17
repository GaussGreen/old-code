//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolSV.cpp
//
//   Description : 
//
//   Date        : 25 April 03
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_VOLSV_CPP
#include "edginc/VolSV.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/VolProcessed.hpp"
#include "edginc/VolRequestRaw.hpp"
#include "edginc/VolRequestTime.hpp"
#include "edginc/FlatVol.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/VolSVJ.hpp"
#include "edginc/VolSVJJ.hpp"
#include "edginc/VolSVCJ.hpp"


const double ACC = 0.0001;

DRLIB_BEGIN_NAMESPACE

DEFINE_TEMPLATE_TYPE(VolSVArray);


VolSV::VolSV(const string& name,
             double initialVol,
             double meanVol,
             double meanReversRate,
             double volVol,
             double correlation):
VolBaseParam(TYPE, name), initialVol(initialVol), meanVol(meanVol), meanReversRate(meanReversRate),
volVol(volVol), correlation(correlation) {
    validatePop2Object();
}


VolSVJSP VolSV::convert(VolSVJ* p) const {
    VolSVJSP volSVJ(new VolSVJ(
        getName(), initialVol, meanVol, meanReversRate, volVol, correlation, 0.0, 0.0, 0.0, 0.0));
    return volSVJ;
}


VolSVJJSP VolSV::convert(VolSVJJ* p) const {
    VolSVJJSP volSVJJ(new VolSVJJ(
        getName(), initialVol, meanVol, meanReversRate, volVol, correlation, 
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
    return volSVJJ;
}


VolSVCJSP VolSV::convert(VolSVCJ* p) const {
    VolSVCJSP volSVCJ(new VolSVCJ(
        getName(), initialVol, correlation, volVol, meanVol, meanReversRate,
        0.0, 0.0, 0.0, 0.0, 0.0, false));
    return volSVCJ;
}


template<> string nameForType<VolSV_SpotDSType>(VolSV_SpotDSType*){
    return "VolSV::MCParams::SpotDiscreteSchemeType";
}
template<> string VolSV_SpotDSTypeHelper::names[VolSV_SpotDSTypeHelper::EnumList::NB_ENUMS] = {
    "EXACT",
    "EULER"
};

template<> string nameForType<VolSV_VarDSType>(VolSV_VarDSType*){
    return "VolSV::MCParams::VarDiscreteSchemeType";
}
template<> string VolSV_VarDSTypeHelper::names[VolSV_VarDSTypeHelper::EnumList::NB_ENUMS] = {
    "EULER",
    "VAR_TRANSFORM_EULER"
};

void VolSV::SVVolParam::ComputeImpVol(const CVolBase*          vol,
                                      const CLatticeDouble&    strikes,
                                      const DateTimeArray&     maturities,
                                      CLatticeDouble&          impV) const{
    // turn the vol into what we must have
    const VolSV* myVol = static_cast<const VolSV *>(vol);
    // then just pass through the parameterised vol
    myVol->ComputeImpVol(strikes, maturities, impV);
}

/*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
VolSurface* VolSV::SVVolParam::spotVolSurfaceFromStrikes(
    const CVolBase*       vol,
    const CDoubleArray&   strikes) const{
    // turn the vol into what we must have
    const VolSV* myVol = 
        static_cast<const VolSV *>(vol);
    // then just pass through the parameterised vol
    return myVol->spotVolSurfaceFromStrikes(strikes);
}

void VolSV::validatePop2Object(){
    try {
        Calibrator::IAdjustable::checkRange(this);
        // build heston
        update();
    }
    catch(exception& e){
        throw ModelException(e, "VolSV::validatePop2Object");
    }
}

void VolSV::ComputeImpVol(const CLatticeDouble&      strikes,
                          const DateTimeArray&       maturities,
                          CLatticeDouble&            impV) const {
    static const string routine("VolSV::ComputeImpVol");
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

/*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
VolSurface* VolSV::spotVolSurfaceFromStrikes(
    const CDoubleArray&   strikes) const{
    static const string routine("VolSV::spotVolSurfaceFromStrikes");
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

/** method that builds a CVolParam. */
CVolParam* VolSV::createVolParam() const{
    return new SVVolParam();
}

IObject* VolSV::defaultCtor(){
    return new VolSV();
}

/** Invoked when Class is 'loaded' */
void VolSV::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(VolSV, clazz);
    SUPERCLASS(VolBaseParam);
//    IMPLEMENTS(IVolatilityBS);
//    IMPLEMENTS(IVolatilityDVF);
    IMPLEMENTS(Calibrator::IAdjustable);
    IMPLEMENTS(VolAJDSuper::ISuperposable);
    IMPLEMENTS(IDynamicsParameter);
    IMPLEMENTS(MarketDataConvert::IConvert<VolSVJ>);
    IMPLEMENTS(MarketDataConvert::IConvert<VolSVJJ>);
    IMPLEMENTS(MarketDataConvert::IConvert<VolSVCJ>);
    EMPTY_SHELL_METHOD(defaultCtor);
    FIELD(initialVol, "initialVol");
    FIELD(meanVol, "meanVol");
    FIELD(meanReversRate, "meanReversRate");
    FIELD(volVol, "volVol");
    FIELD(correlation, "correlation");

    // transient        
    FIELD(baseDate, "");
    FIELD_MAKE_TRANSIENT(baseDate);
    FIELD(heston, "");
    FIELD_MAKE_TRANSIENT(heston);
    // add our fields and their ranges to central list
    Calibrator::IAdjustable::registerField(
        clazz, "initialVol",
        new Range(Heston::RangeDef::initialVol));
    Calibrator::IAdjustable::registerField(
        clazz, "meanVol", 
        new Range(Heston::RangeDef::meanVol));
    Calibrator::IAdjustable::registerField(
        clazz, "meanReversRate",
        new Range(Heston::RangeDef::meanReversRate));
    Calibrator::IAdjustable::registerField(
        clazz, "volVol", 
        new Range(Heston::RangeDef::volVol));
    Calibrator::IAdjustable::registerField(
        clazz, "correlation", 
        new Range(Heston::RangeDef::correlation));
}

VolSV::VolSV(): 
VolBaseParam(TYPE),
initialVol(Heston::DefaultVal::initialVol),
meanVol(Heston::DefaultVal::meanVol),
meanReversRate(Heston::DefaultVal::meanReversRate),
volVol(Heston::DefaultVal::volVol),
correlation(Heston::DefaultVal::correlation){}


CVolProcessed* VolSV::getProcessedVol(const CVolRequest* volRequest,
                                       const CAsset*      /*asset*/) const{
    if (VolRequestRaw::TYPE->isInstance(volRequest) || 
        VolRequestTime::TYPE->isInstance(volRequest)){
        // it's ours or can just use this
        return const_cast<VolSV*>(this);
    }
    else if (ATMVolRequest::TYPE->isInstance(volRequest) ||
             LinearStrikeVolRequest::TYPE->isInstance(volRequest)) {
        // FWD_AT_MAT request for protected assets will ask for this
        // Looking just for something that doesn't break
        HolidaySP noHols(Holiday::noHolidays());
        TimeMetricSP tm(new TimeMetric(1.0, noHols.get()));
        //double flatFXVol = !compVol.empty()? compVol[0] : spotVol[0];
        double flatFXVol = initialVol;
        FlatVolSP flatVol(new FlatVol(this->getName(), 
                                      baseDate, // baseDate
                                      tm.get(),
                                      flatFXVol));
        return flatVol->getProcessedVol(volRequest, 0);
    }
    throw ModelException("VolSV:getProcessedVol", 
                         "Request of type "+
                         volRequest->getClass()->getName()+
                         " not supported for " +this->getName());
}


/*** Build the parameterised vol and cache any values **/
/*** Build the parameterised vol and cache any values **/
void VolSV::buildCache() {
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
void VolSV::calcCumulantComponents(const DateTime&     fromDate,
                                   const DateTime&     toDate,
                                   const ComplexArray& inu,
                                   double              weight,
                                   Complex&            alpha,
                                   ComplexArray&       betas) const{
    static const string method("VolSV::calcCumulantComponents");
    try{
        betas[0] = inu[0];
        double tau = timeMetric->yearFrac(fromDate,
                                          toDate);
        Complex u1 = inu[0] * weight;   // log spot
        Complex u2 = inu[1];            // instantaneous variance
        Complex u3 = (0.5 * weight * (1.0 - weight)) * inu[0];      // integrated variance
        calcJointLapAlphaBeta(tau,
                              u1,
                              u2,
                              u3,
                              alpha,
                              betas[1]);
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

/** Returns the initial values Y_0,..., X^i_0,...
    Naturally, Y_0 = 0. */
double VolSV::getInitialValue(int i) const{
    if (i == 1){
        return (initialVol * initialVol);
    }
    return 0.0;
}

/** Returns the nber of factors. The first factor is understood to be 
    the log spot. */
int VolSV::getNbFactors() const{
    return 2;
}

void VolSV::getMarket(const IModel*     model, 
                       const MarketData* market, 
                       const string&     name){
    VolBaseParam::getMarket(model, market, name);
}

/* Started Log Return */
Complex VolSV::scalelessCumulant(const StFourierProcessLogRtn& process,
                                  const StFourierProductLogRtn& product, 
                                  const Complex& z, 
                                  const DateTime& matDate) const{
    static const string method = "VolSV::scalelessCumulant";
    try{
        
        double tau = timeMetric->yearFrac(baseDate,
                                          matDate);

        Complex alpha, beta;
        calcJointLapAlphaBeta(tau,
                              z,
                              0.0,
                              0.0,
                              alpha,
                              beta);
        return (alpha + beta * Maths::square(initialVol));
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

/* Fwd starting Log Return  */
Complex VolSV::scalelessCumulant(const FwdStFourierProcessLogRtn& process,
                                  const FwdStFourierProductLogRtn& product, 
                                  const Complex& z, 
                                  const DateTime& matDate) const{
    static const string method = "VolSV::scalelessCumulant";
    try{
        
        /* Time fractions from today till start of option and then from
           start of option till maturity */
        double t = timeMetric->yearFrac(baseDate,
                                        product.getStartDate());
        double tau = timeMetric->yearFrac(product.getStartDate(),
                                          matDate);

        /* Calculate alpha, beta component from start date till maturity */
        Complex alpha_tT, beta_tT;
        calcJointLapAlphaBeta(tau,
                              z,
                              0.0,
                              0.0,
                              alpha_tT,
                              beta_tT);
        /* Then, from today till start date */
        Complex alpha_t, beta_t;
        calcJointLapAlphaBeta(t,
                              0.0,
                              beta_tT,
                              0.0,
                              alpha_t,
                              beta_t);
        return (alpha_t + beta_t * Maths::square(initialVol) + alpha_tT);
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

/**ARNAUD**////////////////////////////////////////////////////
/* Fwd Starting Expected Quadratic Variation */

Complex VolSV::cumulant(const FwdStFourierProcessExpQuadVar& process,
                        const FwdStFourierProductExpQuadVar& product, 
                        const Complex& z, 
                        const DateTime& matDate) const{
    static const string method = "VolSV::cumulant";
    try{
        /* Time fractions from today till start of option and then from
           start of option till maturity */
        double t = timeMetric->yearFrac(baseDate,
                                        product.getStartDate());
        double tau = timeMetric->yearFrac(product.getStartDate(),
                                          matDate);

        /* Calculate alpha, beta component from start date till maturity */
        Complex alpha_tT, beta_tT;
		double meanVar = Maths::square(meanVol);
		double phi = exp(-meanReversRate * tau);
		alpha_tT = z * meanVar / meanReversRate * (meanReversRate * tau - 1.0 + phi); 
		beta_tT = z / meanReversRate * (1.0 - phi);
        /* Then, from today till start date */
        Complex alpha_t, beta_t;
        calcJointLapAlphaBeta(t,
                              0.0,
                              beta_tT,
                              0.0,
                              alpha_t,
                              beta_t);
        return (alpha_t + beta_t * Maths::square(initialVol) + alpha_tT);
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}


///////////////////////////////////////////////////////////////////////

/* Started integrated variance */
Complex VolSV::cumulant(const StFourierProcessIntVar& process,
                         const StFourierProductIntVar& product, 
                         const Complex z, 
                         const DateTime& matDate) const {    
    static const string method = "VolSV::cumulant";
    try{

        double tau = timeMetric->yearFrac(baseDate,
                                       matDate);


        Complex alpha, beta;
        calcJointLapAlphaBeta(tau,
                              0.0,
                              0.0,
                              z,
                              alpha,
                              beta);

        return (alpha + beta * Maths::square(initialVol));
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}
    
/* Forward Starting integrated variance */
Complex VolSV::cumulant(const FwdStFourierProcessIntVar& process,
                         const FwdStFourierProductIntVar& product, 
                         const Complex z, 
                         const DateTime& matDate) const {        
    static const string method = "VolSV::cumulant";
    try{

        /* Time fractions from today till start of option and then from
           start of option till maturity */
        double t = timeMetric->yearFrac(baseDate,
                                        product.getStartDate());
        double tau = timeMetric->yearFrac(product.getStartDate(),
                                          matDate);


        /* Calculate alpha, beta component from start date till maturity */
        Complex alpha_tT, beta_tT;          
        calcJointLapAlphaBeta(tau,
                              0.0,
                              0.0,
                              z,
                              alpha_tT,
                              beta_tT);

        /* Then, from today till start date */
        Complex alpha_t, beta_t;
        calcJointLapAlphaBeta(t,
                              0.0,
                              beta_tT,
                              0.0,
                              alpha_t,
                              beta_t);

        return (alpha_t + beta_t * Maths::square(initialVol) + alpha_tT);

    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}
   
/* Started Quadratic Variation */
Complex VolSV::cumulant(const StFourierProcessQuadVar& process,
                         const StFourierProductQuadVar& product, 
                         const Complex z, 
                         const DateTime& matDate) const {    
    static const string method = "VolSV::cumulant";
    try{

        double tau = timeMetric->yearFrac(baseDate,
                                   matDate);

        Complex alpha, beta;
        calcJointLapAlphaBeta(tau,
                              0.0,
                              0.0,
                              z,
                              alpha,
                              beta);
        return (alpha + beta * Maths::square(initialVol));
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

/* Forward Starting quadratic variance */
Complex VolSV::cumulant(const FwdStFourierProcessQuadVar& process,
                         const FwdStFourierProductQuadVar& product, 
                         const Complex z, 
                         const DateTime& matDate) const {        
    static const string method = "VolSV::cumulant";
    try{

        /* Time fractions from today till start of option and then from
           start of option till maturity */
        double t = timeMetric->yearFrac(baseDate,
                                        product.getStartDate());
        double tau = timeMetric->yearFrac(product.getStartDate(),
                                          matDate);


        /* Calculate alpha, beta component from start date till maturity */
        Complex alpha_tT, beta_tT;          
        calcJointLapAlphaBeta(tau,
                              0.0,
                              0.0,
                              z,
                              alpha_tT,
                              beta_tT);

        /* Then, from today till start date */
        Complex alpha_t, beta_t;
        calcJointLapAlphaBeta(t,
                              0.0,
                              beta_tT,
                              0.0,
                              alpha_t,
                              beta_t);

        return (alpha_t + beta_t * Maths::square(initialVol) + alpha_tT);

    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

double VolSV::expectation(const StFourierProcessQuadVar& process,
                          const StFourierProductQuadVar& product, 
                          const DateTime& matDate) const {    
    static const string method = "VolSV::cumulant";
    try{

        double mat = timeMetric->yearFrac(baseDate,
                                          matDate);

        // needs to be delegated to Heston class
        // will do, when I extend this code to other vol classes (Regis)
        double mVar = meanVol * meanVol;
        double tau = meanReversRate * mat;
        //double mrr = meanReversRate;
        double decay = 1.0;
        if (Maths::isPositive(tau)){
            decay = (1.0 - exp(-tau)) / tau;
        }
        return (mVar + (initialVol * initialVol - mVar) * decay);
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

double VolSV::expectation(const FwdStFourierProcessQuadVar& process,
                          const FwdStFourierProductQuadVar& product, 
                          const DateTime& matDate) const {    
    static const string method = "VolSV::cumulant";
    try{
        /* Time fractions from today till start of option and then from
           start of option till maturity */
        double t = timeMetric->yearFrac(baseDate,
                                        product.getStartDate());
        double mat = timeMetric->yearFrac(product.getStartDate(),
                                          matDate);

        // needs to be delegated to Heston class
        // will do, when I extend this code to other vol classes (Regis)
        double mVar = meanVol * meanVol;
        double tau = meanReversRate * mat;
        //double mrr = meanReversRate;
        double decay = 1.0;
        if (Maths::isPositive(tau)){
            decay = (1.0 - exp(-tau)) / tau * exp(-meanReversRate * t);
        }
        return (mVar + (initialVol * initialVol - mVar) * decay);
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

// VarSwap in the Heston Model
// Returns time 0 value for VarSwap between time zero and maturity
double VolSV::varSwap(const DateTime& valueDate, const DateTime& maturity) const {
    static const string method = "VolSV::varSwap";
    try {
        double initialVar = initialVol * initialVol;
        double T = timeMetric->yearFrac(valueDate,maturity);
        if (Maths::isZero(meanReversRate) || Maths::isZero(T)) {
            return initialVar;
        } else {
            double meanVar = meanVol*meanVol;
            double temp =  meanVar + (initialVar - meanVar) * (1.0-exp(-meanReversRate*T)) / (meanReversRate*T);
            return temp;
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// ImpliedVol in the Heston Model (approximation)
// Returns time 0 value for implied ATM vol between maturity1 and maturity2
double VolSV::impliedSqVol(const DateTime& valueDate, 
                           const DateTime& maturity1, 
                           const DateTime& maturity2) const {
    static const string method = "VolSV::impliedSqVol";
    try {
        // some variables needed for both cases, ie meanReversRate = 0 and meanReversRate > 0
        double initialVar = initialVol * initialVol;
        double meanVar = meanVol*meanVol;
        double T = timeMetric->yearFrac(valueDate,maturity1);
        double delta = timeMetric->yearFrac(maturity1,maturity2);
        // degenerated case
        if (Maths::isZero(delta)) {
            return (initialVar - meanVar) * exp(-meanReversRate*T) + meanVar;
        }
        // 
        double addFct1 = 0.0;
        double addFct3 = 0.0;
        double addFct4 = 0.0;
        double VarSwapHestonFwd = 0.0;

        // now, two different cases ... kappa = 0 and kappa > 0 
        if (meanReversRate<ACC) {
            double delta3 = delta*delta*delta;
            addFct1 = correlation * delta3 * initialVar / 2.0;
            addFct3 = delta3 * initialVar / 3.0;
            addFct4 = correlation*correlation*delta3 * initialVar / 6.0;
            VarSwapHestonFwd = initialVar;
        } else {
            // time zero expectation of V_T
            double expInitialVar = (initialVar - meanVar) * exp(-meanReversRate*T) + meanVar; 
            
            // time zero expectation of fwd varswap between T and T+delta (=F_c)
            VarSwapHestonFwd = meanVar
                + (expInitialVar - meanVar) * (1-exp(-meanReversRate*delta)) / (meanReversRate*delta);
            
            // some more auxiliary variables:
            double meanReversRate2 = meanReversRate * meanReversRate;
            double meanReversRate3 = meanReversRate2 * meanReversRate;
            double kappaDelta = meanReversRate*delta;
            
            // now F_1, F_3, F_4
            addFct1 = (-1.0) * correlation / meanReversRate2 * exp(-kappaDelta) * (
                    meanVar * (2.0 * (exp(kappaDelta)-1.0) - kappaDelta*(1+exp(kappaDelta)))
                        + expInitialVar * (1.0+kappaDelta-exp(kappaDelta)) );
            addFct3 = (-1.0) / (4.0*meanReversRate3) * (
                    meanVar * (-2.0*kappaDelta + 5.0 - 4.0*exp(-kappaDelta)*(1.0+kappaDelta) - exp(-2.0*kappaDelta))
                        + expInitialVar * (-2.0 + 4.0*exp(-kappaDelta)*kappaDelta + 2.0*exp(-2.0*kappaDelta)) );
            addFct4 = (-1.0) * correlation*correlation / (2.0*meanReversRate3) * exp(-kappaDelta) * (
                    meanVar * (-6.0 + exp(kappaDelta) * (6.0-2.0*kappaDelta) - 4.0*kappaDelta - kappaDelta*kappaDelta)
                        + expInitialVar * (2.0 * (1.0+kappaDelta) - 2.0*exp(kappaDelta) + kappaDelta*kappaDelta) );
        }
        // put them all together
        double term1 = VarSwapHestonFwd;
        double term2 = addFct1/(2.0*delta);
        double term3 = - addFct3 / (2.0 * VarSwapHestonFwd * delta*delta)*(1.0+VarSwapHestonFwd*delta/4.0)
            - addFct4 / (VarSwapHestonFwd * delta*delta) * (1.0-VarSwapHestonFwd*delta/4.0)
            + addFct1*addFct1 / (VarSwapHestonFwd*VarSwapHestonFwd * delta*delta*delta) * (3.0/4.0 + VarSwapHestonFwd * delta / 16.0);
        return term1 + volVol * term2 + volVol*volVol * term3;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Calculates the components alpha and beta that appear in 
    the exponent of the time-t joint Bilateral Laplace 
    of X_T = (Y_T, V_T, I_T) where Y_T = ln(S_T / F(0, T)) is the dimension-less
    log spot at time T, V_T is the instantaneous variance at time T
    and I_T is the instantaneous variance from 0 to time T
    The Bilateral Laplace transform is of the form
        exp(alpha(tau, u1, u2, u3) 
            + u1 * Y_t 
            + beta(tau, u1, u2, u3) * V_t 
            + u3 * I_t)
    where tau = T - t and the complex numbers u1, u2, u3 are the frequencies 
    wrt Y_T, V_T and I_T, respectively. */
void VolSV::calcJointLapAlphaBeta(double         tau,    // tau = T - t (in years)
                                   const Complex& u1,
                                   const Complex& u2,
                                   const Complex& u3,
                                   Complex&       alpha,
                                   Complex&       beta) const{
    static const string method = "VolSV::calcJointLapAlphaBeta";

    try{    // in case Complex operations fail for some reason
        /* Calculate Heston contribution */
        heston->calcAlphaBeta(tau,
                              u1,
                              u2,
                              u3,
                              alpha,
                              beta);
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

CClassConstSP const VolSV::TYPE =
CClass::registerClassLoadMethod("VolSV", typeid(VolSV), load);

CClassConstSP const VolSV::SVVolParam::TYPE =
CClass::registerClassLoadMethod("VolSV::SVVolParam",
                                typeid(SVVolParam), load);

/** Needed for IAdjustable interface. Returns market data name for vol */
string VolSV::getName() const{
    return CVolBaseParamSurface::getName();
}

/** Called after adjustments have been made to fields (eg calibrator) */
void VolSV::fieldsUpdated(const CFieldArray& fields){
    update();
}

/** Called after adjustments have been made to fields */
void VolSV::update(){
    // refresh heston
    heston = HestonSP(new Heston(initialVol,
                                 meanVol,
                                 meanReversRate,
                                 volVol,
                                 correlation));
}

/** Given an array of gaussian deviates, populate the arrays 'instVars' 
    and 'integratedVars' with simulated values using an Euler
    discretization scheme */
void VolSV::MCGenerateVarPaths(Int2Type<MCParams::VarDSType::EULER>, 
                               const DoubleArray& tradYears,
                               const double*      deviates,
                               DoubleArray&       instVars,
                               DoubleArray&       integratedVars) const{ 


#ifdef DEBUG
    ASSERT(instVars.size() == tradYears.size());
    ASSERT(instVars.size() == integratedVars.size());
#endif
    int nbSteps = tradYears.size();
    double meanVar = Maths::square(meanVol);
    instVars[0] = Maths::square(initialVol);
    integratedVars[0] = 0.0;
    for (int iStep = 1, iLastStep = 0; iStep < nbSteps; 
         iLastStep = iStep++) {

                double dt = tradYears[iStep] - tradYears[iLastStep];
                // drift contribution
                double drift = meanReversRate * (meanVar - instVars[iLastStep]) * dt;
                // diffusion contribution
                double instVarPlus = Maths::max(0.0, instVars[iLastStep]);
                double diffusion = volVol * sqrt(instVarPlus * dt) * deviates[iLastStep];
                // compute inst variance
                instVars[iStep] = instVars[iLastStep] + drift + diffusion;
                // compute integrated var using an Euler scheme too
                integratedVars[iStep] = integratedVars[iLastStep] + instVarPlus * dt;

                
            
            
    }
}
  
#if 1
/** Given an array of gaussian deviates, populate the arrays 'instVars' 
    and 'integratedVars' with simulated values using a variable transform
    scheme together with an Euler discretization scheme */
void VolSV::MCGenerateVarPaths(Int2Type<MCParams::VarDSType::VAR_TRANSFORM_EULER>,
                               const DoubleArray& tradYears,
                               const double*      deviates,
                               DoubleArray&       instVars,
                               DoubleArray&       integratedVars) const{


#ifdef DEBUG
    ASSERT(instVars.size() == tradYears.size());
    ASSERT(instVars.size() == integratedVars.size());
#endif
    int nbSteps = tradYears.size();
    double alpha = meanReversRate * Maths::square(meanVol) 
                   - Maths::square(volVol) / 4.0;
    double vol = initialVol;
    instVars[0] = Maths::square(vol);
    integratedVars[0] = 0.0;
    for (int iStep = 1, iLastStep = 0; iStep < nbSteps; 
         iLastStep = iStep++) {

                double dt = tradYears[iStep] - tradYears[iLastStep];
                // integrated var is function of last step values only
                integratedVars[iStep] = integratedVars[iLastStep] + instVars[iLastStep] * dt;

                // Euler approx for the vol
                vol = sqrt(instVars[iLastStep]);
                ASSERT(!Maths::isZero(vol));
                double drift = 0.5 * (alpha / vol - meanReversRate * vol) * dt;
                double diffusion = 0.5 * volVol * sqrt(dt) * deviates[iLastStep];
                vol += drift + diffusion;
                instVars[iStep] = Maths::square(vol);

                
                
    }
}
#else
/** Given an array of gaussian deviates, populate the arrays 'instVars' 
    and 'integratedVars' with simulated values using a variable transform
    scheme together with an Euler discretization scheme */
void VolSV::MCGenerateVarPaths(Int2Type<MCParams::VarDSType::VAR_TRANSFORM_EULER>,
                               const DoubleArray& tradYears,
                               const double*      deviates,
                               DoubleArray&       instVars,
                               DoubleArray&       integratedVars) const{


#ifdef DEBUG
    ASSERT(instVars.size() == tradYears.size());
    ASSERT(instVars.size() == integratedVars.size());
#endif
    int nbSteps = tradYears.size();
    double alpha = meanReversRate * Maths::square(meanVol) 
                   - Maths::square(volVol) / 4.0;
    ASSERT(Maths::isZero(alpha));
    double vol = initialVol;
    instVars[0] = Maths::square(vol);
    integratedVars[0] = 0.0;
    for (int iStep = 1, iLastStep = 0; iStep < nbSteps; 
         iLastStep = iStep++) {
             
                double dt = tradYears[iStep] - tradYears[iLastStep];
                // integrated var is function of last step values only
                integratedVars[iStep] = integratedVars[iLastStep] + instVars[iLastStep] * dt;

                // vol is an OU process, and therefore can be simulated exactly
                vol = sqrt(instVars[iLastStep]);
                double expEffectiveTime = exp(-0.5 * meanReversRate * dt);
                double sqrtVar = 0.5 * volVol
                         * sqrt((1.0 - Maths::square(expEffectiveTime))
                                / meanReversRate);
                vol = expEffectiveTime * vol + sqrtVar * deviates[iLastStep];
                instVars[iStep] = Maths::square(vol);

                
                
    }
}
#endif

template <>
const CClassConstSP VolSVConvert::TYPE = 
CClass::registerInterfaceLoadMethod("MarketDataConvert::IConvert<VolSV>", typeid(VolSVConvert), load);

/* external symbol to allow class to be forced to be linked in */
bool VolSVLinkIn(){
    return (VolSV::TYPE != 0) && (VolSVConvert::TYPE != 0);  
}

DRLIB_END_NAMESPACE
