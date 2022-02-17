//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolAJDSuper.cpp
//
//   Description : Superposition of Affine Jump-Diffusions
//
//   Date        : 10 April 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolAJDSuper.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

void VolAJDSuper::ISuperposable::load(CClassSP& clazz){
    REGISTER_INTERFACE(VolAJDSuper::ISuperposable, clazz);
    EXTENDS(IObject);
}

CClassConstSP const VolAJDSuper::ISuperposable::TYPE = 
CClass::registerInterfaceLoadMethod(
    "VolAJDSuper::ISuperposable", typeid(VolAJDSuper::ISuperposable), load);

void VolAJDSuper::AJDSuperVolParam::ComputeImpVol(const CVolBase*          vol,
                                                  const CLatticeDouble&    strikes,
                                                  const DateTimeArray&     maturities,
                                                  CLatticeDouble&          impV) const{
    // turn the vol into what we must have
    const VolAJDSuper* myVol = static_cast<const VolAJDSuper *>(vol);
    // then just pass through the parameterised vol
    myVol->ComputeImpVol(strikes, maturities, impV);
}

// work around for msvc 7
typedef VolAJDSuper::SuperposableArray VolAJDSuperSuperposableArray;
DEFINE_TEMPLATE_TYPE_WITH_NAME("VolAJDSuper::ISuperposableArray", VolAJDSuperSuperposableArray);

/*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
VolSurface* VolAJDSuper::AJDSuperVolParam::spotVolSurfaceFromStrikes(
        const CVolBase*       vol,
        const CDoubleArray&   strikes) const{
    // turn the vol into what we must have
    const VolAJDSuper* myVol = 
        static_cast<const VolAJDSuper *>(vol);
    // then just pass through the parameterised vol
    return myVol->spotVolSurfaceFromStrikes(strikes);
}

void VolAJDSuper::validatePop2Object(){
    static const string method("VolAJDSuper::validatePop2Object");
    try {
        int nbVols = vols.size();
        // if names present, must be same size as vols
        if (names.size() != nbVols){
            throw ModelException(method,
                                 "#names should be equal to #vols; got "
                                 + Format::toString(names.size())
                                 + " and "
                                 + Format::toString(nbVols)
                                 +", respectively");
        }
        // for each vol..
        for (int iVol = 0; iVol < nbVols; ++iVol){            
            if (names[iVol].empty()){
                throw ModelException(method,
                                     Format::toString(iVol + 1)
                                     + "-th name is missing");
            }
            CClassConstSP clazz(CClass::forName(vols[iVol]));
            // make sure it's a superposable vol
            if (!ISuperposable::TYPE->isAssignableFrom(clazz)){
                throw ModelException(method,
                                     Format::toString(iVol + 1)
                                     + "-th vol (" 
                                     + clazz->getName()
                                     + ") is not a superposable AJD vol");
            }
        }
        // #weights should #vols - 1. Last weight will be calculated so 
        // the sum is 1.0 (that will save the day in terms of calibration when 
        // we have 2 vols. If #vols > 2, we're screwed - we need an optimizer
        // with constraints (ouch!))       
        if (sqWeights.size() != nbVols - 1){
            throw ModelException(method,
                                 "#sqWeights should be equal to #vols - 1; got "
                                 + Format::toString(sqWeights.size())
                                 + " and "
                                 + Format::toString(nbVols)
                                 +", respectively");
        }
        // check weights are in range [0, 1]
        Calibrator::IAdjustable::checkRange(this);
        // size 'used' weights array and fill them in
        usedSqWeights.resize(nbVols);
        usedWeights.resize(nbVols);
        update();
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

void VolAJDSuper::getMarket(const IModel*     model, 
                            const MarketData* market){
    static const string method("VolAJDSuper::getMarket");
    try{
        int nbVols = vols.size();
        supervols.resize(nbVols);
        // for each vol...
        for (int iVol = 0; iVol < nbVols; ++iVol){
            // turn the string into the actual object
            CClassConstSP clazz(CClass::forName(vols[iVol]));
            MarketObjectSP mo(market->GetData(names[iVol], clazz));
            // cast vol into a supervol
            supervols[iVol] = SuperposableSP(SuperposableSP::dynamicCast(mo));
            // must get the market of the superposable vol, but in order to 
            // do that the superposable vol must be provided with the actual name
            // of the vol that it is supperposed into
            supervols[iVol]->getMarket(model, market, getName());
        }
        // get market for base class
        VolBaseParam::getMarket(model, market);
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

void VolAJDSuper::ComputeImpVol(const CLatticeDouble&      strikes,
                                const DateTimeArray&       maturities,
                                CLatticeDouble&            impV) const {
    static const string routine("VolAJDSuper::ComputeImpVol");
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
VolSurface* VolAJDSuper::spotVolSurfaceFromStrikes(
    const CDoubleArray&   strikes) const{
    static const string routine("VolAJDSuper::spotVolSurfaceFromStrikes");
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
CVolParam* VolAJDSuper::createVolParam() const{
    return new AJDSuperVolParam();
}

IObject* VolAJDSuper::defaultCtor(){
    return new VolAJDSuper();
}

/** Invoked when Class is 'loaded' */
void VolAJDSuper::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(VolAJDSuper, clazz);
    SUPERCLASS(VolBaseParam);
//    IMPLEMENTS(IVolatilityBS);
//    IMPLEMENTS(IVolatilityDVF);
    IMPLEMENTS(Calibrator::IAdjustable);
    EMPTY_SHELL_METHOD(defaultCtor);
    FIELD(sqWeights, "Array of 'n-1' squared weights");
    FIELD(vols, "Array of 'n' superposable AJD vols");
    FIELD(names, "Array of 'n' market names");

    // transient        
    FIELD(baseDate, "");
    FIELD_MAKE_TRANSIENT(baseDate);
    FIELD(usedSqWeights, "");
    FIELD_MAKE_TRANSIENT(usedSqWeights);
    FIELD(usedWeights, "");
    FIELD_MAKE_TRANSIENT(usedWeights);
    FIELD(supervols, "");
    FIELD_MAKE_OPTIONAL(supervols); // need a TRANSIENT_BUT_ADJUSTABLE macro

    // add our fields and their ranges to central list
    Calibrator::IAdjustable::registerField(
        clazz, "sqWeights",
        new Range(ClosedBoundary(0.0), ClosedBoundary(1.0)));
}

VolAJDSuper::VolAJDSuper(): VolBaseParam(TYPE){}

/*** Build the parameterised vol and cache any values **/
/*** Build the parameterised vol and cache any values **/
void VolAJDSuper::buildCache() {
    const VolSurface* backbone = getBackboneSurface();
    baseDate = backbone->getBaseDate();
}

/** Returns the number of factors, the first one of which
    being X^1 = Y the log-spot. */
int VolAJDSuper::getNbFactors() const{
    static const string method("VolAJDSuper::getNbFactors");
    try{
        // start off with 1 factor for the log spot
        int nbFactors = 1;
        for (int iVol = 0; iVol < supervols.size(); ++iVol){
            nbFactors += supervols[iVol]->getNbFactors() - 1;   // remove the log-spot factor
        }
        return nbFactors;
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

/** Returns the initial values Y_0,..., X^i_0,...
    Naturally, Y_0 = 0. */
void VolAJDSuper::getInitialValues(DoubleArray& initVals) const{
    static const string method("VolAJDSuper::getIntialValues");
    try{
        initVals[0] = 0.0;  // log spot
        int index = 1;
        for (int iVol = 0; iVol < supervols.size(); ++iVol){
            const ISuperposable& super = *supervols[iVol];     // for ease
            int nbFactors = super.getNbFactors();             
            for (int iFactor = 1; iFactor < nbFactors; ++iFactor, ++index){
                initVals[index] = super.getInitialValue(iFactor);
            }
        }
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

void VolAJDSuper::calcCumulantComponents(const DateTime&     fromDate,  // time t
                                         const DateTime&     toDate,    // time T
                                         const ComplexArray& inu,         // frequencies
                                         Complex&            inalpha,
                                         ComplexArray&       inbeta) const{
    static const string method("VolAJDSuper::calcCumulantComponents");
    try{
        int index = 1, index2 = 1;
        inalpha = 0.0;
        inbeta[0] = inu[0];     // log spot
        for (int iVol = 0; iVol < supervols.size(); ++iVol){
            const ISuperposable& super = *supervols[iVol];     // for ease
            int nbFactors = super.getNbFactors(); 
            ComplexArray u(nbFactors);
            Complex alpha;
            ComplexArray betas(nbFactors);
            u[0] = inu[0];
            int iFactor = 1;
            for (; iFactor < nbFactors; ++iFactor, ++index){
                u[iFactor] = inu[index];
            }
            double weight = usedWeights[iVol];
            super.calcCumulantComponents(fromDate,
                                         toDate,
                                         u,
                                         weight,
                                         alpha,
                                         betas);
            inalpha += alpha;
            iFactor = 1;
            for (; iFactor < nbFactors; ++iFactor, ++index2){
                inbeta[index2] = betas[iFactor];
            }
        }
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

//////////// class FourierProcessAJDSuper ///////////////////

/* Started Log Return */
Complex VolAJDSuper::scalelessCumulant(const StFourierProcessLogRtn& process,
                                       const StFourierProductLogRtn& product, 
                                       const Complex& z, 
                                       const DateTime& matDate) const{
    static const string method = "VolAJDSuper::scalelessCumulant";
    try{
        int nbFactors = getNbFactors();
        ComplexArray u(nbFactors);
        u[0] = z;
        int iFactor = 1;
        for (; iFactor < nbFactors; ++iFactor){
            u[iFactor] = 0.0;
        }
        Complex alpha;
        ComplexArray betas(nbFactors);
        calcCumulantComponents(baseDate,
                               matDate,
                               u,
                               alpha,
                               betas);
        DoubleArray initVals(nbFactors);
        getInitialValues(initVals);
        Complex cumulant = alpha;
        for (iFactor = 0; iFactor < nbFactors; ++iFactor){
            cumulant += betas[iFactor] * initVals[iFactor];
        }
        return cumulant;
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

/* Fwd starting Log Return  */
Complex VolAJDSuper::scalelessCumulant(const FwdStFourierProcessLogRtn& process,
                                       const FwdStFourierProductLogRtn& product, 
                                       const Complex& z, 
                                       const DateTime& matDate) const{
    static const string method = "VolAJDSuper::scalelessCumulant";
    try{
        int nbFactors = getNbFactors();
        ComplexArray u(nbFactors);
        u[0] = z;
        int iFactor = 1;
        for (; iFactor < nbFactors; ++iFactor){
            u[iFactor] = 0.0;
        }
        Complex alpha;
        ComplexArray betas(nbFactors);
        calcCumulantComponents(product.getStartDate(),
                               matDate,
                               u,
                               alpha,
                               betas);
        Complex cumulant = alpha;
        u[0] = 0.0;        
        for (iFactor = 1; iFactor < nbFactors; ++iFactor){
            u[iFactor] = betas[iFactor];
        }
        calcCumulantComponents(baseDate,
                               product.getStartDate(),
                               u,
                               alpha,
                               betas);
        cumulant += alpha;
        DoubleArray initVals(nbFactors);
        getInitialValues(initVals);
        for (iFactor = 0; iFactor < nbFactors; ++iFactor){
            cumulant += betas[iFactor] * initVals[iFactor];
        }
        return cumulant;
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

#if 0
/* Started integrated variance */
Complex VolAJDSuper::cumulant(const StFourierProcessIntVar& process,
                         const StFourierProductIntVar& product, 
                         const Complex z, 
                         const DateTime& matDate) const {    
    static const string method = "VolAJDSuper::cumulant";
    try{

        double tau = timeMetric->yearFrac(baseDate,
                                       matDate);


        Complex alpha, betas;
        calcJointLapAlphaBeta(tau,
                              0.0,
                              0.0,
                              z,
                              alpha,
                              betas);

        return (alpha + betas * Maths::square(initialVol));
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}
    
/* Forward Starting integrated variance */
Complex VolAJDSuper::cumulant(const FwdStFourierProcessIntVar& process,
                         const FwdStFourierProductIntVar& product, 
                         const Complex z, 
                         const DateTime& matDate) const {        
    static const string method = "VolAJDSuper::cumulant";
    try{

        /* Time fractions from today till start of option and then from
           start of option till maturity */
        double t = timeMetric->yearFrac(baseDate,
                                        product.getStartDate());
        double tau = timeMetric->yearFrac(product.getStartDate(),
                                          matDate);


        /* Calculate alpha, betas component from start date till maturity */
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
Complex VolAJDSuper::cumulant(const StFourierProcessQuadVar& process,
                         const StFourierProductQuadVar& product, 
                         const Complex z, 
                         const DateTime& matDate) const {    
    static const string method = "VolAJDSuper::cumulant";
    try{

        double tau = timeMetric->yearFrac(baseDate,
                                   matDate);

        Complex alpha, betas;
        calcJointLapAlphaBeta(tau,
                              0.0,
                              0.0,
                              z,
                              alpha,
                              betas);

        // Merton
        Complex merton;
        merton = -crashRate*tau*(1. -1./sqrt(1.-2.* z *crashSizeUncertainty*crashSizeUncertainty)*
            exp(z*(pow((crashSizeMean - .5*crashSizeUncertainty*crashSizeUncertainty), 2.)/
            (1.-2.*z*crashSizeUncertainty*crashSizeUncertainty))));

        alpha += merton;        

        return (alpha + betas * Maths::square(initialVol));
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

double VolAJDSuper::expectation(const StFourierProcessQuadVar& process,
                           const StFourierProductQuadVar& product, 
                           const DateTime& matDate) const {    
    static const string method = "VolAJDSuper::cumulant";
    try{

        double mat = timeMetric->yearFrac(baseDate,
                                          matDate);

        /*the expected change in log(S) under a jump */
        // needs to be delegated to Heston + Merton classes
        // will do, when I extend this code to other vol classes (Regis)
        double sqCrashSizeUncertainty = 
            crashSizeUncertainty * crashSizeUncertainty;
        double jgt = crashSizeMean - 0.5 * sqCrashSizeUncertainty;
        double mVar = meanVol * meanVol;
        double tau = meanReversRate * mat;
        double mrr = meanReversRate;

        AdjustRiskPremium(volRiskPrice, mVar, mrr);

        double result = crashRate * (jgt * jgt + sqCrashSizeUncertainty);

        result += mVar + (initialVol * initialVol - mVar) * (1.0 - exp(-tau)) / tau;
        return result;
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}
#endif

CClassConstSP const VolAJDSuper::TYPE =
CClass::registerClassLoadMethod("VolAJDSuper", typeid(VolAJDSuper), load);

CClassConstSP const VolAJDSuper::AJDSuperVolParam::TYPE =
CClass::registerClassLoadMethod("VolAJDSuper::AJDSuperVolParam",
                                typeid(AJDSuperVolParam), load);

/** Needed for IAdjustable interface. Returns market data name for vol */
string VolAJDSuper::getName() const{
    return CVolBaseParamSurface::getName();
}

/** Called after adjustments have been made to fields */
void VolAJDSuper::fieldsUpdated(const CFieldArray& fields){
    update();
}

/** Called after adjustments have been made to fields */
void VolAJDSuper::update(){
    static const string method("VolAJDSuper::update");
    try{
        double sqWeightSum = 0.0;
        double lastSqWeight = 1.0;
        int iWeight = 0;
        for (; iWeight < sqWeights.size(); ++iWeight){
            double sqWeight = sqWeights[iWeight];
            sqWeightSum += sqWeight;
            lastSqWeight = 1.0 - sqWeightSum;
            if (Maths::isNegative(lastSqWeight)){
                throw ModelException(method,
                                     "sqWeight 1 to sqWeight "
                                     + Format::toString(iWeight + 1)
                                     + " sum up to "
                                     + Format::toString(sqWeightSum)
                                     + ", which is greater than 1.0");
            }
            usedSqWeights[iWeight] = sqWeight;
            usedWeights[iWeight] = sqrt(sqWeight);
        }
        usedSqWeights[iWeight] = lastSqWeight;
        usedWeights[iWeight] = sqrt(lastSqWeight);
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}
 
/* external symbol to allow class to be forced to be linked in */
bool VolAJDSuperLinkIn(){
    return (VolAJDSuper::TYPE != 0);
}

DRLIB_END_NAMESPACE
