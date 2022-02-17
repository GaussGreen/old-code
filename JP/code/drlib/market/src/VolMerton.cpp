//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolMerton.cpp
//
//   Description : 
//
//   Author      : Oliver Brockhaus
//
//   Date        : 03 April 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolMerton.hpp"
#include "edginc/Maths.hpp"
#include "edginc/mathlib.hpp"

#define SIM_CRASH_TIME 1

DRLIB_BEGIN_NAMESPACE

/** Returns a processed vol - which combines the vol market data with the
    instrument data in the volRequest */
CVolProcessed* VolMerton::getProcessedVol(
    const CVolRequest* volRequest,
    const CAsset*      asset) const{
    static const string  routine("VolMerton::getProcessedVol");
    try{
        // intercept VolRequestMerton
        if (VolRequestMerton::TYPE->isInstance(volRequest)){
            VolMertonConstSP vol(VolMertonConstSP::attachToRef(this));
            return new VolProcessedMerton(vol);
        }
        // otherwise delegate to base class
        return VolBaseParam::getProcessedVol(volRequest, asset);

    } catch (exception& e){
        throw ModelException(e, routine, "Exception "+ getName());
    }
}

CVolProcessed* VolMerton::getProcessedVol(
    const CVolRequest* volRequest,
    const CAsset*      eqAsset,
    const FXAsset*     fxAsset,
    const Correlation* eqFXCorr) const{
    // intercept VolRequestMerton
    if (VolRequestMerton::TYPE->isInstance(volRequest)){
        throw ModelException("VolMerton::getProcessedVol",
                             "ccy struck not supported");
    }
    // otherwise delegate to base class
    return VolBaseParam::getProcessedVol(volRequest, 
                                         eqAsset, 
                                         fxAsset, 
                                         eqFXCorr);
}

void VolMerton::MertonVolParam::ComputeImpVol(  const CVolBase*          vol,
                                                const CLatticeDouble&    strikes,
                                                const DateTimeArray&     maturities,
                                                CLatticeDouble&          impV) const{
    // turn the vol into what we must have
    const VolMerton* myVol = static_cast<const VolMerton *>(vol);
    // then just pass through the parameterised vol
    myVol->ComputeImpVol(strikes, maturities, impV);
}

/*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
VolSurface* VolMerton::MertonVolParam::spotVolSurfaceFromStrikes(
    const CVolBase*       vol,
    const CDoubleArray&   strikes) const{
    // turn the vol into what we must have
    const VolMerton* myVol = 
        static_cast<const VolMerton *>(vol);
    // then just pass through the parameterised vol
    return myVol->spotVolSurfaceFromStrikes(strikes);
}

void VolMerton::validatePop2Object(){
    try {
        Calibrator::IAdjustable::checkRange(this);
        // build heston + mertonCrash
        update();
    }
    catch(exception& e){
        throw ModelException(e, "VolMerton::validatePop2Object");
    }
}

void VolMerton::ComputeImpVol(  const CLatticeDouble&      strikes,
                                const DateTimeArray&       maturities,
                                CLatticeDouble&            impV) const {
    static const string routine("VolMerton::ComputeImpVol");
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

/** adjust mean vol and reversion rate for risk premium */
void VolMerton::AdjustRiskPremium(  double volRiskPrice, 
                                    double& mVar, 
                                    double& mReversion) const{
    mVar *= mReversion;
    mReversion += volRiskPrice;
    // should we in fact trap mReversion < 0 cases ? which gives "mean
    // explosion"
    if(Maths::isZero(mReversion))
        throw ModelException("VolMerton::AdjustRiskPremium", "divide by zero");

    mVar /= mReversion;
}

/** calculate a simple expected variance */
double VolMerton::FutureVariance(double mat) const{
    /*the expected change in log(S) under a jump */
    double sqCrashSizeUncertainty 
               = crashSizeUncertainty * crashSizeUncertainty;
    double jgt = crashSizeMean - 0.5 * sqCrashSizeUncertainty;
    double mVar= meanVol * meanVol;
    double tau = meanReversRate * mat;
    double mrr = meanReversRate;

    AdjustRiskPremium(volRiskPrice, mVar, mrr);

    double result = crashRate * (jgt * jgt + sqCrashSizeUncertainty);

    result += mVar + (initialVol * initialVol - mVar) * (1.0 - exp(-tau)) / tau;
    return result;
}

    
/*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
VolSurface* VolMerton::spotVolSurfaceFromStrikes(
    const CDoubleArray&   strikes) const{
    static const string routine("VolMerton::spotVolSurfaceFromStrikes");
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
CVolParam* VolMerton::createVolParam() const{
    return new MertonVolParam();
}

IObject* VolMerton::defaultCtor(){
    return new VolMerton();
}

/** Invoked when Class is 'loaded' */
void VolMerton::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(VolMerton, clazz);
    SUPERCLASS(VolBaseParam);
//    IMPLEMENTS(IVolatilityBS);
//    IMPLEMENTS(IVolatilityDVF);
    IMPLEMENTS(IDynamicsParameter);
    IMPLEMENTS(Calibrator::IAdjustable);
    EMPTY_SHELL_METHOD(defaultCtor);
    FIELD(initialVol, "initialVol");
    FIELD(meanVol, "meanVol");
    FIELD(meanReversRate, "meanReversRate");
    FIELD(crashRate, "crashRate");
    FIELD(crashSizeMean, "crashSizeMean");
    FIELD(crashSizeUncertainty, "crashSizeUncertainty");
    FIELD(volRiskPrice,"coefficient of market price of risk for stochastic vol");
    FIELD_MAKE_OPTIONAL(volRiskPrice);
    FIELD(beta, "beta");
    FIELD_MAKE_OPTIONAL(beta);

    // transient
    FIELD(baseDate, "");
    FIELD_MAKE_TRANSIENT(baseDate);
    FIELD(heston, "");
    FIELD_MAKE_TRANSIENT(heston);
    FIELD(mertonCrash, "");
    FIELD_MAKE_TRANSIENT(mertonCrash);
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
        clazz, "crashRate",
        new Range(MertonCrash::RangeDef::crashRate));
    Calibrator::IAdjustable::registerField(
        clazz, "crashSizeMean", 
        new Range(MertonCrash::RangeDef::crashSizeMean));
    Calibrator::IAdjustable::registerField(
        clazz, "crashSizeUncertainty", 
        new Range(MertonCrash::RangeDef::crashSizeUncertainty));
}

VolMerton::VolMerton(): 
    VolBaseParam(TYPE),
    initialVol(Heston::DefaultVal::initialVol),
    meanVol(Heston::DefaultVal::meanVol),
    meanReversRate(Heston::DefaultVal::meanReversRate),
    crashRate(MertonCrash::DefaultVal::crashRate),
    crashSizeMean(MertonCrash::DefaultVal::crashSizeMean),
    crashSizeUncertainty(MertonCrash::DefaultVal::crashSizeUncertainty),
    volRiskPrice(0.0),
    beta(0.0) {}

/*** Build the parameterised vol and cache any values **/
/*** Build the parameterised vol and cache any values **/
void VolMerton::buildCache() {
    const VolSurface* backbone = getBackboneSurface();
    baseDate = backbone->getBaseDate();
}

///////////////////////////////////////////////////////////////////////////////////
// FourierProcessMerton
///////////////////////////////////////////////////////////////////////////////////

/* Started Log Return */
Complex VolMerton::scalelessCumulant(const StFourierProcessLogRtn& process,
                                  const StFourierProductLogRtn& product, 
                                  const Complex& z, 
                                  const DateTime& matDate) const{
    static const string method = "VolMerton::scalelessCumulant";
    try{
        static DateTime lastMatDate(0, 0);
        static double tau = 0.0;

        if (matDate != lastMatDate) {
            tau = timeMetric->yearFrac(baseDate,
                                       matDate);
            lastMatDate = matDate;
        }

        Complex alpha, beta;
        calcJointLapAlphaBeta(tau,
                              z,
                              0.0,
                              alpha,
                              beta);
        return alpha + beta * Maths::square(initialVol);
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

/* Fwd starting Log Return  */
Complex VolMerton::scalelessCumulant(const FwdStFourierProcessLogRtn& process,
                                  const FwdStFourierProductLogRtn& product, 
                                  const Complex& z, 
                                  const DateTime& matDate) const{
    static const string method = "VolMerton::scalelessCumulant";
    try{
        static DateTime lastMatDate(0, 0);
        static double tau = 0.0;
        static double t = 0.0;

        if (matDate != lastMatDate) {
            /* Time fractions from today till start of option and then from
               start of option till maturity */
            t = timeMetric->yearFrac(baseDate,
                                     product.getStartDate());
            tau = timeMetric->yearFrac(product.getStartDate(),
                                       matDate);
            lastMatDate = matDate;
        }

        /* Calculate alpha, beta component from start date till maturity */
        Complex alpha_tT, beta_tT;
        calcJointLapAlphaBeta(tau,
                              z,
                              0.0,
                              alpha_tT,
                              beta_tT);
        /* Then, from today till start date */
        Complex alpha_t, beta_t;
        calcJointLapAlphaBeta(t,
                              0.0,
                              beta_tT,
                              alpha_t,
                              beta_t);
        return alpha_t + beta_t * Maths::square(initialVol) + alpha_tT;
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

/* Started integrated variance */
Complex VolMerton::cumulant(const StFourierProcessIntVar& process,
                            const StFourierProductIntVar& product, 
                            const Complex z, 
                            const DateTime& matDate) const {    
    static const string method = "VolMerton::cumulant";
    try{
        static DateTime lastMatDate(0, 0);
        static double tau = 0.0;

        if (matDate != lastMatDate) {
            tau = timeMetric->yearFrac(baseDate,
                                       matDate);
            lastMatDate = matDate;
        }

        Complex alpha, beta;
        /*
        calcIntVarLapAlphaBeta(tau,
                               0.0,
                               z,
                               alpha,
                               beta);
        */
        return alpha + beta * Maths::square(initialVol);
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

/* Started integrated variance */
Complex VolMerton::cumulant(const StFourierProcessQuadVar& process,
                         const StFourierProductQuadVar& product, 
                         const Complex z, 
                         const DateTime& matDate) const {    
    static const string method = "VolMerton::cumulant";
    try{
        static DateTime lastMatDate(0, 0);
        static double tau = 0.0;

        if (matDate != lastMatDate) {
            tau = timeMetric->yearFrac(baseDate,
                                       matDate);
            lastMatDate = matDate;
        }

        Complex alpha, beta;
        /*
        calcIntVarLapAlphaBeta(tau,
                               0.0,
                               z,
                               alpha,
                               beta);
        */

        // Merton
        Complex merton;
        merton = -crashRate*tau*(1. -1./sqrt(1.-2.* z *crashSizeUncertainty*crashSizeUncertainty)*
            exp(z*(pow((crashSizeMean - .5*crashSizeUncertainty*crashSizeUncertainty), 2.)/
            (1.-2.*z*crashSizeUncertainty*crashSizeUncertainty))));

        alpha += merton;        

        return alpha + beta * Maths::square(initialVol);
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}
    
/* Forward Starting integrated variance */
Complex VolMerton::cumulant(const FwdStFourierProcessIntVar& process,
                         const FwdStFourierProductIntVar& product, 
                         const Complex z, 
                         const DateTime& matDate) const {        
    static const string method = "VolMerton::cumulant";
    try{
        static DateTime lastMatDate(0, 0);
        static double tau = 0.0;
        static double t = 0.0;

        if (matDate != lastMatDate) {
            /* Time fractions from today till start of option and then from
               start of option till maturity */
            t = timeMetric->yearFrac(baseDate,
                                     product.getStartDate());
            tau = timeMetric->yearFrac(product.getStartDate(),
                                       matDate);
            lastMatDate = matDate;
        }

        /* Calculate alpha, beta component from start date till maturity */
        Complex alpha_tT;
        Complex alpha_t, beta_t;
        /*
        calcIntVarLapAlphaBeta(tau,
                               0.0,
                               z,
                               alpha_tT,
                               beta_tT);

        // Then, from today till start date
        calcIntVarLapAlphaBeta(t,
                               beta_tT,
                               0.0,
                               alpha_t,
                               beta_t);
        */
        
        return alpha_t + beta_t * Maths::square(initialVol) + alpha_tT;

    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}
   
double VolMerton::expectation(const StFourierProcessQuadVar& process,
                           const StFourierProductQuadVar& product, 
                           const DateTime& matDate) const {    
    static const string method = "VolMerton::cumulant";
    try{
        static DateTime lastMatDate(0, 0);
        static double mat = 0.0;

        if (matDate != lastMatDate) {
            mat = timeMetric->yearFrac(baseDate,
                                       matDate);
            lastMatDate = matDate;
        }

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

/* Calculate the components alpha and beta that appear in the exponent 
   of the time-t joint Laplace transform of X_T = (Y_T, V_T) in the 
   Heston Merton model where Y_T = ln(S_T / F(0, T))
   The Bilateral Laplace transform is of the form
        exp(alpha(tau, u1, u2) + u1 * Y_t + beta(tau, u1, u2) * V_t)
   where tau = T - t and the complex u1, u2 are the frequencies wrt Y_T and 
   V_T, respectively. */
void VolMerton::calcJointLapAlphaBeta(  double         tau,    // tau = T - t (in years)
                                        const Complex& u1,
                                        const Complex& u2,
                                        Complex&       alpha,
                                        Complex&       beta) const{
    static const string method = "VolMerton::calcJointLapAlphaBeta";

    try{    // in case Complex operations fail for some reason
        /* Calculate Heston contribution */
        heston->calcAlphaBeta(tau,
                              u1,
                              u2,
                              0.0,
                              alpha,
                              beta);

        /* Calculate Merton jump contribution to alpha */
        alpha += mertonCrash->calcAlpha(tau, u1);
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

CClassConstSP const VolMerton::TYPE =
CClass::registerClassLoadMethod("VolMerton", typeid(VolMerton), load);

CClassConstSP const VolMerton::MertonVolParam::TYPE =
CClass::registerClassLoadMethod("VolMerton::MertonVolParam",
                                typeid(MertonVolParam), load);

/** Needed for IAdjustable interface. Returns market data name for vol */
string VolMerton::getName() const{
    return CVolBaseParamSurface::getName();
}

/** Needed for IAdjustable interface. Returns market data name for vol */
DateTime VolMerton::getBaseDate() const{
    return baseDate;
}

/** Needed for IAdjustable interface. Returns market data name for vol */
//DateTime VolMerton::getBaseDate() const{
//    return CVolBaseParamSurface::getBaseDate();
//}

/** Called after adjustments have been made to fields (eg calibrator) */
void VolMerton::fieldsUpdated(const CFieldArray& fields){
    update();
}

/** Called after adjustments have been made to fields */
void VolMerton::update(){
    // refresh heston
    heston = HestonSP(new Heston(initialVol,
                                 meanVol,
                                 meanReversRate,
                                 0 /*volVol*/,
                                 0 /*correlation*/));
    // refresh mertonCrash
    mertonCrash = MertonCrashSP(new MertonCrash(crashRate,
                                                crashSizeMean,
                                                crashSizeUncertainty));
}
 

/* external symbol to allow class to be forced to be linked in */
bool VolMertonLinkIn(){
    return (VolMerton::TYPE != 0);
}

///////////////////////////////////////////////////////////////////////////////////
// VolProcessedMerton
///////////////////////////////////////////////////////////////////////////////////

// constructor
VolProcessedMerton::VolProcessedMerton(const VolMertonConstSP& vol):
CObject(TYPE),
volMerton(vol){
#if 0
    initialVol          = vol->initialVol;
    meanVol             = vol->meanVol;
    meanReversRate      = vol->meanReversRate;
    crashRate           = vol->crashRate;
    crashSizeMean       = vol->crashSizeMean;
    crashSizeUncertainty= vol->crashSizeUncertainty;
    timeMetric          = vol->timeMetric;
    beta                = vol->beta;
    name                = vol->getName();
    baseDate            = vol->getBaseDate();
#endif
    lnOnePlusMeanMinusHalfVar
                        = log(1.0 + volMerton->crashSizeMean) 
                          - 0.5*Maths::square(volMerton->crashSizeUncertainty);

}

/** calculates the trading time between two dates */
double VolProcessedMerton::calcTradingTime( const DateTime &date1, 
                                            const DateTime &date2) const 
{
    return volMerton->timeMetric->yearFrac(date1, date2);
}

/** identifies the market data name of the volatility */
string VolProcessedMerton::getName() const 
{
    return name;
}

/** returns the base date of the vol */
DateTime VolProcessedMerton::getBaseDate() const{
    return volMerton->baseDate;
}

/** retrieve time measure for the vol */
TimeMetricConstSP VolProcessedMerton::GetTimeMetric() const
{
    return volMerton->timeMetric;
}

/** Calculates the number of jumps n such that 
    the probability of more than n jumps between dates date1 and date2 
    is less than epsilon. */
void VolProcessedMerton::Quantile(
    const DateTime &date1,
    const DateTime &date2,
    double         epsilon,
    int            maxJumps,
    double*        quantile,
    int*           nbJumps ) const
{

    double dt       = volMerton->timeMetric->yearFrac(date1,date2);
    double proba    = exp(-dt*volMerton->crashRate);
    *quantile = 0;
    int    iStep;

    for( iStep=1; iStep<=maxJumps; iStep++ )
    {
        *quantile += proba;
        proba *= volMerton->crashRate * dt / (double)iStep;
        if( *quantile>1-epsilon ) break;
    }
    *nbJumps = iStep;
}

/** Calculates drift factor between 2 dates */
double VolProcessedMerton::CalcDrift(
    const DateTime &date1,
    const DateTime &date2) const
{
    double dt   = volMerton->timeMetric->yearFrac(date1,date2);
    double drift= exp( - volMerton->crashSizeMean * volMerton->crashRate * dt 
                       - 0.5 * this->CalcVar(date1,date2) );
    return drift;
}

/** Calculates variance between 2 dates */
double VolProcessedMerton::CalcVar( 
    const DateTime &date1,
    const DateTime &date2) const
{
    double dt1  = volMerton->timeMetric->yearFrac(volMerton->baseDate,date1);
    double var1 = volMerton->meanReversRate > 0 ? 
                  volMerton->meanVol*volMerton->meanVol*dt1
                  + ( volMerton->initialVol*volMerton->initialVol - volMerton->meanVol*volMerton->meanVol ) 
                  * ( 1 - exp(-volMerton->meanReversRate*dt1) ) 
                  / volMerton->meanReversRate :
                  volMerton->initialVol*volMerton->initialVol*dt1;
    double dt2  = volMerton->timeMetric->yearFrac(volMerton->baseDate,date2);
    double var2 = volMerton->meanReversRate > 0 ? 
                  volMerton->meanVol*volMerton->meanVol*dt2
                  + ( volMerton->initialVol*volMerton->initialVol - volMerton->meanVol*volMerton->meanVol ) 
                  * ( 1 - exp(-volMerton->meanReversRate*dt2) ) 
                  / volMerton->meanReversRate :
                  volMerton->initialVol*volMerton->initialVol*dt2;
    return var2-var1;
}

/** Calculates variance at t=time (in years)  */
double VolProcessedMerton::CalcVar( 
    const double time) const
{
    return volMerton->meanVol*volMerton->meanVol
           + ( volMerton->initialVol*volMerton->initialVol - volMerton->meanVol*volMerton->meanVol ) 
           * exp(-volMerton->meanReversRate*time);
}

/** Calculates volatility between 2 dates */
double VolProcessedMerton::CalcVol( 
    const DateTime &date1,
    const DateTime &date2) const
{
    double dt   = volMerton->timeMetric->yearFrac(date1,date2);
    double var  = CalcVar(date1,date2);
    return sqrt(var/dt);
}

/** Calculates jump factor. */
double VolProcessedMerton::CalcJump(double noise,
                                    int    numJumps) const
{
    double expJump = 1;
    if( numJumps>0 ) {
        expJump *= exp( (double)numJumps * ( 
                        lnOnePlusMeanMinusHalfVar
                        + volMerton->crashSizeUncertainty * noise / sqrt( (double)numJumps ) 
                   ) );
    }
    return expJump;
}

/** Calculates variance between a series of dates. If the dateList
        has n dates in it, n-1 vars will be calculated. */
void VolProcessedMerton::CalcVar(const DateTimeArray& dateList,
                                 TCalcType            calcType, 
                                 CDoubleArray&        vars) const
{
    int numDates = dateList.size();
    static char routine[] = "VolProcessedMerton::CalcVar";
    if (numDates < 2 ){
        throw ModelException(routine, "Must supply at least two dates");
    }
    if (vars.size() < numDates-1){
        throw ModelException(routine, "Double array too short");
    }
    for (int i = 0; i < numDates-1; i++){
        vars[i] = CalcVar(dateList[i], dateList[i+1]);
    }
    switch (calcType){
        int i; // MSVC is broken - doesn't support separate variables in loop
    case forward:
        // do nothing
        break;
    case fromFirst:
        for (i = 1; i < numDates-1; i++){
            vars[i] += vars[i-1];
        }
        break;
    case toLast:
        for (i = numDates-2; i >= 0; i--){
            vars[i] += vars[i+1];
        }
        break;
    default:
        throw ModelException(routine, "Unknown calculate type");
    }
    return;
}

/** Calculates variance beginning at dateFrom. If the dateList
    has n dates in it, n variances will be calculated. */
void VolProcessedMerton::CalcVar(const DateTime&      dateFrom,
                                 const DateTimeArray& datesTo,
                                 TCalcType            calcType, 
                                 CDoubleArray&        vars) const
{
    int numDates = datesTo.size();
    DateTimeArray dateSeries(datesTo.size()+1);

    dateSeries[0] = dateFrom;
    for(int i=0;i<numDates;++i)
    {
        dateSeries[i+1] = datesTo[i];
    }
    CalcVar(dateSeries, calcType, vars);
}


/** Calculates vols between a series of dates */
void VolProcessedMerton::CalcVol(const DateTimeArray& dateList, 
                                 TCalcType            calcType, 
                                 CDoubleArray&        vols) const
{
    int numDates = dateList.size();
    static char routine[] = "VolProcessedMerton::CalcVol";
    if (numDates < 2 && vols.size() > 0 ) {
        throw ModelException(routine, "Must supply at least two dates");
    }
    if (vols.size() < numDates-1){
        throw ModelException(routine, "Double array too short");
    }
    switch (calcType){
        int i; // MSVC is broken - doesn't support separate variables in loop
    case forward:
        for (i = 0; i < numDates-1; i++){
            vols[i] = CalcVol(dateList[i], dateList[i+1]);
        }
        break;
    case fromFirst:
        // probably not the most efficient
        for (i = 1; i < numDates; i++){
            vols[i-1] = CalcVol(dateList[0], dateList[i]);
        }
        break;
    case toLast:
        // probably not the most efficient
        for (i = 0; i < numDates-1; i++){
            vols[i] = CalcVol(dateList[i], dateList[numDates-1]);
        }
        break;
    default:
        throw ModelException(routine, "Unknown calculate type");
    }
    return;
}

/** Calculates vols beginning at dateFrom. If the dateList
    has n dates in it, n variances will be calculated. */
void VolProcessedMerton::CalcVol(const DateTime&      dateFrom,
                                 const DateTimeArray& datesTo,
                                 TCalcType            calcType, 
                                 CDoubleArray&        vols) const
{
    int numDates = datesTo.size();
    DateTimeArray dateSeries(datesTo.size()+1);

    dateSeries[0] = dateFrom;
    for(int i=0;i<numDates;++i)
    {
        dateSeries[i+1] = datesTo[i];
    }
    CalcVol(dateSeries, calcType, vols);
}

/** Calculates drifts between a series of dates */
void VolProcessedMerton::CalcDrift(const DateTimeArray& dateList, 
                                   TCalcType            calcType, 
                                   CDoubleArray&        drifts) const
{
    int numDates = dateList.size();
    static char routine[] = "VolProcessedMerton::CalcDrift";
    if (numDates < 2 && drifts.size() > 0 ) {
        throw ModelException(routine, "Must supply at least two dates");
    }
    if (drifts.size() < numDates-1){
        throw ModelException(routine, "Double array too short");
    }
    switch (calcType){
        int i; // MSVC is broken - doesn't support separate variables in loop
    case forward:
        for (i = 0; i < numDates-1; i++){
            drifts[i] = CalcDrift(dateList[i], dateList[i+1]);
        }
        break;
    case fromFirst:
        // probably not the most efficient
        for (i = 1; i < numDates; i++){
            drifts[i-1] = CalcVol(dateList[0], dateList[i]);
        }
        break;
    case toLast:
        // probably not the most efficient
        for (i = 0; i < numDates-1; i++){
            drifts[i] = CalcVol(dateList[i], dateList[numDates-1]);
        }
        break;
    default:
        throw ModelException(routine, "Unknown calculate type");
    }
    return;
}

/** Calculates drifts beginning at dateFrom. If the dateList
    has n dates in it, n drifts will be calculated. */
void VolProcessedMerton::CalcDrift(const DateTime&      dateFrom,
                                   const DateTimeArray& datesTo,
                                   TCalcType            calcType, 
                                   CDoubleArray&        drifts) const
{
    int numDates = datesTo.size();
    DateTimeArray dateSeries(datesTo.size()+1);

    dateSeries[0] = dateFrom;
    for(int i=0;i<numDates;++i)
    {
        dateSeries[i+1] = datesTo[i];
    }
    CalcDrift(dateSeries, calcType, drifts);
}

/** Preprocesses fields required for generateCumulatives() */
void VolProcessedMerton::setupCumulatives(
    const DateTime&      dateFrom,
    const DateTimeArray& datesTo,
    const double*        driverSpotVars)
{
    try {
        int iStep,jStep;
        double quantile;

        // counter bounds
        numSteps        = datesTo.size();
        Quantile(dateFrom,datesTo[numSteps-1],0.0001,100,&quantile,&maxMaxNumJumps);

        // allocate
        factor          = CDoubleMatrixSP(new CDoubleMatrix(numSteps,maxMaxNumJumps));
        proba           = CDoubleMatrixSP(new CDoubleMatrix(numSteps,maxMaxNumJumps));
        volCont         = CDoubleArraySP(new CDoubleArray(numSteps));
        maxNumJumps     = CIntArraySP(new CIntArray(numSteps));
        volJump         = CDoubleArraySP(new CDoubleArray(maxMaxNumJumps));
        meanJump        = CDoubleArraySP(new CDoubleArray(maxMaxNumJumps));
        driverSpotVol   = CDoubleArraySP(new CDoubleArray(numSteps));
        driverFwdVol    = CDoubleArraySP(new CDoubleArray(numSteps));

        // calculate maxNumJumps, factor and proba
        for( iStep=0; iStep<maxMaxNumJumps; iStep++ ) {
            (*volJump)[iStep] = volMerton->crashSizeUncertainty * sqrt((double)iStep);
            (*meanJump)[iStep] = (double)iStep * volMerton->crashSizeMean;
        }

        // local
        DoubleArray fwdVar(numSteps);
        DoubleArray spotVar(numSteps);

        // calculate maxNumJumps, factor and proba
        CalcVar(dateFrom,datesTo,forward,fwdVar);
        CalcVar(dateFrom,datesTo,fromFirst,spotVar);

        for( iStep=0; iStep<numSteps; iStep++ ) {
/*
       for( iStep=0; iStep<numSteps; iStep++ ) {
            (*volCont)[iStep]=sqrt(fwdVar[iStep]);
            Quantile(dateFrom,datesTo[iStep],0.0001,100,&quantile,&(*maxNumJumps)[iStep]);
            (*factor)[iStep][0] = 1.0 / sqrt( spotVar[iStep] );
            double dt = volMerton->timeMetric->yearFrac(dateFrom,datesTo[iStep]);
            double thisProba = exp(-dt*volMerton->crashRate) / quantile;
            (*proba)[iStep][0]  = thisProba;
            for( jStep=1; jStep<(*maxNumJumps)[iStep]; jStep++ ) {
                (*factor)[iStep][jStep] = 1.0 / sqrt( spotVar[iStep] 
                    + (double)jStep * volMerton->crashSizeUncertainty * volMerton->crashSizeUncertainty );
                thisProba *= volMerton->crashRate * dt / (double)jStep;
                (*proba)[iStep][jStep] = thisProba;
            }
            // driver vols
            (*driverSpotVol)[iStep] = sqrt(driverSpotVars[iStep]);
            (*driverFwdVol)[iStep] 
                = iStep>0 ? 
                    sqrt(driverSpotVars[iStep]-driverSpotVars[iStep-1])
                    : (*driverSpotVol)[iStep];
        }
*/
            bool temporalGauss = true;
            // double startVar = temporalGauss ? fwdVar[iStep] : spotVar[iStep];
            DateTime startDate
                            = (temporalGauss && (iStep>0)) ? datesTo[iStep-1] : dateFrom;
            double dt       = volMerton->timeMetric->yearFrac( startDate, datesTo[iStep]);
            double var      = temporalGauss ? fwdVar[iStep] : spotVar[iStep];

            (*volCont)[iStep]=sqrt(fwdVar[iStep]);
            Quantile(startDate,datesTo[iStep],0.0001,100,&quantile,&(*maxNumJumps)[iStep]);
            (*factor)[iStep][0] = 1.0 / sqrt(var);

            double thisProba = exp(-dt*volMerton->crashRate) / quantile;
            (*proba)[iStep][0]  = thisProba;
            for( jStep=1; jStep<(*maxNumJumps)[iStep]; jStep++ ) {
                (*factor)[iStep][jStep] = 1.0 / sqrt( var
                    + (double)jStep * volMerton->crashSizeUncertainty * volMerton->crashSizeUncertainty );
                thisProba *= volMerton->crashRate * dt / (double)jStep;
                (*proba)[iStep][jStep] = thisProba;
            }

            // driver vols
            (*driverSpotVol)[iStep] = sqrt(driverSpotVars[iStep]);
            (*driverFwdVol)[iStep] 
                = iStep>0 ? 
                    sqrt(driverSpotVars[iStep]-driverSpotVars[iStep-1])
                    : (*driverSpotVol)[iStep];
        }
    }
    catch(exception& e){
        throw ModelException(e, "VolProcessedMerton::SetupCumulatives");
    }
}

/** Samples from copula of VolProcessedMerton */
void VolProcessedMerton::generateCumulatives(
        const double*       random,
        const double*       randomJump,
        const double*       numJumps,
        double*             path) const
{
    int iStep;
    double thisPath = 0;
    double thisGaussianPath = 0;
    path[0] = 0;
    for (iStep = 0; iStep < numSteps; iStep++)
    {
        // Merton path
        double lastPath = thisPath;
        thisPath += (*volCont)[iStep]  * random[iStep];
        if( numJumps[iStep]>0 )
        {
            if( numJumps[iStep]<maxMaxNumJumps ) {
                thisPath += (*volJump)[(int)numJumps[iStep]] * randomJump[iStep] 
                               + (*meanJump)[(int)numJumps[iStep]];
            }
            else {
                thisPath += sqrt(numJumps[iStep]) * volMerton->crashSizeUncertainty * randomJump[iStep] 
                               + numJumps[iStep] * volMerton->crashSizeMean;
            }
        }

        // Map to uniform margin
        double* thisFactor = (*factor)[iStep];
        double* thisProba = (*proba)[iStep];
        if( true /*temporalGauss*/ ) {
            // F^{X_{t_k}-X_{t_{k-1}}}(X_{t_k}-X_{t_{k-1}})
            double thisGaussianIncrement = 0;
            for(int jStep = 0; jStep < maxMaxNumJumps/*(*maxNumJumps)[iStep]*/; jStep++ ) {
                thisGaussianIncrement += thisProba[jStep] 
                    * N1( (thisPath - lastPath - (*meanJump)[jStep] ) * thisFactor[jStep] );
            }
            thisGaussianPath += (*driverFwdVol)[iStep] * N1Inverse(thisGaussianIncrement);
            path[iStep] = N1( thisGaussianPath / (*driverSpotVol)[iStep] );
        }
        else {
            path[iStep] = 0;
            // F^{X_{t_k}}(X_{t_k})
            for(int jStep = 0; jStep < maxMaxNumJumps/*(*maxNumJumps)[iStep]*/; jStep++ ) {
                path[iStep] += thisProba[jStep] 
                    * N1( (thisPath - (*meanJump)[jStep] ) * thisFactor[jStep] );
            }
        }
    }
    // Manos' special!
    for (iStep = 0; iStep < numSteps; iStep++)
    {
        path[iStep] = 1-path[iStep];
    }
}

static void volProcessedMertonLoad(CClassSP& clazz){
    REGISTER(VolProcessedMerton, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IVolProcessed);
}

VolProcessedMerton::VolProcessedMerton(const CClassConstSP& clazz):
    CObject(clazz) {}

CClassConstSP const VolProcessedMerton::TYPE = CClass::registerClassLoadMethod(
    "VolProcessedMerton", typeid(VolProcessedMerton), volProcessedMertonLoad);


DRLIB_END_NAMESPACE
