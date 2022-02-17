//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolSVJ.cpp
//
//   Description : 
//
//   Date        : 22 April 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_VOLSVJ_CPP
#include "edginc/VolSVJ.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/VolProcessed.hpp"
#include "edginc/VolRequestRaw.hpp"
#include "edginc/VolRequestTime.hpp"
#include "edginc/FlatVol.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/VolSV.hpp"
#include "edginc/VolSVCJ.hpp"
#include "edginc/VolSVJJ.hpp"



DRLIB_BEGIN_NAMESPACE

DEFINE_TEMPLATE_TYPE(VolSVJArray);

VolSVJ::VolSVJ(const string& name,
               double initialVol,
               double meanVol,
               double meanReversRate,
               double volVol,
               double correlation,
               double crashRate,
               double crashSizeMean,
               double crashSizeUncertainty,
               double volRiskPrice):
VolBaseParam(TYPE, name), initialVol(initialVol), meanVol(meanVol), meanReversRate(meanReversRate),
volVol(volVol), correlation(correlation), crashRate(crashRate), crashSizeMean(crashSizeMean),
crashSizeUncertainty(crashSizeUncertainty), volRiskPrice(volRiskPrice) {
    validatePop2Object();
}


VolSVSP VolSVJ::convert(VolSV* p) const {
    VolSVSP volSV;
    if(Maths::isZero(volRiskPrice)) {
        if( Maths::isZero(crashRate) || 
            (Maths::isZero(crashSizeMean) && Maths::isZero(crashSizeUncertainty) )) {
                // Convert if no jumps
                volSV = VolSVSP(new VolSV(getName(), initialVol, meanVol, meanReversRate, volVol, correlation));
        }
    }
    
    return volSV;
}
    

VolSVJJSP VolSVJ::convert(VolSVJJ* p) const {
    VolSVJJSP volSVJJ;
    if(Maths::isZero(volRiskPrice)) {
        volSVJJ = VolSVJJSP(new VolSVJJ(
            getName(), initialVol, meanVol, meanReversRate, volVol, correlation, crashRate, crashSizeMean, crashSizeUncertainty,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
    }

    return volSVJJ;
}


VolSVCJSP VolSVJ::convert(VolSVCJ* p) const {
    VolSVCJSP volSVCJ;
    if(Maths::isZero(volRiskPrice)) {
        volSVCJ = VolSVCJSP(new VolSVCJ(
            getName(), initialVol, correlation, volVol, meanVol, meanReversRate,
            crashRate, crashSizeMean, crashSizeUncertainty, 0.0, 0.0, false));
    }
    return volSVCJ;
}


template<> string nameForType<VolSVJ_SpotDSType>(VolSVJ_SpotDSType*){
    return "VolSVJ::MCParams::SpotDiscreteSchemeType";
}
template<> string VolSVJ_SpotDSTypeHelper::names[VolSVJ_SpotDSTypeHelper::EnumList::NB_ENUMS] = {
    "EXACT",
    "EULER"
};

template<> string nameForType<VolSVJ_VarDSType>(VolSVJ_VarDSType*){
    return "VolSVJ::MCParams::VarDiscreteSchemeType";
}
template<> string VolSVJ_VarDSTypeHelper::names[VolSVJ_VarDSTypeHelper::EnumList::NB_ENUMS] = {
    "EULER",
    "VAR_TRANSFORM_EULER"
};


void VolSVJ::SVJVolParam::ComputeImpVol(const CVolBase*          vol,
                                        const CLatticeDouble&    strikes,
                                        const DateTimeArray&     maturities,
                                        CLatticeDouble&          impV) const{
    // turn the vol into what we must have
    const VolSVJ* myVol = static_cast<const VolSVJ *>(vol);
    // then just pass through the parameterised vol
    myVol->ComputeImpVol(strikes, maturities, impV);
}

/*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
VolSurface* VolSVJ::SVJVolParam::spotVolSurfaceFromStrikes(
    const CVolBase*       vol,
    const CDoubleArray&   strikes) const{
    // turn the vol into what we must have
    const VolSVJ* myVol = 
        static_cast<const VolSVJ *>(vol);
    // then just pass through the parameterised vol
    return myVol->spotVolSurfaceFromStrikes(strikes);
}

void VolSVJ::validatePop2Object(){
    try {
        if (!Maths::isZero(volRiskPrice)){
            throw ModelException("volRiskPrice is retired -- its value should read zero");
        }
        Calibrator::IAdjustable::checkRange(this);
        // build heston + mertonCrash
        update();
    }
    catch(exception& e){
        throw ModelException(e, "VolSVJ::validatePop2Object");
    }
}

void VolSVJ::ComputeImpVol(const CLatticeDouble&      strikes,
                           const DateTimeArray&       maturities,
                           CLatticeDouble&            impV) const {
    static const string routine("VolSVJ::ComputeImpVol");
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
VolSurface* VolSVJ::spotVolSurfaceFromStrikes(
    const CDoubleArray&   strikes) const{
    static const string routine("VolSVJ::spotVolSurfaceFromStrikes");
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
CVolParam* VolSVJ::createVolParam() const{
    return new SVJVolParam();
}

IObject* VolSVJ::defaultCtor(){
    return new VolSVJ();
}

/** Invoked when Class is 'loaded' */
void VolSVJ::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(VolSVJ, clazz);
    SUPERCLASS(VolBaseParam);
//    IMPLEMENTS(IVolatilityBS);
//    IMPLEMENTS(IVolatilityDVF);
    IMPLEMENTS(Calibrator::IAdjustable);
    IMPLEMENTS(VolAJDSuper::ISuperposable);
    IMPLEMENTS(IDynamicsParameter);
    IMPLEMENTS(MarketDataConvert::IConvert<VolSV>);
    IMPLEMENTS(MarketDataConvert::IConvert<VolSVJJ>);
    IMPLEMENTS(MarketDataConvert::IConvert<VolSVCJ>);
    EMPTY_SHELL_METHOD(defaultCtor);
    FIELD(initialVol, "initialVol");
    FIELD(meanVol, "meanVol");
    FIELD(meanReversRate, "meanReversRate");
    FIELD(volVol, "volVol");
    FIELD(correlation, "correlation");
    FIELD(crashRate, "crashRate");
    FIELD(crashSizeMean, "crashSizeMean");
    FIELD(crashSizeUncertainty, "crashSizeUncertainty");
    FIELD(volRiskPrice,
                 "coefficient of market price of risk for stochastic vol");
    FIELD_MAKE_OPTIONAL(volRiskPrice);

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
        clazz, "volVol", 
        new Range(Heston::RangeDef::volVol));
    Calibrator::IAdjustable::registerField(
        clazz, "correlation", 
        new Range(Heston::RangeDef::correlation));
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

VolSVJ::VolSVJ(): 
VolBaseParam(TYPE),
initialVol(Heston::DefaultVal::initialVol),
meanVol(Heston::DefaultVal::meanVol),
meanReversRate(Heston::DefaultVal::meanReversRate),
volVol(Heston::DefaultVal::volVol),
correlation(Heston::DefaultVal::correlation),
crashRate(MertonCrash::DefaultVal::crashRate),
crashSizeMean(MertonCrash::DefaultVal::crashSizeMean),
crashSizeUncertainty(MertonCrash::DefaultVal::crashSizeUncertainty),
volRiskPrice(0.0) {}

//specific method implemented in order to be used in MCPathConfigSVJ for volRequest (vol + trading time) for CCYP options
//Taken form SRM Should be careful that it cannot be called for other products
CVolProcessed* VolSVJ::getProcessedVol(const CVolRequest* volRequest,
                                       const CAsset*      /*asset*/) const{
    if (VolRequestRaw::TYPE->isInstance(volRequest) || 
        VolRequestTime::TYPE->isInstance(volRequest)){
        // it's ours or can just use this
        return const_cast<VolSVJ*>(this);
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
    throw ModelException("VolSVJ:getProcessedVol", 
                         "Request of type "+
                         volRequest->getClass()->getName()+
                         " not supported for " +this->getName());
}

/*** Build the parameterised vol and cache any values **/
void VolSVJ::buildCache() {
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
void VolSVJ::calcCumulantComponents(const DateTime&     fromDate,
                                    const DateTime&     toDate,
                                    const ComplexArray& inu,
                                    double              weight,
                                    Complex&            alpha,
                                    ComplexArray&       betas) const{
    static const string method("VolSVJ::calcCumulantComponents");
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
        // correct for drift
        alpha += inu[0] * (mertonCrash->calcAlpha(tau, 1.0) 
                            - mertonCrash->calcAlpha(tau, weight));
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

/** Returns the initial values Y_0,..., X^i_0,...
    Naturally, Y_0 = 0. */
double VolSVJ::getInitialValue(int i) const{
    if (i == 1){
        return (initialVol * initialVol);
    }
    return 0.0;
}


double VolSVJ::calcPricePutCliquet(double				strike,
								   double				coupon,
								   DateTime				valueDate,
								   DateTimeArray		maturities,
								   DateTimeArray		qMaturities,
								   YieldCurveWrapper	discount,
								   CAssetWrapper		asset) const{
	int nbMat = maturities.size();
	int nbqMat = qMaturities.size();
	
	DoubleArray eta(nbMat);
	DoubleArray q(nbMat);
    DoubleArray p(nbMat);
	DoubleArray r(nbMat);
	DoubleArray b(nbMat);
	DoubleArray fwd(nbMat);
	DoubleArray sigmaSq(nbMat);

	// discount factor for coupon maturities
	DoubleArray qb(nbqMat);
	// proba that payoff is not made before coupon maturity i
	DoubleArray qp(nbqMat);

	asset->fwdValue(maturities,fwd);

	double logCrashSizeMean = log(1+crashSizeMean)-0.5*Maths::square(crashSizeUncertainty);
   
    int i, j, l;
	for ( i=1; i<nbMat; i++){
		sigmaSq[i] = initialVol+(meanVol-initialVol)*exp(-meanReversRate*maturities[0].yearFrac(maturities[i]));
		sigmaSq[i] = Maths::square(sigmaSq[i])*maturities[i-1].yearFrac(maturities[i]);
	}

    for ( i=1; i<nbMat; i++){
		eta[i] = fwd[i]/fwd[i-1]-1-0.5*sigmaSq[i]-crashRate*maturities[i-1].yearFrac(maturities[i])*crashSizeMean;
    }

	for ( i=1; i<nbMat; i++){
		q[i]=0;
	    int f = 1;
		for ( l=0; l<10; l++){

			if (l > 0)
				f *=l;

			double sq = sqrt(l*Maths::square(crashSizeUncertainty)+sigmaSq[i]);
			if (sq==0){
				q[i] += 0;
			} else {
				q[i] += N1((log(strike)-eta[i]-l*logCrashSizeMean)/sqrt(l*Maths::square(crashSizeUncertainty)+sigmaSq[i]))
					*pow(crashRate*maturities[i-1].yearFrac(maturities[i]),l)/f*exp(-crashRate*maturities[i-1].yearFrac(maturities[i]));
			}
		}
	}
	p[1]=q[1];	
	for ( i=2; i<nbMat; i++){
		p[i]=p[i-1]*q[i]*(1-q[i-1])/q[i-1];
	}

	for ( i=1; i<nbMat; i++){
        r[i] = 0;
		int f = 1;
		for (l=0; l<10; l++){
			if (l > 0)
				f *=l;
			double sq = sqrt(l*Maths::square(crashSizeUncertainty)+sigmaSq[i]);
			if (sq==0){
                r[i] += 0;
			} else {
				r[i] += exp(eta[i]+l*logCrashSizeMean+0.5*(l*Maths::square(crashSizeUncertainty)+sigmaSq[i]))
					*N1((log(strike)-eta[i]-l*logCrashSizeMean)/
						sqrt(l*Maths::square(crashSizeUncertainty)+sigmaSq[i])-sqrt(l*Maths::square(crashSizeUncertainty)+sigmaSq[i]))
					*pow(crashRate*maturities[i-1].yearFrac(maturities[i]),l)/f*exp(-crashRate*maturities[i-1].yearFrac(maturities[i]));
			}
		}
	}

	for ( i=1; i<nbMat; i++){
		b[i] = discount->pv(valueDate,maturities[i]);
	}

	double varLeg =  0;
	for ( i=1; i<nbMat; i++){
		varLeg += strike*b[i]*p[i]-b[i]*p[i]/q[i]*r[i];
	}
	
	j = 1;
	for ( i=1; i<nbqMat; i++){
		if (i == 1){
			qp[1] = 1;
		} else{
			qp[i] = qp[i-1];
		}
		while ((maturities[j]<qMaturities[i])&&(j<nbMat)){
			qp[i] *= 1-q[j];
			j ++;
		}
	}
	for ( i=1; i<nbqMat; i++){
		qb[i] = discount->pv(valueDate,qMaturities[i]);
	}

	double fixLeg =  0;
	for ( i=1; i<nbqMat; i++){
		fixLeg += qb[i]*coupon*qMaturities[i-1].yearFrac(qMaturities[i])*qp[i];	
	}

	int iqLast;
	for ( i=1; i<nbMat; i++){
		iqLast = 0;
		while ((qMaturities[iqLast]<=maturities[i])&&(iqLast<nbqMat)){
				iqLast++;
		}
		iqLast--; 
		fixLeg += b[i]*coupon*qMaturities[iqLast].yearFrac(maturities[i])*p[i];
	}

	return varLeg-fixLeg;
}
/** Returns the nber of factors. The first factor is understood to be 
    the log spot. */
int VolSVJ::getNbFactors() const{
    return 2;
}

void VolSVJ::getMarket(const IModel*     model, 
                       const MarketData* market, 
                       const string&     name){
    VolBaseParam::getMarket(model, market, name);
}

//////////// class FourierProcessSVJ ///////////////////

/* Started Log Return */
Complex VolSVJ::scalelessCumulant(const StFourierProcessLogRtn& process,
                                  const StFourierProductLogRtn& product, 
                                  const Complex& z, 
                                  const DateTime& matDate) const{
    static const string method = "VolSVJ::scalelessCumulant";
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
Complex VolSVJ::scalelessCumulant(const FwdStFourierProcessLogRtn& process,
                                  const FwdStFourierProductLogRtn& product, 
                                  const Complex& z, 
                                  const DateTime& matDate) const{
    static const string method = "VolSVJ::scalelessCumulant";
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

/*ARNAUD*/
///Expected Quad Var///

Complex VolSVJ::cumulant(const FwdStFourierProcessExpQuadVar& process,
                        const FwdStFourierProductExpQuadVar& product, 
                        const Complex& z, 
                        const DateTime& matDate) const{
    static const string method = "VolSVJ::cumulant";
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
		double gammabar = log(1 + crashSizeMean) - 0.5 * Maths::square(crashSizeUncertainty);
		alpha_tT = z * meanVar / meanReversRate * (meanReversRate * tau - 1.0 + phi) + z * crashRate * tau * (Maths::square(gammabar) + Maths::square(crashSizeUncertainty)); 
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




////////////////////////////////////////////////////////////////////////////

/* Started integrated variance */
Complex VolSVJ::cumulant(const StFourierProcessIntVar& process,
                         const StFourierProductIntVar& product, 
                         const Complex& z, 
                         const DateTime& matDate) const {    
    static const string method = "VolSVJ::cumulant";
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
Complex VolSVJ::cumulant(const FwdStFourierProcessIntVar& process,
                         const FwdStFourierProductIntVar& product, 
                         const Complex& z, 
                         const DateTime& matDate) const {        
    static const string method = "VolSVJ::cumulant";
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
Complex VolSVJ::cumulant(const StFourierProcessQuadVar& process,
                         const StFourierProductQuadVar& product, 
                         const Complex& z, 
                         const DateTime& matDate) const {    
    static const string method = "VolSVJ::cumulant";
    try{
        double tau = timeMetric->yearFrac(baseDate, matDate);

        Complex alpha, beta;
        calcJointLapAlphaBeta(tau,
                              0.0,
                              0.0,
                              z,
                              alpha,
                              beta);

        // Merton
        alpha += mertonCrash->calcQuadVarCumulant(tau, z);

        return (alpha + beta * Maths::square(initialVol));
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

/* Forward Starting integrated variance */
Complex VolSVJ::cumulant(const FwdStFourierProcessQuadVar& process,
                         const FwdStFourierProductQuadVar& product, 
                         const Complex& z, 
                         const DateTime& matDate) const {        
    static const string method = "VolSVJ::cumulant";
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

        // Merton
        Complex alpha;
        alpha = mertonCrash->calcQuadVarCumulant(tau, z);

        return (alpha_t + beta_t * Maths::square(initialVol) + alpha_tT + alpha);

    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

double VolSVJ::expectation(const StFourierProcessQuadVar& process,
                           const StFourierProductQuadVar& product, 
                           const DateTime& matDate) const {    
    static const string method = "VolSVJ::cumulant";
    try{
        double tau = timeMetric->yearFrac(baseDate, matDate);
        double jumpTerm = mertonCrash->calcAnnualQuadVar(tau);
        double volTerm = heston->calcAnnualQuadVar(tau);
        return (jumpTerm + volTerm);
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

double VolSVJ::expectation(const FwdStFourierProcessQuadVar& process,
                           const FwdStFourierProductQuadVar& product, 
                           const DateTime& matDate) const {    
    static const string method = "VolSVJ::expectation";
    try{
        /* Time fractions from today till start of option and then from
           start of option till maturity */
        double t = timeMetric->yearFrac(baseDate,
                                        product.getStartDate());
        double tau = timeMetric->yearFrac(product.getStartDate(),
                                          matDate);
        double jumpTerm = mertonCrash->calcAnnualQuadVar(tau);
        double volTerm = heston->calcAnnualQuadVar(t,tau);
        return (jumpTerm + volTerm);
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

/** Calculates the components alpha and beta that appear in 
    the exponent of the time-t joint Bilateral Laplace 
    of X_T = (Y_T, V_T, I_T) where Y_T = ln(S_T / F(0, T)) is the dimension-less
    log spot at time T, V_T is the instantaneous variance at time T
    and I_T is the integrated variance from 0 to time T
    The Bilateral Laplace transform is of the form
        exp(alpha(tau, u1, u2, u3) 
            + u1 * Y_t 
            + beta(tau, u1, u2, u3) * V_t 
            + u3 * I_t)
    where tau = T - t and the complex numbers u1, u2, u3 are the frequencies 
    wrt Y_T, V_T and I_T, respectively. */
void VolSVJ::calcJointLapAlphaBeta(double         tau,    // tau = T - t (in years)
                                   const Complex& u1,
                                   const Complex& u2,
                                   const Complex& u3,
                                   Complex&       alpha,
                                   Complex&       beta) const{
    static const string method = "VolSVJ::calcJointLapAlphaBeta";

    try{    // in case Complex operations fail for some reason
        /* Calculate Heston contribution */
        heston->calcAlphaBeta(tau,
                              u1,
                              u2,
                              u3,
                              alpha,
                              beta);

        /* Calculate Merton jump contribution to alpha */
        alpha += mertonCrash->calcAlpha(tau, u1);
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

CClassConstSP const VolSVJ::TYPE =
CClass::registerClassLoadMethod("VolSVJ", typeid(VolSVJ), load);

CClassConstSP const VolSVJ::SVJVolParam::TYPE =
CClass::registerClassLoadMethod("VolSVJ::SVJVolParam",
                                typeid(SVJVolParam), load);

/** Needed for IAdjustable interface. Returns market data name for vol */
string VolSVJ::getName() const{
    return CVolBaseParamSurface::getName();
}

/** Called after adjustments have been made to fields (eg calibrator) */
void VolSVJ::fieldsUpdated(const CFieldArray& fields){
    update();
}

/** Called after adjustments have been made to fields */
void VolSVJ::update(){
    // refresh heston
    heston = HestonSP(new Heston(initialVol,
                                 meanVol,
                                 meanReversRate,
                                 volVol,
                                 correlation));
    // refresh mertonCrash
    mertonCrash = MertonCrashSP(new MertonCrash(crashRate,
                                                crashSizeMean,
                                                crashSizeUncertainty));
}

VolSVJ::IPathGen* VolSVJ::createPathGen(const IPathGenMaker& maker) const{
    return maker.make(initialVol,
                      meanVol,
                      meanReversRate,
                      volVol,
                      correlation,
                      crashRate,
                      crashSizeMean,
                      crashSizeUncertainty,
                      volRiskPrice);
}

/** Given an array of gaussian deviates, populate the arrays 'instVars' 
    and 'integratedVars' with simulated values using an Euler
    discretization scheme */
void VolSVJ::MCGenerateVarPaths(Int2Type<MCParams::VarDSType::EULER>, 
                               const DoubleArray&   tradYears,
                               const double*        deviates,
                               DoubleArray&         instVars,
                               DoubleArray&         integratedVars) const{
#ifdef DEBUG
    ASSERT(instVars.size() == tradYears.size());
    ASSERT(instVars.size() == integratedVars.size());
#endif
    int nbSteps = tradYears.size();
    double meanVar = Maths::square(meanVol);
    instVars[0] = Maths::square(initialVol);
    integratedVars[0] = 0.0;
    for (int iStep = 1, iLastStep = 0; iStep < nbSteps; iLastStep = iStep++){
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
void VolSVJ::MCGenerateVarPaths(Int2Type<MCParams::VarDSType::VAR_TRANSFORM_EULER>,
                               const DoubleArray&   tradYears,
                               const double*        deviates,
                               DoubleArray&         instVars,
                               DoubleArray&         integratedVars) const{
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
void VolSVJ::MCGenerateVarPaths(Int2Type<MCParams::VarDSType::VAR_TRANSFORM_EULER>,
                               const DoubleArray& tradYears,
                               const double*        deviates,
                               DoubleArray&         instVars,
                               DoubleArray&         integratedVars) const{
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
const CClassConstSP VolSVJConvert::TYPE = 
CClass::registerInterfaceLoadMethod("MarketDataConvert::IConvert<VolSVJ>", typeid(VolSVJConvert), load);

/* external symbol to allow class to be forced to be linked in */
bool VolSVJLinkIn(){
    return (VolSVJ::TYPE != 0) && ( VolSVJConvert::TYPE != 0);
}

DRLIB_END_NAMESPACE
