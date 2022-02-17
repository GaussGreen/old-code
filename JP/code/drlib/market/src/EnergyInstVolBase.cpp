//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyInstVolBase.cpp
//
//   Description : Energy 2 factor vol surface. Based on drcommodityvolsurfacei.cpp
//                    and drcommodityInstaneousvolsurface.cpp in FXLIB.
//                 
//
//   Author      : Sean Chen
//
//   Date        : June 01, 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Addin.hpp"
#include "edginc/EnergyInstVolBase.hpp"

DRLIB_BEGIN_NAMESPACE

EnergyInstVolBase::EnergyInstVolBase(const CClassConstSP& clazz):EnergyVolBase(clazz)
{
}

EnergyInstVolBase::~EnergyInstVolBase()
{
}

EnergyInstVolBase::EnergyInstVolBase():EnergyVolBase(TYPE)
{
}

void EnergyInstVolBase::validatePop2Object()
{
    static const string method = "EnergyInstVolBase::validatePop2Object()";

    // do some checking...
    
    if (alpha<=0.0)
        ModelException(method,"Alpha is not strictly positive");
    fAlphas.push_back(alpha);

    fMaxPricingDate = baseDate;
}


const DoubleArray& EnergyInstVolBase::getAlphas() const
{
    return fAlphas;
}

const DateTime& EnergyInstVolBase::getMaxPricingDate() const
{
    return fMaxPricingDate;
}

const DateTimeArray& EnergyInstVolBase::getFutureMaturityDates() const
{
    return fFutureMaturityDates;
}

const DateTimeArray& EnergyInstVolBase::getOptionExpiryDates() const
{
    return fOptionExpiryDates;
}
        
const DoubleMatrix& EnergyInstVolBase::getSigmas() const
{
    return fSigmas;
}

const DoubleArray& EnergyInstVolBase::getRatios() const
{
    return fRatios;
}

// For Capital T model only
void EnergyInstVolBase::calculateLongTermMarketVol(
    double& /* marketVol */, 
    int&     /* contractID */, 
    int&     /*numCalibratableBenchmark */) const
{
    static const string method = "EnergyInstVolBase::calculateLongTermMarketVol()";
    ModelException(method,"calculateLongTermMarketVol not implemented for volModel");
}

// For Capital T model only
void EnergyInstVolBase::calculateLongTermModelVol(
    double& /* marketVol  */,
    int&     /* contractID */, 
    int&     /*numCalibratableBenchmark */) const
{
    static const string method = "EnergyInstVolBase::calculateLongTerModelVol()";
    ModelException(method,"calculateLongTermModelVol not implemented for volModel");
}

double EnergyInstVolBase::getATMVol(
                                 const DateTime& expiry, 
                                 const DateTime& futuresExpiry) const
{
    static const string method = "EnergyInstVolBase::getATMVol()";

    if (expiry>futuresExpiry)
        ModelException(method,"he option expiry must be less or equal to the futures expiry");

    DateTime today = getBaseDate();
    if (expiry <= today)
        return 0.0;
    else if (  expiry > getMaxPricingDate() )
        ModelException(method,"Expiry date greater than maximum pricing date");
    else
    {
        double alpha, beta;
        double sigma1, sigma2, sigma1Bar, sigma2Bar;
        double futuresExpiryTimePrev, futuresExpiryTimeNow, expiryTime;
        double xPrev, yPrev, zPrev, xNow, yNow, zNow;
        double sum1_1, sum2_2, sum1Bar_1Bar, sum2Bar_2Bar, sum1_1Bar, sum2_1Bar, sum1_2;
        const DateTimeArray& futuresExpiries = getFutureMaturityDates();
        
        const DoubleMatrix& locSigmas = getSigmas();
        const DoubleArray& Alphas = getAlphas();
        alpha = Alphas[0];
        beta = Alphas[1];

        futuresExpiryTimePrev = 0.0;
        xPrev = 1.0;
        yPrev = 1.0;
        zPrev = 1.0;

        int i=0;
        double x_T = 0.0;
        sum1_1 = 0.0;
        sum2_2 = 0.0;
        sum1Bar_1Bar = 0.0;
        sum2Bar_2Bar = 0.0;
        sum1_1Bar = 0.0;
        sum2_1Bar = 0.0;
        sum1_2 = 0.0;

        // Compute X(T) from sigmas
        int maturityIndex=-1;    
        int numFuturesExpiries = futuresExpiries.size();    

        // Determine which contract we are dealing with
        for(i=0; i<numFuturesExpiries; i++) 
        {
            if(futuresExpiries[i] == futuresExpiry)
            {
                maturityIndex = i;
            }
        }
        if(maturityIndex == -1 || maturityIndex >= locSigmas.numRows())
        {
            ModelException(method,"EnergyInstVolSeasonal::ATMVol: FuturesExpiry must be one of the sigma dates!");
        }
        sigma1 = locSigmas[maturityIndex][0];
        sigma1Bar = locSigmas[maturityIndex][1];

        if(Maths::isZero(sigma1Bar))
        {
            ModelException(method,"EnergyInstVolSeasonal::ATMVol: Sigma1Bar is zero!");
        }
        x_T = sigma1/sigma1Bar; // X(T) can be computed from sigma1 and sigma1Bar.
        
        i=0;
        do
        {
            sigma1 = locSigmas[i][0];
            sigma1Bar = locSigmas[i][1];
            sigma2Bar = locSigmas[i][2];
            sigma2 = locSigmas[i][3];

            if(futuresExpiries[i]<=today)
            {
                continue;
            }

            futuresExpiryTimeNow = ((expiry < futuresExpiries[i]) ? expiry : futuresExpiries[i]).daysDiff(today) / 365.0;
            xNow = exp(alpha*futuresExpiryTimeNow);
            yNow = exp(beta*futuresExpiryTimeNow);
            zNow = exp((alpha+beta)*futuresExpiryTimeNow);

            sum1_1 += sigma1Bar*sigma1Bar * (xNow*xNow - xPrev*xPrev);
            sum2_2 += sigma2*sigma2 * (yNow*yNow - yPrev*yPrev);
            sum1Bar_1Bar += sigma1Bar*sigma1Bar * (futuresExpiryTimeNow - futuresExpiryTimePrev);
            sum2Bar_2Bar += sigma2Bar*sigma2Bar * (futuresExpiryTimeNow - futuresExpiryTimePrev);
            sum1_1Bar += sigma1Bar*sigma1Bar * (xNow - xPrev);
            sum2_1Bar += sigma2*sigma1Bar * (yNow - yPrev);
            sum1_2 += sigma1Bar*sigma2 * (zNow - zPrev);

            futuresExpiryTimePrev = futuresExpiryTimeNow;
            xPrev = xNow;
            yPrev = yNow;
            zPrev = zNow;
        }
        while (futuresExpiries[i++] <= expiry);

        expiryTime = expiry.daysDiff(today) / 365.0;
        futuresExpiryTimeNow = futuresExpiry.daysDiff(today) / 365.0;
        xNow = exp(alpha*futuresExpiryTimeNow);
        yNow = exp(beta*futuresExpiryTimeNow);
        zNow = exp((alpha+beta)*futuresExpiryTimeNow);

        return sqrt( (sum1_1*x_T*x_T/(2.0*alpha)/(xNow*xNow) \
                      + sum2_2/(2.0*beta)/(yNow*yNow) \
                      + 2.0*sum1_1Bar*x_T/alpha/xNow \
                      + 2.0*sum2_1Bar/beta/yNow \
                      + 2.0*sum1_2*x_T/(alpha+beta)/zNow \
                      + sum1Bar_1Bar + sum2Bar_2Bar ) \
                      / expiryTime );
    }
	return 0.0;
}

double EnergyInstVolBase::getSmileVolByStrike(const DateTime& expiry, double strike) const
{
    static const string method = "EnergyInstVolBase::getSmileVolByStrike()";
    ModelException(method, "This is not implemented yet.");
    return 0;
}

double EnergyInstVolBase::getSmileVolByDelta(const DateTime& expiry, double delta) const
{
    static const string method = "EnergyInstVolBase::getSmileVolByDelta()";
    ModelException(method, "This is not implemented yet.");
    return 0;
}

void EnergyInstVolBase::analyticApproximation(
    double* fwdPositive_,
    double* fwdNegative_,
    double* volPositive_,
    double* volNegative_,
    double* correlation_,
    const DateTime& expiry,
    const DateTimeArray& fixingDates,
    const DoubleArray& weights,
    const DateTimeArray& maturities,
    const DoubleArray& prices) const
{
    static const string method = "EnergyInstVolBase::analyticApproximation()";

    int i, j;

    double fwdPositive = 0.0, fwdNegative = 0.0, volPositive = 0.0, volNegative = 0.0, correlation = 0.0;
    const DateTime& today = this->getBaseDate();
    const DateTime& maxFixingDate = fixingDates.back(); // here assume fixing dates are ordered (will throw an error if vector is empty)

    if (maxFixingDate>this->getMaxPricingDate())
    {
        ModelException(method,"Last fixing date is greater than the maximum permissible pricing date");
    }
    else if (today<expiry)
    {
        // Extract vol surface information and initialise vol surface variables
        DoubleArray alphas(this->getAlphas());
        double alpha = alphas[0], beta = alphas[1];
        double boundaryTime;
        const DateTimeArray& boundaryDates = this->getFutureMaturityDates();
        DoubleArray boundaryTimes, expAlpha_S, expBeta_S;
        boundaryTimes.push_back(0.0);
        expAlpha_S.push_back(1.0);
        expBeta_S.push_back(1.0);
        for (i=0; boundaryDates[i]<maxFixingDate; ++i)
        {
            boundaryTime = boundaryDates[i].daysDiff(today)/365.0;
            boundaryTimes.push_back(boundaryTime);
            expAlpha_S.push_back(exp(alpha*boundaryTime));
            expBeta_S.push_back(exp(beta*boundaryTime));
        }

        // Extract fixing information and initialise fixing variables
        DateTime fixingDatePrevious = today;
        DateTime maturityDate, fixingDate;
        double x_T, sigma1, sigma1Bar, maturityTime, fixingTime, fixingPrice, fixingVol, weight;
        DoubleArray posX_T, negX_T;
        DoubleArray posFixingTimes, negFixingTimes, posFixingPrices, negFixingPrices, posFixingVols, negFixingVols;
        DoubleArray posExpAlpha_T, negExpAlpha_T, posExpBeta_T, negExpBeta_T, posExpAlpha_t, negExpAlpha_t, posExpBeta_t, negExpBeta_t;
        int intervalEndpoint = 0;
        int maturityIndex;
        int numFuturesExpiries = boundaryDates.size();
        IntArray posIntervalEndpoints, negIntervalEndpoints;
        const DoubleMatrix& locSigmas = this->getSigmas();

        for (i=0; i<fixingDates.size(); ++i)
        {
            fixingDate = fixingDates[i];
            weight = weights[i];
            maturityDate = maturities[i];

            if (fixingDate<fixingDatePrevious)
                ModelException(method,"Fixing dates are non-increasing or historic");
            if (maturityDate<fixingDate)
                ModelException(method,"Fixing date is after contract maturity date");
            
            maturityTime = maturityDate.daysDiff(today)/365.0;
            fixingTime = fixingDate.daysDiff(today)/365.0;
            fixingPrice = weight * prices[i];
            fixingVol = this->getATMVol(fixingDate, maturityDate);
            fixingDatePrevious = fixingDate;

            while (boundaryDates[intervalEndpoint]<fixingDate)
                ++intervalEndpoint; // interval endpoints are by definition increasing

            // Compute X(T) from sigmas
            maturityIndex=-1;    
            for(j=0; j<numFuturesExpiries; j++) // Determine which contract we are dealing with
            {
                if(boundaryDates[j] == maturityDate)
                {
                    maturityIndex = j;
                }
            }
            if(maturityIndex == -1 || maturityIndex >= locSigmas.numRows())
            {
                ModelException(method," FuturesExpiry must be one of the sigma dates!");
            }
            sigma1 = locSigmas[maturityIndex][0];
            sigma1Bar = locSigmas[maturityIndex][1];

            if(Maths::isZero(sigma1Bar))
            {
                ModelException(method," Sigma1Bar is zero!");
            }
            x_T = sigma1/sigma1Bar; // X(T) can be computed from sigma1 and sigma1Bar.

            if (weight>0.0)
            {
                posX_T.push_back(x_T);
                posFixingTimes.push_back(fixingTime);
                posFixingPrices.push_back(fixingPrice);
                posFixingVols.push_back(fixingVol);

                posExpAlpha_T.push_back(exp(alpha*maturityTime));
                posExpBeta_T.push_back(exp(beta*maturityTime));
                posExpAlpha_t.push_back(exp(alpha*fixingTime));
                posExpBeta_t.push_back(exp(beta*fixingTime));

                posIntervalEndpoints.push_back(intervalEndpoint+1);
            }
            else if (weight<0.0)
            {
                negX_T.push_back(x_T);
                negFixingTimes.push_back(fixingTime);
                negFixingPrices.push_back(fixingPrice);
                negFixingVols.push_back(fixingVol);

                negExpAlpha_T.push_back(exp(alpha*maturityTime));
                negExpBeta_T.push_back(exp(beta*maturityTime));
                negExpAlpha_t.push_back(exp(alpha*fixingTime));
                negExpBeta_t.push_back(exp(beta*fixingTime));

                negIntervalEndpoints.push_back(intervalEndpoint+1);
            }
        }

        // Calculate fwds and vols
        double expiryTime = expiry.daysDiff(today)/365.0;
        double variance, corr;
        this->AnalyticVariance(
            &fwdPositive, 
            &variance, 
            posX_T,
            posFixingTimes, 
            posFixingPrices, 
            posFixingVols, 
            posExpAlpha_T, 
            posExpBeta_T, 
            posExpAlpha_t, 
            posExpBeta_t, 
            posIntervalEndpoints,
            boundaryTimes,
            expAlpha_S,
            expBeta_S);
        volPositive = (variance==0.0) ? 0.0 : sqrt(log(variance)/expiryTime);
        this->AnalyticVariance(
            &fwdNegative, 
            &variance, 
            negX_T,
            negFixingTimes, 
            negFixingPrices, 
            negFixingVols, 
            negExpAlpha_T, 
            negExpBeta_T, 
            negExpAlpha_t, 
            negExpBeta_t, 
            negIntervalEndpoints,
            boundaryTimes,
            expAlpha_S,
            expBeta_S);
        volNegative = (variance==0.0) ? 0.0 : sqrt(log(variance)/expiryTime);
        this->AnalyticCorrelation(
            &corr, 
            fwdPositive,
            volPositive,
            posX_T,
            posFixingTimes,
            posFixingPrices,
            posExpAlpha_T, 
            posExpBeta_T, 
            posExpAlpha_t, 
            posExpBeta_t, 
            posIntervalEndpoints,
            fwdNegative,
            volNegative,
            negX_T,
            negFixingTimes,
            negFixingPrices,
            negExpAlpha_T, 
            negExpBeta_T, 
            negExpAlpha_t, 
            negExpBeta_t, 
            negIntervalEndpoints,
            boundaryTimes,
            expAlpha_S,
            expBeta_S);
        correlation = corr/expiryTime;
    }

    if (fwdPositive_)
        *fwdPositive_ = fwdPositive;
    if (fwdNegative_)
        *fwdNegative_ = fwdNegative;
    if (volPositive_)
        *volPositive_ = volPositive;
    if (volNegative_)
        *volNegative_ = volNegative;
    if (correlation_)
        *correlation_ = correlation;
}

void EnergyInstVolBase::AnalyticVariance(
    double* fwd,
    double* var,
    const DoubleArray& x_T,
    const DoubleArray& fixingTimes,
    const DoubleArray& fixingPrices,
    const DoubleArray& fixingVols,
    const DoubleArray& expAlpha_T,
    const DoubleArray& expBeta_T,
    const DoubleArray& expAlpha_t,
    const DoubleArray& expBeta_t,
    const IntArray& intervalEndpoints,
    const DoubleArray& boundaryTimes,
    const DoubleArray& expAlpha_S,
    const DoubleArray& expBeta_S) const
{
    int i, j;

    if (fixingTimes.empty())
    {
        *fwd = 0.0;
        *var = 0.0;
    }
    else
    {
        double expectation, xT, expAlphaT, expBetaT;
        DoubleArray expectations;
        int numFixings = fixingTimes.size();
        double fixingTime = fixingTimes[numFixings-1];
        double fixingPrice = fixingPrices[numFixings-1];
        double fixingVol = fixingVols[numFixings-1];
        double compositeFuture = fixingPrice;
        double diagonalVar = fixingPrice*fixingPrice*exp(fixingVol*fixingVol*fixingTime), nonDiagonalVar = 0.0;
        for (i=0; i<numFixings-1; ++i)
        {
            fixingTime = fixingTimes[i];
            fixingPrice = fixingPrices[i];
            fixingVol = fixingVols[i];

            xT =x_T[i];
            expAlphaT = expAlpha_T[i];
            expBetaT = expBeta_T[i];
            this->AnalyticExpectations(&expectations, i, fixingTime, expAlpha_t[i], expBeta_t[i], intervalEndpoints, boundaryTimes, expAlpha_S, expBeta_S);
            compositeFuture += fixingPrice;

            for (j=i+1; j<numFixings; ++j)
            {
                expectation = exp(
                    expectations[0]*xT*x_T[j]/(expAlphaT*expAlpha_T[j]) +
                    expectations[1]/(expBetaT*expBeta_T[j]) +
                    expectations[2] +
                    expectations[3]*xT/(expAlphaT*expBeta_T[j]) + expectations[3]*x_T[j]/(expBetaT*expAlpha_T[j]) +
                    expectations[4]*xT/expAlphaT + expectations[4]*x_T[j]/expAlpha_T[j] +
                    expectations[5]/expBetaT + expectations[5]/expBeta_T[j]);

                nonDiagonalVar += fixingPrice*fixingPrices[j]*expectation;
            }

            diagonalVar += fixingPrice*fixingPrice*exp(fixingVol*fixingVol*fixingTime);
        }

        *fwd = compositeFuture;
        *var = (diagonalVar+2*nonDiagonalVar)/(compositeFuture*compositeFuture);
    }
}
        
void EnergyInstVolBase::AnalyticCorrelation(
    double* corr,
    double fwd1,
    double vol1,
    const DoubleArray& x_T1,
    const DoubleArray& fixingTimes1,
    const DoubleArray& fixingPrices1,
    const DoubleArray& expAlpha_T1,
    const DoubleArray& expBeta_T1,
    const DoubleArray& expAlpha_t1,
    const DoubleArray& expBeta_t1,
    const IntArray& intervalEndpoints1,
    double fwd2,
    double vol2,
    const DoubleArray& x_T2,
    const DoubleArray& fixingTimes2,
    const DoubleArray& fixingPrices2,
    const DoubleArray& expAlpha_T2,
    const DoubleArray& expBeta_T2,
    const DoubleArray& expAlpha_t2,
    const DoubleArray& expBeta_t2,
    const IntArray& intervalEndpoints2,
    const DoubleArray& boundaryTimes,
    const DoubleArray& expAlpha_S,
    const DoubleArray& expBeta_S) const
{
    int i, j;

    if (fixingTimes1.empty() || fixingTimes2.empty())
    {
        *corr = 0.0;
    }
    else
    {
        double expectation;
        DoubleArray expectations;
        int numFixings1 = fixingTimes1.size(), numFixings2 = fixingTimes2.size();
        double invertedVar = 0.0;
        for (i=0; i<numFixings1; ++i)
        {
            for (j=0; j<numFixings2; ++j)
            {
                if (fixingTimes1[i]<fixingTimes2[j])
                    this->AnalyticExpectations(&expectations, i, fixingTimes1[i], expAlpha_t1[i], expBeta_t1[i], intervalEndpoints1, boundaryTimes, expAlpha_S, expBeta_S);
                else
                    this->AnalyticExpectations(&expectations, j, fixingTimes2[j], expAlpha_t2[j], expBeta_t2[j], intervalEndpoints2, boundaryTimes, expAlpha_S, expBeta_S);
                
                expectation = exp(
                    expectations[0]*x_T1[i]*x_T2[j]/(expAlpha_T1[i]*expAlpha_T2[j]) +
                    expectations[1]/(expBeta_T1[i]*expBeta_T2[j]) +
                    expectations[2] +
                    expectations[3]*x_T1[i]/(expAlpha_T1[i]*expBeta_T2[j]) + expectations[3]*x_T2[j]/(expBeta_T1[i]*expAlpha_T2[j]) +
                    expectations[4]*x_T1[i]/expAlpha_T1[i] + expectations[4]*x_T2[j]/expAlpha_T2[j] +
                    expectations[5]/expBeta_T1[i] + expectations[5]/expBeta_T2[j]);

                invertedVar += fixingPrices1[i]*fixingPrices2[j]*expectation;
            }
        }

        *corr = log(invertedVar/(fwd1*fwd2))/(vol1*vol2);
    }
}

void EnergyInstVolBase::AnalyticExpectations(
    DoubleArray* expectations,
    int fixingIndex,
    double t,
    double expAlpha_t,
    double expBeta_t,
    const IntArray& intervalEndpoints,
    const DoubleArray& boundaryTimes,
    const DoubleArray& expAlpha_S,
    const DoubleArray& expBeta_S) const
{
    int intervalEndpoint = intervalEndpoints[fixingIndex];
    double sigma1Bar, sigma2, sigma2Bar;
    double alphaAlpha = 0.0, betaBeta = 0.0, gammaGamma = 0.0, alphaBeta = 0.0, alphaGamma = 0.0, betaGamma = 0.0;
    
    const DoubleMatrix& locSigmas = this->getSigmas();
    
    for (int i=1; i<intervalEndpoint; ++i)
    {
        sigma1Bar = locSigmas[i-1][1];
        sigma2Bar = locSigmas[i-1][2];
        sigma2 = locSigmas[i-1][3];

        alphaAlpha += sigma1Bar*sigma1Bar*(expAlpha_S[i]*expAlpha_S[i] - expAlpha_S[i-1]*expAlpha_S[i-1]);
        betaBeta += sigma2*sigma2*(expBeta_S[i]*expBeta_S[i] - expBeta_S[i-1]*expBeta_S[i-1]);
        gammaGamma += (sigma1Bar*sigma1Bar + sigma2Bar*sigma2Bar)*(boundaryTimes[i]-boundaryTimes[i-1]);
        alphaBeta += sigma1Bar*sigma2*(expAlpha_S[i]*expBeta_S[i] - expAlpha_S[i-1]*expBeta_S[i-1]);
        alphaGamma += sigma1Bar*sigma1Bar*(expAlpha_S[i] - expAlpha_S[i-1]);
        betaGamma += sigma2*sigma1Bar*(expBeta_S[i] - expBeta_S[i-1]);
    }

    sigma1Bar = locSigmas[intervalEndpoint-1][1];
    sigma2Bar = locSigmas[intervalEndpoint-1][2];
    sigma2 = locSigmas[intervalEndpoint-1][3];

    alphaAlpha += sigma1Bar*sigma1Bar*(expAlpha_t*expAlpha_t - expAlpha_S[intervalEndpoint-1]*expAlpha_S[intervalEndpoint-1]);
    betaBeta += sigma2*sigma2*(expBeta_t*expBeta_t - expBeta_S[intervalEndpoint-1]*expBeta_S[intervalEndpoint-1]);
    gammaGamma += (sigma1Bar*sigma1Bar + sigma2Bar*sigma2Bar)*(t - boundaryTimes[intervalEndpoint-1]);
    alphaBeta += sigma1Bar*sigma2*(expAlpha_t*expBeta_t - expAlpha_S[intervalEndpoint-1]*expBeta_S[intervalEndpoint-1]);
    alphaGamma += sigma1Bar*sigma1Bar*(expAlpha_t - expAlpha_S[intervalEndpoint-1]);
    betaGamma += sigma2*sigma1Bar*(expBeta_t - expBeta_S[intervalEndpoint-1]);
    
    const DoubleArray& Alphas = this->getAlphas();
    double alpha = Alphas[0];
    double beta = Alphas[1];

    alphaAlpha /= (2*alpha);
    betaBeta /= (2*beta);
    // gammaGamma /= 1.0;
    alphaBeta /= (alpha+beta);
    alphaGamma /= alpha;
    betaGamma /= beta;

    expectations->clear();
    expectations->push_back(alphaAlpha);
    expectations->push_back(betaBeta);
    expectations->push_back(gammaGamma);
    expectations->push_back(alphaBeta);
    expectations->push_back(alphaGamma);
    expectations->push_back(betaGamma);
}

double EnergyInstVolBase::getInstVolRatio(
    const double /* alpha */,
    const double /* beta */,
    const double v1Bar,
    const double v1,
    const double v2Bar,
    const double v2,
    const double alphaFactor1,
    const double betaFactor1,
    const double alphaFactor2,
    const double betaFactor2) const
{

    double ratio1, ratio2, ratio3;
        double NormV1, NormV2, NormV2Bar;

    NormV1 = v1/v1Bar;
    NormV2 = v2/v1Bar;
    NormV2Bar = v2Bar/v1Bar;
    
    ratio1 = (1 + NormV1 / alphaFactor1 + NormV2 / betaFactor1) * (1 + NormV1 / alphaFactor1 + NormV2 / betaFactor1) + NormV2Bar * NormV2Bar;

    ratio2 = (1 + NormV1 / alphaFactor2 + NormV2 / betaFactor2) * (1 + NormV1 / alphaFactor2 + NormV2 / betaFactor2) + NormV2Bar * NormV2Bar;
  
    ratio3 = sqrt(ratio1/ratio2); 

    return(ratio3);
}


class EnergyInstVolBaseHelper
{ 
public:
	
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(EnergyInstVolBase, clazz);
        SUPERCLASS(EnergyVolBase);
        FIELD(alpha, "Alpha");
        FIELD(fFutureMaturityDates, "Future maturity date array");
        FIELD(fOptionExpiryDates, "Option expiry date array");

        FIELD(fAlphas, "Alpha and Beta array");
        FIELD_MAKE_TRANSIENT(fAlphas);
        FIELD(fMaxPricingDate, "Max Pricing Date");
        FIELD_MAKE_TRANSIENT(fMaxPricingDate);
        FIELD(fRatios, "Ratios array");
        FIELD_MAKE_TRANSIENT(fRatios);
        FIELD(fSigmas, "Sigmas Matrix");
        FIELD_MAKE_TRANSIENT(fSigmas);
    }
};

CClassConstSP const EnergyInstVolBase::TYPE =CClass::registerClassLoadMethod("EnergyInstVolBase", 
                  typeid(EnergyInstVolBase), EnergyInstVolBaseHelper::load);

// definition of TYPE for MarketWrapper template class
DEFINE_TEMPLATE_TYPE(EnergyInstVolBaseWrapper);

DRLIB_END_NAMESPACE
