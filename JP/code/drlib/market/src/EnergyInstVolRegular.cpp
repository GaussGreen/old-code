//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyInstVolRegular.cpp
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
#include "edginc/EnergyInstVolRegular.hpp"


DRLIB_BEGIN_NAMESPACE

EnergyInstVolRegular::EnergyInstVolRegular():EnergyInstVolCalibrated(TYPE)
{
}

EnergyInstVolRegular::~EnergyInstVolRegular()
{
}

void EnergyInstVolRegular::validatePop2Object()
{
    static const string method = "EnergyInstVolRegular::validatePop2Object()";

    EnergyInstVolCalibrated::validatePop2Object();

    // do some checking...
    
    // Check appropriateness of the parameters
    if (fFirstContract<=0 || fSecondContract<=fFirstContract)
        ModelException(method,"Inconsistent choice of calibration contracts");
    if (fInstVolRatio<=0.0 || fInstVolRatio>=1.0)
        ModelException(method,"Instantaneous volatility ratio is not between 0 and 1");
    if (fInstCorr<=-1.0 || fInstCorr>=1.0)
        ModelException(method,"Instantaneous correlation is not between -1 and 1");

}

DoubleMatrix EnergyInstVolRegular::calibSigmas(
	const DateTimeArray & futureMaturities,
	const DateTimeArray & optionExpiries,
	const DoubleArray & atmVols,
	double x,
	double z) 
	const
{
	static const string method = "EnergyInstVolRegular::calibSigmas()";
	return deriveSigmas(futureMaturities, optionExpiries, atmVols, x, z);
}

double EnergyInstVolRegular::calibX() const
{
	double alpha = fAlphas[0];
	return getNormalizedSigma1(alpha,fSecondContract,fFirstContract,fInstVolRatio,fInstCorr);
}

double EnergyInstVolRegular::calibZ(double x) const
{
	double alpha = fAlphas[0];
	return getNormalizedSigma2bar(alpha,x,fSecondContract,fFirstContract,fInstVolRatio);
}

int EnergyInstVolRegular::deriveSigmas(
                           const int numBenchmarks,
                           const DoubleArray& vols) 
{
    static const string method = "EnergyInstVolRegular::deriveSigmas()";    

    // using local variables
    DateTimeArray sigmasDates=fFutureMaturityDates;
    DateTime baseDate = getBaseDate();
    double alpha = fAlphas[0];
    double beta = fAlphas[1];
    if (alpha<=0.0)
        ModelException(method,"Alpha is not strictly positive");
    if (beta<=0.0)
        ModelException(method,"Beta is not strictly positive");

    // Create placeholders for calibration outputs and intermediate calculations
    double v1, v1Bar, v2Bar, v2;
    double v1_0, v1Bar_0, v2Bar_0, v2_0;
    double x1, x2, x3, x4, x5, x6;

    int i;
    double vol;
    double expiryPrev, expiryNow, maturityPrev, maturityNow;
    double expiryPrevFactor, expiryNowFactor, maturityPrevFactor, maturityNowFactor;
    double expiryPrevFactorS, expiryNowFactorS, maturityPrevFactorS, maturityNowFactorS;
    double sum1, sum2, sum3, sum4, sum5, sum6;
    double ratioInst, ratioFactor, ratioFactorS;

    int expiryMonth;
    double NSigma1, NSigma2Bar, NSigma2;
    double VarTotal, VarPrev;

    fSigmas = DoubleMatrix(numBenchmarks, 4);

    // Define constants
    NSigma1 = getNormalizedSigma1(alpha,fSecondContract,fFirstContract,fInstVolRatio,fInstCorr);
    NSigma2Bar = getNormalizedSigma2bar(alpha,NSigma1,fSecondContract,fFirstContract,fInstVolRatio);

    // Check expiry dates so calibration is still possible even if the first option contract has expired
    DateTimeArray optionExpiries = fOptionExpiryDates;
    if (numBenchmarks>0 && optionExpiries[0]<=baseDate) 
        optionExpiries[0] = sigmasDates[0];

    // Check maturity dates so calibration is still possible even if the first future contract is expiring
    int isExpiryDate = (sigmasDates[0]==baseDate) ? 1 : 0;
    if (isExpiryDate)
    {
        fSigmas[0][0] = 0.0;
        fSigmas[0][1] = 0.0;
        fSigmas[0][2] = 0.0;
        fSigmas[0][3] = 0.0;

        fRatios.push_back(1.0);
    }
    
    expiryNow = 0.0;
    expiryNowFactor = 1.0;
    expiryNowFactorS = 1.0;
    maturityNow = 0.0;
    maturityNowFactor = 1.0;
    maturityNowFactorS = 1.0;
    sum1 = 0.0;
    sum2 = 0.0;
    sum3 = 0.0;
    sum4 = 0.0;
    sum5 = 0.0;
    sum6 = 0.0;
    ratioFactor = exp(alpha*sigmasDates[isExpiryDate].daysDiff(baseDate)/365.0);
    ratioFactorS = exp(beta*sigmasDates[isExpiryDate].daysDiff(baseDate)/365.0);

    // Loop through all points on the curve iteratively
    for (i=isExpiryDate; i<numBenchmarks; ++i)
    {
        expiryPrev = expiryNow;
        expiryPrevFactor = expiryNowFactor;
        expiryPrevFactorS = expiryNowFactorS;
        expiryNow = optionExpiries[i].daysDiff(baseDate)/365.0;
        expiryNowFactor = exp(alpha*expiryNow);
        expiryNowFactorS = exp(beta*expiryNow);
        maturityPrev = maturityNow;
        maturityPrevFactor = maturityNowFactor;
        maturityPrevFactorS = maturityNowFactorS;
        maturityNow = sigmasDates[i].daysDiff(baseDate)/365.0;
        maturityNowFactor = exp(alpha*maturityNow);
        maturityNowFactorS = exp(beta*maturityNow);
        expiryMonth = optionExpiries[i].toMDY().month - 1;
        //NSigma2 = fNSigma2s[expiryMonth];
		NSigma2 = 0;

        vol = vols[i];
        VarTotal = vol * vol * expiryNow;

        if (expiryNow<=0.0 || maturityNow<expiryNow || expiryNow<=maturityPrev)
            ModelException(method,"Unexpected expiry or maturity dates");
        if (vol<=0.0)
            ModelException(method,"Unexpected volatility level");

        VarPrev = sum1 / (maturityNowFactor * maturityNowFactor) 
                + sum2 / maturityNowFactor 
                + sum3 
                + sum4 / (maturityNowFactorS * maturityNowFactorS) 
                + sum5 / maturityNowFactorS 
                + sum6 / (maturityNowFactor * maturityNowFactorS); 

        // integrated calib terms
        x1 = NSigma1*NSigma1 * (expiryNowFactor*expiryNowFactor - maturityPrevFactor*maturityPrevFactor) 
                / (2.0 * alpha * maturityNowFactor*maturityNowFactor);
        x2 = 2.0 * NSigma1 * (expiryNowFactor - maturityPrevFactor) 
                / (alpha * maturityNowFactor);
        x3 = (1 + NSigma2Bar * NSigma2Bar) * (expiryNow - maturityPrev);

        if (NSigma2 == 0)
        {
            x4 = 0;
            x5 = 0;
            x6 = 0;
        }
        else  // seasonal calib terms
        {
                x4 = NSigma2*NSigma2 * (expiryNowFactorS*expiryNowFactorS - maturityPrevFactorS*maturityPrevFactorS) 
                    / (2.0 * beta * maturityNowFactorS*maturityNowFactorS);
            x5 = 2.0 * NSigma2 * (expiryNowFactorS - maturityPrevFactorS) 
                    / (beta * maturityNowFactorS);
            x6 = 2.0 * NSigma1 * NSigma2 * (expiryNowFactor * expiryNowFactorS - maturityPrevFactor * maturityPrevFactorS) 
                    / ((alpha + beta) * maturityNowFactor*maturityNowFactorS);
        }
        
        v1Bar = sqrt((VarTotal- VarPrev) / (x1 + x2 + x3 + x4 + x5 + x6));
        v1    = NSigma1 * v1Bar;
        v2Bar = fabs(NSigma2Bar * v1Bar);
        v2    = NSigma2 * v1Bar;    // 
    
        // Add some constraints so it is consistent with what we have in production
        if(v1 < 0)
        {
            v1 *= -1;
            v1Bar *= -1;
            v2 *= -1;
        }

        if (i == 0) 
        {
            v1_0 = v1;
            v1Bar_0 = v1Bar;
            v2Bar_0 = v2Bar;
            v2_0 = v2;
        }
        
        // inst vol ratio to first contract
        ratioInst = getInstVolRatio(alpha,
                                    beta,
                                    v1Bar_0,
                                    v1_0,
                                    v2Bar_0,
                                    v2_0,
                                    maturityNowFactor,
                                    maturityNowFactorS,
                                    ratioFactor,
                                    ratioFactorS);
        
        // integrated previous variance terms
        sum1 += v1*v1 * (maturityNowFactor*maturityNowFactor - maturityPrevFactor*maturityPrevFactor) 
                / (2.0 * alpha);
        sum2 += 2.0 * v1 * v1Bar * (maturityNowFactor - maturityPrevFactor) 
                / alpha;
        sum3 += (v1Bar*v1Bar + v2Bar*v2Bar) * (maturityNow - maturityPrev);
        if (NSigma2 != 0)
        {
                sum4 += v2*v2 * (maturityNowFactorS*maturityNowFactorS - maturityPrevFactorS*maturityPrevFactorS) 
                    / (2.0 * beta);
            sum5 += 2.0 * v2 * v1Bar * (maturityNowFactorS - maturityPrevFactorS) 
                    / beta;
            sum6 += 2.0 * v1 * v2 * (maturityNowFactor * maturityNowFactorS - maturityPrevFactor * maturityPrevFactorS) 
                    / (alpha + beta);
        }

        fSigmas[i][0] = v1;
        fSigmas[i][1] = v1Bar;
        fSigmas[i][2] = v2Bar;
        fSigmas[i][3] = v2;

        fRatios.push_back(ratioInst);
    }
    fMaxPricingDate = (i==0) ? baseDate : fFutureMaturityDates[i-1];

    return i; // number of benchmark points successfully calibrated
}


DoubleMatrix EnergyInstVolRegular::deriveSigmas(
	const DateTimeArray & fFutureMaturityDates, // future maturities
	const DateTimeArray & fOptionExpiryDates, // option maturities
	const DoubleArray & vols, // atm vols,
	double x,
	double z
	) const
{
	static const string method = "EnergyInstVolRegular::deriveSigmas()";

	// using local variables
	DateTimeArray sigmasDates=fFutureMaturityDates;
	DateTime baseDate = getBaseDate();
	double alpha = fAlphas[0];
	double beta = fAlphas[1];
	if (alpha<=0.0)
		ModelException(method,"Alpha is not strictly positive");
	if (beta<=0.0)
		ModelException(method,"Beta is not strictly positive");

	// Create placeholders for calibration outputs and intermediate calculations
	double v1, v1Bar, v2Bar, v2;
	double v1_0, v1Bar_0, v2Bar_0, v2_0;
	double x1, x2, x3, x4, x5, x6;

	int i;
	double vol;
	double expiryPrev, expiryNow, maturityPrev, maturityNow;
	double expiryPrevFactor, expiryNowFactor, maturityPrevFactor, maturityNowFactor;
	double expiryPrevFactorS, expiryNowFactorS, maturityPrevFactorS, maturityNowFactorS;
	double sum1, sum2, sum3, sum4, sum5, sum6;
	double ratioInst, ratioFactor, ratioFactorS;

	int expiryMonth;
	double NSigma1, NSigma2Bar, NSigma2;
	double VarTotal, VarPrev;

	int numBenchmarks = fOptionExpiryDates.size();
	DoubleMatrix fSigmasOut = DoubleMatrix(numBenchmarks, 4);

	// Define constants
	/*
	NSigma1 = getNormalizedSigma1(alpha,fSecondContract,fFirstContract,fInstVolRatio,fInstCorr);
	NSigma2Bar = getNormalizedSigma2bar(alpha,NSigma1,fSecondContract,fFirstContract,fInstVolRatio);
	*/
	NSigma1 = x;
	NSigma2Bar = z;

	// Check expiry dates so calibration is still possible even if the first option contract has expired
	DateTimeArray optionExpiries = fOptionExpiryDates;
	if (numBenchmarks>0 && optionExpiries[0]<=baseDate) 
		optionExpiries[0] = sigmasDates[0];

	// Check maturity dates so calibration is still possible even if the first future contract is expiring
	int isExpiryDate = (sigmasDates[0]==baseDate) ? 1 : 0;
	if (isExpiryDate)
	{
		fSigmasOut[0][0] = 0.0;
		fSigmasOut[0][1] = 0.0;
		fSigmasOut[0][2] = 0.0;
		fSigmasOut[0][3] = 0.0;

		//fRatios.push_back(1.0);
	}

	expiryNow = 0.0;
	expiryNowFactor = 1.0;
	expiryNowFactorS = 1.0;
	maturityNow = 0.0;
	maturityNowFactor = 1.0;
	maturityNowFactorS = 1.0;
	sum1 = 0.0;
	sum2 = 0.0;
	sum3 = 0.0;
	sum4 = 0.0;
	sum5 = 0.0;
	sum6 = 0.0;
	ratioFactor = exp(alpha*sigmasDates[isExpiryDate].daysDiff(baseDate)/365.0);
	ratioFactorS = exp(beta*sigmasDates[isExpiryDate].daysDiff(baseDate)/365.0);

	// Loop through all points on the curve iteratively
	for (i=isExpiryDate; i<numBenchmarks; ++i)
	{
		expiryPrev = expiryNow;
		expiryPrevFactor = expiryNowFactor;
		expiryPrevFactorS = expiryNowFactorS;
		expiryNow = optionExpiries[i].daysDiff(baseDate)/365.0;
		expiryNowFactor = exp(alpha*expiryNow);
		expiryNowFactorS = exp(beta*expiryNow);
		maturityPrev = maturityNow;
		maturityPrevFactor = maturityNowFactor;
		maturityPrevFactorS = maturityNowFactorS;
		maturityNow = sigmasDates[i].daysDiff(baseDate)/365.0;
		maturityNowFactor = exp(alpha*maturityNow);
		maturityNowFactorS = exp(beta*maturityNow);
		expiryMonth = optionExpiries[i].toMDY().month - 1;
		//NSigma2 = fNSigma2s[expiryMonth];
		NSigma2 = 0;

		vol = vols[i];
		VarTotal = vol * vol * expiryNow;

		if (expiryNow<=0.0 || maturityNow<expiryNow || expiryNow<=maturityPrev)
			ModelException(method,"Unexpected expiry or maturity dates");
		if (vol<=0.0)
			ModelException(method,"Unexpected volatility level");

		VarPrev = sum1 / (maturityNowFactor * maturityNowFactor) 
			+ sum2 / maturityNowFactor 
			+ sum3 
			+ sum4 / (maturityNowFactorS * maturityNowFactorS) 
			+ sum5 / maturityNowFactorS 
			+ sum6 / (maturityNowFactor * maturityNowFactorS); 

		// integrated calib terms
		x1 = NSigma1*NSigma1 * (expiryNowFactor*expiryNowFactor - maturityPrevFactor*maturityPrevFactor) 
			/ (2.0 * alpha * maturityNowFactor*maturityNowFactor);
		x2 = 2.0 * NSigma1 * (expiryNowFactor - maturityPrevFactor) 
			/ (alpha * maturityNowFactor);
		x3 = (1 + NSigma2Bar * NSigma2Bar) * (expiryNow - maturityPrev);

		if (NSigma2 == 0)
		{
			x4 = 0;
			x5 = 0;
			x6 = 0;
		}
		else  // seasonal calib terms
		{
			x4 = NSigma2*NSigma2 * (expiryNowFactorS*expiryNowFactorS - maturityPrevFactorS*maturityPrevFactorS) 
				/ (2.0 * beta * maturityNowFactorS*maturityNowFactorS);
			x5 = 2.0 * NSigma2 * (expiryNowFactorS - maturityPrevFactorS) 
				/ (beta * maturityNowFactorS);
			x6 = 2.0 * NSigma1 * NSigma2 * (expiryNowFactor * expiryNowFactorS - maturityPrevFactor * maturityPrevFactorS) 
				/ ((alpha + beta) * maturityNowFactor*maturityNowFactorS);
		}

		v1Bar = sqrt((VarTotal- VarPrev) / (x1 + x2 + x3 + x4 + x5 + x6));
		v1    = NSigma1 * v1Bar;
		v2Bar = fabs(NSigma2Bar * v1Bar);
		v2    = NSigma2 * v1Bar;    // 

		// Add some constraints so it is consistent with what we have in production
		if(v1 < 0)
		{
			v1 *= -1;
			v1Bar *= -1;
			v2 *= -1;
		}

		if (i == 0) 
		{
			v1_0 = v1;
			v1Bar_0 = v1Bar;
			v2Bar_0 = v2Bar;
			v2_0 = v2;
		}

		// inst vol ratio to first contract
		ratioInst = getInstVolRatio(alpha,
			beta,
			v1Bar_0,
			v1_0,
			v2Bar_0,
			v2_0,
			maturityNowFactor,
			maturityNowFactorS,
			ratioFactor,
			ratioFactorS);

		// integrated previous variance terms
		sum1 += v1*v1 * (maturityNowFactor*maturityNowFactor - maturityPrevFactor*maturityPrevFactor) 
			/ (2.0 * alpha);
		sum2 += 2.0 * v1 * v1Bar * (maturityNowFactor - maturityPrevFactor) 
			/ alpha;
		sum3 += (v1Bar*v1Bar + v2Bar*v2Bar) * (maturityNow - maturityPrev);
		if (NSigma2 != 0)
		{
			sum4 += v2*v2 * (maturityNowFactorS*maturityNowFactorS - maturityPrevFactorS*maturityPrevFactorS) 
				/ (2.0 * beta);
			sum5 += 2.0 * v2 * v1Bar * (maturityNowFactorS - maturityPrevFactorS) 
				/ beta;
			sum6 += 2.0 * v1 * v2 * (maturityNowFactor * maturityNowFactorS - maturityPrevFactor * maturityPrevFactorS) 
				/ (alpha + beta);
		}

		fSigmasOut[i][0] = v1;
		fSigmasOut[i][1] = v1Bar;
		fSigmasOut[i][2] = v2Bar;
		fSigmasOut[i][3] = v2;

		//fRatios.push_back(ratioInst);
	}
	//fMaxPricingDate = (i==0) ? baseDate : fFutureMaturityDates[i-1];

	return fSigmasOut; // calibrated model vols
}

class EnergyInstVolRegularHelper
{

public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(EnergyInstVolRegular, clazz);
        clazz->setPublic();
        SUPERCLASS(EnergyInstVolCalibrated);
        EMPTY_SHELL_METHOD(defaultInstVolRegular);

        FIELD(fFirstContract, "First Contract (days)");
        FIELD(fSecondContract, "Second Contract (days)");
        FIELD(fInstVolRatio, "Inst Vol Ratio");
        FIELD(fInstCorr, "Inst Corr");
    }

    static IObject* defaultInstVolRegular()
    {
        return new EnergyInstVolRegular();
    }

};

CClassConstSP const EnergyInstVolRegular::TYPE =CClass::registerClassLoadMethod("EnergyInstVolRegular", 
                  typeid(EnergyInstVolRegular), EnergyInstVolRegularHelper::load);

// definition of TYPE for MarketWrapper template class
DEFINE_TEMPLATE_TYPE(EnergyInstVolRegularWrapper);


DRLIB_END_NAMESPACE


