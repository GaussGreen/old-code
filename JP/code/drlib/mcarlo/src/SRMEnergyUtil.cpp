//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : SRMEnergyUtil.cpp
//
//   Description : Helper for SRM Energy - used for holding intermediate data
//
//   Author      : 
//
//   Date        : 
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SRMEnergyUtil.hpp"
#include "edginc/EnergyInstVolRegular.hpp"
#include "edginc/EnergyInstVolSeasonal.hpp"
#include "edginc/EnergyInstVolExplicitRegular.hpp"
#include <cassert>

#if 0
#include <fstream> // for logging to file
#endif

DRLIB_BEGIN_NAMESPACE

///// constructs and populates SRMEnergyUtil object
SRMEnergyUtilBase::SRMEnergyUtilBase(
	int	_nbFactors,
    const DateTime & _today,
    EnergyFuturesCurveConstSP _futureCurve,
	const string & _corrInstrStart,
	const string & _corrInstrMaturity
	): // eg 6M
		nbFactors(_nbFactors),
        today(_today), 
		futureCurve(_futureCurve),
		corrInstrStart(_corrInstrStart),
		corrInstrMaturity(_corrInstrMaturity),
		underlyer(_futureCurve->getEnergyUnderlyer()),
		instVolBase(_futureCurve->getEnergyInstVolBase()),
        beta(-999.999)
{
	const DoubleArray & alphas = instVolBase->getAlphas();
	alpha = alphas[0];
}

// populates sqrtYearFracs with sqrt of year fracs between sim dates.
// todayIdx is the index of today in the diffusionDates array.
void SRMEnergyUtilBase::computeSqrtYearFrac(
	int todayIdx, 
	const DateTimeArray & diffusionDates,
	vector<double>& sqrtYearFracs
	)
{
	size_t nDiffusionSteps = diffusionDates.size() - todayIdx - 1;
	ASSERT(sqrtYearFracs.size() >= nDiffusionSteps);

	DateTime prevDate = diffusionDates[todayIdx];
	for (size_t i = 0, dateIdx = todayIdx + 1; 
		i < nDiffusionSteps; 
		++i, ++dateIdx) 
	{
		sqrtYearFracs[i] = sqrt( prevDate.yearFrac(diffusionDates[dateIdx]) );
		prevDate = diffusionDates[dateIdx];
	}	
}

/** Same as above except we assume the first diffusion date is the 
	base date: diffusionDates[0] == baseDate.
	sqrtYearFracs should be at least of size diffusionDates.size() - 1 */
void SRMEnergyUtilBase::computeSqrtYearFrac(
								const DateTimeArray & diffusionDates,
								vector<double>& sqrtYearFracs)
{
	ASSERT(sqrtYearFracs.size()+1 >= (size_t) diffusionDates.size());
	for (int i = 0; i+1 < diffusionDates.size(); ++i) 
		sqrtYearFracs[i] = sqrt( diffusionDates[i].yearFrac(diffusionDates[i+1]) );
}

// compute step-wise k factors (from diffusion time ti to ti+1): 
// used in states diffusion
void SRMEnergyUtilBase::computeKFactors(vector<double> & _vK, 
								  double _alpha, 
								  const vector<double> & _dtSqrts)
{
	ASSERT(_dtSqrts.size() == _vK.size());

	size_t i;
	for (i = 0; i < _dtSqrts.size(); ++i) {
		_vK[i] = exp(-_alpha * _dtSqrts[i]*_dtSqrts[i]);
	}
	return;
}

// compute cumulative (from 0 to Ti) k factors: 
// used in computing lnFwd after diffusion is performed
void SRMEnergyUtilBase::computeCumKFactors(vector<double> & _vKT, 
								   double _alpha, 
								   const DateTime & _today,
								   const DateTimeArray & _fwdMaturities)
{
	ASSERT(_fwdMaturities.size() == _vKT.size());

	double yrFrac;
	for (int i = 0; i < _fwdMaturities.size(); ++i) {
		yrFrac = _today.yearFrac(_fwdMaturities[i]);
		_vKT[i] = exp(-_alpha * yrFrac);
	}
	return;
}

/////// Energy model dependent util starts here: ////////////////////////

SRMEnergyUtilOil::SRMEnergyUtilOil(
	int	_nbFactors,
	const DateTime & _today,
	EnergyFuturesCurveConstSP _futureCurve,
	const string & _corrInstrStart,
	const string & _corrInstrMaturity
	)
	: 
	SRMEnergyUtilBase(
		_nbFactors,
		_today, 
		_futureCurve, 
		_corrInstrStart, 
		_corrInstrMaturity
		)
{
    const string method = "SRMEnergyUtilOil::SRMEnergyUtilOil";
	static CClassConstSP instVolRegularClazz = CClass::forName("EnergyInstVolRegular");
    static CClassConstSP instVolExplicitRegularClazz = CClass::forName("EnergyInstVolExplicitRegular");

	if (instVolRegularClazz->isAssignableFrom(instVolBase->getClass())) 
	{
		EnergyInstVolRegularConstSP instVolRegular = 
			EnergyInstVolRegularConstSP::dynamicCast(instVolBase);

        // calibrate x and z
        x = SRMEnergyUtilOil::calibX(
            instVolRegular->getAlphas(), 
            instVolRegular->getFirstContract(), 
            instVolRegular->getSecondContract(), 
            instVolRegular->getVolRatio(), 
            instVolRegular->getInstCorr());

        z = SRMEnergyUtilOil::calibZ(
            x, 
            instVolRegular->getAlphas(), 
            instVolRegular->getFirstContract(), 
            instVolRegular->getSecondContract(), 
            instVolRegular->getVolRatio());

        y = SRMEnergyUtilOil::calibY();

        // calibrate model vols
        sigmas = SRMEnergyUtilOil::calibSigmas(
            instVolRegular->getBaseDate(), 
            instVolRegular->getFutureMaturityDates(), 
            instVolRegular->getOptionExpiryDates(), 
            instVolRegular->getVols(), 
            instVolRegular->getAlphas(), 
            x, 
            z);

#if 0
		// show calibrated model vols
		ofstream file("c:/debugEnergySRM3.txt", ios_base::app);
		if (file.is_open()) {
			file << "Calibrated Model Vols: " << endl;
			const DateTimeArray & bmDates = instVolRegular->getOptionExpiryDates();
			for (size_t i = 0; i < bmDates.size(); ++i) {
				file << bmDates[i].getDate() << "\t" 
					<< sigmas[i][0] << "\t"
					<< sigmas[i][1] << "\t"
					<< sigmas[i][2] << "\t"
					<< sigmas[i][3] << endl;
			}
			file << endl;
		}
#endif
	}
	else if (instVolExplicitRegularClazz->isAssignableFrom(instVolBase->getClass()))
    { 
        // otherwise just grab sigmas without calibration
        EnergyInstVolExplicitRegularConstSP instVolExplicitRegular = 
            EnergyInstVolExplicitRegularConstSP::dynamicCast(instVolBase);

        x = instVolExplicitRegular->getX();
        y = instVolExplicitRegular->getY();
        z = instVolExplicitRegular->getZ();
        sigmas = instVolExplicitRegular->getSigmas();
	}
    else
        throw ModelException(method, "Unsupported energy volatility class!");
}

/** extend calibrated model vols from benchmark dates to timeline.
	note: today is assumed to be in the timeline, and we only
	consider extending the model vols to diffusion dates that
	are strictly > today */
void SRMEnergyUtilOil::setTimeLine(DateTimeArrayConstSP timeline)
{
	const string method = "SRMEnergyUtilOil::setTimeLine()";

    // first get future dates from timeline (not include today)
    DateTimeArray futureTimeline = today.getFutureDates((*timeline));

    // step 1: compute total vol rate at each sim date
    DoubleArray totalVolRates = SRMEnergyUtilOil::calcTotalVolRatesAtSimDates(
        instVolBase->getOptionExpiryDates(),
        futureTimeline,
        futureTimeline);

    // step 2: pass the total vol rates and sim dates to recalibrate
    extSigmas = SRMEnergyUtilOil::calibSigmas(
        instVolBase->getBaseDate(), 
        futureTimeline, 
        futureTimeline, 
        totalVolRates, 
        instVolBase->getAlphas(), 
        x, 
        z);

#if 0
    // show extended model vols
    ofstream file("c:/debugEnergySRM3.txt", ios_base::app);
    if (file.is_open()) {
        file << "Extended Model Vols: " << endl;
        for (size_t i = 0; i < futureTimeline.size(); ++i) {
            file << futureTimeline[i].getDate() << "\t" 
                << extSigmas[i][0] << "\t"
                << extSigmas[i][1] << "\t"
                << extSigmas[i][2] << "\t"
                << extSigmas[i][3] << endl;
        }
        file << endl;
    }
#endif
}

// equivalent to EnergyInstVolCalibrated::getNormalizedSigma1
double SRMEnergyUtilOil::calibX( 
    const DoubleArray & alphas, 
    int firstContract, // should pass secondContract here (due to some old bug in the energy market data code)
    int secondContract, // should pass firstContract here (due to some old bug in the energy market data code)
    double instVolRatio,
    double instCorr) const
{
    double alpha = alphas[0];
    double NormalizedSigma1;
    double firstContractFactor = exp(-alpha*firstContract/365.0);
    double secondContractFactor = exp(-alpha*secondContract/365.0);
    //double contractTimeDiff = secondContractFactor-firstContractFactor;
    double instVolRatioSquared = instVolRatio*instVolRatio;
    double x1,x2,y1,y2;

    //x1 = contractTimeDiff * (1 - 2 * InstVolRatio * InstCorr + instVolRatioSquared);
    //y1 = firstContractFactor * firstContractFactor *  (1 - InstVolRatio * InstCorr);
    //y2 = - firstContractFactor * secondContractFactor * (1 - instVolRatioSquared);
    //y3 = secondContractFactor * secondContractFactor * InstVolRatio * ( InstCorr - InstVolRatio);
    //x2 = y1 + y2 +y3;
    x1 = 1 - (2 * instVolRatio * instCorr) + instVolRatioSquared;
    y1 = instVolRatio * instCorr * (secondContractFactor + firstContractFactor);
    y2 = -((instVolRatioSquared * secondContractFactor) + firstContractFactor);
    x2 = y1 + y2;
    NormalizedSigma1 = x1/x2;

    return NormalizedSigma1;
}

// equivalent to EnergyInstVolCalibrated::getNormalizedSigma2bar
double SRMEnergyUtilOil::calibZ(
    double x,
    const DoubleArray & alphas, 
    int firstContract, // should pass secondContract here (due to some old bug in the energy market data code)
    int secondContract, // should pass firstContract here (due to some old bug in the energy market data code)
    double instVolRatio) const
{
    // uses x as input to simplify calculation

    double alpha = alphas[0];
    double NormalizedSigma2bar;
    double firstContractFactor = exp(-alpha*firstContract/365.0);
    double secondContractFactor = exp(-alpha*secondContract/365.0);
    double instVolRatioSquared = instVolRatio*instVolRatio;
    double x1,x2;

    x1 = (1 + x * firstContractFactor);
    x2 = (1 + x * secondContractFactor);

    NormalizedSigma2bar = sqrt ((x1*x1 - instVolRatioSquared *x2*x2) / (instVolRatioSquared -1));
    return NormalizedSigma2bar;
}

// equivalent to EnergyInstVolBase::getInstVolRatio
double SRMEnergyUtilOil::getInstVolRatio(    
    const double alpha,
    const double beta,
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

// equivalent to EnergyInstVolRegular::deriveSigmas
DoubleMatrix SRMEnergyUtilOil::calibSigmas(
    const DateTime & baseDate,
    const DateTimeArray & fFutureMaturityDates, 
    const DateTimeArray & fOptionExpiryDates, 
    const DoubleArray & vols, // atm vols,
    const DoubleArray & alphas, 
    double x,
    double z) const
{
    static const string method = "SRMEnergyUtilOil::calibSigmas()";

    // using local variables
    DateTimeArray sigmasDates=fFutureMaturityDates;
    double alpha = alphas[0];
    //double beta = alphas[1];
    if (alpha<=0.0)
        ModelException(method,"Alpha is not strictly positive");
    //if (beta<=0.0)
    //    ModelException(method,"Beta is not strictly positive");

    // Create placeholders for calibration outputs and intermediate calculations
    double v1, v1Bar, v2Bar, v2;
    double v1_0, v1Bar_0, v2Bar_0, v2_0;
    double x1, x2, x3/*, x4, x5, x6*/;

    int i;
    double vol;
    double expiryPrev, expiryNow, maturityPrev, maturityNow;
    double expiryPrevFactor, expiryNowFactor, maturityPrevFactor, maturityNowFactor;
    //double expiryPrevFactorS, expiryNowFactorS, maturityPrevFactorS, maturityNowFactorS;
    double sum1, sum2, sum3/*, sum4, sum5, sum6*/;
    double ratioInst, ratioFactor/*, ratioFactorS*/;

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
    //expiryNowFactorS = 1.0;
    maturityNow = 0.0;
    maturityNowFactor = 1.0;
    //maturityNowFactorS = 1.0;
    sum1 = 0.0;
    sum2 = 0.0;
    sum3 = 0.0;
    //sum4 = 0.0;
    //sum5 = 0.0;
    //sum6 = 0.0;
    ratioFactor = exp(alpha*sigmasDates[isExpiryDate].daysDiff(baseDate)/365.0);
    //ratioFactorS = exp(beta*sigmasDates[isExpiryDate].daysDiff(baseDate)/365.0);

    // Loop through all points on the curve iteratively
    for (i=isExpiryDate; i<numBenchmarks; ++i)
    {
        expiryPrev = expiryNow;
        expiryPrevFactor = expiryNowFactor;
        //expiryPrevFactorS = expiryNowFactorS;
        expiryNow = optionExpiries[i].daysDiff(baseDate)/365.0;
        expiryNowFactor = exp(alpha*expiryNow);
        //expiryNowFactorS = exp(beta*expiryNow);
        maturityPrev = maturityNow;
        maturityPrevFactor = maturityNowFactor;
        //maturityPrevFactorS = maturityNowFactorS;
        maturityNow = sigmasDates[i].daysDiff(baseDate)/365.0;
        maturityNowFactor = exp(alpha*maturityNow);
        //maturityNowFactorS = exp(beta*maturityNow);
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
            /*+ sum4 / (maturityNowFactorS * maturityNowFactorS) 
            + sum5 / maturityNowFactorS 
            + sum6 / (maturityNowFactor * maturityNowFactorS)*/; 

            // integrated calib terms
            x1 = NSigma1*NSigma1 * (expiryNowFactor*expiryNowFactor - maturityPrevFactor*maturityPrevFactor) 
            / (2.0 * alpha * maturityNowFactor*maturityNowFactor);
        x2 = 2.0 * NSigma1 * (expiryNowFactor - maturityPrevFactor) 
            / (alpha * maturityNowFactor);
        x3 = (1 + NSigma2Bar * NSigma2Bar) * (expiryNow - maturityPrev);

        //if (NSigma2 == 0)
        //{
        //x4 = 0;
        //x5 = 0;
        //x6 = 0;
        //}
        //else  // seasonal calib terms
        //{
        //    x4 = NSigma2*NSigma2 * (expiryNowFactorS*expiryNowFactorS - maturityPrevFactorS*maturityPrevFactorS) 
        //        / (2.0 * beta * maturityNowFactorS*maturityNowFactorS);
        //    x5 = 2.0 * NSigma2 * (expiryNowFactorS - maturityPrevFactorS) 
        //        / (beta * maturityNowFactorS);
        //    x6 = 2.0 * NSigma1 * NSigma2 * (expiryNowFactor * expiryNowFactorS - maturityPrevFactor * maturityPrevFactorS) 
        //        / ((alpha + beta) * maturityNowFactor*maturityNowFactorS);
        //}

        v1Bar = sqrt((VarTotal- VarPrev) / (x1 + x2 + x3 /*+ x4 + x5 + x6*/));
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
        ratioInst = SRMEnergyUtilOil::getInstVolRatio(alpha,
            0 /*beta*/,
            v1Bar_0,
            v1_0,
            v2Bar_0,
            v2_0,
            maturityNowFactor,
            1.0 /*maturityNowFactorS*/,
            ratioFactor,
            1.0 /*ratioFactorS*/);

        // integrated previous variance terms
        sum1 += v1*v1 * (maturityNowFactor*maturityNowFactor - maturityPrevFactor*maturityPrevFactor) 
            / (2.0 * alpha);
        sum2 += 2.0 * v1 * v1Bar * (maturityNowFactor - maturityPrevFactor) 
            / alpha;
        sum3 += (v1Bar*v1Bar + v2Bar*v2Bar) * (maturityNow - maturityPrev);
        //if (NSigma2 != 0)
        //{
        //    sum4 += v2*v2 * (maturityNowFactorS*maturityNowFactorS - maturityPrevFactorS*maturityPrevFactorS) 
        //        / (2.0 * beta);
        //    sum5 += 2.0 * v2 * v1Bar * (maturityNowFactorS - maturityPrevFactorS) 
        //        / beta;
        //    sum6 += 2.0 * v1 * v2 * (maturityNowFactor * maturityNowFactorS - maturityPrevFactor * maturityPrevFactorS) 
        //        / (alpha + beta);
        //}

        fSigmasOut[i][0] = v1;
        fSigmasOut[i][1] = v1Bar;
        fSigmasOut[i][2] = v2Bar;
        fSigmasOut[i][3] = v2;

        //fRatios.push_back(ratioInst);
    }
    //fMaxPricingDate = (i==0) ? baseDate : fFutureMaturityDates[i-1];

    return fSigmasOut; // calibrated model vols
}


vector<double> SRMEnergyUtilOil::getVectorGamma() const 
{
	vector<double> gamma(nbFactors);
	DateTime maturity = underlyer->expiryDate(corrInstrMaturity);
	double length = today.yearFrac(maturity);

	double alpha1 = alpha;
	double alpha2 = 0.0;

	switch (nbFactors) {
	case 1:
		gamma[0] = y + x * exp(-alpha1 * length);
		break;
	case 2:
		gamma[0] = y + x * exp(-alpha1 * length);
		gamma[1] = z * exp(-alpha2 * length);
		break;
	}

	return gamma;
}

/*
** Do a PCA to find a set of orthogonal basis for energy.
** Output the sqrt of the inter-energy corr matrix, e.g. the A as in 
** Dickinson's paper.
*/
DoubleMatrix SRMEnergyUtilOil::getMatrixAByPCA() const
{
	const string & method = "SRMEnergyUtilOil::getMatrixAByPCA";

#define  Nbt2Mat   25
	// hard-coded fwd maturities (in months) for doing energy PCA
	double    t2Mat[Nbt2Mat] = {0.5,1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 
								9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 
								17.5, 18.5, 19.5, 20.5, 21.5, 22.5, 23.5, 24.5};
	bool areBetasAllDiff = true;
	long i, j, k, l;
	
	/*
	double  **rCovMtx = NULL;   // [1..Nbt2Mat][1..Nbt2Mat]       
	double  **A       = NULL;   // [1..NbFactor][1..NbFactor]     
	double   *D       = NULL;   // [1..Nbt2Mat]                   
	double   *e       = NULL;   // [1..Nbt2Mat]                   
	double  **b       = NULL;   // [1..NbFactor][1..1]            
	double  **RHO     = NULL;   // [0..NbFactor-1][0..NbFactor-1]  
	double  **Alpha   = NULL;   // [0..NbFactor][0..Nbt2Mat ]     
	*/
	DoubleMatrix rCovMtx;
	DoubleMatrix A;
	vector<double> D;
	vector<double> e;
	DoubleMatrix b;
	DoubleMatrix RHO;
	DoubleMatrix Alpha;

	DoubleMatrix OMtx = DoubleMatrix(nbFactors, nbFactors);
	vector<double> beta(nbFactors); beta[0] = alpha; beta[1] = 0.0;
	vector<double> rho(3, 0.0); // assume independence among energy factors


	if (nbFactors == 1L)	{
		OMtx[0][0] = 1.0;
		return OMtx;
	}

	/* create the RHO matrix */
	/*RHO = (double **) DR_Matrix(DOUBLE,
		0, NbFactor - 1,
		0, NbFactor - 1);
	if (RHO == NULL) goto RETURN;

	Alpha = (double **) DR_Matrix(DOUBLE,
		0, NbFactor - 1,
		0, Nbt2Mat - 1);
	if (Alpha == NULL) goto RETURN;
	*/
    RHO = DoubleMatrix(nbFactors, nbFactors);
	Alpha = DoubleMatrix(nbFactors, Nbt2Mat);

	switch (nbFactors)
	{
	case 2L:
		RHO[0][0] = RHO[1][1] = 1.0;
		RHO[0][1] = RHO[1][0] = rho[0];
		break;
	default:
		throw ModelException(method, "Invalid number of energy factor!");
	}

	/* find out if betas are all different */
	switch (nbFactors)
	{
	case 2L:
		areBetasAllDiff = (!Maths::equals(beta[0],beta[1]));
		break;
	default:
		throw ModelException(method, "Invalid number of energy factor!");
	}

	if (!areBetasAllDiff) {
		/* Some betas are the same, therefore degeneration. */
		/* Use the lower triangular decomp instead          */
		throw ModelException(method, "Betas levels are the same. "
			"Cannot complete decomposition!");
	}

	/* create the rate covariance matrix */
	/* rCovMtx[1..Nbt2Mat][1..Nbt2Mat]   */

	/*
	rCovMtx = (double **) DR_Matrix(DOUBLE,
		1, Nbt2Mat,
		1, Nbt2Mat);
	if (rCovMtx == NULL) goto RETURN;
	*/

	for (i = 0; i < Nbt2Mat; i++)
	{
		Alpha[0][i] = y + x * exp(-beta[0] * t2Mat[i] / 12.0);
		Alpha[1][i] = z * exp(-beta[1] * t2Mat[i] / 12.0);
	}

	/*
	for (i=1; i<=Nbt2Mat; i++)
	{
		for (l=1; l<=Nbt2Mat; l++)
		{
			rCovMtx[i][l] = 0. ; 
			for (j = 0; j < NbFactor; j++)
			{
				for(k = 0; k < NbFactor; k++)
				{
					rCovMtx[i][l] += 
						Alpha[k][l - 1] * Alpha[j][i - 1] * RHO[j][k];
				} // k 
			} // j 
		} // l 
	} // i 
	*/

	/* create the rate covariance matrix */
	/* rCovMtx[1..Nbt2Mat][1..Nbt2Mat]   */
	rCovMtx = DoubleMatrix(Nbt2Mat, Nbt2Mat);
	for (i = 0; i < Nbt2Mat; i++) {
		for (l = 0; l < Nbt2Mat; l++) {
			rCovMtx[i][l] = 0.0;
			for (j = 0; j < nbFactors; j++) {
				for (k = 0; k < nbFactors; k++) {
					rCovMtx[i][l] += 
						Alpha[k][l] * Alpha[j][i] * RHO[j][k];
				} // k 
			} // j 
		} // l 
	} // i 

	/* eigenise the rate covariance matrix              */
	/* on completion, rCovMtx contains the eigenvectors */
	/* and D contains the eigenvalues                   */
	/*
	D = (double *)DR_Array(DOUBLE, 1, Nbt2Mat);
	e = (double *)DR_Array(DOUBLE, 1, Nbt2Mat);
	if ((D == NULL) || (e == NULL)) goto RETURN;

	tred2(rCovMtx,Nbt2Mat,D,e);
	if (tqli(D,e,Nbt2Mat,rCovMtx) == FAILURE) goto RETURN;
	eigsrt(D,rCovMtx,Nbt2Mat);

	// create the inverse of A //
	A = (double **) DR_Matrix(DOUBLE,
		1, NbFactor,
		1, NbFactor);
	if (A == NULL) goto RETURN;

	for (i=1; i<=NbFactor; i++)
	{
		for (l=1; l<=NbFactor; l++)
		{
			A[i][l] = 0.0;
			for (k=1; k<=Nbt2Mat; k++)
			{
				A[i][l] += 1.0/sqrt(D[i]) *
					rCovMtx[k][i] *
					Alpha[l-1][k-1] ;
			} // k 
		} // l 
	} // i  
	*/
	EigenVectorAnalysisSP eigenVecData(rCovMtx.computeEigenVectors());
	/* create the inverse of A */
	A = DoubleMatrix(nbFactors, nbFactors);
	for (i = 0; i < nbFactors; i++) {
		// sqrt the eigenValues
		double& eigenValue = (*eigenVecData->eigenValues)[i];
		if (!Maths::isPositive(eigenValue)){
			// we do 1/sqrt so a small +ve is bad news too
			throw ModelException(method, "Zero eigenvalue!");
		}
		eigenValue = sqrt(eigenValue);
		for (l = 0; l < nbFactors; l++) {
			A[i][l] = 0.0;
			for (k = 0; k < Nbt2Mat; k++) {
				/*
				A[i][l] += 1.0/(*eigenVecData->eigenValues)[i] *
				(*eigenVecData->eigenVectors)[i][k] *
				(Alpha[l] * exp(-beta[l] * t2Mat[k] / 12.0));
				*/
				A[i][l] += 1.0/(*eigenVecData->eigenValues)[i] *
					(*eigenVecData->eigenVectors)[i][k] *
					Alpha[l][k];
			} // k 
		} // l 
	} // i 

	// find A by inverting inverse of A //
	/*
	b = (double **) DR_Matrix(DOUBLE,
		1, NbFactor,
		1, 1);
	for (i=1; i<=NbFactor; i++) b[i][1] = 0.0;

	if (gaussjORI(A, NbFactor, b, 1) == FAILURE) goto RETURN;

	// finally, place answer in OMtx //
	for (i=0; i<NbFactor; i++)
	{
		for(j=0; j<NbFactor; j++)
		{
			OMtx[i][j] = A[i+1][j+1];
		}
	}
	*/
	DoubleMatrix corrDecomposition = A.computeInverse(); 

#undef  Nbt2Mat

	return corrDecomposition;
}

/** computes total volatility rates at simulation dates from calibrated
	model vols by assuming option expires = underlying future maturities.

	The volatility rate at simulation dates s_i will be like the option 
	implied vol rate, with both the option maturity and underlying future 
	maturity equal s_i (equivalent to market quotes). 
**/
DoubleArray SRMEnergyUtilOil::calcTotalVolRatesAtSimDates(
	const DateTimeArray & bmDates, // vol benchmark dates
	const DateTimeArray & timeline, // simulation timeline 
    const DateTimeArray & maturities // future maturity associated with each timeline date
	) const
{
	DoubleArray impVolRates(timeline.size());

	double G1, G2, G3, G4; // see document for definitions of these G's
	double dG1, dG2, dG3, dG4;
	vector<int> sim2BmIndices = DateTime::getCeilingProjection(timeline, bmDates);
	size_t startIdx;
	int prevSim2BmIdx, currSim2BmIdx;
	DateTime prevSimDate, currSimDate;
    DateTime currMaturity;
	double sigma1Bar, sigma2Bar, sigma1, currTotalVar;
	double expAlphaT, prevExpAlphaT; // e^(alpha*time) terms
	double t, T, dt;

	G1 = 0.0;
	G2 = 0.0;
	G3 = 0.0;
	G4 = 0.0;
	startIdx = 0;
	prevSim2BmIdx = 0;
	prevSimDate = today;
	prevExpAlphaT = 1.0;

	if (timeline[0] == today) { // check if 1st sim date == today
		impVolRates[0] = 0;
		++startIdx;
	}

	// loop thru each date on the timeline:
	for (int i = startIdx; i < timeline.size(); ++i) 
	{
		currSimDate = timeline[i];
		currSim2BmIdx = sim2BmIndices[i];
        currMaturity = maturities[i];
		t = today.yearFrac(currSimDate);
        T = today.yearFrac(currMaturity);

		if (prevSim2BmIdx == currSim2BmIdx) 
		{
			dt = prevSimDate.yearFrac(currSimDate);
			expAlphaT = exp(alpha*t);
			sigma1 = sigmas[currSim2BmIdx][0];
			sigma1Bar = sigmas[currSim2BmIdx][1];
			sigma2Bar = sigmas[currSim2BmIdx][2];

			dG1 = sigma1Bar * sigma1Bar * dt;
			dG2 = sigma2Bar * sigma2Bar * dt;
			dG4 = sigma1Bar * sigma1 * (expAlphaT - prevExpAlphaT)/alpha;
			dG3 = sigma1 * sigma1 * (expAlphaT*expAlphaT - prevExpAlphaT*prevExpAlphaT)*0.5/alpha;
		}
		else {
			// take care of first stub:
			dt = prevSimDate.yearFrac(bmDates[prevSim2BmIdx]);
			expAlphaT = exp(alpha * today.yearFrac(bmDates[prevSim2BmIdx]));
			sigma1 = sigmas[prevSim2BmIdx][0];
			sigma1Bar = sigmas[prevSim2BmIdx][1];
			sigma2Bar = sigmas[prevSim2BmIdx][2];
			
			dG1 = sigma1Bar * sigma1Bar * dt;
			dG2 = sigma2Bar * sigma2Bar * dt;
			dG4 = sigma1Bar * sigma1 * (expAlphaT - prevExpAlphaT)/alpha;
			dG3 = sigma1 * sigma1 * (expAlphaT*expAlphaT - prevExpAlphaT*prevExpAlphaT)*0.5/alpha;
			
			prevExpAlphaT = expAlphaT;

			// then do the middle intervals:
			for (int j = prevSim2BmIdx + 1; j <= currSim2BmIdx - 1; ++j) {
				dt = bmDates[j - 1].yearFrac(bmDates[j]);
				expAlphaT = exp(alpha * today.yearFrac(bmDates[j]));
				sigma1 = sigmas[j][0];
				sigma1Bar = sigmas[j][1];
				sigma2Bar = sigmas[j][2];
				
				dG1 += sigma1Bar * sigma1Bar * dt;
				dG2 += sigma2Bar * sigma2Bar * dt;
				dG4 += sigma1Bar * sigma1 * (expAlphaT - prevExpAlphaT)/alpha;
				dG3 += sigma1 * sigma1 * (expAlphaT*expAlphaT - prevExpAlphaT*prevExpAlphaT)*0.5/alpha;

				prevExpAlphaT = expAlphaT;
			}

			// then take care of the end stub:
			dt = bmDates[currSim2BmIdx - 1].yearFrac(currSimDate);
			expAlphaT = exp(alpha*t);
			sigma1 = sigmas[currSim2BmIdx][0];
			sigma1Bar = sigmas[currSim2BmIdx][1];			
			sigma2Bar = sigmas[currSim2BmIdx][2];

			dG1 += sigma1Bar * sigma1Bar * dt;
			dG2 += sigma2Bar * sigma2Bar * dt;
			dG4 += sigma1Bar * sigma1 * (expAlphaT - prevExpAlphaT)/alpha;
			dG3 += sigma1 * sigma1 * (expAlphaT*expAlphaT - prevExpAlphaT*prevExpAlphaT)*0.5/alpha;
		}

		// update G's
		G1 += dG1;
		G2 += dG2;
		G3 += dG3;
		G4 += dG4;

		// compute total variance at current sim date from G's
		//currTotalVar = G1 + G2 + exp(-2*alpha*t)*G3 + 2*exp(-alpha*t)*G4;
        currTotalVar = G1 + G2 + exp(-2*alpha*T)*G3 + 2*exp(-alpha*T)*G4;

		// compute total vol rate at current sim date
		impVolRates[i] = sqrt(currTotalVar / t);

		// update book keeping variables:
		prevSimDate = currSimDate;
		prevSim2BmIdx = currSim2BmIdx;
		prevExpAlphaT = expAlphaT;
	}

	return impVolRates;
}


SRMEnergyUtilTier2::SRMEnergyUtilTier2(
	int	_nbFactors,
	const DateTime & _today,
	EnergyFuturesCurveConstSP _futureCurve,
	const string & _corrInstrStart,
	const string & _corrInstrMaturity
	)
	: 
	SRMEnergyUtilBase(
		_nbFactors,
		_today, 
		_futureCurve, 
		_corrInstrStart, 
		_corrInstrMaturity
		)
    {
        // get model instantaneous vols (calibration is not needed for spreads):
        sigmas = instVolBase->getSigmas();

#if 0
		// show calibrated model vols
		ofstream file("c:/debugEnergySRM3.txt", ios_base::app);
		if (file.is_open()) {
			file << "Calibrated Spread Model Vols: " << endl;
			const DateTimeArray & sigmaDates = instVolBase->getOptionExpiryDates();
			for (size_t i = 0; i < sigmaDates.size(); ++i) {
				file << sigmaDates[i].getDate() << "\t" << sigmas[i][0] << endl;
			}
			file << endl;
		}
#endif
    }

vector<double> SRMEnergyUtilTier2::getVectorGamma() const
{
    const string & method = "SRMEnergyUtilTier2::getVectorGamma";
    vector<double> gamma;
    if (nbFactors == 1) 
    {
        DateTime maturity = underlyer->expiryDate(corrInstrMaturity);
	    double length = today.yearFrac(maturity);
        gamma = vector<double>(1);
        gamma[0] = exp(-alpha * length);
    }
    else 
        throw ModelException(method, "Multi-factor energy spread diffusion not support!");

    return gamma;
}

DoubleMatrix SRMEnergyUtilTier2::getMatrixAByPCA() const
{
    const string & method = "SRMEnergyUtilTier2::getMatrixAByPCA";
    DoubleMatrix A;
    if (nbFactors == 1) 
    {
        A = DoubleMatrix(1, 1);
        A[0][0] = 1.0;
    }
    else
        throw ModelException(method, "Multi-factor energy spread diffusion not support!");

    return A;
}

// Extend input instantaneous vols to the timeline
void SRMEnergyUtilTier2::setTimeLine(DateTimeArrayConstSP timeline)
{
    extSigmas = DoubleMatrix(timeline->size(), 1);
    const DateTimeArray & sigmaDates = instVolBase->getOptionExpiryDates();
    vector<int> sim2SigmaDateIdxs = DateTime::getCeilingProjection((*timeline), sigmaDates);
    
    DateTime currSimDate;
    DateTime prevSimDate;
    int currSim2SigmaDateIdx;
    int prevSim2SigmaDateIdx;
    double A, sumA, sumSigmaSqrA;
    double exp2AlphaT, prevExp2AlphaT;

    prevSimDate = today;
    prevSim2SigmaDateIdx = 0;
    prevExp2AlphaT = 1.0;

    for (int i = 0; i < timeline->size(); ++i) 
    {
        currSim2SigmaDateIdx = sim2SigmaDateIdxs[i];
        currSimDate = (*timeline)[i];

        if (currSim2SigmaDateIdx == prevSim2SigmaDateIdx) // just copy the sigmas:
        {
            exp2AlphaT = exp(2.0 * alpha * today.yearFrac(currSimDate));
            extSigmas[i][0] = sigmas[currSim2SigmaDateIdx][0];
        }
        else { 
            // do the beginning stub:
            exp2AlphaT = exp(2.0 * alpha * today.yearFrac(sigmaDates[prevSim2SigmaDateIdx]));
            sumA = exp2AlphaT - prevExp2AlphaT;
            sumSigmaSqrA = sigmas[prevSim2SigmaDateIdx][0] * sigmas[prevSim2SigmaDateIdx][0] * sumA;

            prevExp2AlphaT = exp2AlphaT;

            // then do the middle intervals:
            for (int j = prevSim2SigmaDateIdx + 1; j <= currSim2SigmaDateIdx - 1; ++j)
            {
                exp2AlphaT = exp(2.0 * alpha * today.yearFrac(sigmaDates[j]));
                A = exp2AlphaT - prevExp2AlphaT;
                sumA += A;
                sumSigmaSqrA += sigmas[j][0] * sigmas[j][0] * A;                

                prevExp2AlphaT = exp2AlphaT;
            }

            // lastly do the end stub:
            exp2AlphaT = exp(2.0 * alpha * today.yearFrac(currSimDate));
            A = exp2AlphaT - prevExp2AlphaT;
            sumA += A;
            sumSigmaSqrA += sigmas[currSim2SigmaDateIdx][0] * sigmas[currSim2SigmaDateIdx][0] * A;

            // compute the instantaneous vol for the current sim interval:
            extSigmas[i][0] = sqrt(sumSigmaSqrA / sumA);
        }

        // update book keeping variables:
        prevSimDate = currSimDate;
        prevSim2SigmaDateIdx = currSim2SigmaDateIdx;
        prevExp2AlphaT = exp2AlphaT;
    }

#if 0
	// show extended model vols
	ofstream file("c:/debugEnergySRM3.txt", ios_base::app);
	if (file.is_open()) {
		file << "Extended Tier2 Spread Model Vols: " << endl;
		for (size_t i = 0; i < timeline->size(); ++i) {
			file << (*timeline)[i].getDate() << "\t" << extSigmas[i][0] << endl;
		}
		file << endl;
	}
#endif
}



/** Create appropriate instance of energy util class */
SRMEnergyUtilBaseSP SRMEnergyUtilCreate(
	string				        _modelType,
	int							_nbFactors,
	const DateTime &			_today,
	EnergyFuturesCurveConstSP	_futureCurve,
	const string &				_corrInstrStart,
	const string &				_corrInstrMaturity
	)
{
	const string & method = "SRMEnergyUtilCreate";
	SRMEnergyUtilBaseSP enrgUtilBase;

    if (CString::equalsIgnoreCase(_modelType, EnrgOilStr) ||
        CString::equalsIgnoreCase(_modelType, EnrgSampras)) {
		enrgUtilBase = SRMEnergyUtilBaseSP(
			new SRMEnergyUtilOil(
			    _nbFactors,
			    _today,
			    _futureCurve,
			    _corrInstrStart,
			    _corrInstrMaturity 
            )
	    );
    }
    else if (CString::equalsIgnoreCase(_modelType, EnrgGasStr)) {
		throw ModelException(method, "Energy gas model not available yet!");
    }
    else if (CString::equalsIgnoreCase(_modelType, EnrgPowStr)) {
		throw ModelException(method, "Energy power model not available yet!");
    }
    else 
        throw ModelException(method, "Unknown energy model type: " + _modelType);

	return enrgUtilBase;
}


DRLIB_END_NAMESPACE
