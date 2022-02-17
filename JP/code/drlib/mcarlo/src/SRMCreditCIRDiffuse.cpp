//----------------------------------------------------------------------------
//
//   Group       : QR Cross Asset
//
//   Filename    : SRMCreditCIRDiffuse.cpp
//
//   Description : CIR credit path generation
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/QMCCreditDiffuse.hpp"
#include "edginc/SRMCreditCIRDiffuse.hpp"
#include "edginc/SRMUtil.hpp"

#include "edginc/MemoryProfiler.hpp"
#include "edginc/Format.hpp"
#include "edginc/IQMCRNGManager.hpp"

#include <algorithm>

DRLIB_BEGIN_NAMESPACE

/** constructor */
SRMCreditCIRDiffuse::SRMCreditCIRDiffuse(IQMCDiffusibleInterestRateSP srmRatesDiffuse, bool isFullMC) :
    SRMCreditDiffuse(srmRatesDiffuse),
    srmCreditCIRUtil(NULL),
    c(-999.),
    d(-999.),
    h(-999.),
    beta(-999.),
    fullMC(isFullMC)
{}

/** destructor */
SRMCreditCIRDiffuse::~SRMCreditCIRDiffuse() {
}

/** initialization */
void SRMCreditCIRDiffuse::setSRMCreditCIRDiffuse(
    int                    _randomIndex,
    const DateTime&        _today,
    SRMCreditCIRUtilSP     _srmCreditCIRUtil,
    const double           _crFxCorr,
    const vector<double>&  _prob)    // historic dates
{
    static const string method("SRMCreditCIRDiffuse::setSRMCreditCIRDiffuse");

    randomIndex=_randomIndex;
    crFxCorr=_crFxCorr;
    today=_today;
    probStart=-1; // FIXME probStart : see QMCRatesDiffuse how to calc it
    prob=_prob;
    setRecoveryRate(_srmCreditCIRUtil->getCdsCurve()->getRecovery());

	const int numDaysForNormalIntensity = 1;
    DateTime startDateForNormalIntensity(MaturityPeriod::toDate(matYearsForNormalIntensity, "A", today));
    DateTime endDateForNormalIntensity(MaturityPeriod::toDate(numDaysForNormalIntensity, "D", startDateForNormalIntensity));
    double yearFracForNormalIntensity = SRMYearFrac(startDateForNormalIntensity, endDateForNormalIntensity);
    double logOfSPRatioForNormalIntensity = _srmCreditCIRUtil->logSurvProbRatio(startDateForNormalIntensity, endDateForNormalIntensity);
    if (Maths::isZero(yearFracForNormalIntensity)){
        throw ModelException(method, "Division by zero in calculation of normalIntensity");
    }
    double normalIntensity = - logOfSPRatioForNormalIntensity / yearFracForNormalIntensity;
    if (Maths::isNegative(normalIntensity)) {
        throw ModelException(method, "normalIntensity cannot be negative");
    }

    // sigma
	double vol = _srmCreditCIRUtil->getSpotVol();
    double qLeft=_srmCreditCIRUtil->getQLeft();
    double qRight=_srmCreditCIRUtil->getQRight();

    // check q
    if (!Maths::equals(qLeft,qRight)) { // force qLeft = qRight
        throw ModelException(method, "qLeft and qRight must be equal");
    }
    if (Maths::isNegative(qLeft) || Maths::isPositive(qLeft - 1)) {
        throw ModelException(method, "qLeft cannot be < 0 or > 1");
    }

    double _c = qLeft*vol*vol*normalIntensity;
    double _d = (1.0 - qLeft)*vol*vol*normalIntensity*normalIntensity;
	double _beta = _srmCreditCIRUtil->getBeta();

	setModelParametersExceptTheta(_c,_d,_beta);

    DefaultRatesSP defRates(_srmCreditCIRUtil->getCdsCurve()->defaultRates());
    const DateTimeArray& dates = defRates->getDates();
    vector<double> logProbs(dates.size());
    for (int i=0; i<dates.size(); ++i)
        logProbs[i] = log(defRates->calcDefaultPV(today, dates[i]));
	calibrateTheta(dates,logProbs);
}



void SRMCreditCIRDiffuse::setModelParametersExceptTheta(double _c, double _d, double _beta)
{
    static const string method("SRMCreditCIRDiffuse::setModelParametersExceptTheta");
	c = _c;
	d = _d;
	beta = _beta;

	cIsZero = Maths::isZero(c);
    // set h, hIsZero
    double hSquared = beta*beta + 2.0*c;
    if (Maths::isNegative(hSquared)) {
        throw ModelException(method, "Cannot take square root to calculate h");
    }
    h = sqrt(hSquared);
    hIsZero = Maths::isZero(h);
    // set minIntensity
    minIntensity = cIsZero ? 0.0 : -d/c;
}


void SRMCreditCIRDiffuse::setThetaDates(const DateTimeArray& dates)
{
    static const string method("SRMCreditCIRDiffuse::setThetaDates");

	numThetas = dates.size();
    QLIB_VERIFY((numThetas > 0), "No Thetas were generated");
    theta.resize(numThetas);

    yearFracsForTheta.resize(numThetas);
    yearFracsForTheta = SRMUtil::computeYearFrac(today, dates);

	// set expHYearFracsForTheta
    expHYearFracsForTheta.resize(numThetas);
    for (int i = 0; i < numThetas; i++) 
        expHYearFracsForTheta[i] = exp(h*yearFracsForTheta[i]);
}

void SRMCreditCIRDiffuse::calibrateTheta(const DateTimeArray& dates, const vector<double> &logSurvProb)
{
    static const string method("SRMCreditCIRDiffuse::CalibrateTheta");

	setThetaDates(dates);

    // set initial intensity
    double firstLogSurvProb = logSurvProb[0];

	double firstYearFrac = yearFracsForTheta[0];
    if (Maths::isZero(firstYearFrac)){
        throw ModelException(method, "Division by zero in calculation of initial intensity");
    }
    initialIntensity = - firstLogSurvProb / firstYearFrac;
	// check bound on initialIntensity
	if (Maths::isNegative(c*initialIntensity + d)) {
		throw ModelException(method, "Square-root CIR process cannot start "
			"with this initial intensity");
	}



    // calibrate and store thetas
    double prevYearFrac = 0.0;
    double prevExpHYearFrac = 1.0;
    double prevLogSurvProb = 0.0;
    for (int i = 0; i < numThetas; ++i) {
        // update
        double thisYearFrac = yearFracsForTheta[i];
        double thisExpHYearFrac = expHYearFracsForTheta[i];
        double thisLogSurvProb = logSurvProb[i];

        // calculate and store theta
        double denom = g2(prevYearFrac, thisYearFrac, prevExpHYearFrac, thisExpHYearFrac);
        if (Maths::isZero(denom)) {
            throw ModelException(method, "Division by zero in calibration "
                                         "of mean reversion levels");
        }
        double sum = 0.0;
        double prevYearFracJ = 0.0;
        double prevExpHYearFracJ = 1.0;
        for (int j = 0; j < i; ++j) {
            double thisYearFracJ = yearFracsForTheta[j];
            double thisExpHYearFracJ = expHYearFracsForTheta[j];
            sum += theta[j]*(
                             (
                              g2(prevYearFracJ, thisYearFrac, prevExpHYearFracJ, thisExpHYearFrac)
                             -g2(prevYearFracJ, prevYearFrac, prevExpHYearFracJ, prevExpHYearFrac)
                             )
                            -
                             (
                              g2(thisYearFracJ, thisYearFrac, thisExpHYearFracJ, thisExpHYearFrac)
                             -g2(thisYearFracJ, prevYearFrac, thisExpHYearFracJ, prevExpHYearFrac)
                             )
                            );
            prevYearFracJ = thisYearFracJ;
            prevExpHYearFracJ = thisExpHYearFracJ;
        }
        double deltaLogSurvProb = thisLogSurvProb
                                 -prevLogSurvProb;
        double deltaF = f2(0.0, thisYearFrac, 1.0, thisExpHYearFrac)
                       -f2(0.0, prevYearFrac, 1.0, prevExpHYearFrac);
        double deltaB = b2(0.0, thisYearFrac, 1.0, thisExpHYearFrac)
                       -b2(0.0, prevYearFrac, 1.0, prevExpHYearFrac);
        theta[i] = ( - deltaLogSurvProb
                     - deltaF
                     - initialIntensity*deltaB
                     - sum
                    )
                   / denom;
        // update
        prevYearFrac = thisYearFrac;
        prevExpHYearFrac = thisExpHYearFrac;
        prevLogSurvProb = thisLogSurvProb;
    }
}

void SRMCreditCIRDiffuse::setTheta(const DateTimeArray& dates, double _initialIntensity, const vector<double> &_theta)
{
    static const string method("SRMCreditCIRDiffuse::setTheta");
	QLIB_VERIFY((_initialIntensity > minIntensity), "Initial Intensity too low");
	setThetaDates(dates);
	QLIB_VERIFY(_theta.size() == numThetas, "wrong number of theta given");
	initialIntensity = _initialIntensity;
	theta = _theta;
}

/** populates simulationTheta by looking up in piecewise constant theta */
void SRMCreditCIRDiffuse::storeSimulationTheta() {
    static const string method("SRMCreditCIRDiffuse::storeSimulationTheta");
    const int numDates = sqrtYearFrac.size();
    simulationTheta.resize(numDates);
    int index = 0;
    double currentYearFracFromToday = 0.0;
    for (int i = 0; i < numDates; i++) {
        currentYearFracFromToday += sqrtYearFrac[i]*sqrtYearFrac[i];
        while ( (index < numThetas) && (!Maths::isNegative(currentYearFracFromToday - yearFracsForTheta[index])) ) {
            ++index;
        }
        // index == numThetas || currentYearFracFromToday < yearFracsForTheta[index]
        if (index == numThetas)
            simulationTheta[i] = theta.back(); // use flat extrapolation beyond the end of the CDS curve
        else
            simulationTheta[i] = theta[index];
    }
}

/** populates esdfForwardYearFracs and expHesdfForwardYearFracs */
void SRMCreditCIRDiffuse::storeMeasurementTAndExpHT(const DateTimeArray& dates) {
    static const string method("SRMCreditCIRDiffuse::storeMeasurementTAndExpHT");
    int numDates = dates.size();
    expHesdfForwardYearFracs.resize(numDates);
    esdfForwardYearFracs.resize(numDates);
    for (int i = 0; i < numDates; i++) {
        double yearFrac = SRMYearFrac(today, dates[i]);
        esdfForwardYearFracs[i] = yearFrac;
        expHesdfForwardYearFracs[i] = exp(h*yearFrac);
    }
}

/** returns B(x,y) - see model documentation */
double SRMCreditCIRDiffuse::b2(double x, double y, double expHx, double expHy) {
    static const string method("SRMCreditCIRDiffuse::b2");
    if (hIsZero) {
        double denom = 2 + beta*(y - x);
        if (Maths::isZero(denom)) {
            throw ModelException(method, "Division by zero");
        }
        return 2.0*(y - x)/denom;
    }
    else {
        double denom = (h + beta)*expHy + (h - beta)*expHx;
        if (Maths::isZero(denom)) {
            throw ModelException(method, "Division by zero");
        }
        return 2.0*(expHy - expHx)/denom;
    }
}

    /** core function used in a2 - see model documentation */
double SRMCreditCIRDiffuse::g2(double x, double y, double expHx, double expHy) {
    static const string method("SRMCreditCIRDiffuse::g2");
    if (cIsZero) {
        if (hIsZero) {
            return 0.5*(y - x)*(y - x);
        }
        else {
            if (Maths::isZero(beta)) {
                throw ModelException(method, "Division by zero");
            }
            return -(b2(x, y, expHx, expHy) - (y - x))/beta;
        }
    }
    else
    {
        if (hIsZero) {
            double denom = 2.0 + beta*(y - x);
            if (Maths::isZero(denom)) {
                throw ModelException(method, "Division by zero");
            }
            return -1.0/c*(beta*(y - x) + 2*log(2.0/denom));
        }
        else {
            if (Maths::isZero(expHx)) {
                throw ModelException(method, "Division by zero");
            }
            double denom = 2.0*h + (beta + h)*(expHy/expHx - 1.0);
            if (Maths::isZero(denom)) {
                throw ModelException(method, "Division by zero");
            }
            return -1.0/c*((beta + h)*(y - x) + 2*log(2.0*h/denom));
        }
    }
}

    /** core function used in a2 - see model documentation */
double SRMCreditCIRDiffuse::f2(double x, double y, double expHx, double expHy) {
    static const string method("SRMCreditCIRDiffuse::f2");
    if (cIsZero) {
        if (hIsZero) {
            return -d/6.0*(y - x)*(y - x)*(y - x);
        }
        else {
            double b = b2(x, y, expHx, expHy);
            if (Maths::isZero(beta)) {
                throw ModelException(method, "Division by zero");
            }
            return d/(2.0*beta*beta)*(b - (y - x) + 0.5*beta*b*b);
        }
    }
    else {
        return d/c*(b2(x, y, expHx, expHy) - (y - x) + beta*g2(x, y, expHx, expHy));
    }
}

/** returns A(x,y) - see model documentation */
double SRMCreditCIRDiffuse::a2(double x, double y, double expHx, double expHy) {
    static const string method("SRMCreditCIRDiffuse::a2");
    double a = f2(x, y, expHx, expHy);
    double prevG = g2(x, y, expHx, expHy);
    for (int k = 0; k < numThetas; ++k) {
        double thisG = g2(x < yearFracsForTheta[k] ? yearFracsForTheta[k] : x,
                          y < yearFracsForTheta[k] ? yearFracsForTheta[k] : y,
                          x < yearFracsForTheta[k] ? expHYearFracsForTheta[k] : expHx,
                          y < yearFracsForTheta[k] ? expHYearFracsForTheta[k] : expHy
                          );
        a += theta[k]*(prevG - thisG);
        prevG = thisG;
    }
    return a;
}

/** returns B(x,y) - see model documentation */
double SRMCreditCIRDiffuse::b2frac(double x, double y, double frac, double hfrac) {
	static const string method("SRMCreditCIRDiffuse::b2");
	if (Maths::isZero(hfrac)) {
		double denom = 2.0 + beta*(y - x);
		if (Maths::isZero(denom)) {
			throw ModelException(method, "Division by zero");
		}
		return 2.0*(y - x)/denom;
	}
	else 
	{
		double exph = exp(hfrac*(y-x));
		double denom = (hfrac + beta)*exph + (hfrac - beta);
		if (Maths::isZero(denom)) throw ModelException(method, "Division by zero");
		return 2.0*(exph - 1.0)/denom;
	}
}

/** core function used in a2 - see model documentation */
double SRMCreditCIRDiffuse::g2frac(double x, double y, double frac, double hfrac) {
	static const string method("SRMCreditCIRDiffuse::g2");
	if (cIsZero) {
		if (Maths::isZero(hfrac)) {
			return 0.5*(y - x)*(y - x);
		}
		else {
			if (Maths::isZero(beta)) {
				throw ModelException(method, "Division by zero");
			}
			return -(b2frac(x, y, frac, hfrac) - (y - x))/beta;
		}
	}
	else
	{
		if (Maths::isZero(hfrac)) {
			double denom = 2.0 + beta*(y - x);
			if (Maths::isZero(denom)) {
				throw ModelException(method, "Division by zero");
			}
			return -1.0/(frac*c)*(beta*(y - x) + 2.0*log(2.0/denom));
		}
		else {
			double denom = 2.0*hfrac + (beta + hfrac)*(exp(hfrac*(y-x)) - 1.0);
			if (Maths::isZero(denom)) {
				throw ModelException(method, "Division by zero");
			}
			return -1.0/(c*frac)*((beta + hfrac)*(y - x) + 2.0*log(2.0*hfrac/denom));
		}
	}
}

/** core function used in a2 - see model documentation */
double SRMCreditCIRDiffuse::f2frac(double x, double y, double frac, double hfrac) {
	static const string method("SRMCreditCIRDiffuse::f2");

	if (cIsZero) {
		if (Maths::isZero(hfrac)) {
			return -d*frac*frac/6.0*(y - x)*(y - x)*(y - x);
		}
		else {
			double b = b2frac(x, y, frac, hfrac);
			if (Maths::isZero(beta)) {
				throw ModelException(method, "Division by zero");
			}
			return d*frac*frac/(2.0*beta*beta)*(b - (y - x) + 0.5*beta*b*b);
		}
	}
	else {
		return d*frac/c*(b2frac(x, y, frac, hfrac) - (y - x) + beta*g2frac(x, y, frac, hfrac));
	}
}

/** returns A(x,y) - see model documentation */
double SRMCreditCIRDiffuse::a2frac(double x, double y, double frac, double hfrac) {
	static const string method("SRMCreditCIRDiffuse::a2");
	double a = f2frac(x, y, frac, hfrac);
	double prevG = g2frac(x, y, frac, hfrac);
	for (int k = 0; k < numThetas; ++k) {
		double thisG = g2frac(x < yearFracsForTheta[k] ? yearFracsForTheta[k] : x,
			y < yearFracsForTheta[k] ? yearFracsForTheta[k] : y, 
			frac, 
			hfrac);
		a += frac*theta[k]*(prevG - thisG);
		prevG = thisG;
	}
	return a;
}


/** finalizes the timeline, allocates memory */
void SRMCreditCIRDiffuse::finalizePathGenerator(DateTimeArrayConstSP allDatesSP)
{
    static const string method("SRMCreditCIRDiffuse::finalize");
    try {

        const DateTimeArray& simDates = SRMUtil::calcDiffusionDates(today, *allDatesSP);
        const int numSimDates = simDates.size();
        assert(numSimDates == SRMUtil::getNumSimDates(today, *allDatesSP));

        calcFirstAndLastDiffusionIdx(simDates, *allDatesSP);

        if (fullMC)
        {
            processAllDates(allDatesSP);
            calcRemapToIRAssetIdx(getForwardForwardDates());

            expProb.resize(getNumESDFDates());
        }

        // Query asset what we saw in the generators?
        const DateTime maxDiffusionDate = getDiffusionBound()->getMaxDiffDate();
        const DateTime maxCurveMat      = getDiffusionBound()->getMaxCurveMat();

        if (maxDiffusionDate < today)
            return; // this asset is not really required for anything, even if it was created (example: common markets with weight 0)

        if (isWholeTimelineSurvivalRequested()) // we need this no matter fullMC or fastMC is running
            wholeTimelineLogSurvProb.resize(lastDiffusionIdx-todayIdx+1,0.0);


        if (srmCreditCIRUtil.get() && srmCreditCIRUtil->getMomentMatchingFlag())
        {// info for moment matching
            srmCreditCIRUtil->computeLogFwdProbSimple(getSpotDates() /*sdfRequestedDates*/, originalProbs);
            srmCreditCIRUtil->computeLogFwdProbSimple(getForwardForwardDates() /*esdfForwardDates*/, originalExpProbs);
        }


        sqrtYearFrac = SRMUtil::computeSqrtYearFrac(simDates);

        storeMeasurementTAndExpHT(getForwardForwardDates() /*esdfForwardDates*/);

        storeSimulationTheta();
  
        if (fullMC)
            trimToDiffusion();
        
        if (!fullMC) // we are in the fastMC domain -- need to compute everything right now and clean up some vectors
        {
            computeLogDetermSurvivalProb(getSpotDates(), prob);
            for(size_t i=0; i<prob.size(); ++i)
                prob[i] = exp(-prob[i]);
            computeLogDetermSurvivalProb(getForwardForwardDates(), originalExpProbs);
            for(size_t i=0; i<originalExpProbs.size(); ++i)
                originalExpProbs[i] = -originalExpProbs[i];

            if (isWholeTimelineSurvivalRequested()) 
                computeLogDetermSurvivalProb(simDates, wholeTimelineLogSurvProb);

            // releasing some memory:
            expHesdfForwardYearFracs.clear();
            esdfForwardYearFracs.clear();
            simulationTheta.clear();
        }

    } catch (exception& e) {
        throw ModelException(e, method);
    }
}


void SRMCreditCIRDiffuse::computeLogDetermSurvivalProb(
    const DateTimeArray& dates, 
    vector<double> &logSurvProbs)
{
    logSurvProbs.resize(dates.size());
    for(size_t j=0; j<logSurvProbs.size(); ++j)
    {
        double tJ = SRMYearFrac(today, dates[j]);
        double expHtJ = exp(h*tJ);
        double a = a2(0.0, tJ, 1.0, expHtJ);
        double b = b2(0.0, tJ, 1.0, expHtJ);
        logSurvProbs[j] = a + b*initialIntensity;
    }
}



#define DUMPSIZE(x) cerr << #x << " p=" << (&(x)) << " " << typeid(*this).name() <<  " size= " << (x.size()) << " * " << sizeof(x[0]) <<  endl
/** generates path across all dates */
void SRMCreditCIRDiffuse::generatePathAndSurvivalRates(
                                       IQMCRNGManagerSP rngMgr)
{
    if (!fullMC)
        return; // do nothing

    assert(randomIndex >=0); // truly an assertion here shall fail if evaluates to false
    const double* randoms = rngMgr->getCorrelatedRandoms(randomIndex); // for ease

    if(expProbIndexes.empty()) // means that the asset was not "finalized" -- i.e. was deemed unneeded at all
        return;

    int probDatePos = probStart; // position in probIndexes/prob.
    /* set up position in expProbIndexes/expSvob. If expProbIndexes[0] is 0 then
       we just store 0 for the state variables and start at expDFIndexes[1] */
    // to do: review use of probStart and subtracting it in expDatePos formula
    int expProbDatePos = expProbIndexes.front() == 0? 1: 0;
    int probDateIdx = probIndexes[probDatePos]; // date when we save survival probability
    int expProbDateIdx = expProbIndexes[expProbDatePos]; // date when we save expected survival probability
    int stopIdx = Maths::min(probDateIdx, expProbDateIdx); // when to do something

    double intensity = initialIntensity; // state variable
    double accumulatedIntensity = 0.0; // int_0^t intensity(u) du , ie - ln( survival probability )
    
#if defined(MEMPROF)    
        DUMPSIZE(yearFracsForTheta);
        DUMPSIZE(probIndexes);
        DUMPSIZE(prob);
        DUMPSIZE(originalProbs);
        DUMPSIZE(expProb);
        DUMPSIZE(originalExpProbs);
        DUMPSIZE(esdfForwardYearFracs);
        DUMPSIZE(expHesdfForwardYearFracs);
        DUMPSIZE(expHYearFracsForTheta);
        DUMPSIZE(simulationTheta);
        DUMPSIZE(sqrtYearFrac);
        DUMPSIZE(theta);
        cerr << "UtilSP= " << srmCreditCIRUtil.get() << endl;
        getTimeLogic()->debug();
#endif
    // loop over all dates
    for (int i = 0; i < lastDiffusionIdx - todayIdx /*numDates*/; /* increment in body */) {
        
        double W3 = randoms[i];

        double rootDelT = sqrtYearFrac[i];

        // update accumulatedIntensity
        accumulatedIntensity += intensity*rootDelT*rootDelT;

        assert(c*intensity + d >= 0.0); 
        // if this assert fails - something is corrupted, as flooring below should be taking care of everything

         // update intensity
        double intensityVol = sqrt(c*intensity + d);
        intensity += (simulationTheta[i] - beta*intensity)*rootDelT*rootDelT
                     + intensityVol*rootDelT*W3;
        // add cups component, if necessary
        if (sigmaFX) {
            double fXVol = (*sigmaFX)[i];
            intensity -= crFxCorr*fXVol*intensityVol*rootDelT*rootDelT;
        }

        if (!cIsZero && Maths::isNegative(intensity - minIntensity)) {
            intensity = minIntensity; // floor to ensure existence of square-root
        }

        i++; /* up the loop counter ...
             this comes back to the fact with eg 5 dates and 4 diffusion steps; it's easier to identify the index
             of the date you want rather than the index of the following one */

        if (isWholeTimelineSurvivalRequested())
            wholeTimelineLogSurvProb[i] = accumulatedIntensity;

        if (i + todayIdx == stopIdx) {
            // hit an "event", ie record sthg
            if (stopIdx == probDateIdx) { // save survival probability
	        prob[probDatePos] = accumulatedIntensity; // exponentiate later
                probDatePos++;
                probDateIdx = probIndexes[probDatePos];
            }
            if (stopIdx == expProbDateIdx) {
	        expProb[expProbDatePos].intensity = intensity;
                expProbDatePos++;
                expProbDateIdx = expProbIndexes[expProbDatePos];
            }
            stopIdx = Maths::min(probDateIdx, expProbDateIdx); // refresh
        }
    }
    // then do exp on prob vector
    for (unsigned int j = probStart; j < prob.size(); j++) {
        prob[j] = exp(-prob[j]);
    }
}

/** accesses the expected value ExpSDF(md, fd) where md is a
    simulated measurement date and fd is some future date after the
    measurement is made. */
double SRMCreditCIRDiffuse::getExpectedSurvivalDiscFactor(FwdIdx i, FwdIdx j)
{
    // to avoid one virtual call specify the needed class explicitly
    return exp(SRMCreditCIRDiffuse::getLnExpectedSurvivalDiscFactor(i, j));
}

/** accesses the natural log of the expected value ExpSDF(md, fd)
    where md is a simulated measurement date and fd is some future
    date after the measurement is made. */
double SRMCreditCIRDiffuse::getLnExpectedSurvivalDiscFactor(FwdIdx i, FwdIdx j)
{
    if (!fullMC)
        return getOriginalLnExpectedSurvivalDiscFactor(i,j);
    double tI = esdfForwardYearFracs[i];
    double tJ = esdfForwardYearFracs[j];
    double expHtI = expHesdfForwardYearFracs[i];
    double expHtJ = expHesdfForwardYearFracs[j];
    double a = a2(tI, tJ, expHtI, expHtJ);
    double b = b2(tI, tJ, expHtI, expHtJ);
    int iEDF = getTimeLogic()->getReqEDFIdx(i);
    double intensity = expProb[iEDF].intensity;
    return - a - b*intensity;
}
// To save memory we resize some of the arrays to the size that is actually needed (i.e. to lastIdx)
// The semantics should not change whether this function is called or not.

void SRMCreditCIRDiffuse::trimToDiffusion()
{
    int trimSize = lastDiffusionIdx + 1;
    simulationTheta.resize(trimSize);
    sqrtYearFrac.resize(trimSize);
    
    // remove extra/temp data inside getTimeLogic()
    DateTime maxDiff = getDiffusionBound()->getMaxDiffDate();
    DateTime maxCurveMat = getDiffusionBound()->getMaxCurveMat();
    getTimeLogic()->trim(maxDiff, maxCurveMat); // CR specific: we use aggregated SVs for CR
    getTimeLogic()->cleanup();
}

/** Accessing the expected value E[ SDF(md, fd)^fraction] where md is a
    simulated measurement date and fd is some future date after the
    measurement is made, and a fraction is a power. */
double SRMCreditCIRDiffuse::getExpectedFractionalSurvivalDiscFactor(
									        FwdIdx measurementDateIdx,
											FwdIdx futureDateIdx,
											double fraction)
{
    return exp(getLnExpectedFractionalSurvivalDiscFactor(measurementDateIdx, futureDateIdx, fraction));
}

    /** Accessing the natural log of the expected value E[ SDF(md, fd)^fraction]
        where md is a simulated measurement date and fd is some future
        date after the measurement is made. */
double SRMCreditCIRDiffuse::getLnExpectedFractionalSurvivalDiscFactor(
                                        FwdIdx measurementDateIdx,
                                        FwdIdx futureDateIdx,
                                        double fraction)
{
	/** to compute E(exp(-fraction \int \lambda_u du ) ), we observe that fraction \lambda_u is still a CIR process,
	with parameters beta = beta, c = fraction*c, d = fraction*fraction*d, 
	theta = fraction*theta, lambda0 = fraction*lambda0.*/

    if (!fullMC)
        return getOriginalLnExpectedSurvivalDiscFactor(measurementDateIdx,futureDateIdx)*fraction;

	static const string method("SRMCreditCIRDiffuse::getLnExpectedFractionalSurvivalDiscFactor");

	double tI = esdfForwardYearFracs[measurementDateIdx];
	double tJ = esdfForwardYearFracs[futureDateIdx];
	double hSquared = beta*beta + 2.0*c*fraction;
	QLIB_VERIFY(!Maths::isNegative(hSquared), "Cannot take square root to calculate h_frac");
	double hfrac = sqrt(hSquared);

	double a = a2frac(tI, tJ, fraction, hfrac);
	double b = b2frac(tI, tJ, fraction, hfrac);
	int iEDF = getTimeLogic()->getReqEDFIdx(measurementDateIdx);
	double intensity = expProb[iEDF].intensity*fraction;
	return - a - b*intensity;
}



DRLIB_END_NAMESPACE
