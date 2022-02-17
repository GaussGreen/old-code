//----------------------------------------------------------------------------
//
//   Group       : Convertibles Research
//
//   Filename    : DDEModule.cpp
//
//   Description : util functions for DDE, such as risk eq vol calib, sprd vol creation
//
//   Author      : Qing Hou
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolSurface.hpp"
#include "edginc/DDEModule.hpp"
#include "edginc/AssetDDE.hpp"
#include "edginc/Black.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/Actual365F.hpp"
#include "edginc/FDUtils.hpp"
#include "edginc/FDSolver1FGeneric.hpp"
#include "edginc/CDSHelper.hpp"
#include "edginc/Addin.hpp"

DRLIB_BEGIN_NAMESPACE

const double DDEModule::DEF_IMPVAR_EPSILON = 1.e-5;
const string DDECalib::DEF_TEMPLATE_FREQ    = "3M";     // quarterly
const int    DDECalib::DEF_NB_STEPPERYEAR   = 50;
const int    DDECalib::DEF_NB_STEPPERBENCHMK = 20;
const double DDECalib::DEF_MAX_STEP_SIZE    = 0.5;		// semi-annual
const double DDECalib::DEF_TOLERANCE_SPRD   = 0.00005;
const double DDECalib::DEF_TOLERANCE_VOL    = 0.0005;
const string DDECalib::DEF_MAX_HORIZON      = "30Y";
const string DDECalib::DEF_MIN_HORIZON      = "1M";
const int    DDECalibMC::DEF_NB_PATH        = 10000;
const int    DDECalibFD::DEF_NB_STOCKSTEP   = 200;
const int    DDECalibFD::DEF_NB_MAX_ITER    = 10;
const double DDECalibFD::DEF_TRUNCATION     = 5;
const double DDECalibFD::DEF_CREDIT_STRK    = 4;
const int    DDECalibFD::DEF_NB_SEGMENT_ENDS = 5;


////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  FD1FDDE quick pricer. 
//  Tailored to the FD1F calibration and to quick pricing for DDE output
//
////////////////////////////////////////////////////////////////////////////////////////////////////

class DDEModuleFD1FPricerBase {
public:
	DDEModuleFD1FPricerBase();
	virtual ~DDEModuleFD1FPricerBase();

	// initialize. set or allocate memories
	// subsequent pricing must have date belonging to input "dates" array
	virtual void init(const DDEModule *ddeModule, const DateTimeArray *dates);
	
protected:
    // functions for backward pricing
    void doPrice(const DateTime &date, double *prices, double *deltas=0, double *gammas=0)
    {
		doPrice(getTimeStep(date), prices, deltas, gammas);
    }
    void doPrice(int curTimeStep, double *prices, double *deltas=0, double *gammas=0);
	virtual void payoffAtMat() = 0;

    // function for forward pricing of vanilla call option
    void doPriceVanCallForward(const DateTime &date, int numStrike, const double *strikes, 
            bool firstStepSmooth, bool useDensity, double *prices, double *deltas=0, double *gammas=0, double *zeroBond=0)
    {
		doPriceVanCallForward(getTimeStep(date), numStrike, strikes, firstStepSmooth, useDensity, prices, deltas, gammas, zeroBond);
    }
    void doPriceVanCallForward(int curTimeStep, int numStrike, const double *strikes, 
            bool firstStepSmooth, bool useDensity, double *prices, double *deltas=0, double *gammas=0, double *zeroBond=0);

    // common functions
	void setupPriceArray(int maxNbPrices);
    int getSegIdx(int curTimeStep) const;
    int getSegIdx(double T) const;
    int getTimeStep(const DateTime &date) const;

    virtual void updateStockArray(int nStockStep, const double *stockArray){};
    virtual double getVol(int curTimeStep) = 0;
    virtual void getSpreads(int curTimeStep) = 0;
	virtual bool setupFDTermsOnce() const = 0;

private:
    // function for backward calculation
    void setupFDTermsPerStep(int curTimeStep);

    // function for forward calculation
    // need to get spreads per time step and vol per time step and cache
    void preSetupFDTermsForward(int curTimeStep);
    // need to calc drift term per finer step with in each time step
    void setupFDTermsForward(int curTimeStep, const double *priceArray);

    void setupFDTermsForwardDensity(int curTimeStep);
private:
	void clear();
	void clearPriceArray();

protected:
	int			numPriceArrays;
	int			nTimeSteps;
	int			nStockSteps;

	double	   *Ts;
	double     *stockArray;		// just a pointer to one of the stock array
	double     **priceArray;
    double     *sprdArray;

    DateTimeArray timePts;

private:
	int			settedStep;		// keep track of time step that has been set
    int         nSegments;      // ie. # of stock arrays
    int        *segEndIdx;          // time step index of seg ends

	double		stockNow;
	double	   *dtFDs;
	double	   *irs;
	double	   *divys;

	double     **priceOldArray;
    double     *priceFwdArray;  // for forward calculation

    double     *segDtFDs;       // only used for forward calculation for now
    double     *segEndTs;       // only used for forward calculation for now
	FDSolver1FGeneric *fdSolvers;   // one solver per segment, contain relevant stock array

	DoubleArraySP *driftArrays;		// these provide access to the arrays, the actual mem is managed by FDTerm
	DoubleArraySP *diffusionArrays;
	DoubleArraySP *discountArrays;
	FDTermStructureSP *driftTerms;
	FDTermStructureSP *diffusionTerms;
	FDTermStructureSP *discountTerms;
	FDTermStructureSP couponTerm;
};


DDEModuleFD1FPricerBase::DDEModuleFD1FPricerBase()
{
	numPriceArrays = 0;
	nSegments = 0;
	segEndIdx = 0;
	sprdArray = 0;
    segEndTs = 0;
    segDtFDs = 0;
    fdSolvers = 0;
	priceFwdArray = 0;
	priceArray = 0;
	priceOldArray = 0;
	dtFDs = 0;
	Ts = 0;
	irs = 0;
	divys = 0;
	discountArrays = 0;
	driftArrays = 0;
	diffusionArrays = 0;
	discountTerms = 0;
	driftTerms = 0;
	diffusionTerms = 0;
}

DDEModuleFD1FPricerBase::~DDEModuleFD1FPricerBase()
{
	clear();
}

/** clean up */
void DDEModuleFD1FPricerBase::clear()
{
	if (segEndIdx != 0){
		delete [] segEndIdx;
		segEndIdx = 0;
	}
	if (fdSolvers != 0){
		delete [] fdSolvers;
		fdSolvers = 0;
	}
	nSegments = 0;
	if (dtFDs != 0){
		delete [] dtFDs;
		dtFDs = 0;
	}
	if (segDtFDs != 0){
		delete [] segDtFDs;
		segDtFDs = 0;
	}
	if (segEndTs != 0){
		delete [] segEndTs;
		segEndTs = 0;
	}
	if (Ts != 0){
		delete [] Ts;
		Ts = 0;
	}
	if (irs != 0){
		delete [] irs;
		irs = 0;
	}
	if (divys != 0){
		delete [] divys;
		divys = 0;
	}
	if (discountArrays != 0){
		delete [] discountArrays;
		discountArrays = 0;
	}
	if (driftArrays != 0){
		delete [] driftArrays;
		driftArrays = 0;
	}
	if (diffusionArrays != 0){
		delete [] diffusionArrays;
		diffusionArrays = 0;
	}
	if (discountTerms != 0){
		delete [] discountTerms;
		discountTerms = 0;
	}
	if (driftTerms != 0){
		delete [] driftTerms;
		driftTerms = 0;
	}
	if (diffusionTerms != 0){
		delete [] diffusionTerms;
		diffusionTerms = 0;
	}
	clearPriceArray();
}

/* initialize. set or allocate memories */
void DDEModuleFD1FPricerBase::init(const DDEModule *ddeModule, const DateTimeArray *outDates)
{
	static const string method = "DDEModuleFD1FPricerBase::init";

	try {
		// combine DDEModule time template with requested output dates
		const DateTimeArray &timeTemplate = ddeModule->calibSchd->getDateArray();
		if ( !outDates || outDates->size() == 0 )
			timePts = timeTemplate;
		else
			timePts = DateTime::merge(*outDates, timeTemplate);

		// remove historical dates.
		DateTime valueDate = timeTemplate[0];
		for (vector<DateTime>::iterator iter(timePts.begin()); iter != timePts.end(); /* inc in loop */)
		{
			if ( valueDate.isGreater(*iter) )
				iter = timePts.erase(iter);
			else
				++iter;
		}
		if( timePts.size() == 0 )
			throw ModelException(method, "Output dates are all < valueDate!");
		
		// if DDEModule has DDECalibFD, use FD parameters from this calib
		// otherwise we use default value for FD grid set up
		const AssetDDE *asset = ddeModule->asset;
		const TimeMetric *timeMetric = asset->getTimeMetric().get();
		DDECalibConstSP calibParams = ddeModule->calibParams;
		DDECalibFDConstSP calibFD;
		if( DDECalibFD::TYPE->isInstance(calibParams.get()) )
			calibFD = DDECalibFDConstSP::dynamicCast(calibParams);
		else
        {
            // create backward FD for pricing as default
            DDECalibFDSP calibFDTmp = DDECalibFDSP(new DDECalibFD(false));
            calibFDTmp->validatePop2Object();
            calibFD = calibFDTmp;
        }

		nStockSteps = calibFD->nStockSteps;
		double TruncationStd = calibFD->TruncationStd;
		int gridType = calibFD->gridType;
		MaturityPeriodArraySP segExpiries = calibFD->segExpiries;
        IntArraySP segStepPerYears = calibFD->segStepPerYears;

		// calculate some term structure. remember that 1st point of timePts is valueDate
		int i, nTimeSteps = timePts.size() - 1;
		DoubleArray pvs(nTimeSteps+1), fwds(nTimeSteps+1);
		asset->getPv(valueDate, timePts, pvs);
		asset->getFwd(timePts, fwds);
		stockNow = fwds[0];

		// set up the segments. 
		for(nSegments = 0; nSegments < segExpiries->size(); nSegments++)
		{
			if( timePts.back() <= (*segExpiries)[nSegments]->toDate(valueDate) ) break;
		}
		nSegments++;
        fdSolvers = new FDSolver1FGeneric[nSegments];
		segEndIdx = new int[nSegments];
		segEndTs = new double[nSegments];
        segDtFDs = new double[nSegments];

        sprdArray = new double [nStockSteps+1];
        priceFwdArray = new double [nStockSteps+1];

		// find BS var for the grid set up
		int curSegIdx, k=0;
		LinearStrikeVolRequest volRequest(	stockNow, 
											valueDate,
											timePts.back(),
											false /* not fwd starting */);
		CVolProcessedBSSP volBS(asset->getUnderlyerProcessedVol(&volRequest));
		for(curSegIdx=0; curSegIdx<nSegments; curSegIdx++)
		{
			DateTime segEndDate;
			if( curSegIdx == (nSegments-1) )
			{
				// last segment needs to end at last timePts
				// otherwise, may create sudden P&L when deal maturity cross a segment date
				segEndDate = timePts.back();
			}
			else
			{
				segEndDate = (*segExpiries)[curSegIdx]->toDate(valueDate);
			}
			while( k <= nTimeSteps && timePts[k] <= segEndDate )
				k++;
			segEndIdx[curSegIdx] = k-1; // hold the timeStep that is <= segEndDate
            segEndTs[curSegIdx] = timeMetric->yearFrac(timePts[0], segEndDate);
            segDtFDs[curSegIdx] = 1./(*segStepPerYears)[Maths::min(curSegIdx, segStepPerYears->size()-1)];

			// set up stock grids
			double varMax = volBS->CalcVar(valueDate, segEndDate);
            varMax = TruncationStd*sqrt(varMax);
            // don't let the grid get too extreme
            double DEF_MAX_STDEV = 14;
            if( varMax > DEF_MAX_STDEV ) varMax = DEF_MAX_STDEV;     
			double stockMid = stockNow;
			double stockMax = stockMid*exp(varMax);
			double stockMin = stockMid*exp(-varMax);
		    // allocate memory for the fdSolver and create grid
            if (fdSolvers[curSegIdx].allocateMem(nStockSteps) == FAILURE)
			    throw ModelException(method, "fdSolver.allocateMem failure");
			if( fdSolvers[curSegIdx].createGrid(stockMin,stockMax,nStockSteps,gridType,stockMid) != SUCCESS)
				throw ModelException(method, "fdSolver.createGrid failure");
		}

		// allocate memory for the term structures
		settedStep = 0;
		dtFDs = new double[nTimeSteps];
		Ts = new double[nTimeSteps+1];
        Ts[0] = 0;
		irs	= new double[nTimeSteps];
		divys = new double[nTimeSteps];
		discountArrays	= new DoubleArraySP[nTimeSteps];
		driftArrays		= new DoubleArraySP[nTimeSteps];
		diffusionArrays = new DoubleArraySP[nTimeSteps];
		discountTerms	= new FDTermStructureSP[nTimeSteps];
		driftTerms		= new FDTermStructureSP[nTimeSteps];
		diffusionTerms	= new FDTermStructureSP[nTimeSteps];
		for (i=0; i<nTimeSteps; i++)
		{
			Ts[i+1] = timeMetric->yearFrac(timePts[0], timePts[i+1]);
            double dt = Ts[i+1] - Ts[i];
			// max step size for FD calculation. later disregarded if larger then template time step size
			if(calibParams->useEqualStep )
				dtFDs[i] = 1./calibParams->nStepPerYear;
			else
			{
				dtFDs[i] = dt/calibParams->nStepPerBenchMk;
				if( dtFDs[i] > calibParams->maxStepSize ) 
					dtFDs[i] = calibParams->maxStepSize;
			}

			irs[i] = -log(pvs[i+1]/pvs[i])/dt;
			divys[i] = irs[i] - log(fwds[i+1]/fwds[i])/dt;

			discountArrays[i]	= DoubleArraySP(new DoubleArray(nStockSteps+1));
			driftArrays[i]		= DoubleArraySP(new DoubleArray(nStockSteps+1));
			diffusionArrays[i]	= DoubleArraySP(new DoubleArray(nStockSteps+1));
			discountTerms[i]	= FDTermStructureSP(new FDTermStructure(discountArrays[i]));
			driftTerms[i]		= FDTermStructureSP(new FDTermStructure(driftArrays[i]));
			diffusionTerms[i]	= FDTermStructureSP(new FDTermStructure(diffusionArrays[i]));

			// if needed, setupFDTerms for each step
			if( setupFDTermsOnce() ) 
			{
				curSegIdx = getSegIdx(i+1);
				stockArray = fdSolvers[curSegIdx].getSpotsPtr();
				setupFDTermsPerStep(i+1);
				settedStep = i+1;
			}
		}

		couponTerm = FDTermStructureSP(new FDTermStructure(0.0));

	} catch( exception &e) {
		throw ModelException(e, method, "Failed");
	}
}

// given time step, return segment index
int DDEModuleFD1FPricerBase::getSegIdx(int curTimeStep) const
{
	if( curTimeStep >  segEndIdx[nSegments-1] )
		throw ModelException("DDEModuleFD1FPricerBase::getSegIdx", "Internal error. curTimeStep is larger than last segEndIdx");

	int curSegIdx=0;
	while( curSegIdx < nSegments && segEndIdx[curSegIdx] < curTimeStep )
	{
		curSegIdx++;
	}
	return curSegIdx;
}

// given time, return segment index
int DDEModuleFD1FPricerBase::getSegIdx(double T) const
{
	int curSegIdx=0;
    while( curSegIdx < nSegments && Maths::isNegative(segEndTs[curSegIdx] - T) )
	{
		curSegIdx++;
	}
	if( curSegIdx == nSegments )
		throw ModelException("DDEModuleFD1FPricerBase::getSegIdx", "Internal error. T is larger than last segEndTs");
	
    return curSegIdx;
}

// clear memory for price array
void DDEModuleFD1FPricerBase::clearPriceArray()
{
	int i;
	if (sprdArray != 0){
		delete [] sprdArray;
		sprdArray = 0;
	}
	if (priceFwdArray != 0){
		delete [] priceFwdArray;
		priceFwdArray = 0;
	}
	if (priceArray != 0){
		for (i=0; i<numPriceArrays; i++)
			delete [] priceArray[i];
		delete [] priceArray;
		priceArray = 0;  
	}
	if (priceOldArray != 0){
		for (i=0; i<numPriceArrays; i++)
			delete [] priceOldArray[i];
		delete [] priceOldArray;
		priceOldArray = 0;  
	}
	numPriceArrays = 0;
}

// alloate or reallocate memory for price array
void DDEModuleFD1FPricerBase::setupPriceArray(int nbPrices)
{
	if( numPriceArrays == nbPrices ) return;

	if( numPriceArrays )
		clearPriceArray();

	numPriceArrays = nbPrices;
	priceArray = new double*[numPriceArrays];
	priceOldArray = new double*[numPriceArrays];
	for (int i=0; i<numPriceArrays; i++)
	{
		priceArray[i] = new double[nStockSteps+1];
		priceOldArray[i] = new double [nStockSteps+1];
	}
}

void DDEModuleFD1FPricerBase::setupFDTermsPerStep(int curTimeStep)
{
    // get pricer specific diffusion vol and spreads
    updateStockArray(nStockSteps, stockArray);
    double vol = getVol(curTimeStep);
    getSpreads(curTimeStep);

    double *discounts = &discountArrays[curTimeStep-1]->front();
    double *drifts = &driftArrays[curTimeStep-1]->front();
    double *diffusions = &diffusionArrays[curTimeStep-1]->front();

    for (int j=0; j<=nStockSteps ; ++j) {
        discounts[j] = irs[curTimeStep-1] + sprdArray[j];
        drifts[j] = (discounts[j] - divys[curTimeStep-1]) * stockArray[j];
        diffusions[j] = 0.5 * vol * vol * stockArray[j] * stockArray[j];
    }
}

void DDEModuleFD1FPricerBase::setupFDTermsForwardDensity(int curTimeStep)
{
    double *drifts = &driftArrays[curTimeStep-1]->front();
    for (int j=0; j<=nStockSteps ; ++j) {
        drifts[j] = (divys[curTimeStep-1] -irs[curTimeStep-1] - sprdArray[j]) * stockArray[j];
    }
}

void DDEModuleFD1FPricerBase::preSetupFDTermsForward(int curTimeStep)
{
    // get pricer specific diffusion vol and spreads
    updateStockArray(nStockSteps, stockArray);
    double vol = getVol(curTimeStep);
    getSpreads(curTimeStep);

    DoubleArray &diffusions = *(diffusionArrays[curTimeStep-1]);
    for (int i=0; i<=nStockSteps ; ++i) {
        diffusions[i] = .5 * vol * vol * stockArray[i] * stockArray[i];
    }
}

void DDEModuleFD1FPricerBase::setupFDTermsForward(int curTimeStep, const double *C)
{
    int j, j1, j2;
    double qr = divys[curTimeStep-1] - irs[curTimeStep-1]; 
    DoubleArray &driftArray = *(driftArrays[curTimeStep-1]);

    double sum=0, dCdK, dCdKOld=0, sprdOld=0;
    for(j=nStockSteps; j>=0; j--)
    {
        j1 = (j?(j-1):0);
        j2 = (j<nStockSteps?(j+1):j);
        dCdK = (C[j2]-C[j1])/(stockArray[j2]-stockArray[j1]);
        // if dCdK is zero, dCdK old must be zero as well, use the spread directly
        // otherwise, accumulate sum (integral) and compute the effective spread
        // enforce convexity. otherwise, may amplify fluctuations and failure to price
        if( dCdK > dCdKOld ) dCdK = dCdKOld;
        if( Maths::isZero(dCdK) )
            driftArray[j] = stockArray[j] * qr;
        else
        {
            sum += 0.5*(dCdK - dCdKOld) * (sprdArray[j] + sprdOld);
            driftArray[j] = stockArray[j] * (qr - sum / dCdK);
            dCdKOld = dCdK;
        }
        sprdOld = sprdArray[j];
    }
}


int DDEModuleFD1FPricerBase::getTimeStep(const DateTime &date) const
{
		int curTimeStep;

		// find time step to start pricing
		curTimeStep = -1;
		for (int i=1; i<timePts.size(); i++)
		{
			if( date.equals(timePts[i]) )
			{
				curTimeStep = i;
				break;
			}
		}
		if( curTimeStep == -1 )
			throw ModelException("DDEModuleFD1FPricerBase::doPrice", "Internal error. Required date not in timePts");

        return curTimeStep;
}

void DDEModuleFD1FPricerBase::doPrice(int curTimeStep, double *prices, double *deltas, double *gammas)
{
	static const string method = "DDEModuleFD1FPricerBase::doPrice";

	try{
		int i, j, k;

		// get the stock array
		int curSegIdx = getSegIdx(curTimeStep);
        stockArray = fdSolvers[curSegIdx].getSpotsPtr();

        // we no longer raise error if segment is too big. ideally should try to change segment 
        // within a time step if the step is too big, but lack of time to implement now
        // if( curSegIdx > 0 && segEndIdx[curSegIdx] == segEndIdx[curSegIdx-1] )
        //    throw ModelException(method, "Calibration time step too big. Does not allow time step cover an entire segment for backward pricing");


		// if needed, setupFDTerms for step
		if( !setupFDTermsOnce() ) 
		{
			if( curTimeStep < settedStep || curTimeStep > (settedStep + 1) )
				throw ModelException(method, "Internal error. Can not jump around when setupFDTermsPerStep");
			settedStep = curTimeStep;
			setupFDTermsPerStep(curTimeStep);
		}

		// initialize terminal payoff
        payoffAtMat();

		double upBarrier = -1, upPayout = 0, upPayoutDelta = 0;
		double downBarrier = -1, downPayout = 0, downPayoutDelta = 0;

		// do the loop calculations
		for (i = curTimeStep-1; i>=0; i--){
			
            // FD time step may be smaller than DDE template step
			double dtFD = dtFDs[i];
            double dt = Ts[i+1] - Ts[i];
            int nStepFD = 1 + int(dt/dtFD - 1e-10);
            if( nStepFD < 1 ) nStepFD = 1;
            double ddt = dt - (nStepFD - 1) * dtFD;  // initial dt. stub at end
            for(j=0; j<nStepFD; j++)
            {

			for (k=0; k<numPriceArrays; k++) {
				double *tmpPointer = priceOldArray[k];
				priceOldArray[k] = priceArray[k];
				priceArray[k] = tmpPointer;

				fdSolvers[curSegIdx].forceNonNegative(true);

                if (fdSolvers[curSegIdx].solveBarrier(nStockSteps, driftTerms[i], diffusionTerms[i], discountTerms[i], couponTerm, ddt, 
                                          priceArray[k],priceOldArray[k],upBarrier,upPayout,upPayoutDelta,
                                          upBarrier,upPayout,upPayoutDelta,downBarrier,downPayout,
                                          downPayoutDelta,downBarrier,downPayout,downPayoutDelta,0,0) == FAILURE) 
                {
                    throw ModelException(method, "fdSolver.solveBarrier failure");
                }
            }

            ddt = dtFD;
            } // end of j loop

            // may need to switch stock array
            if (curSegIdx > 0 && i == segEndIdx[curSegIdx-1]) {
                // continue to decrease segment in case time step cross segment
                while( curSegIdx > 0 && i == segEndIdx[curSegIdx-1] )
                    curSegIdx--;

                double *stockArrayOld = stockArray;
                stockArray = fdSolvers[curSegIdx].getSpotsPtr();

                for (k=0; k<numPriceArrays; k++) {
                    double *tmpPointer = priceOldArray[k];
                    priceOldArray[k] = priceArray[k];
                    priceArray[k] = tmpPointer;

                    if (FDCubicSplineInterp(nStockSteps+1,
                                            stockArrayOld,
                                            priceOldArray[k],
                                            nStockSteps+1,
                                            stockArray,
                                            priceArray[k]) == FAILURE) 
                    {
                        throw ModelException(method, "FDCubicSplineInterp failure");
                    }
                }

            } // end of checking segment end
		}

		// interp to get final price
		for (k=0; k<numPriceArrays; k++) {
			double price, delta, gamma;
			if (FDInterpolationD(nStockSteps+1,stockArray,priceArray[k],1,&stockNow,&price,&delta,&gamma) == FAILURE)
                throw ModelException(method, "FDCubicSplineInterpD failure");
			if( prices ) prices[k] = price;
			if( deltas ) deltas[k] = delta;
			if( gammas ) gammas[k] = gamma;
		}

	} catch ( exception &e ) {
		throw ModelException(e, method, "Failed for date " + timePts[curTimeStep].toString() );
	}
}

// forward pricing used for quick pricing of vanilla call and for calibration
void DDEModuleFD1FPricerBase::doPriceVanCallForward(int curTimeStep, int numStrike, const double *strikes,
    bool firstStepSmooth, bool useDensity, double *prices, double *deltas, double *gammas, double *zeroBond)
{
    static const string method = "DDEModuleFD1FPricerBase::doPriceVanCallForward";

    try{
        int i, j, k;

        k=0; // there is only 1 price array for forward calculation

        // get the stock array corresponding to the start of the time step
        int curSegIdx = getSegIdx(Ts[curTimeStep-1]);
        stockArray = fdSolvers[curSegIdx].getSpotsPtr();
        for (i=0; i<numStrike; i++) 
        {
            if( strikes[i] < stockArray[0] || strikes[i] > stockArray[nStockSteps] )
                throw ModelException(method, "Strike " + Format::toString(strikes[i]) + " outside range. Need to increase TruncationStd");
        }

        // if first step, initialize payoff at start of fwd evolution
        // if repeat calc at a time step, set price with cached values at start of current step
        // if start next step, cache price from end of last step. may need to interp first if 
        // last step is a segment end date
        if( curTimeStep == 1 )
        {
            for(i=0; i<=nStockSteps; i++)
            {
                if( useDensity )
                    priceFwdArray[i] = Maths::isZero(stockNow - stockArray[i])?0.5:((stockNow > stockArray[i])?1.0:0.0);
                else
                    priceFwdArray[i] = Maths::max(0.0, 1.0 - stockArray[i]/stockNow);
            }
        } 
        else if ( curTimeStep > settedStep )
        {
            for(i=0; i<=nStockSteps; i++)
                priceFwdArray[i] = priceArray[0][i];
        }
        
        for(i=0; i<=nStockSteps; i++)
            priceArray[0][i] = priceFwdArray[i];

        if( curTimeStep < settedStep || curTimeStep > (settedStep + 1) )
            throw ModelException(method, "Internal error. Can not jump around when setupFDTermsPerStep");
        settedStep = curTimeStep;

        // pre setup to compute price independent terms: diffusion, and spreads
        // no need to set up discount term, since we scale the function by S0 * exp(-qT)
        FDTermStructureSP discountTerm(new FDTermStructure(0.0));
        preSetupFDTermsForward(curTimeStep);
        if( useDensity )
            setupFDTermsForwardDensity(curTimeStep);

        double upBarrier = -1, upPayout = 0, upPayoutDelta = 0;
        double downBarrier = -1, downPayout = 0, downPayoutDelta = 0;

        // do fwd evolution up to this step
        i = curTimeStep-1;

        // FD time step may be smaller than DDE template step
        double T = Ts[i];
        int endSegIdx = getSegIdx(Ts[i+1]);
        if( endSegIdx >= nSegments )
		    throw ModelException(method, "Internal error. Step end time is larger than last segEndTs");

        for( ; curSegIdx <= endSegIdx; curSegIdx++)
        {
            double endT = Maths::min(Ts[i+1], segEndTs[curSegIdx]);

            double dtFD = segDtFDs[curSegIdx];
            double dt = endT - T;
            int nStepFD = 1 + int(dt/dtFD - 1e-10);
            if( nStepFD < 1 ) nStepFD = 1;
            if( Maths::isZero(dt) ) nStepFD = 0; // no need to evolve if step is 0
            double ddt = dt - (nStepFD - 1) * dtFD;  // initial dt. stub at end
            for(j=0; j<nStepFD; j++)
            {
                // calculate drift term which depends on price
                if( !useDensity )
                    setupFDTermsForward(curTimeStep, priceArray[k]);
            
                double *tmpPointer = priceOldArray[k];
                priceOldArray[k] = priceArray[k];
                priceArray[k] = tmpPointer;
            
                fdSolvers[curSegIdx].forceNonNegative(true); 
            
                if (fdSolvers[curSegIdx].solveBarrier(nStockSteps, driftTerms[i], diffusionTerms[i], discountTerm, couponTerm, ddt, 
                    priceArray[k],priceOldArray[k],upBarrier,upPayout,upPayoutDelta,
                    upBarrier,upPayout,upPayoutDelta,downBarrier,downPayout,
                    downPayoutDelta,downBarrier,downPayout,downPayoutDelta,0,0) == FAILURE) 
                {
                    throw ModelException(method, "fdSolver.solveBarrier failure");
                }

                // do smoothing for the first step, otherwise numerical noise can cause wrong greek
                if( firstStepSmooth && curTimeStep == 1 && curSegIdx == 0 && j == 0 )
                {
                    double pv, fwd, vol;
                    pv = exp(- ddt * (irs[0] + sprdArray[nStockSteps/2])); // our grid midpoint is today sprd
                    fwd = stockNow * exp(- ddt * divys[0]) / pv;
                    vol = getVol(curTimeStep) * sqrt(ddt);
                    for(int ii=0; ii<=nStockSteps; ii++)
                    {
                        if( useDensity )
                            priceArray[k][ii] = N1(0.5 * vol + log(fwd/stockArray[ii])/vol);
                        else
                            priceArray[k][ii] = Black::price(true, 1.0, stockArray[ii]/fwd, pv, vol*vol);
                    }
                }

                ddt = dtFD;
        
            } // end of j loop

            T = endT;

            // if beyond segment end, need to switch to next segment
            if( Maths::isPositive(Ts[i+1] - segEndTs[curSegIdx]) )
            {
                double *stockArrayOld = stockArray;
                stockArray = fdSolvers[curSegIdx+1].getSpotsPtr();

                double *tmpPointer = priceOldArray[k];
                priceOldArray[k] = priceArray[k];
                priceArray[k] = tmpPointer;
                
                // spline interpolation results are not reliable outside range of old stock array
                // it may results in negative option price for high strike, or diff from intrinsic for low strike
                // currently ok since we use forceNonNegative(). however, should be better resolved.  
                // setting the end point of all segment to be the same doesn't work for high strike
                if (FDInterpolation(nStockSteps+1,
                    stockArrayOld,
                    priceOldArray[k],
                    nStockSteps+1,
                    stockArray,
                    priceArray[k]) == FAILURE) 
                {
                    throw ModelException(method, "FDCubicSplineInterp failure");
                }                     

                // force function to be positive. if useDensity, force to be < 1
                int ii;
                for(ii=nStockSteps; ii>=0 && priceArray[k][ii]<0; ii--) priceArray[k][ii] = 0;
                if( useDensity )
                    for(ii=0; ii<=nStockSteps && priceArray[k][ii]>1; ii++) priceArray[k][ii] = 1;
                                             
                preSetupFDTermsForward(curTimeStep);
                if( useDensity ) 
                    setupFDTermsForwardDensity(curTimeStep);

            }
        } // end of segment loop

        // calculate zero bond if required
        // also, if useDensity, need to integrate to convert density to option prices
        // use priceOldArray as temp memory to save the density, no need to restore
        double *priceArrayForVanCall = priceArray[k];
        if( useDensity )
        {
            priceArrayForVanCall = priceOldArray[k];
            
            double sum=0, D, DOld=0, kInv, kInvOld=0;
            for(j=nStockSteps; j>=0; j--)
            {
                D = priceArray[k][j];
                kInv = 1./stockArray[j];
                sum += 0.5*(D - DOld) * (kInv + kInvOld);
                priceArrayForVanCall[j] = (D - stockArray[j] * sum);
                DOld = D;
                kInvOld = kInv;
            }
            sum += (1 - DOld) * kInvOld;

            if( zeroBond )
                zeroBond[0] = sum;
        }
        else
            if( zeroBond )
                zeroBond[0] = (priceArray[k][0]-priceArray[k][1])/(stockArray[1]-stockArray[0]);

        // interp to get final price
        for (i=0; i<numStrike; i++) {
            double price, delta, gamma, strike = strikes[i];
            if (FDInterpolationD(nStockSteps+1,stockArray,priceArrayForVanCall,1,&strike,&price,&delta,&gamma) == FAILURE)
                throw ModelException(method, "FDCubicSplineInterpD failure");
            if( prices ) prices[i] = price;
            if( deltas ) deltas[i] = delta;
            if( gammas ) gammas[i] = gamma;
        }

    } catch ( exception &e ) {
        throw ModelException(e, method, "Failed for date " + timePts[curTimeStep].toString() );
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  FD1FDDE quick pricer.
//  geared towards DDE sensitivity output. all FD terms are calculated in one go. 
//  later can price different payoffs using the same grid
//
////////////////////////////////////////////////////////////////////////////////////////////////////

class DDEModuleFD1FPricer : public DDEModuleFD1FPricerBase {
public:
	DDEModuleFD1FPricer() : csFunc(0){};
	~DDEModuleFD1FPricer(){};

	void init(const DDEModule *ddeModule, const DateTimeArray *outDates);

	// price zero bonds with maturity at a particular time step
	void priceZero(const DateTime &date, double &price, double &delta)
	{
		isZero = true;
		doPrice(date, &price, &delta);
	}

	// price european vanilla. input strikes
	void priceOpt(const DateTime &date, double strk, double &price, double &delta)
	{
		isZero = false;
		strike = strk;
		doPrice(date, &price, &delta);
	}

private:
	bool setupFDTermsOnce() const { return true; }
	double getVol(int curTimeStep);
	void getSpreads(int curTimeStep);
	void payoffAtMat();

private:
	bool isZero;
	double strike;
	const SpreadEquityFunc *csFunc;
	CVolProcessedBSSP volBS;
};


void DDEModuleFD1FPricer::init(const DDEModule *ddeModule, const DateTimeArray *outDates)
{
	// get credit spread function
	csFunc = ddeModule->getSpreadFunc();

	// get LN vols. this simple calc if fine for now since risky vol currently is a curve, not surface
	LinearStrikeVolRequest volRequest(	ddeModule->asset->getSpot(),
										ddeModule->valueDate,
										ddeModule->lastDate,
										false /* not fwd starting */);
    volBS = CVolProcessedBSSP::dynamicCast(
							CVolProcessedSP((ddeModule->getProcessedVol(&volRequest))));

	DDEModuleFD1FPricerBase::init(ddeModule, outDates);
	setupPriceArray(1);
}

double DDEModuleFD1FPricer::getVol(int curTimeStep)
{
	if( !volBS )
		throw ModelException("DDEModuleFD1FPricer::getVol", "Internal error. csFunc and/or volBS not set yet");

    return volBS->CalcVol(timePts[curTimeStep-1], timePts[curTimeStep]);
}

void DDEModuleFD1FPricer::getSpreads(int curTimeStep)
{
	if( !csFunc )
		throw ModelException("DDEModuleFD1FPricer::getSpreads", "Internal error. csFunc and/or volBS not set yet");
	csFunc->getSpreadCC(timePts[curTimeStep-1], timePts[curTimeStep], nStockSteps + 1, stockArray, sprdArray, Ts[curTimeStep]-Ts[curTimeStep-1]);
}

void DDEModuleFD1FPricer::payoffAtMat()
{
	// initialize terminal zero bond and vanilla call price payoff
	for (int i=0; i<=nStockSteps ; ++i) 
		priceArray[0][i] = isZero?1:Maths::max(0.0, stockArray[i] - strike);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  DDEModule base class functions
//
////////////////////////////////////////////////////////////////////////////////////////////////////

DDEModule::DDEModule()
: CObject(TYPE), isOneFactor(true), lastGoodCalibDate(DateTime(0,0))
{}

DDEModule::DDEModule(CClassConstSP clazz, bool isOneFactor) 
: CObject(clazz), isOneFactor(isOneFactor), lastGoodCalibDate(DateTime(0,0))
{}

void DDEModule::update(const AssetDDE *asset, const IDDEInitiator *initiator)
{
	this->initiator = initiator;
	this->asset = asset;
	smartCalib();
}

void DDEModule::smartCalib()
{
	bool calibCredit = asset->isCreditTweaked();

	// if there is no time template, must be initial state, build it
    // if value date differ from 1st template date, must be theta tweaked, rebuild
	if( !calibSchd || calibSchd->length() == 0 || 
        calibSchd->firstDate() != valueDate )
	{
		buildTemplate();
		calibCredit = true;
	}

	// create risky vol once. sensitivity tweaks are applied to it later.
	if( riskyVolOW() && !riskyEqVols )
		riskyEqVols = asset->getRiskyVol(valueDate);

	// if clean sprd is null, must be initial calibation
	if( !cleanSprds )
		calibCredit = true;
		
	if( calibCredit || asset->isEquityTweaked() )
	{
		if( calibCredit ) 
			preprocessCredit();
		else
			restoreCreditOnly();
	
		calib();

		// if clean sprd backup is null, must be initial calibation, store initial calib results
		if( !cleanSprdsBak ) store();
	}
	else
		restore();
}

void DDEModule::preprocessCredit()
{
	const DDEParams *ddeParams = asset->getDDEParams();
	createCleanSpreadCurve();
	createNonEqCleanSpreadCurve();
	cleanSprdVols = createCleanSpreadVolCurve(cleanSprds.get(), ddeParams, asset->getTimeMetric().get());
	corrs = createEquitySpreadCorrCurve(cleanSprds.get(), ddeParams);

	// if use related vol, adjust vol and corr
	if( ddeParams->isRelatedVol )
	{
		DoubleArray &vols = *cleanSprdVols->vols;
		DoubleArray &corrArrs = *corrs->corrs;
		if( vols.size() != corrArrs.size() )
			throw(ModelException("DDEModule::preprocessCredit", "Internal error. Sprd vol and corr dimension unequal"));
		
		int i;
		for(i=0; i<vols.size(); i++)
		{
			vols[i] *= corrArrs[i];
			corrArrs[i] = 1;
		}
	}

	// for 2 factor, need to calibrate the spread backbone
	if( !isOneFactor )
	{
		bbCleanSprds = calibLNCleanSpreadCurve(cleanSprds.get(), cleanSprdVols.get(), calibParams.get());
	}
}


// store the initial calibrated results
void DDEModule::store()
{
	if( !!cleanSprds )		
		cleanSprdsBak = CleanSpreadCurveSP(copy(cleanSprds.get()));
	if( !!nonEqCleanSprds )		
		nonEqCleanSprdsBak = CleanSpreadCurveSP(copy(nonEqCleanSprds.get()));
	if( !!bbCleanSprds )	
		bbCleanSprdsBak = CleanSpreadCurveSP(copy(bbCleanSprds.get()));
	if( !!cleanSprdVols ) 
		cleanSprdVolsBak = CleanSpreadVolCurveSP(copy(cleanSprdVols.get()));
	if( !!corrs )
		corrsBak = EquitySpreadCorrCurveSP(copy(corrs.get()));
	if( !!riskyEqVols )
		riskyEqVolsBak = CVolBaseSP(copy(riskyEqVols.get()));
	if( !!spreadFunc )
		spreadFuncBak = SpreadEquityFuncSP(copy(spreadFunc.get()));
}

// restore the initial calibrated results
void DDEModule::restore()
{
	restoreCreditOnly();

	if( !!riskyEqVolsBak )
		riskyEqVols = CVolBaseSP(copy(riskyEqVolsBak.get()));
	if( !!spreadFuncBak )
		spreadFunc = SpreadEquityFuncSP(copy(spreadFuncBak.get()));
}

// no need to restore equity
void DDEModule::restoreCreditOnly()
{
	if( !!cleanSprdsBak )		
		cleanSprds = CleanSpreadCurveSP(copy(cleanSprdsBak.get()));
	if( !!nonEqCleanSprdsBak )		
		nonEqCleanSprds = CleanSpreadCurveSP(copy(nonEqCleanSprdsBak.get()));
	if( !!bbCleanSprdsBak )	
		bbCleanSprds = CleanSpreadCurveSP(copy(bbCleanSprdsBak.get()));
	if( !!cleanSprdVolsBak ) 
		cleanSprdVols = CleanSpreadVolCurveSP(copy(cleanSprdVolsBak.get()));
	if( !!corrsBak )
		corrs = EquitySpreadCorrCurveSP(copy(corrsBak.get()));
}

bool DDEModule::riskyVolOW() const
{
	return (calibParams->useRiskyVolOW && asset->hasRiskyVol());
}

void DDEModule::buildTemplate()
{
	static const string method = "DDEModule::buildTemplate";
	try {

	// ask for calib time template
	valueDate = asset->getEquity()->getValueDate();
	lastDate = initiator->maxMaturity();
    if( lastDate < calibParams->minHorizon->toDate(valueDate) ) 
        lastDate = calibParams->minHorizon->toDate(valueDate);
    if( lastDate > calibParams->maxHorizon->toDate(valueDate) ) 
        lastDate = calibParams->maxHorizon->toDate(valueDate);
    DateTime valueDate0 = (!calibSchd || calibSchd->length()==0 )?valueDate:calibSchd->firstDate();
	DateTimeArraySP timeTemplate = calibParams->getTemplateDates(valueDate0, lastDate);
	DateTimeArray &tempDates = *timeTemplate;

	// insert critical dates if supplied. also remove dumplicate days w.r.t. trading time
	DateTimeArray calibDates;
	initiator->sensitiveDates(calibDates);
	if( calibDates.size() )
	{
		tempDates = DateTime::merge(tempDates, calibDates);
	}

    const TimeMetric *timeMetric = asset->getTimeMetric().get();
    DateTime prevDate = valueDate;
    vector<DateTime>::iterator iter(tempDates.begin());
	while( iter != tempDates.end() )
	{
		if (valueDate.isGreaterOrEqual(*iter) || lastDate.isLess(*iter) ||
            Maths::isZero(timeMetric->volDays(prevDate, *iter)) )
			iter = tempDates.erase(iter);
		else
		{
			prevDate = *iter;
			++iter;
		}
	}
    // insert value date into 1st position
    tempDates.insert(tempDates.begin(), valueDate);

	// set up the strike template using flat piecewise term structure
	DoubleArray tempStrks(tempDates.size());
	if( !calibParams->useAtm )
		initiator->sensitiveStrikes(tempDates, tempStrks, calibStrkIsPct);
	else
	{
		// if useAtmFwd to calibrate, overwrite the product provided strkIsPct flag
		calibStrkIsPct = true;
		// fill up atm strikes
		for(int j=0; j<tempStrks.size(); j++)
			tempStrks[j] = 1;
	}

	// put time/strkTemplate into the schedule
	calibSchd = ScheduleSP(new Schedule(tempDates, tempStrks, Schedule::INTERP_STAIRS));

	} catch (exception &e) {
         throw ModelException(e, method, "Failed");
	}
}

double DDEModule::getCalibStrike(const DateTime &date) const
{
	double strike = calibSchd->interpolate(date);
	if( calibStrkIsPct )
		strike *= asset->fwdValue(calibParams->atmIsFwd?date:valueDate);

    // if ccy struck, we want to factor in the spotFX 
    double fxFactor = asset->getFXSpot();

    // if support next strike, floor/ceiling strike if it's out of surface
    const CVolBase *obj = asset->getVol().get();
    const IAllStrikes* allStrikes = dynamic_cast<const IAllStrikes*>(obj);
    if (allStrikes && calibParams->strkOffSurfaceAdj) {
        DoubleArraySP allK(allStrikes->getAllStrikes());
        int numK = allK->size();

        if (strike < (*allK)[0] * fxFactor) {
            strike = (*allK)[0] * fxFactor;
        }
        else if (strike > (*allK)[numK-1] * fxFactor) {
            strike = (*allK)[numK-1] * fxFactor;
        }
        else {
            for (int i = 0; i < numK && (*allK)[i] * fxFactor >= strike; i++) {
                strike = (*allK)[i] * fxFactor;
            }
        }
    }

    return strike;
}

bool DDEModule::sensShift(Theta *shift)
{
    valueDate = shift->rollDate(valueDate);
    // calibSchd is rebuild later during calib
	return true;		// some member may need to be tweaked. eg. riskyEqVol overwrite
}

// default implementation. precompute risky zero price and vanilla prices, then calib to these prices
// specific DDEModule calls their own detailed calib() function
// create spread function if is 1 factor, or spread backbone if 2 factor
void DDEModule::calib()
{
	static const string method = "DDEModule::calib";

	try {
		string name = asset->getVolName();
		const TimeMetric *timeMetric = asset->getTimeMetric().get();
		const DateTimeArray &dates = calibSchd->getDateArray();

		int i, nbDates = dates.size(), nbStep = nbDates-1;
		DoubleArray	strike0(1), t(nbDates), ndps(nbDates), calls(nbDates);
		DoubleArray strikes(nbDates), pvs(nbDates), fwds(nbDates), bsVars(nbDates), riskyVars(nbDates);
		ExpiryArray expiries(nbStep);
		DoubleMatrix vols(1, nbStep);

		// get times in yrs, pv/fwd/strk, as well as ndp/call to calibrate to
		asset->getPv(valueDate, dates, pvs);
		asset->getFwd(dates, fwds);

		// compute gross var to match BlackScholes price, then strip credit component
		double strike = getCalibStrike(valueDate);
		strike0[0] = strike; // this strike is not important since vol surface has only 1 strike
		LinearStrikeVolRequest volRequest(	strike, 
											valueDate,
											lastDate,
											false /* not fwd starting */);
		CVolProcessedBSSP volBS = CVolProcessedBSSP(asset->getUnderlyerProcessedVol(&volRequest));
		
        // risky vol overwrite. we know that volBS now contains the risky vol overwrite
        if( riskyVolOW() )
        {
            CVolProcessedBSSP volBSRisky(dynamic_cast<CVolProcessedBS*>(getProcessedVol(&volRequest))); // use riskyEqVol to get processed vol
            for(i=0; i<nbDates; i++)
            {
                t[i] = timeMetric->yearFrac(valueDate, dates[i]);
                ndps[i] = cleanSprds->getDefaultPV(valueDate, dates[i]);
                bsVars[i] = volBS->CalcVar(valueDate, dates[i]);
                riskyVars[i] = volBSRisky->CalcVar(valueDate, dates[i]);
            }

            // specific calibration (RG, LN1F, LN2F etc)
            if( isOneFactor )
            {
                DoubleArray sprdVars(nbDates);
                cleanSprdVols->getSpotSprdVars(dates, sprdVars); // find spot sprd var to calibrate to

                calib1F(dates, t, pvs, ndps, sprdVars, fwds, strikes, bsVars, calls, riskyVars); // strike, vsVar, calls not used

                spreadFunc->setTimeMetric(timeMetric);
            }
            else
                calib2F(dates, t, pvs, ndps, fwds, strikes, bsVars, calls, riskyVars);

            lastGoodCalibDate = lastDate; // we know credit always calibrate by itself
            return;
        } // end if risky vol overwrite portion
        			
		for(i=0; i<nbDates; i++)
		{
			t[i] = timeMetric->yearFrac(valueDate, dates[i]);
			if( i>0 )
				expiries[i-1] = ExpirySP(MaturityPeriod::dateSubtract(dates[i], valueDate)); // create relative date diff for expiry
			ndps[i] = cleanSprds->getDefaultPV(valueDate, dates[i]);
			strikes[i] = getCalibStrike(dates[i]);
			if( !Maths::isZero(strikes[i] - strike ) )	
			{
				strike = strikes[i];
				LinearStrikeVolRequest volRequest(	strike, 
													valueDate,
													lastDate,
													false /* not fwd starting */);
				volBS = CVolProcessedBSSP(asset->getUnderlyerProcessedVol(&volRequest));
			}

			// interpolate the vol using our LN request and calc variance
			bsVars[i] = volBS->CalcVar(valueDate, dates[i]);
			calls[i] = Black::price(true /* isCall */, fwds[i], strike, pvs[i], bsVars[i]);
			riskyVars[i] = 0; // initialize
		}
			
		// make sure option is higher than intrinsic value
		DateTimeArray calibDates = dates;
        for(i = 1; i<nbDates; i++)
        {
            double intrinsic = Maths::max(0.0, (fwds[i] - strikes[i] * ndps[i]) * pvs[i]); 
            if( !Maths::isZero(calls[i]) && !Maths::isPositive(calls[i]-intrinsic) ) break;
        }

		if( i<nbDates )
		{
			string msg = "Call price " + Format::toString(calls[i]) 
				+ " less than risky intrinsic Max(0.0, " + Format::toString((fwds[i] - strikes[i] * ndps[i]) * pvs[i]) 
				+ ") on " + dates[i].toString() + " with strike " + Format::toString(strikes[i]) + " and ndp " + Format::toString(ndps[i]);

			if( calibParams->throwIfFail || i==1 )
				throw ModelException(method, msg);

			// output warning for initial calib only.
			if( !riskyEqVols )
			{
				ModelException modelErr(method + ":WARNING:" + msg);
				modelErr.errorLog();
			}
			// keep expiries and vol dimension since want to maintain 0 fwd vol
			calibDates.resize(i);
		}

		// specific calibration (RG, LN1F, LN2F etc)
		if( isOneFactor )
		{
			DoubleArray sprdVars(nbDates);
			cleanSprdVols->getSpotSprdVars(dates, sprdVars); // find spot sprd var to calibrate to

			calib1F(calibDates, t, pvs, ndps, sprdVars, fwds, strikes, bsVars, calls, riskyVars);

			spreadFunc->setTimeMetric(timeMetric);
		}
		else
			calib2F(calibDates, t, pvs, ndps, fwds, strikes, bsVars, calls, riskyVars);
		
		lastGoodCalibDate = dates[0];
		for(i=1; i<nbDates; i++)
		{

			// make sure forward variance is positive
			if( riskyVars[i] <= riskyVars[i-1] )
			{
				if( !riskyEqVols && i<calibDates.size() ) 
				{
					// output warning message for the initial calibration only. avoid excess output from sensitivity runs, esp implied sprd/vol runs
					string msg = "Negative fwd variance for date " + dates[i].toString() + " at strike " + Format::toString(strikes[i]);
					if( calibParams->throwIfFail ) throw ModelException(method, msg);
					ModelException modelErr(method + ":WARNING:" + msg);
					modelErr.errorLog();
				}
		
				for(; i<nbDates; i++)
					vols[0][i-1] = (i==1?0:((1+1e-5)*vols[0][i-2]*sqrt(t[i-1]/t[i])));
				break;
			}

			// convert var to vol
			vols[0][i-1] = sqrt(riskyVars[i]/t[i]);
			lastGoodCalibDate = dates[i];
		}			
		
		riskyEqVols = CVolBaseSP(new VolSurface(name, timeMetric, strike0, vols, &expiries, valueDate));

	}
	catch (exception& e) {
        throw ModelException(e, method, "Failed to calibrate risky eq vol and spread function");
    }

}

void DDEModule::calib1F(const DateTimeArray &dates, const DoubleArray &t, const DoubleArray &pvs, const DoubleArray &ndps, const DoubleArray &sprdVars,
			const DoubleArray &fwds, const DoubleArray &strikes, const DoubleArray &bsVars, const DoubleArray &calls, DoubleArray &riskyVars)
{
	throw ModelException("DDEModule::calib1F", "Not implemented for base DDEModule class");
}

void DDEModule::calib2F(const DateTimeArray &dates, const DoubleArray &t, const DoubleArray &pvs, const DoubleArray &ndps,
			const DoubleArray &fwds, const DoubleArray &strikes, const DoubleArray &bsVars, const DoubleArray &calls, DoubleArray &riskyVars)
{
	throw ModelException("DDEModule::calib1F", "Not implemented for base DDEModule class");
}


void DDEModule::createCleanSpreadCurve()
{
	try {

	cleanSprds = CleanSpreadCurveSP(copy(asset->getCleanSpreadCurve(valueDate).get()));
	DoubleArray &rates = *cleanSprds->cleanSpreads;

	int i;
	for( i=0; i<rates.size(); i++ ) rates[i] *= asset->riskyness;
	
	} catch ( exception &e ) {
		throw ModelException(e, "DDEModule::createCleanSpreadCurve", "Failed to create clean spread curve");
	}
}

void DDEModule::createNonEqCleanSpreadCurve()
{
	try {

	nonEqCleanSprds = CleanSpreadCurveSP(copy(asset->getCleanSpreadCurve(valueDate).get()));
	DoubleArray &rates = *nonEqCleanSprds->cleanSpreads;

	int i;
	for( i=0; i<rates.size(); i++ ) rates[i] *= (1-asset->riskyness);

	} catch ( exception &e ) {
		throw ModelException(e, "DDEModule::createNonEqCleanSpreadCurve", "Failed to create non equity related clean spread curve");
	}
}

CleanSpreadVolCurveSP DDEModule::createCleanSpreadVolCurve(const CleanSpreadCurve *cleanSprds, const DDEParams *ddeParams, const TimeMetric *timeMetric)
{
	try {

	CleanSpreadVolCurveSP volCurve(new CleanSpreadVolCurve());
	DoubleArray &rates = *cleanSprds->cleanSpreads;
	DoubleArray vols(rates.size());

	// determine mean reversion using the average 2Y spread as representative spread
	int horizonDays = (int)(DDEParams::DEF_E2C_HORIZON * 365.25);
	double sprdForMR = cleanSprds->getCleanSpread(cleanSprds->valueDate.rollDate(horizonDays));
	double mr = ddeParams->getMeanReversion(sprdForMR);

	int i;
	for( i=0; i<rates.size(); i++ ) vols[i] =  ddeParams->spreadVol(rates[i]);
	return CleanSpreadVolCurveSP(new CleanSpreadVolCurve(cleanSprds->valueDate, 
														cleanSprds->expiries.get(),
														&vols,
														mr,
														timeMetric));
	} catch ( exception &e ) {
		throw ModelException(e, "DDEModule::createCleanSpreadVolCurve", "Failed to create clean spread vol curve");
	}
}

EquitySpreadCorrCurveSP DDEModule::createEquitySpreadCorrCurve(const CleanSpreadCurve *cleanSprds, const DDEParams *ddeParams)
{
	CleanSpreadVolCurveSP volCurve(new CleanSpreadVolCurve());
	DoubleArray &rates = *cleanSprds->cleanSpreads;
	DoubleArray corrs(rates.size());

	int i;
	for( i=0; i<rates.size(); i++ ) corrs[i] =  ddeParams->spreadEqCorr(rates[i]);
	return EquitySpreadCorrCurveSP(new EquitySpreadCorrCurve(cleanSprds->valueDate, 
														cleanSprds->expiries.get(),
														&corrs));
}

class DDECalibLNSprdData
{
public:
	static double calibFuncLNSprd(double bbSprd, void *calibDataMem);
	
	DDECalibLNSprdData(double target, const DoubleArray *zFac, const DoubleArray *pvs) 
		: target(target), zFac(zFac), pvs(pvs){};

	double target;
	const DoubleArray *zFac;
	const DoubleArray *pvs;
};


// for root finder for spread backbone
double DDECalibLNSprdData::calibFuncLNSprd(double bbSprd, void *calibDataMem)
{
	const DoubleArray &zFac = *((DDECalibLNSprdData *)calibDataMem)->zFac;
	const DoubleArray &pvs = *((DDECalibLNSprdData *)calibDataMem)->pvs;
	double pv=0, target = ((DDECalibLNSprdData *)calibDataMem)->target;
	int i, nbPath = zFac.size();

	for(i=0; i<nbPath; i++)
	{
		pv += pvs[i] * exp(- bbSprd * zFac[i] );
	}
	return target - pv/nbPath;
}


/** period wise calibration of lognormal clean spread backbone to match risky discount factors
 *  use Monte Carlo due to lack of precision with analytical approx */
CleanSpreadCurveSP DDEModule::calibLNCleanSpreadCurve(
				const CleanSpreadCurve		*cleanSprds, 
				const CleanSpreadVolCurve	*volCurve,
				const DDECalib				*calibParams)
{
	static const string method = "DDEModule::calibLNCleanSpreadCurve";

	try{

		DoubleArray rates0(cleanSprds->cleanSpreads->size());
		const DateTimeArray &dates = cleanSprds->getDates();
		const DoubleArray &vols = *volCurve->vols;
		double mr = volCurve->meanReversion;
		const TimeMetric *timeMetric = volCurve->timeMetric.get();

		DateTime valueDate = cleanSprds->valueDate;

		if( !DDECalibMC::TYPE->isInstance(calibParams) )
			throw(ModelException(method, "Currently only support MonteCarlo calibration"));

		const DDECalibMC *calibMC = dynamic_cast<const DDECalibMC *>(calibParams);

		int i, n, nbPath = calibMC->nbPath;
		double bbSprd = 0;
		DoubleArray z(nbPath, 0), zFac(nbPath, 0), pv(nbPath, 1);

		// get all random numbers in one go
		DoubleArray rand(nbPath, 0);
		IRandomSP randGenerator(copy(calibMC->rand.get()));
		randGenerator->init();

		double t = 0;
		Actual365F  dcc;		// dcc for clean spread and mean reverting behavior
		for( i=0; i<rates0.size(); i++ )
		{
			double target = cleanSprds->getDefaultPV(dates[i]);

			randGenerator->fetch(rand.size(), &rand[0]);
	
			// initial guess from analytical approx, increment z, 
			// then root finding which will also update pvs
			double dt = dcc.yearFraction(valueDate, dates[i]) - t;
			double dtVol = timeMetric->yearFrac(i==0?valueDate:dates[i-1], dates[i]);
			double fwdSprd = (1 - cleanSprds->getDefaultPV(i==0?valueDate:dates[i-1], dates[i]))/dt;
			double expFac = exp(mr * t);
			double volExp = vols[i] * expFac * sqrt(dtVol);
			for( n=0; n<nbPath; n++)
			{
				if( i > 0 ) pv[n] *= exp(- bbSprd * zFac[n] );
				z[n] += rand[n + i*nbPath] * volExp;
				zFac[n] = dt * exp( z[n] / expFac );
			}

			// root searching for d factor, ie. asset/debt ratio, to calibrate spread
			DDECalibLNSprdData calibData(target, &zFac, &pv);

			try {
				bbSprd = zbrentUseful(
						&DDECalibLNSprdData::calibFuncLNSprd,	/* (I) function to find the root of			*/
						&calibData,								/* (I) only 1 other parameter to function	*/
						0.1 * fwdSprd,							/* (I) low value for sprd backbone			*/ 
						fwdSprd,								/* (I) high value for sprd backbone			*/
						calibParams->tolSprd);					/* (I) tolerance. 0.05bp					*/
			} catch (exception& e) {
				// give warning message and flat extend calibrated portion
				bbSprd = (i==0?fwdSprd:rates0[i-1]);
				ModelException modelErr(e, method + " failed for date " + dates[i].toString());
				modelErr.errorLog();
			}

			rates0[i] = bbSprd;

			t += dt;

		}
		return CleanSpreadCurveSP(new CleanSpreadCurve(	cleanSprds->valueDate, 
														cleanSprds->expiries.get(),
														&rates0));
	} catch(exception& e) {
        throw ModelException(e, method, "Failed to calibrate.");

	}
}


const CleanSpreadCurve *DDEModule::getCleanSpreadCurve() const
{
	return cleanSprds.get();
}

const CleanSpreadCurve *DDEModule::getNonEqSpreads() const
{
	return nonEqCleanSprds.get();
}

const SpreadEquityFunc *DDEModule::getSpreadFunc() const
{
	if( !isOneFactor )
		throw ModelException("DDEModule::getSpreadCurves should only be called by 1 factor DDE");

	return spreadFunc.get();
}

void DDEModule::getSpreadCurves(const CleanSpreadCurve *&bbSprds,
								const CleanSpreadVolCurve *&vols,
								const EquitySpreadCorrCurve *&eqSprdCorrs) const
{
	if( isOneFactor )
		throw ModelException("DDEModule::getSpreadCurves should only be called by 2 factor DDE");

	bbSprds = bbCleanSprds.get();
	vols = cleanSprdVols.get();
	eqSprdCorrs = corrs.get();
}

CVolProcessed *DDEModule::getProcessedVol(const CVolRequest* volRequest) const
{
	return riskyEqVols->getProcessedVol(volRequest, asset);
}


//////////////////////////////////////////////////////////////////////////////
//																			//
//	routines to calculate extra output to visualize DDE model				//
//  including class of OutputRequest that can pass in extra info			//						
//																			//
//////////////////////////////////////////////////////////////////////////////

class OutputRequestTS : public OutputRequest
{
public:
    static CClassConstSP const TYPE;

	//// creates deep copy
    IObject* clone() const
	{
		OutputRequestTS* copy = new OutputRequestTS(this->getRequestName());
		if( !!outputDates )
			copy->outputDates = DateTimeArraySP(outputDates.clone());
		return copy;
	}

    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(OutputRequestTS, clazz);
        SUPERCLASS(OutputRequest);
        EMPTY_SHELL_METHOD(defaultOutputRequestTS);
        FIELD(outputDates,				"dates to calc requested output");
		FIELD_MAKE_OPTIONAL(outputDates);
    }

	static IObject* defaultOutputRequestTS(){
		return new OutputRequestTS(OutputRequest::DDE_SPRD_DELTA);
	}

protected:
    OutputRequestTS(const string &name) : OutputRequest(TYPE, name){} 
    OutputRequestTS(const CClassConstSP &clazz, const string &name) : OutputRequest(clazz, name){} 

public:
	DateTimeArraySP	outputDates;
};

CClassConstSP const OutputRequestTS::TYPE = CClass::registerClassLoadMethod(
              "OutputRequestTS", typeid(OutputRequestTS), load);


// local class for DDE delta output
class DDEOutSens : public CObject
{
public:
    static CClassConstSP const TYPE;

    DDEOutSens(const DateTimeArray& dates, 
        const DoubleArray&  levels, 
        const DoubleArray&  senses) : CObject(TYPE), 
                                    dates(new DateTimeArray(dates)), 
                                    levels(new DoubleArray(levels)), 
                                    senses(new DoubleArray(senses)){} 

    static void load(CClassSP& clazz){
        clazz->setPrivate(); // make visible to EAS/spreadsheet
        REGISTER(DDEOutSens, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultDDEOutSens);
        FIELD(dates,                "dates");
        FIELD_MAKE_OPTIONAL(dates);
        FIELD(levels,               "levels");
        FIELD_MAKE_OPTIONAL(levels);
        FIELD(senses,               "senses");
        FIELD_MAKE_OPTIONAL(senses);
    }

    static IObject* defaultDDEOutSens(){
        return new DDEOutSens();
    }

private:
    DDEOutSens() : CObject(TYPE){} 

public:
    DateTimeArrayConstSP    dates;
    DoubleArrayConstSP      levels;
    DoubleArrayConstSP      senses;
};

CClassConstSP const DDEOutSens::TYPE = CClass::registerClassLoadMethod(
              "DDEOutSens", typeid(DDEOutSens), load);

typedef smartPtr<DDEOutSens>      DDEOutSensSP;

// we use weight to cache the base level of sprd/vol to which the delta is applied
DDEOutSensSP DDEModule::calcDDEVisual(const OutputRequest *request) const
{
	static const string method = "DDEModule::calcDDEVisual";

	try {

	// use outDates if provided, or use the timeTemplate. remove any historical 
	const OutputRequestTS *requestLoc = dynamic_cast<const OutputRequestTS *>(request);
	const TimeMetric *timeMetric = asset->getTimeMetric().get();
	DateTimeArray outDates;
	if( requestLoc && !!requestLoc->outputDates && requestLoc->outputDates->size() )
		outDates = *requestLoc->outputDates;
	else
		outDates = calibSchd->getDates();

	for (vector<DateTime>::iterator iter(outDates.begin()); iter != outDates.end(); /* inc in loop */)
	{
		if ( *iter <= valueDate )
			iter = outDates.erase(iter);
		else
			++iter;
	}

	DoubleArray outValues(outDates.size());
	DoubleArray baseValues(outDates.size());
	DoubleArray pvs(outDates.size());
	asset->getPv(valueDate, outDates, pvs);
	double price, delta, t;

	DDEModuleFD1FPricer pricer;
	pricer.init(this, &outDates);
	string name = request->getRequestName();
	if( name == OutputRequest::DDE_SPRD_DELTA )
	{
		DoubleArray ndps(outDates.size());
		for(int i=0; i<outDates.size(); i++)
		{
			pricer.priceZero(outDates[i], price, delta);
			t = timeMetric->yearFrac(valueDate, outDates[i]);
			baseValues[i] = -log(price/pvs[i])/t;
			outValues[i] = -log(1+delta/price)/t;
		}
	}
	else if ( name == OutputRequest::DDE_VOL_DELTA || name == OutputRequest::DDE_ATMVOL_DELTA ||
			  name == OutputRequest::DDE_VOL_SKEW || name == OutputRequest::DDE_ATMVOL_SKEW )
	{
		bool isAtm = ((name == OutputRequest::DDE_ATMVOL_DELTA) || (name == OutputRequest::DDE_ATMVOL_SKEW));
		bool isSkew = ((name == OutputRequest::DDE_VOL_SKEW) || (name == OutputRequest::DDE_ATMVOL_SKEW));
		double spot = asset->getSpot(), strike, bsVar, varUp, varDown;
		double dx = 0.005 * spot; // tweak size
		DoubleArray fwds(outDates.size());
		asset->getFwd(outDates, fwds);
		for(int i=0; i<outDates.size(); i++)
		{
            try {

			strike = isAtm?(calibParams->atmIsFwd?fwds[i]:spot):getCalibStrike(outDates[i]);
			t = timeMetric->yearFrac(valueDate, outDates[i]);
			// get LN vols as guess for impliedVar calculation
			LinearStrikeVolRequest volRequest( strike,
												valueDate,
												outDates[i],
												false /* not fwd starting */);
			CVolProcessedBSSP volBS = CVolProcessedBSSP::dynamicCast(
									CVolProcessedSP(asset->getUnderlyerProcessedVol(&volRequest)));
			bsVar = volBS->CalcVar(valueDate, outDates[i]);
			
			if( isSkew )
			{	
				// get vol diff from different strike
				pricer.priceOpt(outDates[i], strike + dx, price, delta);
				varUp = Black::impliedVariance(true, fwds[i], strike + dx, pvs[i], bsVar, price, bsVar*DEF_IMPVAR_EPSILON);
				pricer.priceOpt(outDates[i], strike - dx, price, delta);
				varDown = Black::impliedVariance(true, fwds[i], strike - dx, pvs[i], bsVar, price, bsVar*DEF_IMPVAR_EPSILON);
				// get the calculated option vol
				pricer.priceOpt(outDates[i], strike, price, delta);
				bsVar = Black::impliedVariance(true, fwds[i], strike, pvs[i], bsVar, price, bsVar*DEF_IMPVAR_EPSILON);
			}
			else
			{
				// get vol diff from different initial spot
				pricer.priceOpt(outDates[i], strike, price, delta);
				varUp = Black::impliedVariance(true, fwds[i]*(1 + dx/spot), strike, pvs[i], bsVar, price + delta * dx, bsVar*DEF_IMPVAR_EPSILON);
				varDown = Black::impliedVariance(true, fwds[i]*(1 - dx/spot), strike, pvs[i], bsVar, price - delta * dx, bsVar*DEF_IMPVAR_EPSILON);
				// get the calculated option vol
				bsVar = Black::impliedVariance(true, fwds[i], strike, pvs[i], bsVar, price, bsVar*DEF_IMPVAR_EPSILON);
			}

			baseValues[i] = sqrt(bsVar/t);
			outValues[i] = 0.5*(sqrt(varUp/t) - sqrt(varDown/t))/dx;
           
            } catch ( exception ) {
                baseValues[i] = -1.;    // set warning value. don't fail.
                outValues[i] = -1.;
            }
		}
	}
	else
		throw ModelException(method, "Internal error. Request not recognized");

	return DDEOutSensSP(new DDEOutSens(outDates, baseValues, outValues));

	} catch ( exception &e ) {
		throw ModelException (e, method, "Failed for request " + request->getRequestName());
	}
}

// output request for implied vol output
class OutputRequestImpVol : public OutputRequestTS
{
public:
    static CClassConstSP const TYPE;

    //// creates deep copy
    IObject* clone() const
    {
        OutputRequestImpVol* copy = new OutputRequestImpVol();
        if( !!outputDates )
            copy->outputDates = DateTimeArraySP(outputDates.clone());
        if( !!outputStrikes )
            copy->outputStrikes = DoubleArraySP(outputStrikes.clone());
        return copy;
    }

    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(OutputRequestImpVol, clazz);
        SUPERCLASS(OutputRequestTS);
        EMPTY_SHELL_METHOD(defaultOutputRequestImpVol);
        FIELD(outputStrikes,    "strikes to calc requested output");
        FIELD_MAKE_OPTIONAL(outputStrikes);
    }

    static IObject* defaultOutputRequestImpVol(){
        return new OutputRequestImpVol();
    }

private:
    OutputRequestImpVol() : OutputRequestTS(TYPE, OutputRequest::DDE_IMPLIED_VOL){} 

public:
    DoubleArraySP   outputStrikes;
};

CClassConstSP const OutputRequestImpVol::TYPE = CClass::registerClassLoadMethod(
              "OutputRequestImpVol", typeid(OutputRequestImpVol), load);

smartPtr<CVolBase> DDEModule::calcDDEImpVol(const OutputRequest *request) const
{
    static const string method = "DDEModule::calcDDEImpVol";

    try {

        const OutputRequestTS *requestLoc = dynamic_cast<const OutputRequestTS *>(request);
        DateTimeArray outDates;
        if( requestLoc && !!requestLoc->outputDates && requestLoc->outputDates->size() )
            outDates = *requestLoc->outputDates;
        else
            outDates = calibSchd->getDates();

        for (vector<DateTime>::iterator iter(outDates.begin()); iter != outDates.end(); /* inc in loop */)
        {
            if ( *iter <= valueDate )
                iter = outDates.erase(iter);
            else
                ++iter;
        }

        const OutputRequestImpVol *requestLoc2 = dynamic_cast<const OutputRequestImpVol *>(request);
        DoubleArray outStrikes;
        if( requestLoc2 && !!requestLoc2->outputStrikes && requestLoc2->outputStrikes->size() )
            outStrikes = *requestLoc2->outputStrikes;
        else 
        {
            const IAllStrikes* allStrikes = dynamic_cast<const IAllStrikes*>(asset->getVol().get());
            if (allStrikes) {
                outStrikes = *allStrikes->getAllStrikes();

                if (outStrikes.empty()) {
                    throw ModelException(method, "Found no strikes");
                }
            }
            else {
                throw ModelException(method, "Strikes input is required if "
                                     "request is not OutputRequestImpVol and "
                                     "asset vol is not VolSurface");
            }
        }

        DoubleMatrix vols(outStrikes.size(), outDates.size());
        DoubleArray fwds(outDates.size()), pvs(outDates.size());
        asset->getPv(valueDate, outDates, pvs);
        asset->getFwd(outDates, fwds);
        const TimeMetric *timeMetric = asset->getTimeMetric().get();
        int i, j;
        double price, delta, t, bsVar;

        DDEModuleFD1FPricer pricer;
        pricer.init(this, &outDates);

        for(j=0; j<outStrikes.size(); j++)
        for(i=0; i<outDates.size(); i++)
        {

            t = timeMetric->yearFrac(valueDate, outDates[i]);
            // get LN vols as guess for impliedVar calculation
            LinearStrikeVolRequest volRequest( outStrikes[j],
                                                valueDate,
                                                outDates[i],
                                                false /* not fwd starting */);
            CVolProcessedBSSP volBS = CVolProcessedBSSP::dynamicCast(
                                    CVolProcessedSP(asset->getUnderlyerProcessedVol(&volRequest)));
            bsVar = volBS->CalcVar(valueDate, outDates[i]);
            
            pricer.priceOpt(outDates[i], outStrikes[j], price, delta);
            try {
                bsVar = Black::impliedVariance(true, fwds[i], outStrikes[j], pvs[i], bsVar, price, bsVar*DEF_IMPVAR_EPSILON);
                vols[j][i] = sqrt(bsVar/t);
            } catch ( exception ) {
                vols[j][i] = -1.; // don't fail
            }

        } // end of loop i
    
        string name = asset->getVolName();
        ExpiryArray expiries(outDates.size());
        for(i=0; i<outDates.size(); i++)
            expiries[i] = ExpirySP(MaturityPeriod::dateSubtract(outDates[i], valueDate));
        return CVolBaseSP(new VolSurface(name, timeMetric, outStrikes, vols, &expiries, valueDate));

    } catch ( exception &e ) {
        throw ModelException (e, method);
    }
}


void DDEModule::recordOutputRequests(Control* control,
									Results* results) const
{
	if( control && control->isPricing() )
	{
		OutputRequest* request;
		request = NULL;
		if (control->requestsOutput(OutputRequest::DDE_EQ_VOL, request) && !!riskyEqVols)
		{
			CVolBaseSP ddeEqVol(copy(riskyEqVols.get()));
			results->storeRequestResult(request, ddeEqVol, OutputNameSP(new OutputName(asset->getTrueName(), OutputRequest::DDE_EQ_VOL)));
		}
		request = NULL;
		if (control->requestsOutput(OutputRequest::DDE_SPRD_BBONE, request) && !!bbCleanSprds)
		{
			CleanSpreadCurveSP bbSprd(copy(bbCleanSprds.get()));
			results->storeRequestResult(request, bbSprd, OutputNameSP(new OutputName(asset->getTrueName(), OutputRequest::DDE_SPRD_BBONE)));
		}
		request = NULL;
		if (control->requestsOutput(OutputRequest::DDE_SPRD_VOL, request) && !!cleanSprdVols)
		{
			CleanSpreadVolCurveSP sprdVol(copy(cleanSprdVols.get()));
			results->storeRequestResult(request, sprdVol, OutputNameSP(new OutputName(asset->getTrueName(), OutputRequest::DDE_SPRD_VOL)));
		}
		request = NULL;
		if (control->requestsOutput(OutputRequest::DDE_EQ_SPRD_COR, request) && !!corrs)
		{
			EquitySpreadCorrCurveSP corr(copy(corrs.get()));
			results->storeRequestResult(request, corr, OutputNameSP(new OutputName(asset->getTrueName(), OutputRequest::DDE_EQ_SPRD_COR)));
		}
		request = NULL;
		if (control->requestsOutput(OutputRequest::DDE_SPRD_FUNC, request) && !!spreadFunc)
		{
			SpreadEquityFuncSP sprdFunc(copy(spreadFunc.get()));
			results->storeRequestResult(request, sprdFunc, OutputNameSP(new OutputName(asset->getTrueName(), OutputRequest::DDE_SPRD_FUNC)));
		}
		request = NULL;
		if (control->requestsOutput(OutputRequest::DDE_LAST_GOOD_CALIB_DATE, request) )
		{
			DateTimeSP lastCalibDt(copy(&lastGoodCalibDate));
			results->storeRequestResult(request, lastCalibDt, OutputNameSP(new OutputName(asset->getTrueName(), OutputRequest::DDE_LAST_GOOD_CALIB_DATE)));
		}
		request = NULL;
		if (control->requestsOutput(OutputRequest::DDE_CALIB_SCHEDULE, request) )
		{
			DateTimeArray calibDts = calibSchd->getDates();
			calibDts.erase(calibDts.begin());
			DoubleArray calibStrks(calibDts.size());
			for(int i=0; i<calibDts.size(); i++)
			{
				calibStrks[i] = getCalibStrike(calibDts[i]);
			}
			ScheduleSP calibSchdOut(new Schedule(calibDts, calibStrks, Schedule::INTERP_STAIRS));
			results->storeRequestResult(request, calibSchdOut, 
				OutputNameSP(new OutputName(asset->getTrueName(), request->getRequestName())));
		}
		request = NULL;
		if (control->requestsOutput(OutputRequest::DDE_SPRD_DELTA, request) )
		{
			DDEOutSensSP visual = this->calcDDEVisual(request);
			results->storeRequestResult(request, visual, 
				OutputNameSP(new OutputName(asset->getTrueName(), request->getRequestName())));
		}
		if (control->requestsOutput(OutputRequest::DDE_VOL_DELTA, request) )
		{
			DDEOutSensSP visual = this->calcDDEVisual(request);
			results->storeRequestResult(request, visual, 
				OutputNameSP(new OutputName(asset->getTrueName(), request->getRequestName())));
		}
		if (control->requestsOutput(OutputRequest::DDE_ATMVOL_DELTA, request) )
		{
			DDEOutSensSP visual = this->calcDDEVisual(request);
			results->storeRequestResult(request, visual, 
				OutputNameSP(new OutputName(asset->getTrueName(), request->getRequestName())));
		}
		if (control->requestsOutput(OutputRequest::DDE_VOL_SKEW, request) )
		{
			DDEOutSensSP visual = this->calcDDEVisual(request);
			results->storeRequestResult(request, visual, 
				OutputNameSP(new OutputName(asset->getTrueName(), request->getRequestName())));
		}
		if (control->requestsOutput(OutputRequest::DDE_ATMVOL_SKEW, request) )
		{
			DDEOutSensSP visual = this->calcDDEVisual(request);
			results->storeRequestResult(request, visual, 
				OutputNameSP(new OutputName(asset->getTrueName(), request->getRequestName())));
		}
		if (control->requestsOutput(OutputRequest::DDE_IMPLIED_VOL, request) )
		{
			CVolBaseSP impVol = this->calcDDEImpVol(request);
			results->storeRequestResult(request, impVol, 
				OutputNameSP(new OutputName(asset->getTrueName(), request->getRequestName())));
		}
	}
}

class DDEModuleHelper{
public:
    /** Invoked when Class is 'loaded' */

    static void load(CClassSP& clazz){ 
		clazz->setPrivate(); // make invisible to EAS/spreadsheet
		REGISTER(DDEModule, clazz);
		SUPERCLASS(CObject);
		EMPTY_SHELL_METHOD(defaultDDEModule);
		IMPLEMENTS(Theta::IShift);
		FIELD(isOneFactor,	"");
		FIELD(valueDate,		"");
		FIELD(lastDate,		"");
		FIELD(calibSchd, 	"");
		FIELD_MAKE_TRANSIENT(calibSchd); 
		FIELD(cleanSprds,		"");
		FIELD_MAKE_TRANSIENT(cleanSprds); 
		FIELD(nonEqCleanSprds,	"");
		FIELD_MAKE_TRANSIENT(nonEqCleanSprds); 
		FIELD(bbCleanSprds, 	"");
		FIELD_MAKE_TRANSIENT(bbCleanSprds); 
		FIELD(cleanSprdVols,	"");
		FIELD_MAKE_TRANSIENT(cleanSprdVols); 
		FIELD(corrs,			"");
		FIELD_MAKE_TRANSIENT(corrs); 
		FIELD(riskyEqVols,		"");
		FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(riskyEqVols); // tweakable to allow greek calculation with risky vol overwrite
		FIELD(spreadFunc,		"");
		FIELD_MAKE_TRANSIENT(spreadFunc); 
		FIELD(calibStrkIsPct,	"");
		FIELD_MAKE_TRANSIENT(calibStrkIsPct); 
		FIELD(calibParams,			"");
		FIELD_MAKE_TRANSIENT(calibParams); 
	}
	static IObject* defaultDDEModule(){
		return new DDEModule();
	}
};

CClassConstSP const DDEModule::TYPE = CClass::registerClassLoadMethod(
              "DDEModule", typeid(DDEModule), DDEModuleHelper::load);


//////////////////////////////////////////////////////////////////////////////
//																			//
//	Spread Equity Function implementation for LN							//
//																			//
//////////////////////////////////////////////////////////////////////////////


class SpreadEquityFuncLN: public SpreadEquityFunc
{
public:
// not done yet
    static CClassConstSP const TYPE;
	friend class DDEModuleLN;

	void getSpreadCC(const DateTime &start,
					const DateTime &end,
					int		nbSpot, 
					const double *spots,
					double *spreads,
					double dt=0) const;

	static void load(CClassSP& clazz);

	static IObject * defaultSpreadEquityFuncLN()
	{
		return new SpreadEquityFuncLN();
	}

private:
	DateTimeArraySP	dates;		// first date should be basedate/valuedate
	CDoubleArraySP	factor;		// backbone factor
	CDoubleArraySP	power;		// power, ie. ratio of sprd vol / eq vol
	CDoubleArraySP	dts;		// time interval between dates.

	SpreadEquityFuncLN() : SpreadEquityFunc(TYPE){};
};


void SpreadEquityFuncLN::getSpreadCC(
						const DateTime &start,
						const DateTime &end,
						int	nbSpot, 
						const double *spots,
						double *spreads,
						double dt) const
{
	// sanity check
	if( end <= start ) 
		throw(ModelException("SpreadEquityFuncLN::getSpreadCC", "Start date >= end date"));
	if( Maths::isZero(dt) ) 
		dt = metric->yearFrac(start, end);	
	if( Maths::isZero(dt) ) // now if dt is still zero. need to throw error
		throw ModelException("SpreadEquityFuncLN::getSpreadCC. Internal error. Time fraction is 0 between start " 
							+ start.toString() + " and end " + end.toString());

	// interpolate the spread/eq backbone and vol
	DateTimeArray &dates = *this->dates;
	DoubleArray &factor	= *this->factor;
	DoubleArray &power = *this->power;

	if( dates[0] > start )
		throw ModelException("SpreadEquityFuncLN::getSpreadCC. start date is before the value date of spread function");

	//internal sanity check
	if( dates.size() < 2 )
		throw ModelException("SpreadEquityFuncLN::getSpreadCC. Internal error. Dimension must be >= 2");
	if( dates.size() != (factor.size() + 1) || dates.size() != (power.size() + 1) )
		throw ModelException("SpreadEquityFuncLN::getSpreadCC. Internal error. Date dimension must be factor/power dimension + 1");

	int i=0, i1=-1, i2=-1;
	while( i < dates.size() && dates[i] <= start )
		i++;

	i2 = i;
	i1 = i2 - 1;

	// reset spreads first
	for( i=0; i<nbSpot; i++) spreads[i] = 0;
	
	while( i2 < dates.size() && dates[i1] <= end )
	{
		// calculate spreads for the spots
		double fac = factor[i1];
		DateTime t1 = (dates[i1]>start?dates[i1]:start);
		DateTime t2 = (dates[i2]<end?dates[i2]:end);
		fac *= metric->yearFrac(t1, t2);

		for( i=0; i<nbSpot; i++)
		{
			spreads[i] += fac * ::pow(spots[i], -power[i1]);
		}

		i1++;
		i2++;
	}

	// use last available info if end is beyond spread function
	if( dates[dates.size()-1] < end )
	{
		i2 = dates.size()-1;
		i1 = i2-1;
		double fac = factor[i1];
		DateTime t1 = (dates[i2]>start?dates[i2]:start);
		DateTime t2 = end;
		fac *= metric->yearFrac(t1, t2);

		for( i=0; i<nbSpot; i++)
		{
			spreads[i] += fac * ::pow(spots[i], -power[i1]);
		}

	}

	for( i=0; i<nbSpot; i++)
	{
		spreads[i] /= dt;
	}
}

void SpreadEquityFuncLN::load(CClassSP& clazz){ 
		clazz->setPrivate(); // make invisible to EAS/spreadsheet
		REGISTER(SpreadEquityFuncLN, clazz);
		SUPERCLASS(SpreadEquityFunc);
		EMPTY_SHELL_METHOD(SpreadEquityFuncLN::defaultSpreadEquityFuncLN);
		FIELD(dates,	"");
		FIELD(factor,	"");
		FIELD(power,	"");
		FIELD(dts,		"");
}


CClassConstSP const SpreadEquityFuncLN::TYPE = CClass::registerClassLoadMethod(
              "SpreadEquityFuncLN", typeid(SpreadEquityFuncLN), SpreadEquityFuncLN::load);



////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  FD1FDDECalib used for DDEModuleLN
//  calibration routines needed to reproduce zero bond and vanilla option prices
//
////////////////////////////////////////////////////////////////////////////////////////////////////

class FD1FDDECalib : public DDEModuleFD1FPricerBase {
public:
	static const double DEF_UNIT_TOL_X;
	static const double DEF_UNIT_TOL_F;

	FD1FDDECalib();
	~FD1FDDECalib();

	// initialize. set or allocate memories 
	void init(const DDEModule *ddeModule, const DoubleArray &t, const DoubleArray &pvs, const DoubleArray &ndps, const DoubleArray &sprdVars,
		const DoubleArray &fwds, const DoubleArray &strikes, const DoubleArray &bsVars, const DoubleArray &calls);

	// calibration for a time step
	void calibViaSimplePrice(int curTimeStep,  double curPower, const double *guess, double *result);

	static void newtPriceWrap(double *x,int n,double *fvec,double **fjac)
	{
		curObj->priceNSens(x, n, fvec, fjac);
	}

private:
	bool setupFDTermsOnce() const { return false; }
	double getVol(int curTimeStep);
	void getSpreads(int curTimeStep);
	void payoffAtMat();
    void updateStockArray(int nStockStep, const double *stockArray);

private:

	// calc price and sensitivity/jacobi
	void priceNSens(double *x,int n,double *fvec,double **fjac);

	// calc price alone
	void simplePrice(int n, double *x, double *fvec);

private:
    const DDECalibFD *calibFD;

	int			curTimeStep; // internal member. the idx of step being calibrated
	int			nbJacobi; 

	const double	*t, *pvs, *ndps, *sprdVars, *fwds, *strikes, *bsVars, *calls;

	double		target[2], errxDivisor[2], errfDivisor[2];
    double      xBound[2];
	const double *curX;  // current value of roots
    double      curPower;
    double	   *powSpotArray;

	bool		isRiskyVolOW;
	double		dRiskyVolOW;

	static FD1FDDECalib *curObj;
};

const double FD1FDDECalib::DEF_UNIT_TOL_X = 0.01;
const double FD1FDDECalib::DEF_UNIT_TOL_F = 0.01;
FD1FDDECalib *FD1FDDECalib::curObj = 0;

FD1FDDECalib::FD1FDDECalib()
{
		nbJacobi = 0;
		powSpotArray = 0;
}

FD1FDDECalib::~FD1FDDECalib()
{
		if (powSpotArray != 0){
			delete [] powSpotArray;
			powSpotArray = 0;
		}
}

/* initialize. set or allocate memories */
void FD1FDDECalib::init(const DDEModule *ddeModule, const DoubleArray &t, const DoubleArray &pvs, const DoubleArray &ndps, const DoubleArray &sprdVars,
	const DoubleArray &fwds, const DoubleArray &strikes, const DoubleArray &bsVars, const DoubleArray &calls)
{
		static const string method = "FD1FDDECalib::init";

		isRiskyVolOW = ddeModule->riskyVolOW();

		DDEModuleFD1FPricerBase::init(ddeModule, 0); // no extra dates
		setupPriceArray(isRiskyVolOW?1:2);	// price zero and option

		calibFD = dynamic_cast<const DDECalibFD *>(ddeModule->calibParams.get());
		if( calibFD == 0 )
			throw ModelException(method, "Interal error. Can not call FD1FDDECalib if calibParams is not type DDECalibFD");

		this->t			= &t[0];
		this->pvs		= &pvs[0];
		this->ndps		= &ndps[0];
		this->sprdVars	= &sprdVars[0];
		this->fwds		= &fwds[0];
		this->strikes	= &strikes[0];
		this->bsVars	= &bsVars[0];
		this->calls		= &calls[0];

		powSpotArray = new double [nStockSteps+1];
}

// called to calibration for each time step
void FD1FDDECalib::calibViaSimplePrice(int curTimeStep, double curPower, const double *guess, double *result)
{
		static const string method = "FD1FDDECalib::calibViaSimplePrice";

		int i, ntrial;
		double errx, errf;

		this->curTimeStep = curTimeStep;
        this->curPower = curPower;

		// find target and error divisor for zero bond and vanilla call
        // need to control the min time span for accuracy measure in case of tiny span
        const double minDt = 1./12.; // one month
		i = curTimeStep;
		target[0] = pvs[i] * ndps[i];
		errxDivisor[0] = (calibFD->tolSprd/DEF_UNIT_TOL_X) * ::pow(fwds[i], curPower);
        errfDivisor[0] = (calibFD->tolSprd/DEF_UNIT_TOL_F) * Maths::max(minDt, t[i]-t[i-1]) * ndps[i] * pvs[i];	// sprd sensitivity
		if( errfDivisor[0] < 1e-6 ) errfDivisor[0] = 1e-6;
		if( isRiskyVolOW )
			dRiskyVolOW = guess[1];
		else
		{
			target[1] = calls[i];
			errxDivisor[1] = (calibFD->tolEqVol/DEF_UNIT_TOL_X);
			errfDivisor[1] = (calibFD->tolEqVol/DEF_UNIT_TOL_F) * Black::vega(true /*iscall*/, fwds[i], strikes[i], pvs[i], Maths::max(minDt, t[i]), sqrt(bsVars[i]/t[i])); // vega sensitivity
            errfDivisor[1] *= guess[1]*Maths::min(t[i], Maths::max(minDt, t[i]-t[i-1]))/sqrt(bsVars[i]*t[i]);// adjust for forward period
            if( errfDivisor[1] < 1e-6 ) errfDivisor[1] = 1e-6;
		}
        // initialize xBounds
        xBound[0] = -1;
        xBound[1] = -1;

		// use Newton to find root. Notice that numerical recipe array index starts from 1
		curObj = this;
		for(i=0; i<numPriceArrays; i++) result[i] = guess[i]/errxDivisor[i];
		try {
			nbJacobi = 0;
			mnewt(calibFD->maxIter, result-1, numPriceArrays, DEF_UNIT_TOL_X, DEF_UNIT_TOL_F, 
				&ntrial, &errx, &errf, newtPriceWrap);
		} catch ( exception &e ) {
			throw ModelException(e, "FD1FDDECalib::calibViaSimplePrice", "Failed in call to mnewt");
		}

		const double RESOLUTION = 0.5 * DEF_UNIT_TOL_X;
		for(i=0; i<numPriceArrays; i++){
			if( Maths::isNegative(result[i]) )
			{
				if( Maths::isNegative(fabs(result[i])-RESOLUTION) )
					result[i] = 0;
				else
					throw ModelException("FD1FDDECalib::calibViaSimplePrice", "Calibrated " + Format::toString(i?"vol ":"spread factor") 
							+ Format::toString(result[i]*errxDivisor[i]) + " is less than 0");
			}
			result[i] *= errxDivisor[i];
		}
}

void FD1FDDECalib::priceNSens(double *x,int n,double *fvec,double **fjac)
{
    // make sure the parameters are positive or 0 
    // otherwise, set to ~0 and set the other input to the optimal bound value
    int j;
    for(j=1; j<=n; j++)
        if( Maths::isNegative(x[j]) )
        {
            x[j] = DBL_EPSILON;
            if( !isRiskyVolOW ) x[n+1-j] = xBound[n-j];
        }

    // calculate price. 
	simplePrice(n, x+1, fvec+1);
		
	// if jacobiOnce, only calculate jacobi once during the 1st sens call
	if( nbJacobi < calibFD->maxNbJacobi )
	{

		// tweak to get sensitivity, ie. jacobi
		double dx;
		vector<double> f(n+1);
		for(int j=1; j<=n; j++)
		{
			dx = 0.5*DEF_UNIT_TOL_X;
			x[j] += dx;
			simplePrice(n, x+1, &f[1]);
			for(int i=1; i<=n; i++)
				fjac[i][j] = (f[i] - fvec[i])/dx;
			x[j] -= dx;
		}

		nbJacobi++;
	}

    // calculate the optimal x value if the other x value is 0
    if( !isRiskyVolOW )
    {
        double fac = fjac[1][1] * fjac[1][2] + fjac[2][1] * fjac[2][2];
        for(j=1; j<=n; j++)
        {
            xBound[j-1] = x[j] + (x[n+1-j] * fac - fvec[1]*fjac[1][j] - fvec[2]*fjac[2][j])/(fjac[1][j]*fjac[1][j] + fjac[2][j]*fjac[2][j]);
            if( !Maths::isPositive( xBound[j-1] ) ) xBound[j-1] = DBL_EPSILON;
        }
    }
}

void FD1FDDECalib::simplePrice(int n, double *x, double *prices)
{
    curX = x;
    if(calibFD->isForward)
    {
        doPriceVanCallForward(curTimeStep, numPriceArrays-1, &strikes[curTimeStep], calibFD->firstStepSmooth, calibFD->useDensity, &prices[1], 0, 0, &prices[0]);       
        for(int k=0; k<numPriceArrays; k++)
            prices[k] *= fwds[curTimeStep] * pvs[curTimeStep];
    }
    else
        doPrice(curTimeStep, prices);
    
    for(int k=0; k<numPriceArrays; k++)
		prices[k] = (prices[k]-target[k])/errfDivisor[k];
}

// use input variable value to update drift/diffusion/discount terms
double FD1FDDECalib::getVol(int curTimeStep)
{
	return isRiskyVolOW?dRiskyVolOW:(curX[1] * errxDivisor[1]);
}

void FD1FDDECalib::updateStockArray(int nStockSteps, const double *stockArray)
{
    // set up S^power to speed up calculation
    int i;
    for (i=0; i<=nStockSteps ; ++i) powSpotArray[i] = ::pow(stockArray[i], -curPower);
}

void FD1FDDECalib::getSpreads(int curTimeStep)
{
    int i;
	double sprdBB = curX[0] * errxDivisor[0];
	for (i=0; i<=nStockSteps ; ++i)
		sprdArray[i] = sprdBB * powSpotArray[i];
}

// initialize terminal zero bond and vanilla call option payoff
void FD1FDDECalib::payoffAtMat()
{
    if( !isRiskyVolOW )
    {
        if( strikes[curTimeStep] < stockArray[0] || strikes[curTimeStep] > stockArray[nStockSteps] )
            throw ModelException("FD1FDDECalib::payoffAtMat", "Strike outside range. Need to increase TruncationStd");
    }
    
	for (int i=0; i<=nStockSteps ; ++i) {
		priceArray[0][i] = 1;
		if( !isRiskyVolOW )
			priceArray[1][i] = Maths::max(0.0, stockArray[i] - strikes[curTimeStep]);
	}
}


//////////////////////////////////////////////////////////////////////////////
//																			//
//			DDEModuleLN1F													//
//																			//
//////////////////////////////////////////////////////////////////////////////


class DDEModuleLN : public DDEModule
{
public:
    static CClassConstSP const TYPE;

	static void load(CClassSP &clazz)
	{
		clazz->setPrivate(); // make visible to EAS/spreadsheet
		REGISTER(DDEModuleLN, clazz);
		SUPERCLASS(DDEModule);
		EMPTY_SHELL_METHOD(defaultDDEModuleLN);
	}
	static IObject* defaultDDEModuleLN(){
		return new DDEModuleLN();
	}

	void calib1F(const DateTimeArray &dates, const DoubleArray &t, const DoubleArray &pvs, const DoubleArray &ndps, const DoubleArray &sprdVars,
		const DoubleArray &fwds, const DoubleArray &strikes, const DoubleArray &bsVars, const DoubleArray &calls, DoubleArray &riskyVars);

	void calib2F(const DateTimeArray &dates, const DoubleArray &t, const DoubleArray &pvs, const DoubleArray &ndps,
		const DoubleArray &fwds, const DoubleArray &strikes, const DoubleArray &bsVars, const DoubleArray &calls, DoubleArray &riskyVars);

	DDEModuleLN(bool isOneFactor=true) : DDEModule(TYPE, isOneFactor){};

private:
};

CClassConstSP const DDEModuleLN::TYPE = CClass::registerClassLoadMethod(
              "DDEModuleLN", typeid(DDEModuleLN), DDEModuleLN::load);



// class used in MC calibration of DDE
class DDELN1FCalibData
{
public:
	// for risky vol
	DDELN1FCalibData(const DoubleArray *spot, const DoubleArray *ndp, double pv, double strike, double call) 
		: pv(pv), strike(strike), call(call), spot(spot), ndp(ndp){};

	static double calibFuncLN1FVol(double voldt, void *calibDataMem)
	{
		const DoubleArray &spot = *((DDELN1FCalibData *)calibDataMem)->spot;
		const DoubleArray &ndp = *((DDELN1FCalibData *)calibDataMem)->ndp;
		double pv = ((DDELN1FCalibData *)calibDataMem)->pv;
		double strike = ((DDELN1FCalibData *)calibDataMem)->strike;
		double call = ((DDELN1FCalibData *)calibDataMem)->call;
		int n, nbPath = spot.size();

		double opt=0, varDt = voldt * voldt;
		for(n=0; n<nbPath; n++)
		{
			opt += Black::price(true /* isCall */, spot[n], strike, 1, varDt) * ndp[n];
		}
		return call - opt * pv/nbPath;
	}

	// for spread backbone
	DDELN1FCalibData(const DoubleArray *ndp, const DoubleArray *powSpot, double zero) 
		: zero(zero), ndp(ndp), powSpot(powSpot){};

	static double calibFuncLN1FSprd(double fac, void *calibDataMem)
	{
		const DoubleArray &ndp = *((DDELN1FCalibData *)calibDataMem)->ndp;
		const DoubleArray &powSpot = *((DDELN1FCalibData *)calibDataMem)->powSpot;
		double zero = ((DDELN1FCalibData *)calibDataMem)->zero;
		int n, nbPath = powSpot.size();

		double pv=0;
		for(n=0; n<nbPath; n++)
		{
			pv += ndp[n] * exp(- fac * powSpot[n]);
		}
		return zero - pv/nbPath;
	}

	double pv, strike, call, zero;
	const DoubleArray *spot;
	const DoubleArray *ndp;
	const DoubleArray *powSpot;
};

void DDEModuleLN::calib1F(const DateTimeArray &dates, const DoubleArray &t, const DoubleArray &pvs, const DoubleArray &ndps, const DoubleArray &sprdVars,
		const DoubleArray &fwds, const DoubleArray &strikes, const DoubleArray &bsVars, const DoubleArray &calls, DoubleArray &riskyVars)
{
	static const string method = "DDEModuleLN::calib1F";

	try{
		// functional form Sprd = c * pow(Eq, p).
		// given equity vol and this form, p is calc to match spread var
		// fitting c and equity vol to match vanilla calls

		int nbDate = dates.size();
		smartPtr<SpreadEquityFuncLN> func(new SpreadEquityFuncLN());
		func->dates = DateTimeArraySP(new DateTimeArray());
		func->factor = CDoubleArraySP(new DoubleArray());
		func->power = CDoubleArraySP(new DoubleArray());		
		func->dts = CDoubleArraySP(new DoubleArray());		
		// insert value date into the spread function dates such that factor/power applies between each date pairs
		func->dates->push_back(valueDate);

		DoubleArray &factor=*func->factor, &power=*func->power;
		
		if( DDECalibFD::TYPE->isInstance(calibParams.get()) )
		{
			double guess[2], result[2], curPower, dRiskyVol, fac;
			FD1FDDECalib rootFinder;
			rootFinder.init(this, t, pvs, ndps, sprdVars, fwds, strikes, bsVars, calls);

			// use risky vol to set power and set initial guess for riskyVol
            // in case strike is very high w.r.t. fwd, call value will be too small for impliedVar()
            // use bsVar instead
			int i = 1;
			if( !riskyVolOW() ) // use RG to estimate the initial risky var
                riskyVars[i] = Maths::isZero(calls[i])?bsVars[i]:Black::impliedVariance(true, fwds[i], strikes[i]*ndps[i], pvs[i], bsVars[i], calls[i], bsVars[i]*DEF_IMPVAR_EPSILON);

			curPower = sqrt(sprdVars[i]/bsVars[i]);
			fac = -log(ndps[i]) * ::pow(fwds[0], curPower)/t[i];	// not including time frac
			dRiskyVol = sqrt(riskyVars[i]/t[i]);					// not including time frac
			for( i=1; i<nbDate; i++ )
			{	
				if( i>1 )
					// find spread power. defined to approximately match spread vol up to end of last time step
					curPower = sqrt(sprdVars[i-1]/bsVars[i-1]);

				if( riskyVolOW() )
					dRiskyVol = sqrt((riskyVars[i] - riskyVars[i-1]) / (t[i]-t[i-1]));

                // calibration for a time step
				guess[0] = fac;
				guess[1] = dRiskyVol;
				try {
					rootFinder.calibViaSimplePrice(i, curPower, guess, result);
				} catch (exception& e) {
					if( calibParams->throwIfFail || i==1)
						throw(ModelException(e, method, "Failed to calibrate at " + Format::toString(t[i]) + "Y"));
					else if( !riskyEqVols ) 
					{
						// output warning message for the initial calibration only. avoid excess output from sensitivity runs, esp implied sprd/vol runs
						ModelException modelErr(e, method + ":WARNING: Failed to calibrate at " + Format::toString(t[i]) + "Y");
						modelErr.errorLog();
					}
					break;
				}

				if( !riskyVolOW() )
				{
				// output the calibrated risky var
				dRiskyVol = result[1];
				riskyVars[i] = riskyVars[i-1] + dRiskyVol * dRiskyVol * (t[i]-t[i-1]);
				}

				// store spread results
				fac = result[0];
				factor.push_back(fac);
				power.push_back(curPower);
				func->dates->push_back(dates[i]);
				func->dts->push_back(t[i]-t[i-1]);
			}

		}
		else if( DDECalibMC::TYPE->isInstance(calibParams.get()) )
		{

		const DDECalibMC *calibMC = dynamic_cast<const DDECalibMC *>(calibParams.get());
		
		int i, n, nbPath = calibMC->nbPath;
		double dRiskyVol, p, fac;
		DoubleArray hazard(nbPath), spot(nbPath), ndp(nbPath), powSpot(nbPath);

		DoubleArray rand(nbPath, 0);
		IRandomSP randGenerator(copy(calibMC->rand.get()));
		randGenerator->init();

		// define the initial p to be ratio of sprd vol and black-schole vol
		power.push_back(sqrt(sprdVars[1]/(riskyVolOW()?riskyVars[1]:bsVars[1])));
		factor.push_back(-log(ndps[1]) * ::pow(fwds[0], power[0])/t[1]);
		func->dates->push_back(dates[1]);
		func->dts->push_back(t[1]);

		// spot and hazard rate array
		for( n=0; n<nbPath; n++)
		{
			hazard[n] = -log(ndps[1]);
			spot[n] = fwds[1]/ndps[1];
			ndp[n] = ndps[1];
		}
		for( i=1; i<nbDate; i++ )
		{	
			double guess;

			if( riskyVolOW() )
				dRiskyVol = sqrt((riskyVars[i] - riskyVars[i-1]) / (t[i]-t[i-1]));
			else
			{

			// calibrate eq vol to match vanilla call
			// root searching for risky vol
			guess = sqrt(bsVars[i] - bsVars[i-1]);
			DDELN1FCalibData calibData(&spot, &ndp, pvs[i], strikes[i], calls[i]);

			try {
				dRiskyVol = zbrentUseful(
						&DDELN1FCalibData::calibFuncLN1FVol,	/* (I) function to find the root of			*/
						&calibData,								/* (I) only 1 other parameter to function	*/
						0.01*guess,								/* (I) low value for risky vol				*/ 
						2*guess,								/* (I) high value for risky vol				*/
						calibParams->tolEqVol);					/* (I) tolerance. 1% of input vol			*/

				if( Maths::isZero(dRiskyVol) ) dRiskyVol = 0;
				if( Maths::isNegative(dRiskyVol) ) 
					throw ModelException(method, "Calibrated risky vol " + Format::toString(dRiskyVol) + " is < 0");

			} catch (exception& e) {
				if( calibParams->throwIfFail || i==1)
					throw(ModelException(e, method, "Failed to calibrate risky vol at " + Format::toString(t[i]) + "Y"));
				else if( !riskyEqVols ) 
				{
					// output warning message for the initial calibration only. avoid excess output from sensitivity runs, esp implied sprd/vol runs
					ModelException modelErr(e, method + ":WARNING: Failed to calibrate risky vol at " + Format::toString(t[i]) + "Y");
					modelErr.errorLog();
				}
				break;
			}

			// output the calibrated risky var
			riskyVars[i] = riskyVars[i-1] + dRiskyVol * dRiskyVol;
			
			} // end of checking riskyVolOW()


			if( i < (nbDate-1) )
			{
				// find spread power. defined to match spread vol
				p = sqrt(sprdVars[i]/bsVars[i]);

				// find spread backbone to match the risky zero bond
				double volDrift=-0.5*dRiskyVol * dRiskyVol;
				randGenerator->fetch(rand.size(), &rand[0]);
				for( n=0; n<nbPath; n++)
				{
					spot[n] *= exp(volDrift + dRiskyVol * rand[n]);
					powSpot[n] = ::pow(spot[n], -p);
				}

				// find spread backbone
				guess = log(ndps[i]/ndps[i+1])*::pow(fwds[i]/ndps[i], p);
				DDELN1FCalibData calibData2(&ndp, &powSpot, ndps[i+1]);

				try {
					fac = zbrentUseful(
							&DDELN1FCalibData::calibFuncLN1FSprd,	/* (I) function to find the root of			*/
							&calibData2,							/* (I) only 1 other parameter to function	*/
							0.01*guess,								/* (I) low value for sprd backbone			*/ 
							5*guess,								/* (I) high value for sprd backbone			*/
							calibParams->tolSprd);					/* (I) tolerance. 1% of input bbone			*/

					if( Maths::isZero(fac) ) fac = 0;
					if( Maths::isNegative(fac) ) 
						throw ModelException(method, "Calibrated sprd func factor " + Format::toString(fac) + " is < 0");

				} catch (exception& e) {
					if( calibParams->throwIfFail || i==1)
						throw(ModelException(e, method, "Failed to calibrate spread backbone at " + Format::toString(t[i]) + "Y"));
					else if( !riskyEqVols ) 
					{
						// output warning message for the initial calibration only. avoid excess output from sensitivity runs, esp implied sprd/vol runs
						ModelException modelErr(e, method + ":WARNING: Failed to calibrate spread backbone at " + Format::toString(t[i]) + "Y");
						modelErr.errorLog();
					}
					break;
				}
				
				double dfwd = fwds[i+1]/fwds[i];
				for( n=0; n<nbPath; n++)
				{
					hazard[n] = fac * powSpot[n];
					spot[n] *= dfwd * exp(hazard[n]);
					ndp[n] *= exp(-hazard[n]);
				}

				// store results
				factor.push_back(fac/(t[i+1]-t[i]));
				power.push_back(p);
				func->dates->push_back(dates[i+1]);
				func->dts->push_back(t[i+1]-t[i]);
			}
		}
		

		} else 			
			throw(ModelException(method, "Currently only support FD1F or MonteCarlo calibration"));
		
		// create spread function
		spreadFunc = func;

	} catch(exception& e) {
        throw ModelException(e, method, "Failed to calibrate.");
	}

}

void DDEModuleLN::calib2F(const DateTimeArray &dates, const DoubleArray &t, const DoubleArray &pvs, const DoubleArray &ndps,
		const DoubleArray &fwds, const DoubleArray &strikes, const DoubleArray &bsVars, const DoubleArray &calls, DoubleArray &riskyVars)
{
	throw ModelException("DDEModuleLN::calib2F not implemented yet");
}



//////////////////////////////////////////////////////////////////////////////
//																			//
//	Spread Equity Function implementation for RG							//
//																			//
//////////////////////////////////////////////////////////////////////////////


class SpreadEquityFuncRG: public SpreadEquityFunc
{
public:
// not done yet
    static CClassConstSP const TYPE;
	friend class DDEModuleRG;

	void getSpreadCC(const DateTime &start,
					const DateTime &end,
					int		nbSpot, 
					const double *spots,
					double *spreads,
					double dt=0) const;

	static void load(CClassSP& clazz);

	static IObject * defaultSpreadEquityFuncRG()
	{
		return new SpreadEquityFuncRG;
	}

private:
	CleanSpreadCurveSP	sprds; // clean spread curve
	
	SpreadEquityFuncRG() : SpreadEquityFunc(TYPE) {};
};

void SpreadEquityFuncRG::getSpreadCC(
						const DateTime &start,
						const DateTime &end0,
						int	nbSpot, 
						const double *spots,
						double *spreads,
						double dt) const
{	
	// sanity check
	DateTime end(end0);
	if( end < start ) 
		throw(ModelException("SpreadEquityFuncRG::getSpreadCC", "Start date > end date"));
	if( end < start.rollDate(1) )
		end = start.rollDate(1);

	// calculate spreads for the spots
	if( Maths::isZero(dt) ) dt = metric->yearFrac(start, end);
	double spread = -log(sprds->getDefaultPV(start, end))/dt;
	for(int i=0; i<nbSpot; i++)
	{
		spreads[i] = spread;
	}

}

void SpreadEquityFuncRG::load(CClassSP& clazz){ 
		clazz->setPrivate(); // make invisible to EAS/spreadsheet
		REGISTER(SpreadEquityFuncRG, clazz);
		SUPERCLASS(SpreadEquityFunc);
		EMPTY_SHELL_METHOD(SpreadEquityFuncRG::defaultSpreadEquityFuncRG);
		FIELD(sprds,	"");
}

CClassConstSP const SpreadEquityFuncRG::TYPE = CClass::registerClassLoadMethod(
              "SpreadEquityFuncRG", typeid(SpreadEquityFuncRG), SpreadEquityFuncRG::load);


//////////////////////////////////////////////////////////////////////////////
//																			//
//			DDEModuleRG class												//
//																			//
//////////////////////////////////////////////////////////////////////////////


class DDEModuleRG : public DDEModule
{
public:
    static CClassConstSP const TYPE;

	static void load(CClassSP &clazz)
	{
		clazz->setPrivate(); // make visible to EAS/spreadsheet
		REGISTER(DDEModuleRG, clazz);
		SUPERCLASS(DDEModule);
		EMPTY_SHELL_METHOD(defaultDDEModuleRG);
	}
	static IObject* defaultDDEModuleRG(){
		return new DDEModuleRG();
	}

	void calib1F(const DateTimeArray &dates, const DoubleArray &t, const DoubleArray &pvs, const DoubleArray &ndps, const DoubleArray &,//spot sprd var is not used
		const DoubleArray &fwds, const DoubleArray &strikes, const DoubleArray &bsVars, const DoubleArray &calls, DoubleArray &riskyVars);

	void calib2F(const DateTimeArray &dates, const DoubleArray &t, const DoubleArray &pvs, const DoubleArray &ndps,
		const DoubleArray &fwds, const DoubleArray &strikes, const DoubleArray &bsVars, const DoubleArray &calls, DoubleArray &riskyVars)
	{}; // empty

	void createSpreadFunc();

	DDEModuleRG() : DDEModule(TYPE, true){};

private:
};

typedef smartConstPtr<DDEModuleRG> DDEModuleRGConstSP;
typedef smartPtr<DDEModuleRG>      DDEModuleRGSP;

CClassConstSP const DDEModuleRG::TYPE = CClass::registerClassLoadMethod(
              "DDEModuleRG", typeid(DDEModuleRG), DDEModuleRG::load);


void DDEModuleRG::calib1F(const DateTimeArray &dates, const DoubleArray &t, const DoubleArray &pvs, const DoubleArray &ndps, const DoubleArray &,
		const DoubleArray &fwds, const DoubleArray &strikes, const DoubleArray &bsVars, const DoubleArray &calls, DoubleArray &riskyVars)
{
	static const string method = "DDEModuleRG::calib";

	try {
		if( !riskyVolOW() )
		{
            int i, nbDates = dates.size();
            for(i=1; i<nbDates; i++)
            {
                double guess = bsVars[i];
                riskyVars[i] = Maths::isZero(calls[i])?bsVars[i]:Black::impliedVariance(true /* isCall */, fwds[i], strikes[i] * ndps[i], pvs[i], guess, calls[i], guess*DEF_IMPVAR_EPSILON);
            }
		}

		createSpreadFunc();

    } catch( exception &e ) {
		throw ModelException(e, method);
	}
}

void DDEModuleRG::createSpreadFunc()
{
	//create spreadFunc
	smartPtr<SpreadEquityFuncRG> func(new SpreadEquityFuncRG());
	func->sprds = CleanSpreadCurveSP(copy(cleanSprds.get()));
	spreadFunc = SpreadEquityFuncSP(func);
}


//////////////////////////////////////////////////////////////////////////////
//																			//
//	DDEModule creator														//
//																			//
//////////////////////////////////////////////////////////////////////////////


/* create DDEModule depending on required type */
DDEModuleSP DDEModule::createModule(const string &ddeType, const DDECalib *calib, const AssetDDE *asset, const IDDEInitiator *initiator)
{
	static string method = "DDEModule::createModule";
	DDEModuleSP ddeModule;

	try {
		// create appropriate types of DDEModule
		if( ddeType == "RG" )
		{
			ddeModule = DDEModuleSP(new DDEModuleRG());
		}
		else if ( ddeType == "LN1F" )
		{
			ddeModule = DDEModuleSP(new DDEModuleLN(true));
		}
		else if ( ddeType == "LN2F" )
		{
			ddeModule = DDEModuleSP(new DDEModuleLN(false));
		}
		else
			throw(ModelException(method, "Only RiskGrow and LN DDE implemented"));
	
		ddeModule->calibParams = DDECalibConstSP(calib);
		ddeModule->update(asset, initiator);

	} catch (exception &e) {
		throw(ModelException(e, method));
	}

	return ddeModule;
}



////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  DDE module calibration addin
//
////////////////////////////////////////////////////////////////////////////////////////////////////

class DDEModuleAddin: public CObject {
    static CClassConstSP const TYPE;

    // input parameters
    CAssetWrapper       asset;
    IObjectSP           initDDE; // dde initiator instance
    DDECalibSP          calibParams;
    MarketDataSP        market;
    string              ycName; // yield curve name
    string              ddeType;
    OutputRequestArraySP outputRequests;

    class DummyModel : public CModel, virtual public Asset::ICanHaveDDE
    {
    public:
        static CClassConstSP const TYPE;
        DummyModel() : CModel(TYPE){}
        
        void Price(CInstrument*  instrument, 
                       CControl*     control, 
                       CResults*     results){}
        bool isDDE() const { return true; }

        IModel::WantsRiskMapping wantsRiskMapping() const {
            return riskMappingIrrelevant;
        }
    };

    static IObjectSP getCalibOutput(DDEModuleAddin* params) {
        static const string routine = "DDEModuleAddin::getCalibOutput";
        try {
            // create asset dde and directly set its credit curve member (normally done in getMarket)
            DummyModel model;
            CAsset::getAssetMarketData( &model, 
                                params->market.get(),
                                "N", //ccyTreatment
                                params->ycName,
                                params->asset);
            AssetDDESP assetDDE = AssetDDESP::dynamicCast(params->asset.getSP());

            // create DDE initiator using input calib dates
            const IDDEInitiator *initDDE = dynamic_cast<const IDDEInitiator *>(params->initDDE.get());
            if( !initDDE )
                throw ModelException(routine, "Member initDDE is not instance of IDDEInitiator");

            // create dde module
            DDEModuleSP ddeModule(DDEModule::createModule(params->ddeType, params->calibParams.get(), assetDDE.get(), initDDE));

            // create control and get output into result
            SensitivityArraySP sens;
            CControl ctrl(sens, params->outputRequests, false, "");
            ResultsSP results(new Results);
            ddeModule->recordOutputRequests(&ctrl, results.get());

            return results;
        } 
        catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    /** for reflection */
    DDEModuleAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        clazz->setPublic();         // make visible to EAS/spreadsheet
        REGISTER(DDEModuleAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultDDEModuleAddin);
        FIELD(ddeType,           "DDE modele type (RG, LN1F etc)");
        FIELD(ycName,            "risk free yield curve");
        FIELD(asset,             "asset");
        FIELD(market,                   "market data");
        FIELD(calibParams,              "calibration parameters");
        FIELD(initDDE,                  "calibration dates and strikes");
        FIELD(outputRequests,           "output requests");

        Addin::registerClassObjectMethod("DDE_CALIB",
                                         Addin::MARKET,
                                         "calibrate asset/credit under DDE",
                                         TYPE,
                                         false,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)getCalibOutput);
    }

    static IObject* defaultDDEModuleAddin(){
        return new DDEModuleAddin();
    }
    
};
   
CClassConstSP const DDEModuleAddin::TYPE = CClass::registerClassLoadMethod(
    "DDEModuleAddin", typeid(DDEModuleAddin), load);

CClassConstSP const DDEModuleAddin::DummyModel::TYPE = CClass::registerInterfaceLoadMethod(
    "DDEModuleAddin::DummyModel", typeid(DDEModuleAddin::DummyModel), 0);


////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  IDDEInitiator and Simple class of DDE initiator
//
////////////////////////////////////////////////////////////////////////////////////////////////////

IDDEInitiator::~IDDEInitiator(){} // empty

void IDDEInitiator::sensitiveDates(DateTimeArray &dates) const
{
	// default implementation
	dates.resize(1);
	dates[0] = maxMaturity();
}

// given dates, get strikes. dates may include model required dates plus the sensitiveDates.
void IDDEInitiator::sensitiveStrikes(const DateTimeArray	dates,
									DoubleArray				&strikes,	// same dimension as dates
									bool					&strikeIsPct) const
{
	if( dates.size() != strikes.size() )
		throw ModelException("IDDEInitiator::sensitiveStrikes", "Dates and strikes dimension mismatch");
		
	// default implementation
	strikeIsPct = true;
	for(int i=0; i<dates.size(); i++) strikes[i] = 1; // ATM strike
}


class SimpleDDEInitiator : public CObject, virtual public IDDEInitiator
{
public:
    static CClassConstSP const TYPE;

    DateTime maxMaturity() const { return calibSchd->lastDate(); }
    void sensitiveDates(  DateTimeArray &dates ) const
    {
        dates = calibSchd->getDates();
    }
    void sensitiveStrikes(  const DateTimeArray     dates,
                            DoubleArray             &strikes,
                            bool                    &strikeIsPct) const
    {
        if( dates.size() != strikes.size() )
            throw ModelException("SimpleDDEInitiator::sensitiveStrikes", "Internal error. Dates and Strikes dimension mismatch");
        
        strikeIsPct = strkIsPct;

        for(int i=dates.size()-1; i>=0; i--)
        {
            try {
                strikes[i] = calibSchd->interpolate(dates[i]);
            } catch (exception &e) {
                if( i==(dates.size()-1) ) 
                    throw ModelException(e, "SimpleDDEInitiator::sensitiveStrikes");

                 // use strike for later date if not in schedule
                strikes[i] = strikes[i+1];
            }
        }
    }

    void validatePop2Object()
    {
        if( !calibSchd || !calibSchd->length() )
            throw ModelException("SimpleDDEInitiator", "Empty calib schedule");
    }

    static void load(CClassSP& clazz){
        clazz->setPublic();         // make visible to EAS/spreadsheet
        REGISTER(SimpleDDEInitiator, clazz);
        EMPTY_SHELL_METHOD(defaultSimpleDDEInitiator);
        SUPERCLASS(CObject);
        FIELD(calibSchd,        "calibration dates and strikes");
        FIELD(strkIsPct, "if strike is percentage relative strk");
        FIELD_MAKE_OPTIONAL(strkIsPct);
    }
    static IObject* defaultSimpleDDEInitiator(){
        return new SimpleDDEInitiator();
    }

private:
    SimpleDDEInitiator() : CObject(TYPE), strkIsPct(false){};

    ScheduleSP		calibSchd;
    bool            strkIsPct; 
};

CClassConstSP const SimpleDDEInitiator::TYPE = CClass::registerClassLoadMethod(
              "SimpleDDEInitiator", typeid(SimpleDDEInitiator), SimpleDDEInitiator::load);



//////////////////////////////////////////////
//                                          //
//          DDECalib                        //
//                                          //
//////////////////////////////////////////////


DDECalib::DDECalib()
:   CObject(TYPE), useProdDateOnly(false), useEqualStep(false), templateFreq(new MaturityPeriod(DEF_TEMPLATE_FREQ)), expiries(DDECalib::createDefExpiries())
, nStepPerBenchMk(DEF_NB_STEPPERBENCHMK), maxStepSize(DEF_MAX_STEP_SIZE)
, nStepPerYear(DEF_NB_STEPPERYEAR), tolSprd(DEF_TOLERANCE_SPRD), tolEqVol(DEF_TOLERANCE_VOL), useAtm(false), atmIsFwd(false)
, strkOffSurfaceAdj(true), useRiskyVolOW(false), throwIfFail(true), maxHorizon(ExpirySP(new MaturityPeriod(DEF_MAX_HORIZON)))
, minHorizon(ExpirySP(new MaturityPeriod(DEF_MIN_HORIZON)))
{}

DDECalib::DDECalib(const CClassConstSP &clazz)
:   CObject(clazz), useProdDateOnly(false), useEqualStep(false), templateFreq(new MaturityPeriod(DEF_TEMPLATE_FREQ)), expiries(DDECalib::createDefExpiries())
, nStepPerBenchMk(DEF_NB_STEPPERBENCHMK), maxStepSize(DEF_MAX_STEP_SIZE)
, nStepPerYear(DEF_NB_STEPPERYEAR), tolSprd(DEF_TOLERANCE_SPRD), tolEqVol(DEF_TOLERANCE_VOL), useAtm(false), atmIsFwd(false)
, strkOffSurfaceAdj(true), useRiskyVolOW(false), throwIfFail(true), maxHorizon(ExpirySP(new MaturityPeriod(DEF_MAX_HORIZON)))
, minHorizon(ExpirySP(new MaturityPeriod(DEF_MIN_HORIZON)))
{}

void DDECalib::validatePop2Object()
{
    if( nStepPerYear == -1 ) nStepPerYear = DEF_NB_STEPPERYEAR;
	if( nStepPerBenchMk == -1 ) nStepPerBenchMk = DEF_NB_STEPPERBENCHMK;
	if( Maths::isZero(maxStepSize + 1) ) maxStepSize = DEF_MAX_STEP_SIZE;
	if( Maths::isZero(tolSprd + 1) ) tolSprd = DEF_TOLERANCE_SPRD;
	if( Maths::isZero(tolEqVol + 1) ) tolEqVol = DEF_TOLERANCE_VOL;

    if( useEqualStep && nStepPerYear < 1 )
        throw ModelException("DDECalib::validatePop2Object", "nStepPerYear must be >= 1 or is -1 (ie. use default)");
    if( !useEqualStep && nStepPerBenchMk < 1 )
        throw ModelException("DDECalib::validatePop2Object", "nStepPerBenchMk must be >= 1 or is -1 (ie. use default)");
    if( !useEqualStep && maxStepSize < 1/365. )
        throw ModelException("DDECalib::validatePop2Object", "maxStepSize must be >= 1/365 or is -1 (ie. use default)");
    if( useEqualStep && !templateFreq )
        throw ModelException("DDECalib::validatePop2Object", "Empty templateFreq when useEqualStep = TRUE");
    if( !useEqualStep && (!expiries || expiries->size() == 0) )
        throw ModelException("DDECalib::validatePop2Object", "Empty expiries when useEqualStep = FALSE");
	if( !Maths::isPositive(tolSprd) || !Maths::isPositive(tolEqVol) )
        throw ModelException("DDECalib::validatePop2Object", "Sprd/Vol tolerance must be > 0 or is -1 (ie. use default)");
}

DateTimeArraySP DDECalib::getTemplateDates( DateTime valueDate,
                                            DateTime lastDate) const
{
    DateTimeArraySP timeTemplate(new DateTimeArray);
    DateTimeArray &tempDates = *timeTemplate;
    tempDates.push_back(valueDate); // first date is value date

    if( !useProdDateOnly ){ // insert dates only if useProdDateOnly is false

    if( useEqualStep )
    {
        DateTime currDate = templateFreq->toDate(valueDate);
        if( valueDate != currDate )
        {
            while( currDate < lastDate )
            {
                tempDates.push_back(currDate);
                currDate = templateFreq->toDate(currDate);
            }
        }
    } 
    else 
    {
        for(int i=0; i<expiries->size(); i++)
        {
            DateTime currDate = (*expiries)[i]->toDate(valueDate);
			if( currDate <= (i==0?valueDate:tempDates[i-1]) )
				throw ModelException("DDECalib::getTemplateDates", "Input calib expiries must be strictly ascending");
            if( currDate >= lastDate ) break;
            tempDates.push_back(currDate);
        }
    }
    } // end of if useProdDateOnly

    tempDates.push_back(lastDate); // make sure last date is part of template

    return timeTemplate;
}

ExpiryArraySP DDECalib::createDefExpiries()
{
    const int DEF_NB_EXPIRY = 11;
    static string DEF_EXPIRIES[DEF_NB_EXPIRY] = {"1Y", "2Y", "3Y", "4Y", "5Y", "6Y", "7Y", "8Y", "10Y", "15Y", "20Y"};

    ExpiryArraySP expiries(new ExpiryArray(DEF_NB_EXPIRY));
    for(int i=0; i<expiries->size(); i++)
    {
        (*expiries)[i] = ExpirySP(new MaturityPeriod(DEF_EXPIRIES[i]));
    }
    return expiries;
}

class DDECalibHelper{
public:
    /** Invoked when Class is 'loaded' */

    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(DDECalib, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultDDECalib);
        FIELD(templateFreq,         "template frequency for calib template. 3M, 1W etc. used if useEqualStep is Y");
        FIELD_MAKE_OPTIONAL(templateFreq);
        FIELD(expiries,             "expiries for calib template. used if useEqualStep is N");
        FIELD_MAKE_OPTIONAL(expiries);
        FIELD(nStepPerYear,  "nb of steps per year for finer step between benchmarks. used if useEqualStep is Y");
        FIELD_MAKE_OPTIONAL(nStepPerYear);
        FIELD(nStepPerBenchMk,  "nb of finer step between benchmarks for calibration. used if useEqualStep is N");
        FIELD_MAKE_OPTIONAL(nStepPerBenchMk);
        FIELD(maxStepSize,	"max size of finer step between benchmarks. used if useEqualStep is N");
        FIELD_MAKE_OPTIONAL(maxStepSize);
        FIELD(tolSprd,       "tolerance for spread");
        FIELD_MAKE_OPTIONAL(tolSprd);
        FIELD(tolEqVol,      "tolerance for equity vol");
        FIELD_MAKE_OPTIONAL(tolEqVol);
        FIELD(useAtm,        "whether always calibrate to atm vol");
        FIELD_MAKE_OPTIONAL(useAtm);
        FIELD(atmIsFwd,      "0 if ATM means atm spot (default), 1 if atm fwd");
        FIELD_MAKE_OPTIONAL(atmIsFwd);
        FIELD(useProdDateOnly,  "1 if don't insert dates using freq/expiries (default), 0 otherwise");
        FIELD_MAKE_OPTIONAL(useProdDateOnly);
        FIELD(useEqualStep,  "0 if use expiries (default), 1 if use template freq");
        FIELD_MAKE_OPTIONAL(useEqualStep);
        FIELD(strkOffSurfaceAdj,  "1 if use lowest strk on surface when calib strk too low (default), 0 otherwise");
        FIELD_MAKE_OPTIONAL(strkOffSurfaceAdj);
        FIELD(throwIfFail,		"1 if use risky vol overwrite when provided, 0 otherwise (default)");
        FIELD_MAKE_OPTIONAL(throwIfFail);
        FIELD(useRiskyVolOW,		"1 if throw error when calib fail (default), 0 otherwise");
        FIELD_MAKE_OPTIONAL(useRiskyVolOW);
        FIELD(maxHorizon,		"max horizon for calibration, in case deal has very long matuirty");
        FIELD_MAKE_OPTIONAL(maxHorizon);
        FIELD(minHorizon,		"min horizon for calibration, in case deal has very short matuirty");
        FIELD_MAKE_OPTIONAL(minHorizon);
    }

    static IObject* defaultDDECalib(){
        return new DDECalib();
    }
};


CClassConstSP const DDECalib::TYPE = CClass::registerClassLoadMethod(
              "DDECalib", typeid(DDECalib), DDECalibHelper::load);


DDECalibMC::DDECalibMC()
:   DDECalib(TYPE), nbPath(DEF_NB_PATH) {}


class DDECalibMCHelper{
public:
    /** Invoked when Class is 'loaded' */

    static void load(CClassSP& clazz){
        clazz->setPublic();         // make visible to EAS/spreadsheet
        REGISTER(DDECalibMC, clazz);
        EMPTY_SHELL_METHOD(defaultDDECalibMC);
        SUPERCLASS(DDECalib);
        FIELD(nbPath,        "nb of MonteCarlo path");
        FIELD_MAKE_OPTIONAL(nbPath);
        FIELD(rand,                 "random number generator");
    }
    static IObject* defaultDDECalibMC(){
        return new DDECalibMC();
    }
};

CClassConstSP const DDECalibMC::TYPE = CClass::registerClassLoadMethod(
              "DDECalibMC", typeid(DDECalibMC), DDECalibMCHelper::load);

// default grid type is exponential cause sinh type causes too few point around atm fwd
// for long term horizon and/or large vol
DDECalibFD::DDECalibFD(bool isForward)
:   DDECalib(TYPE), nStockSteps(DEF_NB_STOCKSTEP), gridType(2), maxIter(DEF_NB_MAX_ITER), 
TruncationStd(DEF_TRUNCATION), maxNbJacobi(DEF_NB_MAX_ITER), isForward(isForward), firstStepSmooth(true),
useDensity(true), creditStrike(DEF_CREDIT_STRK)
{
}

MaturityPeriodArraySP DDECalibFD::createDefSegExpiries()
{
	static string DEF_SEGMENT_ENDS[DEF_NB_SEGMENT_ENDS] = {"1M", "3M", "1Y", "4Y", "16Y"};
	MaturityPeriodArraySP segExpiries = MaturityPeriodArraySP(new MaturityPeriodArray(DEF_NB_SEGMENT_ENDS));
	for(int i = 0; i < DEF_NB_SEGMENT_ENDS; i++)
	{
		(*segExpiries)[i] = MaturityPeriodSP(new MaturityPeriod(DEF_SEGMENT_ENDS[i]));
	}
	return segExpiries;
}

IntArraySP DDECalibFD::createDefSegStepPerYears()
{
	static int DEF_SEGMENT_PPYS[DEF_NB_SEGMENT_ENDS] = {120, 60, 25, 20, 10};
    IntArraySP segStepPerYears = IntArraySP(new IntArray(DEF_NB_SEGMENT_ENDS));
	for(int i = 0; i < DEF_NB_SEGMENT_ENDS; i++)
	{
		(*segStepPerYears)[i] = DEF_SEGMENT_PPYS[i];
	}
	return segStepPerYears;
}


void DDECalibFD::validatePop2Object()
{
    static const string method = "DDECalibFD::validatePop2Object";

    int i;

    DDECalib::validatePop2Object();

    if( nStockSteps == -1 ) nStockSteps = DEF_NB_STOCKSTEP;
    if( maxIter == -1 ) maxIter = DEF_NB_MAX_ITER;
    if( gridType == -1 ) gridType = 0;
    if( Maths::isZero(TruncationStd + 1) ) TruncationStd = DEF_TRUNCATION;
    if( maxNbJacobi == -1 ) maxNbJacobi = DEF_NB_MAX_ITER;
    if( maxNbJacobi > maxIter ) maxNbJacobi = maxIter;
	if( Maths::isZero(creditStrike + 1) ) creditStrike = DEF_CREDIT_STRK;

    if( nStockSteps < 1 )
        throw ModelException(method, "nStockSteps must be >= 1 or is -1 (ie. use default)");
    if( maxIter < 1 )
        throw ModelException(method, "maxIter must be >= 1 or is -1 (ie. use default)");
    if( gridType < 0 || gridType > 2 )
        throw ModelException(method, "gridType must be 0, 1, 2 or is -1 (ie. use default)");
    if( !Maths::isPositive(TruncationStd) )
        throw ModelException(method, "TruncationStd must be > 0 or is -1 (ie. use default)");
    if( maxNbJacobi < 1 )
        throw ModelException(method, "maxNbJacobi must be > 1 or is -1 (ie. use default)");
	if( !Maths::isPositive(creditStrike) )
        throw ModelException(method, "credit strike must be > 0 or is -1 (ie. use default)");
	if( !segExpiries || segExpiries->size() == 0 )
    {
        if( !!segStepPerYears && segStepPerYears->size() > 0 )
            throw ModelException(method, "Can not input segStepPerYears without input segExpiries");

		segExpiries = createDefSegExpiries();
        segStepPerYears = createDefSegStepPerYears();
    }
    else 
    {
        if( !segStepPerYears )
        {
            segStepPerYears = IntArraySP(new IntArray(segExpiries->size()));
	        for(i=0; i<segStepPerYears->size(); i++)
		        (*segStepPerYears)[i] = nStepPerYear;
        }
        else
        {
            if ( segStepPerYears->size() != segExpiries->size() )
                throw ModelException(method, "segment expiries and segment stepPerYears must be same size");

            for(i=0; i<segStepPerYears->size(); i++)
                if( (*segStepPerYears)[i] < 1 )
                    throw ModelException(method, "segment stepPerYears must be >= 1");
        }
    }
    
    // for forward calibration, grid need to have enough low strike points.
    if( isForward )
    {
        if( gridType != 2 )
            throw ModelException(method, "Must use gridType 2 if forward calib");
    }
}

class DDECalibFDHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic();         // make visible to EAS/spreadsheet
        REGISTER(DDECalibFD, clazz);
        EMPTY_SHELL_METHOD(defaultDDECalibFD);
        SUPERCLASS(DDECalib);
        FIELD(nStockSteps,       "nb of stock steps per time slice.");
        FIELD_MAKE_OPTIONAL(nStockSteps);
        FIELD(gridType,          "0 for sin hyperbolic, 1 for linear, 2 for exponential");
        FIELD_MAKE_OPTIONAL(gridType);
        FIELD(maxIter,           "max nb of iteration in Newton root finding");
        FIELD_MAKE_OPTIONAL(maxIter);
        FIELD(TruncationStd,     "num of stdev to truncate");
        FIELD_MAKE_OPTIONAL(TruncationStd);
        FIELD(maxNbJacobi,       "if only calc jacobi for 1st iteration in root finding");
        FIELD_MAKE_OPTIONAL(maxNbJacobi);
        FIELD(segExpiries,              "expiries for FD segments");
        FIELD_MAKE_OPTIONAL(segExpiries);
        FIELD(segStepPerYears,          "num step per year for FD segments (only used for forward calib for now)");
        FIELD_MAKE_OPTIONAL(segStepPerYears);
        FIELD(isForward,         "true if use forward calibration.");
        FIELD_MAKE_OPTIONAL(isForward);
        FIELD(firstStepSmooth,   "true if do smoothing for 1st time step using Black price");
        FIELD_MAKE_OPTIONAL(firstStepSmooth);
        FIELD(useDensity,        "true if forward calc on function (C - K * dC/dK).");
        FIELD_MAKE_OPTIONAL(useDensity);
        FIELD(creditStrike,      "strike relative to ATM. used to calib spread from vanilla call's sensitivity to strike. not relevant if useDensity");
        FIELD_MAKE_OPTIONAL(creditStrike);
    }
    static IObject* defaultDDECalibFD(){
        return new DDECalibFD();
    }
};

CClassConstSP const DDECalibFD::TYPE = CClass::registerClassLoadMethod(
              "DDECalibFD", typeid(DDECalibFD), DDECalibFDHelper::load);


////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Clean spread vol curve addin
//
////////////////////////////////////////////////////////////////////////////////////////////////////

class CleanSpreadVolCurveAddin: public CObject {
    static CClassConstSP const TYPE;

    // input parameters
    YieldCurveSP        yieldCurve;
    CDSParSpreadsSP     parSpreadCurve;
    DateTime            valueDate;
	DDEParamsSP			ddeParams;
	TimeMetricSP		timeMetric;

    static IObjectSP getCleanSpreadVols(CleanSpreadVolCurveAddin* params) {
        static const string routine = "CleanSpreadVolCurveAddin::getCleanSpreads";
        try {

            CleanSpreadCurveSP cleanSpreads = CDSHelper::getCleanSpreadCurve(
                                                    *params->parSpreadCurve,
                                                    params->valueDate,
                                                    true); // returnFwdRates

			CleanSpreadVolCurveSP cleanSpreadVols = CleanSpreadVolCurveSP(
					DDEModule::createCleanSpreadVolCurve(
									cleanSpreads.get(),
									params->ddeParams.get(),
									params->timeMetric.get()));

            DateTime::DateArraySP dateArray(new DateTime::DateArray(0));

			DateTimeArray dateTimeArray(cleanSpreads->getDates());
			for (int i = 0; i < dateTimeArray.size(); i++)
				dateArray->push_back(DateTime::DateSP(new DateTime::Date(dateTimeArray[i].getDate())));

            ObjectArraySP objArray(new ObjectArray(0));

            objArray->push_back(dateArray);
            objArray->push_back(cleanSpreadVols->vols);

            return objArray;
        } 
        catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    /** for reflection */
    CleanSpreadVolCurveAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(CleanSpreadVolCurveAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCleanSpreadVolCurveAddin);
        FIELD(yieldCurve,               "risk free yield curve");
        FIELD(parSpreadCurve,           "CDS par spread curve");
        FIELD(valueDate,         "value date");
        FIELD(ddeParams,				"DDE parameters");
        FIELD(timeMetric,				"time metric");

        Addin::registerClassObjectMethod("CLEAN_SPREAD_VOL_CURVE",
                                         Addin::MARKET,
                                         "derives the vol curve of clean fwd spread given a CDS curve",
                                         TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)getCleanSpreadVols);
    }

    static IObject* defaultCleanSpreadVolCurveAddin(){
        return new CleanSpreadVolCurveAddin();
    }
    
};
   
CClassConstSP const CleanSpreadVolCurveAddin::TYPE = CClass::registerClassLoadMethod(
    "CleanSpreadVolCurveAddin", typeid(CleanSpreadVolCurveAddin), load);



////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Equity spread corr curve addin
//
////////////////////////////////////////////////////////////////////////////////////////////////////

class EquitySpreadCorrCurveAddin: public CObject {
    static CClassConstSP const TYPE;

    // input parameters
    YieldCurveSP        yieldCurve;
    CDSParSpreadsSP     parSpreadCurve;
    DateTime            valueDate;
	DDEParamsSP			ddeParams;

    static IObjectSP getEquitySpreadCorrs(EquitySpreadCorrCurveAddin* params) {
        static const string routine = "EquitySpreadCorrCurveAddin::getCleanSpreads";
        try {

            CleanSpreadCurveSP cleanSpreads = CDSHelper::getCleanSpreadCurve(
                                                    *params->parSpreadCurve,
                                                    params->valueDate,
                                                    true); // returnFwdRates

			EquitySpreadCorrCurveSP eqSpreadCorrs = EquitySpreadCorrCurveSP(
				DDEModule::createEquitySpreadCorrCurve(
									cleanSpreads.get(),
									params->ddeParams.get()));

            DateTime::DateArraySP dateArray(new DateTime::DateArray(0));

			DateTimeArray dateTimeArray(cleanSpreads->getDates());
			for (int i = 0; i < dateTimeArray.size(); i++)
				dateArray->push_back(DateTime::DateSP(new DateTime::Date(dateTimeArray[i].getDate())));

            ObjectArraySP objArray(new ObjectArray(0));

            objArray->push_back(dateArray);
            objArray->push_back(eqSpreadCorrs->corrs);

            return objArray;
        } 
        catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    /** for reflection */
    EquitySpreadCorrCurveAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(EquitySpreadCorrCurveAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultEquitySpreadCorrCurveAddin);
        FIELD(yieldCurve,               "risk free yield curve");
        FIELD(parSpreadCurve,           "CDS par spread curve");
        FIELD(valueDate,         "value date");
        FIELD(ddeParams,				"DDE parameters");

        Addin::registerClassObjectMethod("EQUITY_SPREAD_CORR_CURVE",
                                         Addin::MARKET,
                                         "derives the vol curve of clean fwd spread given a CDS curve",
                                         TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)getEquitySpreadCorrs);
    }

    static IObject* defaultEquitySpreadCorrCurveAddin(){
        return new EquitySpreadCorrCurveAddin();
    }
    
};
   
CClassConstSP const EquitySpreadCorrCurveAddin::TYPE = CClass::registerClassLoadMethod(
    "EquitySpreadCorrCurveAddin", typeid(EquitySpreadCorrCurveAddin), load);



DRLIB_END_NAMESPACE
