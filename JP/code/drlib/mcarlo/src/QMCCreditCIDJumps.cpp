//----------------------------------------------------------------------------
//
//   Group       : QR Cross Asset
//
//   Filename    : QMCCreditCIDJumps.cpp
//
//   Description : Jumpy Asset of SRM CID credit path generation
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/QMCCreditCIDJumps.hpp"
#include "edginc/SRMUtil.hpp"
#include "edginc/IQMCRNGManager.hpp"
#include "edginc/Maths.hpp"
DRLIB_BEGIN_NAMESPACE



void QMCCreditCIDJumps::setQMCCreditCIDJumps(
        DateTime        _today,
        const vector<double>& _jumpImpact,
        const vector<double>& _jumpMeanSize,
        const vector<double>& _jumpMeanReversionSpeed)
{
    // initialize all internal parameters
    // first verify that all the vectors are of the same size
    QLIB_VERIFY( _jumpImpact.size()== _jumpMeanSize.size() && 
        _jumpMeanReversionSpeed.size()==_jumpMeanSize.size(), 
        "inconsistent nb of impacts, mean sizes, and mean reversion speed");

    STOREJUMPS = true;

    today = _today;
    jumpImpact = _jumpImpact;
    jumpMeanSize = _jumpMeanSize;
    jumpMeanReversionSpeed = _jumpMeanReversionSpeed ;
    lambda_tstar.resize(jumpImpact.size());
    if (STOREJUMPS)
    {
        offsets = vector<size_t>(_jumpImpact.size()+1, 0);
        jumpSizes.resize(1); // it will grow in genPath when necessary
    }
}

/** finalize the timelines, allocate necessary memory */
void QMCCreditCIDJumps::finalize(DateTimeArrayConstSP allDates)
{
    // saving the various dates

    DateTimeArray sdfRequestedDates = getTimeLogic()->getDFDates(); //{t}
	sdfRequestedDatesAsDouble.resize(sdfRequestedDates.size());
	for (size_t i=0; i<sdfRequestedDatesAsDouble.size(); i++)
		sdfRequestedDatesAsDouble[i] = today.yearFrac(sdfRequestedDates[i]);

    DateTimeArray esdfForwardDates = getForwardForwardDates(); // {t} + {T}
	esdfForwardDatesAsDouble.resize(esdfForwardDates.size());
    for (int i=0; i<esdfForwardDates.size(); i++)
		esdfForwardDatesAsDouble[i] = today.yearFrac(esdfForwardDates[i]);

	DateTimeArray esdfRequestedDates = getTimeLogic()->getReqEDFDates(); // {t*} 
	esdfRequestedDatesAsDouble.resize(esdfRequestedDates.size());
	for (int i=0; i<esdfRequestedDates.size(); i++)
		esdfRequestedDatesAsDouble[i] = today.yearFrac(esdfRequestedDates[i]);


    for(size_t i=0; i<lambda_tstar.size(); ++i)
        lambda_tstar[i].resize(esdfRequestedDates.size());

    survivalProb.resize(sdfRequestedDates.size()); //the historic data are passed in, resize to accommodate future ones

    if (STOREJUMPS)
    {
        allDatesAsDouble.resize(allDates->size());
        for(int i=0; i<allDates->size(); ++i)
            allDatesAsDouble[i] = today.yearFrac((*allDates)[i]);
    }
}

/** generate path across all dates. */
void QMCCreditCIDJumps::generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/)
{
    
    if (!fullMC)
        return;  // nothing to do -- the results are created on request


	// it means - we are in fullMC
    // jumpDates have already been simulated. The only thing we simulate here are jump Sizes
	// and then we compute lambda_tstar and \int lambda_tstar given these jumpDates and jumpSizes

	// first initialize lambda_tstar
    for (size_t i=0; i<jumpImpact.size(); i++) // loop over the number of Markets
			fill(lambda_tstar[i].begin(),lambda_tstar[i].end(),0.0);


    if (jumpGenerator->getTotalNumberOfJumps()>0)
    {
        fill(survivalProb.begin(),survivalProb.end(),0.0);
		if ((STOREJUMPS)&&(jumpSizes.size()<jumpGenerator->getTotalNumberOfJumps()))
			jumpSizes.resize(jumpGenerator->getTotalNumberOfJumps());  // make sure we have enough space to store the jump sizes
		double actualJumpSize, jumpDate; 
        for (size_t i=0; i<jumpImpact.size(); i++) // loop over the number of Markets
        {
            if (STOREJUMPS) 
    			offsets[i+1]=offsets[i]+jumpGenerator->getNumberOfJumps(i);

            for (size_t j=0; j<jumpGenerator->getNumberOfJumps(i); j++)  // loop over the number of jumps that happened in this Market
            {
                // simulate the size of the jump of our name
                actualJumpSize = getJumpSizeSample(jumpImpact[i], jumpMeanSize[i], rngMgr->getSharedGen()->fetch());
				if (STOREJUMPS) 
                    jumpSizes[j+offsets[i]] = actualJumpSize;
                jumpDate = jumpGenerator->getJumpDateByIdxAsDouble(i,j);

				// calculate and store \lamda_tstar for all {t*} in the lambda_tstar[] array
				size_t firstPostJumpIndex = lower_bound(esdfRequestedDatesAsDouble.begin(),
														esdfRequestedDatesAsDouble.end(),
														jumpDate) - esdfRequestedDatesAsDouble.begin();
				for (size_t k=firstPostJumpIndex; k<esdfRequestedDatesAsDouble.size(); k++)
					lambda_tstar[i][k] += actualJumpSize*exp(-jumpMeanReversionSpeed[i]*(esdfRequestedDatesAsDouble.at(k) - jumpDate));

			    // calculate this jump contribution to survivialProb
				firstPostJumpIndex = lower_bound(sdfRequestedDatesAsDouble.begin(),
														sdfRequestedDatesAsDouble.end(),
														jumpDate) - sdfRequestedDatesAsDouble.begin();
				for (size_t k=firstPostJumpIndex; k<sdfRequestedDatesAsDouble.size(); k++)
						survivalProb[k] += actualJumpSize*SRMUtil::GFACSimple(sdfRequestedDatesAsDouble[k] - jumpDate, jumpMeanReversionSpeed[i]);
            }
        }
		for (size_t k=0; k<sdfRequestedDatesAsDouble.size(); k++)
            survivalProb[k] = exp( - survivalProb[k] );
    } 
    else fill(survivalProb.begin(),survivalProb.end(),1.0);

    //fullMC = true;
}

double QMCCreditCIDJumps::getWholeTimelineLogSurvProb(size_t idx)   
{

	double t=allDatesAsDouble[idx];  
    if (!fullMC)
        return log(calculateConditionalSurvivalDiscFactor(t));


    QLIB_VERIFY(STOREJUMPS, "getLnSurvivalProbWholeTimeline cannot be called unless the jumps are stored (i.e. STOREJUMPS = true)");

	if (jumpGenerator->getTotalNumberOfJumps()==0) return 0.0;

	double logSurvProba = 0.0;
	double actualJumpSize, jumpDate; 
	for (size_t i=0; i<jumpImpact.size(); i++) // loop over the number of Markets
		for (size_t j=0; j<jumpGenerator->getNumberOfJumps(i); j++)  // loop over the number of jumps that happened in this Market
			{
				// retrieve the date of the jump 
				jumpDate = jumpGenerator->getJumpDateByIdxAsDouble(i,j);
				if (jumpDate<t)
				{
					// retrieve the size of the jump 
					actualJumpSize = jumpSizes[j+offsets[i]]; 
			        // calculate this jump contribution to logSurvProba
					logSurvProba -= actualJumpSize*SRMUtil::GFACSimple(t - jumpDate, jumpMeanReversionSpeed[i]);
				}
			}
    return logSurvProba;
}


double QMCCreditCIDJumps::getJumpSizeByIdx(size_t sourceIdx, size_t jumpIdx)
{ 
    QLIB_VERIFY(STOREJUMPS, "cannot be called unless the jumps are stored (i.e. STOREJUMPS = true)");
    QLIB_VERIFY(sourceIdx < offsets.size()-1, "Requested for sourceIdx > numSources");

    return  jumpSizes.at(jumpIdx + offsets.at(sourceIdx)); 
}

/** accessing the diffused SDF, i.e. a non-default probability by a given date */
double QMCCreditCIDJumps::getSurvivalDiscFactor(SpotIdx measurementDateIdx)
{
    if (fullMC)
        return survivalProb[measurementDateIdx];
    else // here we are in fast MC, but not hybrid MC -- will be done later
    {
        double T = sdfRequestedDatesAsDouble[measurementDateIdx];
        return calculateConditionalSurvivalDiscFactor(T);
    }
}

double QMCCreditCIDJumps::calculateConditionalSurvivalDiscFactor(double T)
{
    double sdf = 1.0;
    for(size_t m=0; m<jumpImpact.size(); ++m)
        sdf *= getConditionSurvProba(0.0, T, m, 0.0);
    return sdf;
}



/** Accessing the natural log of the expected value ExpSDF(md, fd)
    where md is a simulated measurement date and fd is some future
    date after the measurement is made. */
double QMCCreditCIDJumps::getLnExpectedSurvivalDiscFactor(
                                    FwdIdx measurementDateIdx,
                                    FwdIdx futureDateIdx)
{
//    QLIB_VERIFY(fullMC, "this function is not available and should not be called in FastMC");



//	size_t index_t = getTimeLogic()->getReqEDFIdx(size_t(measurementDateIdx));
	double tstar = esdfForwardDatesAsDouble.at(measurementDateIdx), T = esdfForwardDatesAsDouble.at(futureDateIdx);
	if (tstar>=T) return 0.0;

    if (!fullMC) // we are in the fast MC
    {
        return log(calculateConditionalSurvivalDiscFactor(T)/calculateConditionalSurvivalDiscFactor(tstar));
    }

	double alpha, beta, lambdat, answer = 0;
	for (size_t i=0; i<jumpImpact.size(); i++)
	{
		lambdat = lambda_tstar[i][getTimeLogic()->getReqEDFIdx(measurementDateIdx)]; // lambda at measurement Date for this given Poisson process
		beta = SRMUtil::GFACSimple(T-tstar, jumpMeanReversionSpeed[i]);
		alpha =  jumpImpact[i] * ( jumpMeanSize[i] * (T-tstar) - log(1 + jumpMeanSize[i]*beta))/ (jumpMeanSize[i] + jumpMeanReversionSpeed[i]) ;
        answer += - alpha * jumpGenerator->getJumpIntensity(i) - lambdat * beta;
    }
    return answer;
}


// DanielNg: need to check the reference date
/** returns log of deterministic surv prob */
double QMCCreditCIDJumps::getOriginalLnSurvivalDiscFactor(SpotIdx measurementDateIdx)
{
    return 0.0;
}

/** returns log of deterministic forward surv prob */
double QMCCreditCIDJumps::getOriginalLnExpectedSurvivalDiscFactor(FwdIdx measurementDateIdx, FwdIdx futureDateIdx)
{
    return 0.0;
}

/** given a uniform (random) variable unif, give the product of a Bernoulli(impact) times Exponential with mean meanSize */
double QMCCreditCIDJumps::getJumpSizeSample(double impact, double meanSize, double unif)
{
    if (unif>impact) return 0;
    return -meanSize*log( unif/impact );
}




double QMCCreditCIDJumps::getConditionSurvProba(
                            double t, 
	                        double bigT, 
	                        size_t market,
	                        double lambda_t)
{
    size_t nbJumps = jumpGenerator->getNumberOfJumps(market);

	double answer=1.0;
	if ( (!Maths::isZero(lambda_t)) && (bigT>t) ) 
	{
		double beta = - SRMUtil::GFACSimple(bigT-t, jumpMeanReversionSpeed[market]);
		answer  = exp(beta*lambda_t);
	}
	for (size_t k=0; k<nbJumps; k++)
	{
		double timeOfJump = jumpGenerator->getJumpDateByIdxAsDouble(market,k);
		if ( (timeOfJump>=t)&&(timeOfJump<bigT) )
		{
			double beta = - SRMUtil::GFACSimple(bigT-timeOfJump, jumpMeanReversionSpeed[market]);
			answer*=1.0-jumpImpact[market]*(1.0 - 1.0/(1.0 - beta*jumpMeanSize[market]));
		}
	}
	return answer;
}



DRLIB_END_NAMESPACE

