//----------------------------------------------------------------------------
//
//   Group       : QR Cross Asset
//
//   Filename    : QMCPureJumps.cpp
//
//   Description : Generator of jumps and retrieval of generated data
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/QMCPureJumps.hpp"
#include "edginc/Format.hpp"
#include "edginc/SRMConstants.hpp"
#include "edginc/IQMCRNGManager.hpp"

#include <numeric>
DRLIB_BEGIN_NAMESPACE

void QMCPureJumps::modifyLowHigh(int _low, int _high)
{
    QLIB_VERIFY( _low>=0 && _high>=-1 && (_high==-1 || _high>_low), 
        "inconsistent low and high parameters of strata");
    if (_low!=low)
    {
        low = _low;
        probaZeroLow = 0.0;
        probaLow = exp( - totalIntensity*lastDate);
        for (int i=1; i<=low; ++i) 
        {
            probaZeroLow += probaLow;
            probaLow *= totalIntensity*lastDate/i;
        }
    }
    if (_high!=high)
    {
        high = _high;
        if (high>low)
        {    
            double proba = probaLow;
            probaZeroHigh = probaZeroLow;
            for (int i=low+1; i<=high; ++i) 
            {
                probaZeroHigh += proba;
                proba *= totalIntensity*lastDate/i;
            }
        }
        else probaZeroHigh = 1.0; // we are in the case high == -1, which stands for high = +\infty
    }
}


double QMCPureJumps::getStrataProbability(QMCStrataConstSP strata)
{
    const QMCStrataCIDJumps* pStrata = dynamic_cast<const QMCStrataCIDJumps*> (strata.get());
    if (pStrata == NULL) return 1.0;
        modifyLowHigh(pStrata->getNJumps().first, pStrata->getNJumps().second);
    return probaZeroHigh-probaZeroLow;
}


/** finalize the timelines, allocate necessary memory */
void QMCPureJumps::finalize(DateTimeArrayConstSP allDates)
{
    // should create jumpDates with decent size, like 3*(Tmax-T0)*lambda_total
    lastDate = today.yearFrac(diffusionBound.getMaxDiffDate());

    if (lastDate > 0) // otherwise -- if the last date is in the past - we do not actually use this pure jumps generator 
        jumpDates.resize(int(3*lastDate*totalIntensity)+1); // padded from the beginning, but it will grow in genPath if necessary

    low = -2;
    high = -2; // will be modified

}

/** generate path across all dates. */
void QMCPureJumps::generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP strata)
{
    bool isStratified = false;
    const QMCStrataCIDJumps* pStrata = dynamic_cast<const QMCStrataCIDJumps*> (strata.get());

    if (pStrata != NULL)
    {
        modifyLowHigh(pStrata->getNJumps().first, pStrata->getNJumps().second);
        isStratified = true;
    }

    if (totalIntensity*lastDate>0)  
    {
        double unif;
        if (!isStratified)
        {
            for (size_t i=0; i<intensities.size(); i++)
            {
                unif = rngMgr->getSharedGen()->fetch();
                // populate offsets[i+1] with Njumps[i], offsets[0] must be always 0
                offsets[i+1] = getPoissonSample(intensities[i]*lastDate, unif); // nbJumps of the i^th Poisson process
            }
        }
        else
        {
            unif = rngMgr->getSharedGen()->fetch();

            // let us first get the sum of the number of jumps over all markets
            size_t nbJumpsToDistribute = getPoissonSampleStratified(totalIntensity*lastDate, unif, low, 
                                                             probaLow, probaZeroHigh - probaZeroLow);

            double freqLeft = totalIntensity;
            for (size_t i=0; i<intensities.size(); ++i)
            {
                double alpha = intensities[i] / freqLeft;
                if (alpha>1.0-SRMConstants::SRM_TINY) 
                {
                    offsets[i+1] = nbJumpsToDistribute; 
                    for (;i<intensities.size(); ++i) offsets[i+1]=0; // no more jumps left to distribute
                }
                else
                {
                    unif = rngMgr->getSharedGen()->fetch();
                    offsets[i+1] = getBinomialSample(alpha, nbJumpsToDistribute, unif);
                }
                nbJumpsToDistribute -= offsets[i+1];
                freqLeft -= intensities[i];
            }

            // do the conditional simulation...
        }
        for(size_t i=1; i<offsets.size(); ++i) 
            offsets[i] += offsets[i-1];   // cumulative
        // grow the jumpDates if necessary
        if (jumpDates.size() < offsets.back()) jumpDates.resize(offsets.back());

        // generate times of jumps for each source for each of Njumps
        for (size_t i=0; i<offsets.back(); i++)
            jumpDates.at(i) = rngMgr->getSharedGen()->fetch()*lastDate;

    }
    else fill(offsets.begin(), offsets.end(), 0);
}



/** after the path was generated, accessing the generated list of jump dates */
size_t      QMCPureJumps::getTotalNumberOfJumps()
{ 
    return offsets.back(); 
}

size_t      QMCPureJumps::getNumberOfJumps(size_t sourceIdx)
{ 
    QLIB_VERIFY(sourceIdx < intensities.size(), "Requested for sourceIdx > numSources");
    return offsets.at(sourceIdx+1) - offsets.at(sourceIdx); 
}

DateTime    QMCPureJumps::getJumpDateByIdx(size_t sourceIdx, size_t jumpIdx)
{ 
    double jumpDateInYears =  getJumpDateByIdxAsDouble(sourceIdx, jumpIdx); 
    return getBaseDate().rollDate(int(0.5 + 365 * jumpDateInYears));
}

double    QMCPureJumps::getJumpDateByIdxAsDouble(size_t sourceIdx, size_t jumpIdx)
{ 
    QLIB_VERIFY(sourceIdx < intensities.size(), "Requested for sourceIdx > numSources");
    return  jumpDates.at(jumpIdx + offsets.at(sourceIdx)); 
}

size_t QMCPureJumps::getPoissonSample(double intensity, double unif)
{
    return getPoissonSampleStratified(intensity, unif, 0, exp(-intensity), 1.0);
}

size_t QMCPureJumps::getPoissonSampleStratified(double intensity, double unif, int _low, 
                                                double _probaLow, double _probaInLowHigh)
{
    if (intensity<SRMConstants::SRM_TINY) return _low;
    
    int nbJumps = _low;
    double proba = _probaLow / _probaInLowHigh;
    double sumProba = proba;
    while ((unif>sumProba)&&(sumProba<1-SRMConstants::SRM_TINY)) // second test to make sure we do not stay in this loop for ever
    {
        nbJumps++;
        proba *= intensity / nbJumps;
        sumProba += proba;
    }; 
    return nbJumps;
}



size_t QMCPureJumps::getBinomialSample(double p, size_t N, double unif)
{
    size_t n=0;
    double proba, Sproba;
    proba = Sproba = pow(1.0-p,int(N));
    while ((unif>Sproba)&&(n<N))
    {
        proba *= p*(N-n)/( (1-p)*(n+1) );
        Sproba += proba;
        ++n;
    }
    return n;
}


void QMCPureJumps::setQMCPureJumps(
        const DateTime&       _baseDate,
        const vector<double>& _intensities)
{
    today=_baseDate;
    intensities=_intensities;
    offsets = vector<size_t>(_intensities.size()+1, 0);
    totalIntensity = accumulate(intensities.begin(), intensities.end(), 0.0);
}

/*void QMCPureJumps::setJumpDate(size_t sourceIdx, size_t jumpIdx, double jumpDateInYears)
{
    QLIB_VERIFY(sourceIdx < intensities.size(), "Requested for sourceIdx > numSources");
    jumpDates.at(jumpIdx + offsets.at(sourceIdx)) = jumpDateInYears;
}*/


DRLIB_END_NAMESPACE

