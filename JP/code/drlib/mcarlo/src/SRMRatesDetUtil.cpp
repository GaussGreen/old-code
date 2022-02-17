#include "edginc/config.hpp"
#include "edginc/SRMRatesDetUtil.hpp"
#include "edginc/SRMUtil.hpp"
#include "edginc/SRMSwaption.hpp"
#include "edginc/IRCalib.hpp"
#include "edginc/SwapTool.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/Actual365F.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/CriticalDateCollector.hpp"
#include "edginc/MaturityTimePeriod.hpp"

DRLIB_BEGIN_NAMESPACE

SRMRatesDetermUtil::SRMRatesDetermUtil(
        const DateTime&      baseDate,
        int                  numFactors,
        const string&        modelParamsKey,
        const string&        smileParamsKey,
        CVolProcessedSP      _processedVol,
        IYieldCurveConstSP   discYC,
        IYieldCurveConstSP   diffYC,
        bool                 _skipFlag,
        double               _flatVolIr,
        const string&        cutoffChoice,
        double               constCutoffValue,
        const string&        corrSwapStart, // eg 1Y  (offset to yc spot date)
        const string&        corrSwapMat,   // eg 10Y (offset to start)
        const string&        corrSwapDCC,   // eg Act/365F
        const string&        corrSwapFreq) // eg 6M
            : 
        SRMRatesUtil(
        baseDate,
        numFactors,
        modelParamsKey,
        smileParamsKey,
        _processedVol,
        discYC,
        diffYC,
        _skipFlag,
        _flatVolIr,
        cutoffChoice,
        constCutoffValue,
        corrSwapStart, // eg 1Y  (offset to yc spot date)
        corrSwapMat,   // eg 10Y (offset to start)
        corrSwapDCC,   // eg Act/365F
        corrSwapFreq)
{


    qLeft = 1.0; 
    qRight = 1.0;

    fwdShift = 0;

}

/** Constructs and populates SRMRatesUtil object suitable using the 'VF' style
calibration. qLeft and qRight are hard coded to 1, as is alpha. */
SRMRatesDetermUtil::SRMRatesDetermUtil(
                const DateTime&      baseDate,
                FXAssetConstSP       fx, /* optional - if specified will 
                                          * add FX smile dates to timeline */
                const string&        modelParamsKey,
                CVolProcessedSP      processedVol,
                IYieldCurveConstSP   discYC,
                IYieldCurveConstSP   diffYC,
                bool                 skipFlag,
                const string&        _corrSwapStart, // eg 1Y  (offset to today)
                const string&        _corrSwapMat,   // eg 10Y (offset to start)
                const string&        _corrSwapDCC,   // eg Act/365F
                const string&        _corrSwapFreq): // eg 6M
SRMRatesUtil(
    baseDate,
    fx, 
    modelParamsKey,
    processedVol,
    discYC,
    diffYC,
    skipFlag,
    _corrSwapStart, // eg 1Y  (offset to today)
    _corrSwapMat,   // eg 10Y (offset to start)
    _corrSwapDCC,   // eg Act/365F
    _corrSwapFreq)
{

}


SRMRatesDetermUtil::~SRMRatesDetermUtil(void)
{
}

void SRMRatesDetermUtil::setTimeLine(DateTimeArrayConstSP simDates)
{
    static const string method("SRMRatesDetermUtil::setTimeLine");

    if (initialized) 
    {
        if (0 && (simDates->size() != dates->size() || ! DateTime::isSubset(*simDates, *dates)))
            throw ModelException(method, "Re-initialized with a different timeline");
        return;
    }

    initialized = true;

    try 
    {
        dates = simDates;
        calcExtendedTimeLine(); // get 'extended' timeline

        spotVol(0.0);
        
        if (!processedVol){
            // this wasn't in the original plan but it seems you can do stuff
            // even if you don't have the swaption vols
            //spotVol(0.0); // populate SpotVol with 1.0 at each point
            return;
        }
        else 
        {
            if (VolProcessedBSIR::TYPE->isAssignableFrom(processedVol->getClass()))
            {
                VolProcessedBSIRSP processedVol = VolProcessedBSIRSP::dynamicCast(this->processedVol);
                swapFrequency = processedVol->getSwapFrequency();
                swapDCC = processedVol->getSwapDCC();
                // get hold of swaptionExpiries, swapStartDates, swapMatDates, vols
                processedVol->getBMDetails(swaptionExpiries, swapStartDates,
                    swapMatDates, swaptionVols);

                // populate LastDate - first find benchmark on/after last sim date
                int bmIdx = (*dates).back().findUpper(swaptionExpiries);
                if (bmIdx == swaptionExpiries.size())
                {
                    bmIdx--;
                }
                LastDate = swapMatDates[bmIdx].max((*dates).back());
                spotVolBM(skipFlag); // calculate SpotVol and tau
            }
            else
            {
                // what to do here?
            }
        }


    } catch (exception& e){
        throw ModelException(e, method, "For currency "+discYC->getCcy());
    }

    //TODO TODO
}

/** this function returns product of spot vols and forward rates
"interest rate basis point vol". The array returned is of the same
length as the supplied array */
DoubleArray SRMRatesDetermUtil::basisPointVol(
                          const DateTimeArray& dates) const
{
    int numDates = dates.size(); // for ease
    DoubleArray rBpVol(numDates);
    for (int n = 0; n < numDates; n++) {
        rBpVol[n] = 0.0;
    }
    return rBpVol;
}


/** vol param */
double SRMRatesDetermUtil::getAlpha(int factorIdx) const
{
    return 0.;
}

void SRMRatesDetermUtil::getAlpha(vector<double>& alphaByFactor) const
{
}

void SRMRatesDetermUtil::getBeta(vector<double>& betaByFactor) const
{
}

/** Returns beta for specified factor */
double SRMRatesDetermUtil::getBeta(int factorIdx) const
{
    return 0.;
    //return model->factors[factorIdx].beta;
}

/** Returns rho[rhoIndex] where rho are model parameters */
double SRMRatesDetermUtil::getRho(int rhoIndex) const
{
    return 0.;
    //return model->rho[rhoIndex];
}

/** Returns the model rho parameters - length corresponds to off diagonal
elements of a symmetric matrix (size of which is number of factors) */
const vector<double>& SRMRatesDetermUtil::getRho() const
{
    static vector<double> dummy(0);

    return (dummy);

}

/** Calls bFactor for the swap used for correlation purposes */
vector<double> SRMRatesDetermUtil::bFactor() const
{
    return bFactor(corrSwapStart,
        corrSwapMat,
        corrSwapDCC,
        MaturityPeriodSP(new MaturityPeriod(corrSwapFreq)));
}

/** Calls bFactor for the swap given by the specified index */
vector<double> SRMRatesDetermUtil::bFactor(int swapIndex) const
{
    return bFactor(swapStartDates[swapIndex], 
        swapMatDates[swapIndex],
        swapDCC, 
        swapFrequency);
}

vector<double> SRMRatesDetermUtil::bFactor(const DateTime&      swapStart,
                                          const DateTime&      swapMat,
                                          DayCountConventionSP swapDayCC,
                                          MaturityPeriodSP     swapFreq) const
{
    static const string method("SRMRatesDetermUtil::bFactor");
    return vector<double>();
}  

DoubleMatrix SRMRatesDetermUtil::triangulation() const
{
    static const string method("SRMRatesDetermUtil::triangulation");
    // for ease
    int Nbfac = numFactors();   /* Number of factors */
    /* initialise matrix (return value) */
    DoubleMatrix TriangMtx(Nbfac, Nbfac);

    return TriangMtx;
}

DoubleMatrix SRMRatesDetermUtil::get3A() const
{
    return DoubleMatrix();
}


// returns the instantaneous vol by factor in the form (factor, time point)
void SRMRatesDetermUtil::instFactorVol(
                                      vector< vector<double> >& vol,				  
                                      const vector<double>& DeltaTime,  // (I)  passed for convenience 
                                      const DateTimeArray& TPDate,      // (I)
                                      const vector<double>& IrFwdRate,  // (I)  pre-computed for performance
                                      int expiryIndex,                  // (I)  index in TPDate
                                      int fwdMatIndex) const            // (I)  in case, the underlying has mat > expiry
{
    static const string method("SRMRatesDetermUtil::instFactorVol");

    try
    {
        // checks 
        if ( ((int)DeltaTime.size() < TPDate.size() - 1) )
        {
            throw ModelException(method, "DeltaTime vector is too short");
        }

        int NbFac = this->numFactors();
        int NbTP = TPDate.size();

        if ((int)vol.size() < NbFac)
        {
            throw ModelException(method, "Factor dimension of vol is too low");
        }

        return;

    }
    catch (exception& e)
    {
        throw ModelException(e, method, "For currency "+discYC->getCcy());
    }
}


// Calculates a 'factor variance' as used by equity/FX vol calibration 
void SRMRatesDetermUtil::irAssetVariance(
                                        vector<double>& irVariance,        // (O)  variance over [T(i-1),  T(i)] for all i <= N
                                        vector<double>& irAssetCovar,      // (O)  covariance in [T(i-1),  T(i)] for all i <= N
                                        const vector<double>& rhoEqIr,            // (I)
                                        const vector<double>& DeltaTime,   // (I)  passed for convenience 
                                        const DateTimeArray& TPDate,       // (I)  needed for simpleFwdCurve
                                        const vector<double>& IrFwdRate,  // (I)  pre-computed for performance
                                        const vector<int>& periodIndex,   // (I)  indices in TPDate
                                        int fwdMatIndex) const             // (I)  indices in TPDate
{
    static const string method("SRMRatesDetermUtil::irAssetVariance");
    
    try
    {
        // checks 
        QLIB_VERIFY(irVariance.size() == irAssetCovar.size(),
            "irVariance and irAssetCovar have different lengths");
        QLIB_VERIFY(irVariance.size() >= periodIndex.size(),
            "irVariance vector is too short");
        QLIB_VERIFY((int) DeltaTime.size() >= TPDate.size() - 1,
            "DeltaTime vector is too short");

        for (size_t k = 0;k<irVariance.size();k++)
        {
            irVariance[k] = 0.0;
            irAssetCovar[k] = 0.0;
        }
    }
    catch (exception& e)
    {
        throw ModelException(e, method, "For currency "+discYC->getCcy());
    }

}



/** Calibration routine : populates SpotVol and tau array */
void SRMRatesDetermUtil::spotVolBM(bool skipFlag) 
{
    //No vol for deterministic rates!
}


DRLIB_END_NAMESPACE
