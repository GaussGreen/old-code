//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : QMCEnergyDiffuse.cpp
//
//   Description : An implementation for a template (new) asset diffusion
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/QMCEnergyDiffuse.hpp"
#include "edginc/SRMEnergyUtil.hpp"
#include "edginc/Atomic.hpp"
#include <cassert>

#if 0
#include <fstream> // for logging to file
#endif

DRLIB_BEGIN_NAMESPACE

/************************************************************************
 * These classes form the interface for all the diffusible assets       *
 * they are convenient for a few generic methods acting on all the      *
 * different assets uniformly and for model-independent access to the   *
 * data the models produce                                              *
 ************************************************************************/

/***************************************************************************
 *  The model-independent interfaces to specific diffusible asset types    *
 ***************************************************************************/

/** The model-independent interface to a diffusible template asset :
    This class should be a base class for any sample template diffusion model */

QMCEnergyDiffuse::QMCEnergyDiffuse(QMCRatesDiffuseSP srmRatesDiffuse) : // NULL must not happen
    ratesAsset(srmRatesDiffuse),
    sigmasFx(srmRatesDiffuse->getSigmaFX()),
    rndIdx(-999),
    todayIdx(-999),
    nbFactors(-999),
    diffBound(srmRatesDiffuse->getDiffusionBound())
{}

void QMCEnergyDiffuse::setQMCEnergyDiffuse(
    EnergyFuturesCurveConstSP _futureCurve, // contains all info of the energy future curve (+ vol info)
    SRMEnergyUtilBaseSP _enrgUtil,
    int _rndIdx,
    const DateTime & _today,
    const vector<double> & _enrgFxCorr, // for doing quanto adjustment
    const vector<double> & _histPrices // historical forward prices
    )
{
    futureCurve = _futureCurve;
    enrgUtil = _enrgUtil;
    rndIdx = _rndIdx;
    today = _today;
    nbFactors = _enrgUtil->numFactors();

    setEnrgFxCorr(_enrgFxCorr);
}


/** finalize the timelines, allocate necessary memory */
void QMCEnergyDiffuse::finalize(DateTimeArrayConstSP allDates)
{
    static const string method("QMCEnergyDiffuse::finalize");

    todayIdx = today.find((*allDates));

    DateTimeArray efpRequestedDates = getForwardDates();
    DateTimeArray efpForwardDates   = getForwardForwardDates();
    assert(!efpRequestedDates.empty());
    assert(!efpForwardDates.empty());

    // stop diffusion at the last measurement date and get # of diffusion steps
    diffusionDates = efpRequestedDates.back().getPastDates((*allDates));
    nDiffusionSteps = diffusionDates.size() - todayIdx - 1; 

    // map expected future measurement date idx to diffusion date idx
    assert(DateTime::isSubset((*allDates), efpRequestedDates));
    efpIndexes = DateTime::getIndexes((*allDates), efpRequestedDates); // [Nefp]
    efpIndexes.push_back(diffusionDates.size() + 1); // make life easier - add request for index off the end

    // Projection of union(t_i, T_j) onto {t_i}:
    // created map (projection)  esdfForwardDates->esdfRequestedDates
    assert(DateTime::isSubset(efpForwardDates, efpRequestedDates));
    fwdIdx2RequestedIdx = DateTime::getProjection(efpForwardDates, efpRequestedDates); 

    // precompute delta t^0.5
    dtSqrts.assign(nDiffusionSteps, 0.0);
    SRMEnergyUtilBase::computeSqrtYearFrac(todayIdx, diffusionDates, dtSqrts);

    // precompute k-factors
    double alpha = enrgUtil->getAlpha();
    ks.assign(nDiffusionSteps, 0.0);
    kts.assign(efpRequestedDates.size(), 0.0);
    kTs.assign(efpForwardDates.size(), 0.0);
    SRMEnergyUtilBase::computeKFactors(ks, alpha, dtSqrts);
    SRMEnergyUtilBase::computeCumKFactors(kts, alpha, today, efpRequestedDates);
    SRMEnergyUtilBase::computeCumKFactors(kTs, alpha, today, efpForwardDates);
}


DRLIB_END_NAMESPACE
