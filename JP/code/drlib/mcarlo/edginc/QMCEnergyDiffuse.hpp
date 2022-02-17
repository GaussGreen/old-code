//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : QMCEnergyDiffuse.hpp
//
//   Description : An implementation for a template (new) asset diffusion
//
//----------------------------------------------------------------------------

#ifndef QMC_ENERGY_DIFFUSE_HPP
#define QMC_ENERGY_DIFFUSE_HPP

#include "edginc/IQMCDiffusibleAsset.hpp"
#include "edginc/EnergyFuturesCurve.hpp"
#include "edginc/QMCRatesDiffuse.hpp"
#include "edginc/QMCHelperBoundedDiffusion.hpp"

DRLIB_BEGIN_NAMESPACE

#define EnrgOilStr  "OIL"
#define EnrgGasStr  "GAS"
#define EnrgPowStr  "POW"
#define EnrgSampras  "SAMPRAS"

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

// forward declare energy util base class and smart pointer
class SRMEnergyUtilBase;
typedef smartPtr<SRMEnergyUtilBase> SRMEnergyUtilBaseSP;
class QMCEnergyDiffuse;
typedef smartPtr<QMCEnergyDiffuse> QMCEnergyDiffuseSP;

class QMCEnergyDiffuse : public IQMCDiffusibleEnergy
{
public:
    // constructor
    QMCEnergyDiffuse(QMCRatesDiffuseSP srmRatesDiffuse); // NULL must not happen

    void setQMCEnergyDiffuse(
        EnergyFuturesCurveConstSP _futureCurve, // contains all info of the energy future curve (+ vol info)
        SRMEnergyUtilBaseSP _enrgUtil, // smart pointer to energy util class
        int _rndIdx,
        const DateTime & _today,
        const vector<double> & _enrgFxCorr, // for doing quanto adjustment
        const vector<double> & _histPrices // historical forward prices
        );

    virtual ~QMCEnergyDiffuse() {}

    /**   Declaration of pure virtual methods to be overridden by actual models */

    /** (1) methods for engine control */
    /** finalize the timelines, allocate necessary memory */
    virtual void finalize(DateTimeArrayConstSP allDates);

    /** get all the dates important for asset */
//  virtual DateTimeArray getAssetDates();

    /** getting the simulation start date */
    DateTime getBaseDate() { return today; }

    /** just a convenient way to get the forward maturity dates */
// FIXME: masks  getForwardDates() from IQMCDiffusibleAssetBase   DateTimeArray  getForwardDates() const { return getForwardForwardDates(); }

    virtual int numFactors() const {return nbFactors; }

    // read off some const quantities handy for doing diffusion
    // should be called by setQMCEnergyDiffuse(...)
    virtual void setEnrgFxCorr(const vector<double> & _enrgFxCorr) = 0;

//------------------------------------------------------------------------------
    /** this is for efficient access to DiscFactor values, but optional */
    //virtual double* getInternalAssetValueArray() { return &assetValues[0]; }
    /** support for MaxDiffDate/MaxCurveMat */
    virtual QMCHelperBoundedDiffusion *    getDiffusionBound() { return & diffBound;}

protected:

    // compute quanto adjustment
    // note: call this only at the end of the model's finalize() method
    virtual void calcQuantoAdjustments() = 0;

    // Phi's are deterministic, so precompute them
    // note: call this only at the end of the model's finalize() method
    virtual void preDiffusePhi() = 0;

    EnergyFuturesCurveConstSP futureCurve;
    SRMEnergyUtilBaseSP enrgUtil;

    QMCRatesDiffuseSP ratesAsset;
    const vector<double> * sigmasFx; // FX vols for quanto adjuetment, obtained from ratesAsset
    vector<double> histPrices; // historical forward prices
    vector<double> initPrices; // [efpForwardDates.size()]: initial fwd price at maturity dates

    DateTimeArray diffusionDates;
    DateTime today;
    size_t nDiffusionSteps;
    int todayIdx;

    vector<int> efpIndexes; // date indexes for when we save expected future prices
    vector<int> fwdIdx2RequestedIdx; // for mapping efpForwardDates->efpRequestedDates

    int rndIdx;
    int nbFactors;

    // for diffusion:
    vector<double> ks;  // k-factors: for doing diffusion [diffusionDates.size() - 1]
    vector<double> kts; // cumulative k-factors: from today to each measurement dates [efpRequestedDates.size()]
    vector<double> kTs; // cumulative k-factors: from today to each future maturities [efpForwardDates.size()]
    vector<double> dtSqrts;

    QMCEnergyBoundedDiffusion  diffBound;


}; // end QMCEnergyDiffuse

DECLARE(QMCEnergyDiffuse);


DRLIB_END_NAMESPACE
#endif //QMC_ENERGY_DIFFUSE_HPP
