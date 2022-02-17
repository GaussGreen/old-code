//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SRMEnergyDiffuse.hpp
//
//   Description : An implementation for a template (new) asset diffusion
//
//----------------------------------------------------------------------------

#ifndef SRM_ENERGY_DIFFUSE_HPP
#define SRM_ENERGY_DIFFUSE_HPP

#include "edginc/IQMCDiffusibleAsset.hpp"
#include "edginc/EnergyFuturesCurve.hpp"
#include "edginc/QMCRatesDiffuse.hpp"
#include "edginc/QMCEnergyDiffuse.hpp"
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

// forward declare energy util base class and smart pointer
class SRMEnergyUtilBase;
typedef smartPtr<SRMEnergyUtilBase> SRMEnergyUtilBaseSP;
typedef smartPtr<QMCEnergyDiffuse> QMCEnergyDiffuseSP;


/**********************************************************/
/******************* Oil Model Stuffs *********************/
/**********************************************************/
class SRMEnergyOilDiffuse : public QMCEnergyDiffuse, public DiffusionDriver
{
public:
    SRMEnergyOilDiffuse(QMCRatesDiffuseSP srmRatesDiffuse);

    /** finalize the timelines, allocate necessary memory */
    virtual void finalize(DateTimeArrayConstSP allDates);

    /** generate path across all dates. */
    virtual void generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/);

    /** Accessing the expected value ExpFP(md, fp) where md is a
        simulated measurement date and fd is some future date after the
        measurement is made. */
    virtual double  getExpectedPrice(   FwdIdx i/*measurementDateIdx*/,
                                        FwdIdx j/*futureDateIdx*/);

    /** Accessing the natural log of the ExpFP(md, fd)
    where md is a simulated measurement date and fd is some future
    date after the measurement is made. */
    virtual double  getLnExpectedPrice( FwdIdx i/*measurementDateIdx*/,
                                        FwdIdx j/*futureDateIdx*/);

    // read off some const quantities handy for doing diffusion
    // should be called by setQMCEnergyDiffuse(...)
    virtual void setEnrgFxCorr(const vector<double> & _enrgFxCorr);

protected:  // Stuffs that can be shared between the Oil model and any of its decendent, e.g. Gas:
    struct DiffusedState{
        double PHI00;
        double PHI01;
        double PHI11;
        double GAMMA0;
        double GAMMA1;
    };

    struct mu{ // mu's are the quanto adjustments for Gamma
        double mu0;
        double mu1;
    };

    // compute quanto adjustment for the 2-factor oil model
    // note: call this only at the end of oil model's finalize() method
    virtual void calcQuantoAdjustments();

    // Phi's are deterministic, so precompute them
    // note: call this only at the end of oil model's finalize() method
    virtual void preDiffusePhi();

    vector<mu>              MUs;    // quanto adjustments for Gamma
    vector<DiffusedState>   states; // diffusion states for each measurement date [efpRequestedDates.size()]

    double rhoFxEnrg0;  // for quanto adjustment
    double rhoFxEnrg1;

private:

}; // end SRMEnergyOilDiffuse

/**********************************************************/
/******************* Gas Model Stuffs *********************/
/**********************************************************/
class SRMEnergyGasDiffuse : public QMCEnergyDiffuse, public DiffusionDriver
{
public:
    SRMEnergyGasDiffuse(QMCRatesDiffuseSP srmRatesDiffuse);

    /** finalize the timelines, allocate necessary memory */
    virtual void finalize(DateTimeArrayConstSP allDates);

    /** generate path across all dates. */
    virtual void generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/) {};

    /** Accessing the expected value ExpFP(md, fp) where md is a
    simulated measurement date and fd is some future date after the
    measurement is made. */
    virtual double  getExpectedPrice(   FwdIdx i/*measurementDateIdx*/,
        FwdIdx j/*futureDateIdx*/);

    /** Accessing the natural log of the ExpFP(md, fd)
    where md is a simulated measurement date and fd is some future
    date after the measurement is made. */
    virtual double  getLnExpectedPrice( FwdIdx i/*measurementDateIdx*/,
        FwdIdx j/*futureDateIdx*/);

    // read off some const quantities handy for doing diffusion
    // should be called by setQMCEnergyDiffuse(...)
    virtual void setEnrgFxCorr(const vector<double> & _enrgFxCorr);

protected:
    struct DiffusedState{
        double PHI11;
        double PHI22;
        double PHI33;
        double PHI21;
        double PHI31;
        double PHI32;
        double GAMMA1;
        double GAMMA2;
        double GAMMA3;
    };

    struct mu{ // mu's are the quanto adjustments for Gamma
        double mu1;
        double mu2;
        double mu3;
    };

    // compute quanto adjustment
    // note: call this only at the end of oil model's finalize() method
    virtual void calcQuantoAdjustments() {};

    // Phi's are deterministic, so precompute them
    // note: call this only at the end of oil model's finalize() method
    virtual void preDiffusePhi();

    // for diffusion:
    vector<double> k2s; // k-factors correspond to beta [diffusionDates.size() - 1]
    vector<double> kt2s; // cumulative k-factors correspond to beta: [efp.size()]
    vector<double> kT2s; // cumulative k-factors correspond to beta: [efpForwardDates.size()]

    vector<mu>              MUs;    // quanto adjustments for Gamma
    vector<DiffusedState>   states; // diffusion states for each measurement date [efpRequestedDates.size()]

    double rhoFxEnrg0;  // for quanto adjustment
    double rhoFxEnrg1;

private:

}; // end SRMEnergyGasDiffuse


////////////////////////////////////////////////////////////////////////
/////////////////////// Tier 2 Energy Stuffs ///////////////////////////
////////////////////////////////////////////////////////////////////////
class SRMEnergyTier2Diffuse : public QMCEnergyDiffuse
{
public:
    // constructor
    //SRMEnergyTier2Diffuse(QMCRatesDiffuseSP srmRatesDiffuse);
    SRMEnergyTier2Diffuse(QMCRatesDiffuseSP srmRatesDiffuse, QMCEnergyDiffuseSP _parent);

    virtual ~SRMEnergyTier2Diffuse() {}

    /** set up pointer to parent energy */
    //virtual void setParent(QMCEnergyDiffuseSP _parent);

    /** finalize the timelines, allocate necessary memory */
    virtual void finalize(DateTimeArrayConstSP allDates);

    /** generate path across all dates. */
    virtual void generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/);

    /** Accessing the expected value ExpFP(md, fp) where md is a
        simulated measurement date and fd is some future date after the
        measurement is made. */
    virtual double  getExpectedPrice(   FwdIdx i/*measurementDateIdx*/,
                                        FwdIdx j/*futureDateIdx*/);

    /** Accessing the natural log of the ExpFP(md, fd)
    where md is a simulated measurement date and fd is some future
    date after the measurement is made. */
    virtual double  getLnExpectedPrice( FwdIdx i/*measurementDateIdx*/,
                                        FwdIdx j/*futureDateIdx*/);

    // read off some const quantities handy for doing diffusion
    // should be called by setQMCEnergyDiffuse(...)
    virtual void setEnrgFxCorr(const vector<double> & _enrgFxCorr);

    virtual void addAggregatedDates(
        const DateTimeArray& spot,
        const DateTimeArray& measurementDates,
        const DateTimeArray& resetDates);

protected:
    friend class QMCEnergyDiffuse;

    // quanto adjustment for energy spread
    virtual void calcQuantoAdjustments();

    // no Phi for energy spread
    virtual void preDiffusePhi() {};


    QMCEnergyDiffuseSP parent;
    vector<double> GAMMAs; // internal state variable
    vector<double> MUs; // quanto adjustments for GAMMA
    double rhoFxEnrg;   // for quanto adjustment

    // for doing interpolation when the tier2 maturities (fwdEDF dates)
    // do not align with its tier1 parent's maturities.
    vector<int> child2ParentLowerIdxs;
    vector<int> child2ParentUpperIdxs;
    vector<double> parentLowerInterpRatios;
    vector<double> parentUpperInterpRatios;
    
    // for mapping b/w tier2 to tier1 measurement (requested) dates
    vector<int> child2ParentReqIdxs;

}; // end SRMEnergyTier2Diffuse



/** Create the appropriate tier 1 energy model class */
QMCEnergyDiffuseSP SRMEnergyDiffuseCreate(
                                     string modelType,
                                     QMCRatesDiffuseSP srmRatesDiffuse
                                     );

DRLIB_END_NAMESPACE
#endif //SRM_ENERGY_DIFFUSE_HPP
