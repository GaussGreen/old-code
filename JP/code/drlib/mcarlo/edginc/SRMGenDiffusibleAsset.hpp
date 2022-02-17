//----------------------------------------------------------------------------
//
//   Group       : QR&D
//
//   Filename    : SRMGenDiffusibleAsset.hpp
//
//   Description : the implementation of generators for creating QMCDiffusibleAsset classes
//
//
//----------------------------------------------------------------------------
#ifndef SRMGenDiffusibleAsset_HPP
#define SRMGenDiffusibleAsset_HPP

#include "edginc/QMCGenDiffusibleAsset.hpp"

#include "edginc/QMCRatesDiffuse.hpp"
// remove SRMHJM* reference as soon as this level has become model-independent
#include "edginc/SRMRatesHJMDiffuse.hpp" 
#include "edginc/SRMRatesLiborDiffuse.hpp"
#include "edginc/SRMRatesDetDiffuse.hpp"
#include "edginc/SRMRatesHJMUtil.hpp" 
#include "edginc/SRMRatesDetUtil.hpp"
#include "edginc/QMCCreditDiffuse.hpp"
// remove SRMCreditHJM* reference as soon as this level has become model-independent 
#include "edginc/SRMCreditHJMDiffuse.hpp"
#include "edginc/SRMCreditHJMUtil.hpp" 
// remove SRMCreditCIR* reference as soon as this level has become model-independent 
#include "edginc/SRMCreditCIRDiffuse.hpp"
#include "edginc/SRMCreditCIRUtil.hpp" 
#include "edginc/QMCCreditCIDJumps.hpp"
#include "edginc/QMCCreditComposite.hpp"
#include "edginc/SRMCreditLiborDiffuse.hpp"
#include "edginc/SRMCreditLiborUtil.hpp" 
#include "edginc/SRMEquityDiffuse.hpp"
#include "edginc/SRMFXDiffuse.hpp"
#include "edginc/SRMEnergyDiffuse.hpp"
#include "edginc/SRMBasisSpreadHJMDiffuse.hpp"
#include "edginc/DECLARE.hpp"
#include "edginc/ISupportPathConfigSRM.hpp"

DRLIB_BEGIN_NAMESPACE


// SRM-specific stuff
class MCARLO_DLL SRMGenDiffusibleFX : public QMCGenDiffusibleFX,
                                      public virtual ISupportPathConfigSRM
{
public:
    SRMGenDiffusibleFX() {}

    // Implement ISupportPathConfigSRM
    virtual void createSRMUtil(MCPathConfigSRM*       mcPathConfig,
                            const DateTime&        today,
                            DateTimeArrayConstSP   simDates);

    // (5.7)
    virtual vector<double> expandAssetAssetCorrs(IQMCGenDiffusibleAsset* asset2, double corr)
    {
        return asset2->expandAssetToFXCorrs(this, corr);
    }
    virtual vector<double> expandAssetToIRCorrs(QMCGenDiffusibleIR* irAsset, double corr);
    virtual vector<double> expandAssetToCRCorrs(QMCGenDiffusibleCredit* crAsset, double corr);
    virtual vector<double> expandAssetToEQCorrs(QMCGenDiffusibleEquity* eqAsset, double corr);
    virtual vector<double> expandAssetToFXCorrs(QMCGenDiffusibleFX* fxAsset, double corr);
    virtual vector<double> expandAssetToBSCorrs(QMCGenDiffusibleBasisSpread* bsAsset, double corr);
    virtual vector<double> expandAssetToENCorrs(QMCGenDiffusibleEnergy* enAsset, double corr);

    virtual void setSRMDiffusibleAsset(
            const DateTime&     today,              // base date of the run
            const MCPathConfigSRM* mcPathConfig,    // to get max/min boundaries, etc
            const IPastValues*  pastValues,         // historic values
            DependenceSP        dependence);        // factorized correlations
    virtual DateTimeArray     getSRMCriticalDates(const MCPathConfigSRM* mcPathConfig,
        const DateTime& start, const DateTime& finish);

    virtual void initializeDiffusibleAssets();

    virtual IQMCDiffusibleFXSP getDiffusibleFX() { return pAssetFX; }
    virtual SRMFXUtilSP       getFXUtil() { return fxUtil; }

    SRMFXUtilSP fxUtil;   // SRMFXUtil for this ccy 
    SRMFXDiffuseSP pAssetFX; // "shell" of a diffusible SRM FX asset
};
DECLARE(SRMGenDiffusibleFX);


class MCARLO_DLL SRMGenDiffusibleIR : public QMCGenDiffusibleIR,
                                      public virtual ISupportPathConfigSRM
{
public:
    ~SRMGenDiffusibleIR() {}
    virtual SRMRatesUtilSP getSRMRatesUtil() = 0;

    // (5.7)
    virtual vector<double> expandAssetAssetCorrs(IQMCGenDiffusibleAsset* asset2, double corr)
    {
        return asset2->expandAssetToIRCorrs(this, corr);
    }
    virtual vector<double> expandAssetToIRCorrs(QMCGenDiffusibleIR* irAsset, double corr);
    virtual vector<double> expandAssetToENCorrs(QMCGenDiffusibleEnergy* enAsset, double corr);
    virtual vector<double> expandAssetToCRCorrs(QMCGenDiffusibleCredit* crAsset, double corr);
    virtual vector<double> expandAssetToEQCorrs(QMCGenDiffusibleEquity* eqAsset, double corr);
    virtual vector<double> expandAssetToFXCorrs(QMCGenDiffusibleFX* fxAsset, double corr);
    virtual vector<double> expandAssetToBSCorrs(QMCGenDiffusibleBasisSpread* bsAsset, double corr);
};
/*
class MCARLO_DLL SRMGenDiffusibleIRHJM : public QMCGenDiffusibleIR,
                                         public virtual ISupportPathConfigSRM
                                         */
class MCARLO_DLL SRMGenDiffusibleIRHJM : public SRMGenDiffusibleIR
{
public:

    ~SRMGenDiffusibleIRHJM() { if (pAssetHJMRates.get()) pAssetHJMRates->detachFX(); }

    // Implement ISupportPathConfigSRM
    virtual void createSRMUtil(MCPathConfigSRM*       mcPathConfig,
                            const DateTime&        today,
                            DateTimeArrayConstSP   simDates);
    virtual void setSRMDiffusibleAsset(
        const DateTime&     today,              // base date of the run
        const MCPathConfigSRM* mcPathConfig,    // to get max/min boundaries, etc
        const IPastValues*  pastValues,         // historic values
        DependenceSP        dependence);        // factorized correlations
    virtual DateTimeArray getSRMCriticalDates(const MCPathConfigSRM* mcPathConfig,
        const DateTime& start, const DateTime& finish);

    virtual void initializeDiffusibleAssets();
    virtual void pushbackICERecord(StateVarDBase& svdbICE, vector<ICE>& vICE);

    virtual IQMCDiffusibleInterestRateSP getDiffusibleIR() { return pAssetHJMRates; }
    virtual QMCRatesUtilSP getRatesUtil() { return ratesHJMUtil; }
    virtual SRMRatesUtilSP getSRMRatesUtil() {return ratesHJMUtil; }

    /** For multi-factor assets, the expansion of the correlations to the multi-factors
    is handles by this method.  The matrix passed in has the correct dimensions
    so just store the results. */
    virtual void expandFactors(SparseDoubleMatrix& outputSparseCorrMatrix ) const;

    // (5.6)
    /** For 1 factor return beta.  For multi-factor, return expanded beta vector 
    of dimension equal to the number of factors */
    virtual vector<double> expandTier2Betas(double beta) const;

    // (5.7)
    /** Just follow the implementations at the SRMGenDiffusibleIR level */
    virtual vector<double> expandAssetAssetCorrs(IQMCGenDiffusibleAsset* asset2, double corr) { return SRMGenDiffusibleIR::expandAssetAssetCorrs(asset2, corr); }
    virtual vector<double> expandAssetToIRCorrs(QMCGenDiffusibleIR* irAsset, double corr) { return SRMGenDiffusibleIR::expandAssetToIRCorrs(irAsset, corr); }
    virtual vector<double> expandAssetToCRCorrs(QMCGenDiffusibleCredit* crAsset, double corr) { return SRMGenDiffusibleIR::expandAssetToCRCorrs(crAsset, corr); }
    virtual vector<double> expandAssetToEQCorrs(QMCGenDiffusibleEquity* eqAsset, double corr) { return SRMGenDiffusibleIR::expandAssetToEQCorrs(eqAsset, corr); }
    virtual vector<double> expandAssetToFXCorrs(QMCGenDiffusibleFX* fxAsset, double corr) { return SRMGenDiffusibleIR::expandAssetToFXCorrs(fxAsset, corr); }
    virtual vector<double> expandAssetToBSCorrs(QMCGenDiffusibleBasisSpread* bsAsset, double corr) { return SRMGenDiffusibleIR::expandAssetToBSCorrs(bsAsset, corr); }
    virtual vector<double> expandAssetToENCorrs(QMCGenDiffusibleEnergy* enAsset, double corr) { return SRMGenDiffusibleIR::expandAssetToENCorrs(enAsset, corr); }

    SRMRatesHJMUtilSP ratesHJMUtil;
    SRMRatesHJMDiffuseSP pAssetHJMRates; // "shell" of a diffusible IR asset
};
DECLARE(SRMGenDiffusibleIRHJM);

/*
class MCARLO_DLL SRMGenDiffusibleIRLibor : public QMCGenDiffusibleIR,
                                           public virtual ISupportPathConfigSRM
                                           */
class MCARLO_DLL SRMGenDiffusibleIRLibor : public SRMGenDiffusibleIR
{
public:

    ~SRMGenDiffusibleIRLibor() { if (pAssetLiborRates.get()) pAssetLiborRates->detachFX(); }

    // Implement ISupportPathConfigSRM
    virtual void createSRMUtil(MCPathConfigSRM*       mcPathConfig,
        const DateTime&        today,
        DateTimeArrayConstSP   simDates);
    virtual void setSRMDiffusibleAsset(
        const DateTime&     today,              // base date of the run
        const MCPathConfigSRM* mcPathConfig,    // to get max/min boundaries, etc
        const IPastValues*  pastValues,         // historic values
        DependenceSP        dependence);        // factorized correlations
    virtual DateTimeArray getSRMCriticalDates(const MCPathConfigSRM* mcPathConfig,
        const DateTime& start, const DateTime& finish);

    virtual void initializeDiffusibleAssets();
    virtual void pushbackICERecord(StateVarDBase& svdbICE, vector<ICE>& vICE);

    virtual IQMCDiffusibleInterestRateSP getDiffusibleIR() { return pAssetLiborRates; }
    virtual QMCRatesUtilSP getRatesUtil() { return ratesLiborUtil; }
    virtual SRMRatesUtilSP getSRMRatesUtil() {return ratesLiborUtil; }

    /** For multi-factor assets, the expansion of the correlations to the multi-factors
    is handles by this method.  The matrix passed in has the correct dimensions
    so just store the results.  */
    virtual void expandFactors(SparseDoubleMatrix& outputSparseCorrMatrix ) const;

    // (5.6)
    /** For 1 factor return beta.  For multi-factor, return expanded beta vector 
    of dimension equal to the number of factors */
    virtual vector<double> expandTier2Betas(double beta) const;

    // (5.7)
    /** Just follow the implementations at the SRMGenDiffusibleIR level */
    virtual vector<double> expandAssetAssetCorrs(IQMCGenDiffusibleAsset* asset2, double corr) { return SRMGenDiffusibleIR::expandAssetAssetCorrs(asset2, corr); }
    virtual vector<double> expandAssetToIRCorrs(QMCGenDiffusibleIR* irAsset, double corr) { return SRMGenDiffusibleIR::expandAssetToIRCorrs(irAsset, corr); }
    virtual vector<double> expandAssetToCRCorrs(QMCGenDiffusibleCredit* crAsset, double corr) { return SRMGenDiffusibleIR::expandAssetToCRCorrs(crAsset, corr); }
    virtual vector<double> expandAssetToEQCorrs(QMCGenDiffusibleEquity* eqAsset, double corr) { return SRMGenDiffusibleIR::expandAssetToEQCorrs(eqAsset, corr); }
    virtual vector<double> expandAssetToFXCorrs(QMCGenDiffusibleFX* fxAsset, double corr) { return SRMGenDiffusibleIR::expandAssetToFXCorrs(fxAsset, corr); }
    virtual vector<double> expandAssetToBSCorrs(QMCGenDiffusibleBasisSpread* bsAsset, double corr) { return SRMGenDiffusibleIR::expandAssetToBSCorrs(bsAsset, corr); }
    virtual vector<double> expandAssetToENCorrs(QMCGenDiffusibleEnergy* enAsset, double corr) { return SRMGenDiffusibleIR::expandAssetToENCorrs(enAsset, corr); }

    SRMRatesLiborUtilSP ratesLiborUtil;
    SRMRatesLiborDiffuseSP pAssetLiborRates; // "shell" of a diffusible IR asset
};
DECLARE(SRMGenDiffusibleIRLibor);


class MCARLO_DLL SRMGenDiffusibleIRDeterm: public QMCGenDiffusibleIR,
                                           public virtual ISupportPathConfigSRM
{
public:

    // Implement ISupportPathConfigSRM
    virtual void createSRMUtil(MCPathConfigSRM*       mcPathConfig,
        const DateTime&        today,
        DateTimeArrayConstSP   simDates);
    /** zero factors so do nothing. */
    virtual void expandFactors(SparseDoubleMatrix& outputSparseCorrMatrix ) const {}

    // (5.6)
    /** Skip since deterministic */
    virtual vector<double> expandTier2Betas(double beta) const;

    // (5.7)
    /** Skip since deterministic */
    virtual vector<double> expandAssetAssetCorrs(IQMCGenDiffusibleAsset* asset2, double corr) { return vector<double>(0); }
    virtual vector<double> expandAssetToIRCorrs(QMCGenDiffusibleIR* irAsset, double corr) { return vector<double>(0); }
    virtual vector<double> expandAssetToCRCorrs(QMCGenDiffusibleCredit* crAsset, double corr) { return vector<double>(0); }
    virtual vector<double> expandAssetToEQCorrs(QMCGenDiffusibleEquity* eqAsset, double corr) { return vector<double>(0); }
    virtual vector<double> expandAssetToFXCorrs(QMCGenDiffusibleFX* fxAsset, double corr) { return vector<double>(0); }
    virtual vector<double> expandAssetToBSCorrs(QMCGenDiffusibleBasisSpread* bsAsset, double corr) { return vector<double>(0); }
    virtual vector<double> expandAssetToENCorrs(QMCGenDiffusibleEnergy* enAsset, double corr) { return vector<double>(0); }

    virtual void setSRMDiffusibleAsset(
        const DateTime&     today,              // base date of the run
        const MCPathConfigSRM* mcPathConfig,    // to get max/min boundaries, etc
        const IPastValues*  pastValues,         // historic values
        DependenceSP        dependence);        // factorized correlations
    virtual DateTimeArray getSRMCriticalDates(const MCPathConfigSRM* mcPathConfig,
        const DateTime& start, const DateTime& finish);

    virtual void initializeDiffusibleAssets();

    virtual IQMCDiffusibleInterestRateSP getDiffusibleIR() { return pAssetDetRates; }
    virtual QMCRatesUtilSP getRatesUtil() { return ratesDetUtil; }

    SRMRatesDetermUtilSP ratesDetUtil;
    SRMRatesDetermDiffuseSP pAssetDetRates; // "shell" of a diffusible IR asset
};
DECLARE(SRMGenDiffusibleIRDeterm);



class MCARLO_DLL SRMGenDiffusibleEquity : public QMCGenDiffusibleEquity,
                                          public virtual ISupportPathConfigSRM 
{
public:

    SRMGenDiffusibleEquity(EquityDiffusionStyle _diffusionStyle) : diffusionStyle(_diffusionStyle) {}
    
    // Implement ISupportPathConfigSRM
    virtual void createSRMUtil(MCPathConfigSRM*       mcPathConfig,
        const DateTime&        today,
        DateTimeArrayConstSP   simDates);

    // (5.7)
    virtual vector<double> expandAssetAssetCorrs(IQMCGenDiffusibleAsset* asset2, double corr)
    {
        return asset2->expandAssetToEQCorrs(this, corr);
    }
    virtual vector<double> expandAssetToIRCorrs(QMCGenDiffusibleIR* irAsset, double corr);
    virtual vector<double> expandAssetToCRCorrs(QMCGenDiffusibleCredit* crAsset, double corr);
    virtual vector<double> expandAssetToEQCorrs(QMCGenDiffusibleEquity* eqAsset, double corr);
    virtual vector<double> expandAssetToFXCorrs(QMCGenDiffusibleFX* fxAsset, double corr);
    virtual vector<double> expandAssetToBSCorrs(QMCGenDiffusibleBasisSpread* bsAsset, double corr);
    virtual vector<double> expandAssetToENCorrs(QMCGenDiffusibleEnergy* enAsset, double corr);

    virtual void setSRMDiffusibleAsset(
        const DateTime&     today,              // base date of the run
        const MCPathConfigSRM* mcPathConfig,    // to get max/min boundaries, etc
        const IPastValues*  pastValues,         // historic values
        DependenceSP        dependence);        // factorized correlations

    virtual DateTimeArray     getSRMCriticalDates(const MCPathConfigSRM* mcPathConfig,
        const DateTime& start, const DateTime& finish);

    virtual void initializeDiffusibleAssets();

    virtual IQMCDiffusibleEQSP getDiffusibleEQ() { return pAssetEquity; }
    virtual SRMEquityUtilSP getEquityUtil() { return equityUtil; }

    virtual void setCurrencyTreatment(CurrencyTreatment ccyTreatment);
private:
    EquityDiffusionStyle diffusionStyle;       // method to use for discretising SDE, e.g. EULER
    
    SRMEquityUtilSP equityUtil;
    SRMEquityDiffuseSP pAssetEquity; // "shell" of a diffusible equity asset
};
DECLARE(SRMGenDiffusibleEquity);


class MCARLO_DLL SRMGenDiffusibleCredit : public QMCGenDiffusibleCredit,
                                          public virtual ISupportPathConfigSRM 
{
public:
    ~SRMGenDiffusibleCredit() {}
    // (5.7)
    virtual vector<double> expandAssetAssetCorrs(IQMCGenDiffusibleAsset* asset2, double corr)
    {
        return asset2->expandAssetToCRCorrs(this, corr);
    }
    virtual vector<double> expandAssetToIRCorrs(QMCGenDiffusibleIR* irAsset, double corr);
    virtual vector<double> expandAssetToCRCorrs(QMCGenDiffusibleCredit* crAsset, double corr);
    virtual vector<double> expandAssetToEQCorrs(QMCGenDiffusibleEquity* eqAsset, double corr);
    virtual vector<double> expandAssetToFXCorrs(QMCGenDiffusibleFX* fxAsset, double corr);
    virtual vector<double> expandAssetToBSCorrs(QMCGenDiffusibleBasisSpread* bsAsset, double corr);
    virtual vector<double> expandAssetToENCorrs(QMCGenDiffusibleEnergy* enAsset, double corr);
};
/*
class MCARLO_DLL SRMGenDiffusibleCreditHJM : public QMCGenDiffusibleCredit,
                                             public virtual ISupportPathConfigSRM 
                                             */
class MCARLO_DLL SRMGenDiffusibleCreditHJM : public SRMGenDiffusibleCredit
{
public:
    // Implement ISupportPathConfigSRM
    virtual void createSRMUtil(MCPathConfigSRM*       mcPathConfig,
        const DateTime&        today,
        DateTimeArrayConstSP   simDates);

    // (5.7)
    /** Just follows implementation of SRMGenDiffusibleCredit */
    virtual vector<double> expandAssetAssetCorrs(IQMCGenDiffusibleAsset* asset2, double corr) { return SRMGenDiffusibleCredit::expandAssetAssetCorrs(asset2, corr); }
    virtual vector<double> expandAssetToIRCorrs(QMCGenDiffusibleIR* irAsset, double corr) { return SRMGenDiffusibleCredit::expandAssetToIRCorrs(irAsset, corr); }
    virtual vector<double> expandAssetToCRCorrs(QMCGenDiffusibleCredit* crAsset, double corr) { return SRMGenDiffusibleCredit::expandAssetToCRCorrs(crAsset, corr); }
    virtual vector<double> expandAssetToEQCorrs(QMCGenDiffusibleEquity* eqAsset, double corr) { return SRMGenDiffusibleCredit::expandAssetToEQCorrs(eqAsset, corr); }
    virtual vector<double> expandAssetToFXCorrs(QMCGenDiffusibleFX* fxAsset, double corr) { return SRMGenDiffusibleCredit::expandAssetToFXCorrs(fxAsset, corr); }
    virtual vector<double> expandAssetToBSCorrs(QMCGenDiffusibleBasisSpread* bsAsset, double corr) { return SRMGenDiffusibleCredit::expandAssetToBSCorrs(bsAsset, corr); }
    virtual vector<double> expandAssetToENCorrs(QMCGenDiffusibleEnergy* enAsset, double corr) { return SRMGenDiffusibleCredit::expandAssetToENCorrs(enAsset, corr); }

    virtual void setSRMDiffusibleAsset(
        const DateTime&     today,              // base date of the run
        const MCPathConfigSRM* mcPathConfig,    // to get max/min boundaries, etc
        const IPastValues*  pastValues,         // historic values
        DependenceSP        dependence);        // factorized correlations

    virtual DateTimeArray     getSRMCriticalDates(const MCPathConfigSRM* mcPathConfig,
        const DateTime& start, const DateTime& finish);

    virtual void initializeDiffusibleAssets();

    virtual IQMCDiffusibleCreditSpreadSP getDiffusibleCredit() { return pAssetCreditHJM; }
    virtual SRMCreditUtilSP getCreditUtil() { return creditHJMUtil; }

    SRMCreditHJMUtilSP creditHJMUtil;
    SRMCreditHJMDiffuseSP pAssetCreditHJM;  // "shell" of a diffusible credit asset
};
DECLARE(SRMGenDiffusibleCreditHJM);

/*
class MCARLO_DLL SRMGenDiffusibleCreditCIR : public QMCGenDiffusibleCredit,
                                             public virtual ISupportPathConfigSRM
                                             */
class MCARLO_DLL SRMGenDiffusibleCreditCIR : public SRMGenDiffusibleCredit
{
public:
    // Implement ISupportPathConfigSRM
    virtual void createSRMUtil(MCPathConfigSRM*       mcPathConfig,
        const DateTime&        today,
        DateTimeArrayConstSP   simDates);

    // (5.7)
    /** Just follows implementation of SRMGenDiffusibleCredit */
    virtual vector<double> expandAssetAssetCorrs(IQMCGenDiffusibleAsset* asset2, double corr) { return SRMGenDiffusibleCredit::expandAssetAssetCorrs(asset2, corr); }
    virtual vector<double> expandAssetToIRCorrs(QMCGenDiffusibleIR* irAsset, double corr) { return SRMGenDiffusibleCredit::expandAssetToIRCorrs(irAsset, corr); }
    virtual vector<double> expandAssetToCRCorrs(QMCGenDiffusibleCredit* crAsset, double corr) { return SRMGenDiffusibleCredit::expandAssetToCRCorrs(crAsset, corr); }
    virtual vector<double> expandAssetToEQCorrs(QMCGenDiffusibleEquity* eqAsset, double corr) { return SRMGenDiffusibleCredit::expandAssetToEQCorrs(eqAsset, corr); }
    virtual vector<double> expandAssetToFXCorrs(QMCGenDiffusibleFX* fxAsset, double corr) { return SRMGenDiffusibleCredit::expandAssetToFXCorrs(fxAsset, corr); }
    virtual vector<double> expandAssetToBSCorrs(QMCGenDiffusibleBasisSpread* bsAsset, double corr) { return SRMGenDiffusibleCredit::expandAssetToBSCorrs(bsAsset, corr); }
    virtual vector<double> expandAssetToENCorrs(QMCGenDiffusibleEnergy* enAsset, double corr) { return SRMGenDiffusibleCredit::expandAssetToENCorrs(enAsset, corr); }

    virtual void setSRMDiffusibleAsset(
        const DateTime&     today,              // base date of the run
        const MCPathConfigSRM* mcPathConfig,    // to get max/min boundaries, etc
        const IPastValues*  pastValues,         // historic values
        DependenceSP        dependence);        // factorized correlations

    virtual DateTimeArray     getSRMCriticalDates(const MCPathConfigSRM* mcPathConfig,
        const DateTime& start, const DateTime& finish);

    virtual void initializeDiffusibleAssets();

    virtual IQMCDiffusibleCreditSpreadSP getDiffusibleCredit() { return pAssetCreditCIR; }
    virtual SRMCreditUtilSP getCreditUtil() { return creditCIRUtil; }

    SRMCreditCIRUtilSP creditCIRUtil;
    SRMCreditCIRDiffuseSP pAssetCreditCIR; // "shell" of a diffusible credit asset
};
DECLARE(SRMGenDiffusibleCreditCIR);

/*
class MCARLO_DLL SRMGenDiffusibleCreditLibor : public QMCGenDiffusibleCredit,
                                               public virtual ISupportPathConfigSRM
                                               */
class MCARLO_DLL SRMGenDiffusibleCreditLibor : public SRMGenDiffusibleCredit
{
public:
    // Implement ISupportPathConfigSRM
    virtual void createSRMUtil(MCPathConfigSRM*       mcPathConfig,
        const DateTime&        today,
        DateTimeArrayConstSP   simDates);

    // (5.7)
    /** Just follows implementation of SRMGenDiffusibleCredit */
    virtual vector<double> expandAssetAssetCorrs(IQMCGenDiffusibleAsset* asset2, double corr) { return SRMGenDiffusibleCredit::expandAssetAssetCorrs(asset2, corr); }
    virtual vector<double> expandAssetToIRCorrs(QMCGenDiffusibleIR* irAsset, double corr) { return SRMGenDiffusibleCredit::expandAssetToIRCorrs(irAsset, corr); }
    virtual vector<double> expandAssetToCRCorrs(QMCGenDiffusibleCredit* crAsset, double corr) { return SRMGenDiffusibleCredit::expandAssetToCRCorrs(crAsset, corr); }
    virtual vector<double> expandAssetToEQCorrs(QMCGenDiffusibleEquity* eqAsset, double corr) { return SRMGenDiffusibleCredit::expandAssetToEQCorrs(eqAsset, corr); }
    virtual vector<double> expandAssetToFXCorrs(QMCGenDiffusibleFX* fxAsset, double corr) { return SRMGenDiffusibleCredit::expandAssetToFXCorrs(fxAsset, corr); }
    virtual vector<double> expandAssetToBSCorrs(QMCGenDiffusibleBasisSpread* bsAsset, double corr) { return SRMGenDiffusibleCredit::expandAssetToBSCorrs(bsAsset, corr); }
    virtual vector<double> expandAssetToENCorrs(QMCGenDiffusibleEnergy* enAsset, double corr) { return SRMGenDiffusibleCredit::expandAssetToENCorrs(enAsset, corr); }

    virtual void setSRMDiffusibleAsset(
        const DateTime&     today,              // base date of the run
        const MCPathConfigSRM* mcPathConfig,    // to get max/min boundaries, etc
        const IPastValues*  pastValues,         // historic values
        DependenceSP        dependence);        // factorized correlations
    virtual DateTimeArray     getSRMCriticalDates(const MCPathConfigSRM* mcPathConfig,
        const DateTime& start, const DateTime& finish);

    virtual void initializeDiffusibleAssets();

    virtual IQMCDiffusibleCreditSpreadSP getDiffusibleCredit() { return pAssetCreditLibor; }
    virtual SRMCreditUtilSP getCreditUtil() { return creditLiborUtil; }

    SRMCreditLiborUtilSP creditLiborUtil;
    SRMCreditLiborDiffuseSP pAssetCreditLibor; // "shell" of a diffusible credit asset
};
DECLARE(SRMGenDiffusibleCreditLibor);


class MCARLO_DLL SRMGenDiffusibleEnergy : public QMCGenDiffusibleEnergy,
                                          public virtual ISupportPathConfigSRM 
{
public:
    // Implement ISupportPathConfigSRM
    virtual void createSRMUtil(MCPathConfigSRM*       mcPathConfig,
        const DateTime&        today,
        DateTimeArrayConstSP   simDates);

    /** Expand correlation matrix for multi-factors. */
    virtual void expandFactors(SparseDoubleMatrix& outputSparseCorrMatrix ) const;

    // (5.6)
    /** For 1 factor return beta.  For multi-factor, return expanded beta vector 
    of dimension equal to the number of factors */
    virtual vector<double> expandTier2Betas(double beta) const;

    // (5.7)
    virtual vector<double> expandAssetAssetCorrs(IQMCGenDiffusibleAsset* asset2, double corr)
    {
        return asset2->expandAssetToENCorrs(this, corr);
    }
    virtual vector<double> expandAssetToIRCorrs(QMCGenDiffusibleIR* irAsset, double corr);
    virtual vector<double> expandAssetToCRCorrs(QMCGenDiffusibleCredit* crAsset, double corr);
    virtual vector<double> expandAssetToEQCorrs(QMCGenDiffusibleEquity* eqAsset, double corr);
    virtual vector<double> expandAssetToFXCorrs(QMCGenDiffusibleFX* fxAsset, double corr);
    virtual vector<double> expandAssetToBSCorrs(QMCGenDiffusibleBasisSpread* bsAsset, double corr);
    virtual vector<double> expandAssetToENCorrs(QMCGenDiffusibleEnergy* enAsset, double corr);

    virtual void setSRMDiffusibleAsset(
        const DateTime&     today,              // base date of the run
        const MCPathConfigSRM* mcPathConfig,    // to get max/min boundaries, etc
        const IPastValues*  pastValues,         // historic values
        DependenceSP        dependence);        // factorized correlations
    virtual DateTimeArray     getSRMCriticalDates(const MCPathConfigSRM* mcPathConfig,
        const DateTime& start, const DateTime& finish);

	virtual void initializeDiffusibleAssets();


	virtual IQMCDiffusibleEnergySP getDiffusibleEnergy() { return pAssetEnrg; }
	virtual SRMEnergyUtilBaseSP getSRMEnrgUtil() { return enrgUtil; }


    QMCEnergyDiffuseSP pAssetEnrg; // shell of a diffusible energy asset
	SRMEnergyUtilBaseSP enrgUtil;
};
DECLARE(SRMGenDiffusibleEnergy);


class MCARLO_DLL SRMGenDiffusibleBasisSpread : public QMCGenDiffusibleBasisSpread,
                                               public virtual ISupportPathConfigSRM 
{
public:
    // Implement ISupportPathConfigSRM
    virtual void createSRMUtil(MCPathConfigSRM*       mcPathConfig,
        const DateTime&        today,
        DateTimeArrayConstSP   simDates);

    // (5.7)
    virtual vector<double> expandAssetAssetCorrs(IQMCGenDiffusibleAsset* asset2, double corr)
    {
        return asset2->expandAssetToBSCorrs(this, corr);
    }
    virtual vector<double> expandAssetToIRCorrs(QMCGenDiffusibleIR* irAsset, double corr);
    virtual vector<double> expandAssetToCRCorrs(QMCGenDiffusibleCredit* crAsset, double corr);
    virtual vector<double> expandAssetToEQCorrs(QMCGenDiffusibleEquity* eqAsset, double corr);
    virtual vector<double> expandAssetToFXCorrs(QMCGenDiffusibleFX* fxAsset, double corr);
    virtual vector<double> expandAssetToBSCorrs(QMCGenDiffusibleBasisSpread* bsAsset, double corr);
    virtual vector<double> expandAssetToENCorrs(QMCGenDiffusibleEnergy* enAsset, double corr);

    virtual void setSRMDiffusibleAsset(
        const DateTime&     today,              // base date of the run
        const MCPathConfigSRM* mcPathConfig,    // to get max/min boundaries, etc
        const IPastValues*  pastValues,         // historic values
        DependenceSP        dependence);        // factorized correlations

    virtual DateTimeArray     getSRMCriticalDates(const MCPathConfigSRM* mcPathConfig,
        const DateTime& start, const DateTime& finish);


	virtual void initializeDiffusibleAssets();


	virtual IQMCDiffusibleBasisIndexSP getDiffusibleBasis() { return pAsset; }
	virtual SRMBasisSpreadUtilSP getSPUtil() { return basisUtil; }


    SRMBasisSpreadHJMDiffuseSP pAsset; // shell of a diffusible asset
	SRMBasisSpreadUtilSP basisUtil;
};
DECLARE(SRMGenDiffusibleBasisSpread);


// CID extension
// TODO: why is QMC stuff here at the SRM level??
class MCARLO_DLL QMCGenDiffusibleCreditCID : public SRMGenDiffusibleCreditCIR,
                                             public virtual ISupportPathConfigSRM 
{
public:
    // Implement ISupportPathConfigSRM
    virtual void createSRMUtil(MCPathConfigSRM*       mcPathConfig,
        const DateTime&        today,
        DateTimeArrayConstSP   simDates);

    // (5.7)
    /** Not implemented yet, dummy implementations at this level.  Future 
    implementations should go to the pathConfig level e.g. SRM or CID */
//    virtual vector<double> expandAssetAssetCorrs(IQMCGenDiffusibleAsset* asset2, double corr) { return vector<double>(0); }
//    virtual vector<double> expandAssetToIRCorrs(QMCGenDiffusibleIR* irAsset, double corr) { return vector<double>(0); }
//    virtual vector<double> expandAssetToCRCorrs(QMCGenDiffusibleCredit* crAsset, double corr) { return vector<double>(0); }
//    virtual vector<double> expandAssetToEQCorrs(QMCGenDiffusibleEquity* eqAsset, double corr) { return vector<double>(0); }
//    virtual vector<double> expandAssetToFXCorrs(QMCGenDiffusibleFX* fxAsset, double corr) { return vector<double>(0); }
//    virtual vector<double> expandAssetToBSCorrs(QMCGenDiffusibleBasisSpread* bsAsset, double corr) { return vector<double>(0); }
//    virtual vector<double> expandAssetToENCorrs(QMCGenDiffusibleEnergy* enAsset, double corr) { return vector<double>(0); }

    virtual void setSRMDiffusibleAsset(
        const DateTime&     today,              // base date of the run
        const MCPathConfigSRM* mcPathConfig,    // to get max/min boundaries, etc
        const IPastValues*  pastValues,         // historic values
        DependenceSP        dependence);        // factorized correlations

    virtual DateTimeArray     getSRMCriticalDates(
        const MCPathConfigSRM* mcPathConfig,
        const DateTime& start, 
        const DateTime& finish) { 
            return SRMGenDiffusibleCreditCIR::getSRMCriticalDates(mcPathConfig,start,finish); }


    virtual void initializeDiffusibleAssets();

    virtual int     numFactors() const // using 1 correlated brownian in Ian's implementation if Full MC
    { 
        return isFullMC ? 1 : 0; 
    }

    QMCGenDiffusibleCreditCID(bool _isCommonFactor) : isCommonFactor(_isCommonFactor), isFullMC(true) {}
    bool    isCommonFactor; // if not - it is idio factor
    bool    isFullMC;

};
DECLARE(QMCGenDiffusibleCreditCID);



DRLIB_END_NAMESPACE
#endif // QMCGenDiffusibleAsset

