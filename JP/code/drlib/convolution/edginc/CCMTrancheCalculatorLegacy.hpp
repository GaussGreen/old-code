//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Author      : Mark A Robson
//
//   Date        : 23 Dec 2005
//
//----------------------------------------------------------------------------

#ifndef QR_CCMTRANCHECALCULATORLEGACY_HPP
#define QR_CCMTRANCHECALCULATORLEGACY_HPP
#include "edginc/CCMConvolution.hpp"
#include "edginc/TrancheLossCalculatorLegacy.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart)
 * pointers and therefore their include files are not required here - they
 * will be required in the .cpp file though)*/
FORWARD_DECLARE(CreditTrancheLossConfig);
FORWARD_DECLARE(CounterPartyCredit);
class CreditMetricsModel; // aim to remove this
class ConvolutionProduct;
class Control; // aim to remove this

// some implementation remains in CDO.cpp
/** Need to scrap this class. */
class CONVOLUTION_DLL CCMTrancheCalculatorLegacy : 
    public CObject,
    public virtual ITrancheLossCalculatorLegacy
{
public:
    static CClassConstSP const TYPE;

    CCMTrancheCalculatorLegacy();

    /** Builds calculator skipping all variant data ie data that changes
        across time */
    CCMTrancheCalculatorLegacy(
        bool                      doCreditMetrics, /* (I) */
        bool                      doRFL,           /* (I) */
        const bool                recoverNotional, /* (I) */
        const ConvolutionProduct* product,         /* (I) */
        double                    lossUnit);

    CCMTrancheCalculatorLegacy(
        bool                           doCreditMetrics, /* (I) */
        bool                           doRFL,           /* (I) */
        const bool                     recoverNotional, /* (I) */
        CreditTrancheLossConfigConstSP tranche,         /* (I) */
        double                         lossUnit,
        CounterPartyCreditConstSP      counterParty);

    CCMTrancheCalculatorLegacy(CCMTrancheCalculatorLegacy* copy);

    virtual ~CCMTrancheCalculatorLegacy();

    /** This setup is the same as the setup in CDO.cpp but
        takes a ConvolutionProduct (the goal is to scrap the one in
        CDO.cpp and break the dependency between it and CDO). The betaOverride
        must be either empty (=> no override) or the size of the other arrays */
    void setupCM(
        const DoubleArray&        survivalProb, /* name survival proba */
        double                    counterPartyProb, /* name survival proba */
        const DoubleArray&        betaOverride); // optional per name override

    /** Same as setupCM but for CCM rather than Credit Metrics. The betaOverride
        must be either empty (=> no override) or the size of the other arrays */
    void setupCCM(
        const DoubleArray&        survivalProb, /* name survival proba */
        const DoubleArray&        floorSurvivalProb, /* from 'senior' curve */
        double                    counterPartyProb, /* name survival proba */
        const DoubleArray&        betaOverride); // optional per name override

    void convolution(const Control*            control,
                     const CreditMetricsModel* model,
                     const int                 timepoint,
                     const int                 cacheSet);

    void convolutionNoDensityCache(const Control*            control,
                                   const CreditMetricsModel* model,
                                   const int                 timepoint,
                                   const int                 cacheSet);

    CCMConvolutionSP convolutionNoCache();

    virtual void loss(
        double        K1,    /* (I) lower strike      */
        double        K2,    /* (I) upper strike      */
        double       &L,     /* (O) tranche loss amt  */
        double       &Lcond) /* (O) tranche loss amt cond on cpty surviving */
        const;
#if 0
    CCMConvolution::NameParamArray& getBasketInfo();
#endif
    pair<DoubleArraySP, DoubleArraySP> getDensities() const;
    
    int hashCode() const;

    bool equal(CCMTrancheCalculatorLegacy* c );

    /** creates basketInfo needed by RFL model */
    static void createBasketInfoRFL(
        CreditTrancheLossConfigConstSP  tranche,         /* (I) */
        CCMConvolution::NameParamArray& basketInfo);     /* (O) */

    /** creates basketInfo neede by non RFL models */
    static void createBasketInfo(
        CreditTrancheLossConfigConstSP  tranche,         /* (I) */
        CCMConvolution::NameParamArray& basketInfo);     /* (O) */

    /** creates basketInfo field for non RFL. */
    static CCMConvolution::NameParamSP createCptyInfo(
        CounterPartyCreditConstSP counterParty);      /* (I) */

    /** create cpty info needed by RFL model */
    static CCMConvolution::NameParamSP createCptyInfoRFL(
        CounterPartyCreditConstSP counterParty);     /* (I) */

    /** populates existing basketInfo parameter with inputs
        needed by CM */
    static void populateBasketInfoCM(
        CreditTrancheLossConfigConstSP  tranche,         /* (I) */
        const bool                      recoverNotional, /* (I) */
        CCMConvolution::NameParamArray& basketInfo);     // (M)  

    /** populates existing basketInfo parameter with additional inputs
        needed by CCM */
    static void populateBasketInfoCCM(
        CreditTrancheLossConfigConstSP  tranche,      /* (I) */
        CCMConvolution::NameParamArray& basketInfo);  // (M)  

    /** populates cpty field. Same as CDO::getCptyInfo but is on this class
        and operates through ConvolutionProduct and is for CreditMetrics -
        see comments above about splitting up class into 2.*/
    static void populateCptyInfoCM(
        CounterPartyCreditConstSP   counterParty, /* (I) */
        CCMConvolution::NameParamSP cpty);        // (M)  

    /** populates cpty field with additional info needed for CCM - needs 
        refactoring */
    static void populateCptyInfoCCM(
        CounterPartyCreditConstSP   counterParty, /* (I) */
        CCMConvolution::NameParam&  cpty);    // (M)  

    /** populates existing basketInfo parameter with additional inputs
        needed by CCM */
    static void populateBasketInfoRFL(
        CreditTrancheLossConfigConstSP  tranche,      /* (I) */
        CCMConvolution::NameParamArray& basketInfo);  // (M)  

    /** populates cpty field with additional info needed for CCM - needs 
        refactoring */
    static void populateCptyInfoRFL(
        CounterPartyCreditConstSP  counterParty, /* (I) */
        CCMConvolution::NameParam& cpty);       // (M)  

    /** Messy utility method for mapping date from one class to another
        Needs tidying up */
    //static CCMConvolution::NameParamSP createCptyInfoCM(
    //    const ConvolutionProduct* product);  // for counterparty

private:
    static void load(CClassSP& clazz);
    static IObject* defaultCCMTrancheCalculatorLegacy();


    /// fields ////
    double                           pastLoss;
    CCMConvolution::NameParamArray   basketInfo;
    CCMConvolution::NameParamSP      cpty;
    DoubleArray                      ntl;
    DoubleArray                      density;
    DoubleArray                      densityCond;
    int                              maxs;
    int                              maxl;
    double                           lossUnit;
}; 

typedef smartPtr<CCMTrancheCalculatorLegacy> CCMTrancheCalculatorLegacySP; 
typedef smartConstPtr<CCMTrancheCalculatorLegacy> CCMTrancheCalculatorLegacyConstSP; 
typedef array<CCMTrancheCalculatorLegacySP,CCMTrancheCalculatorLegacy> CCMTrancheCalculatorLegacyArray; 



DRLIB_END_NAMESPACE
#endif
