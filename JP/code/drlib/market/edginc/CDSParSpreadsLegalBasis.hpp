//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : CDSParSpreadsLegalBasis.hpp
//
//   Description : Legal Basis adjustment for CDSParSpreads
//   
//   Author      : Antoine Gregoire 
//
//   Date        : Jan 2005
//
//
//----------------------------------------------------------------------------

#ifndef EDR_CDSPARSPREADS_LEGAL_BASIS_HPP
#define EDR_CDSPARSPREADS_LEGAL_BASIS_HPP

#include "edginc/CDSParSpreadsAdjustment.hpp"
#include "edginc/TweakableWith.hpp"
#include "edginc/PointwiseTweakableWith.hpp"
#include "edginc/LegalBasisAdditiveParallelTweak.hpp"
#include "edginc/LegalBasisAbsolutePointwise.hpp"
#include "edginc/LegalBasisMultiplierParallelTweak.hpp"
#include "edginc/LegalBasisRelativePointwise.hpp"
#include "edginc/LegalBasisAdditiveRecoveryTweak.hpp"
#include "edginc/LegalBasisMultiplierRecoveryTweak.hpp"
#include "edginc/LegalBasisAdditiveRelativeParallelTweak.hpp"
#include "edginc/QuasiContractualBaseCorrelation.hpp"

DRLIB_BEGIN_NAMESPACE

/** 
 * This class contains the 'Legal Basis' parameters used
 * to adjust CDS par spreads and recovery (a*X+b transformation).
 * */
class MARKET_DLL CDSParSpreadsLegalBasis :  
        public CDSParSpreadsAdjustment,
        virtual public TweakableWith<LegalBasisAdditiveParallelTweak>,
        virtual public ITweakableWithRespectTo<LegalBasisAbsolutePointwise>,
        virtual public TweakableWith<LegalBasisMultiplierParallelTweak>,
        virtual public ITweakableWithRespectTo<LegalBasisRelativePointwise>,
        virtual public TweakableWith<LegalBasisAdditiveRecoveryTweak>,
        virtual public TweakableWith<LegalBasisMultiplierRecoveryTweak>,
        virtual public ITweakableWithRespectTo<LegalBasisAdditiveRelativeParallelTweak>,
        virtual public QuasiContractualBaseCorrelation::IShift
{
    friend class Dummy; // suppress compiler warning
public:
    static CClassConstSP const TYPE;
    
    virtual ~CDSParSpreadsLegalBasis();

    /** Returns the name of this object. This is the name with which
     * it is stored in the market data cache and is the name with
     * which results (eg tweaks) should be reported against */
    virtual string getName() const;
    
    /** checks parameters immediately after object is constructed */
    void validatePop2Object();
    
    /** Checks if the expiries in the CDS par spread curve are consistent
     *  with the ones in this adjustment */
    bool expiriesValid(const ICDSBootstrappableWrapper& cdsToAdjust) const;

    /** Returns the expiries themselves */
    ExpiryArrayConstSP getAdjustmentExpiries() const;

    /** Returns the adjusted recovery */
    double getAdjustedRecovery (const ICDSBootstrappableWrapper& cdsToAdjust) const;

    /** Returns the adjusted par spread curve */
    DoubleArrayConstSP getAdjustedParSpreads (const ICDSBootstrappableWrapper& cdsToAdjust) const;

    /** Returns a single adjusted rate, corresponding to expiry */
    double getAdjustedParSpread(const ICDSBootstrappable* cdsToAdjust,
                                ExpiryConstSP expiry) const; 

    /** Returns the adjusted maximum tweak sizes */
    const DoubleArrayConstSP adjustMaxTweakSizes(const ICDSBootstrappableWrapper& cdsToAdjust,
                                                 DoubleArrayConstSP calculatedMaxTwkSizes,
                                                 IExpiryRiskPropertyConstSP withRespectTo) const;
    // Additive coeficients - Parallel tweak 
    virtual string sensName(LegalBasisAdditiveParallelTweak* shift) const;
    virtual bool sensShift(LegalBasisAdditiveParallelTweak* shift);

    // Additive coeficients - Pointwise tweak
    virtual string sensName(const LegalBasisAbsolutePointwise*) const;
    virtual TweakOutcome sensShift(const PropertyTweak<LegalBasisAbsolutePointwise>& shift);
    virtual ExpiryWindowArrayConstSP sensQualifiers(const LegalBasisAbsolutePointwise*) const;

    // Multiplier coeficients - Parallel tweak
    virtual string sensName(LegalBasisMultiplierParallelTweak* shift) const;
    virtual bool sensShift(LegalBasisMultiplierParallelTweak* shift);

    // Multiplier coeficients - Pointwise tweak
    virtual string sensName(const LegalBasisRelativePointwise*) const;
    virtual TweakOutcome sensShift(const PropertyTweak<LegalBasisRelativePointwise>& shift);
    virtual ExpiryWindowArrayConstSP sensQualifiers(const LegalBasisRelativePointwise*) const;

    // Additive coefficents - parallel relative tweak
    virtual string sensName(const LegalBasisAdditiveRelativeParallelTweak*) const;
    virtual TweakOutcome sensShift(const 
        PropertyTweak<LegalBasisAdditiveRelativeParallelTweak>& tweak);

    // Additive recovery coeficient - Parallel tweak
    virtual string sensName(LegalBasisAdditiveRecoveryTweak* shift) const;
    virtual bool sensShift(LegalBasisAdditiveRecoveryTweak* shift);

    // Multiplier recovery coeficient - Parallel tweak
    virtual string sensName(LegalBasisMultiplierRecoveryTweak* shift) const;
    virtual bool sensShift(LegalBasisMultiplierRecoveryTweak* shift);

    // Remove the legal basis adjustment, ie, set the additive coefficients to
    // zero and the multiplicative coefficients to 1.
    virtual bool sensShift(QuasiContractualBaseCorrelation* shift);
 
private:
    /** Constructor */
    CDSParSpreadsLegalBasis();
    
    /** Used for reflection */
    static IObject* defaultConstructor();
    
    /** Copy constructor : don't use */  
    CDSParSpreadsLegalBasis(
        const CDSParSpreadsLegalBasis& cdsParSpreadsLegalBasis);

    /** Override '=' operator : don't use */
    CDSParSpreadsLegalBasis& operator=(
        const CDSParSpreadsLegalBasis& cdsParSpreadsLegalBasis);

    /** Invoked when Class is 'loaded', used for reflection */
    static void load(CClassSP& clazz);

// FIELDS    
    /** Name of that market object */
    string name;
    
    /** Expiries array, must match those in the par spread curve*/
    ExpiryArraySP expiries;
    
    /** multiplier coefficients to apply to par spread curve */
    DoubleArray parSpreadsMultiplierCoefficients;
    
    /** additive coefficients to apply to par spread curve */
    DoubleArray parSpreadsAdditiveCoefficients;
    
    /** multiplier coefficient to apply to recovery */
    double recoveryMultiplierCoefficient;
    
    /** additive coefficient to apply to recovery */
    double recoveryAdditiveCoefficient;
};

typedef smartConstPtr<CDSParSpreadsLegalBasis> CDSParSpreadsLegalBasisConstSP;
typedef smartPtr<CDSParSpreadsLegalBasis> CDSParSpreadsLegalBasisSP;
#ifndef QLIB_CDSPARSPREADSLEGALBASIS_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<CDSParSpreadsLegalBasis>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<CDSParSpreadsLegalBasis>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<CDSParSpreadsLegalBasis>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<CDSParSpreadsLegalBasis>);
#endif

// support for wrapper class
typedef MarketWrapper<CDSParSpreadsLegalBasis> CDSParSpreadsLegalBasisWrapper;
#ifndef QLIB_CDSPARSPREADSLEGALBASIS_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<CDSParSpreadsLegalBasis>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<CDSParSpreadsLegalBasis>);
#endif

DRLIB_END_NAMESPACE

#endif
