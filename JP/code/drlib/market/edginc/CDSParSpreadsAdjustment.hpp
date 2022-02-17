//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : CDSParSpreadsAdjustment.hpp
//
//   Description : Abstract base class for candidate adjustments to par spread curves
//   
//   Author      : Gordon Stephens
//
//   Date        : Nov 2005
//
//----------------------------------------------------------------------------

#ifndef QLIB_CDSPARSPREADSADJUSTMENT_HPP
#define QLIB_CDSPARSPREADSADJUSTMENT_HPP

#include "edginc/CDSParSpreads.hpp"

DRLIB_BEGIN_NAMESPACE

/** A class that represents an adjustment that can be made to
    an underlying par spreads curve. */
class MARKET_DLL CDSParSpreadsAdjustment : public MarketObject
{
    public:
        static CClassConstSP const TYPE;

        /** Validation */
        virtual bool expiriesValid(const ICDSBootstrappableWrapper& cdsToAdjust) const = 0;
    
        /** Returns the expiries themselves */
        virtual ExpiryArrayConstSP getAdjustmentExpiries() const = 0;

        /** Returns the adjusted recovery */
        virtual double getAdjustedRecovery (const ICDSBootstrappableWrapper& cdsToAdjust) const = 0;

        /** Returns the adjusted par spread curve */
        virtual DoubleArrayConstSP getAdjustedParSpreads(const ICDSBootstrappableWrapper& cdsToAdjust) const = 0;

        /** Returns a single adjusted rate, corresponding to expiry */
        virtual double getAdjustedParSpread(const ICDSBootstrappable* cdsToAdjust,
                                            ExpiryConstSP expiry) const = 0; 

        /** Tweaking support */
        virtual const DoubleArrayConstSP adjustMaxTweakSizes(const ICDSBootstrappableWrapper& cdsToAdjust,
                                                             DoubleArrayConstSP maximumAllowed,
                                                             IExpiryRiskPropertyConstSP wrt) const = 0;
    protected:
        CDSParSpreadsAdjustment(const CClassConstSP& clazz);

    private:
        /** Invoked when Class is 'loaded', used for reflection */
        static void load(CClassSP& clazz);
};

typedef smartPtr<CDSParSpreadsAdjustment> CDSParSpreadsAdjustmentSP;
typedef smartConstPtr<CDSParSpreadsAdjustment> CDSParSpreadsAdjustmentConstSP;
typedef MarketWrapper<CDSParSpreadsAdjustment> CDSParSpreadsAdjustmentWrapper;

DRLIB_END_NAMESPACE

#endif
