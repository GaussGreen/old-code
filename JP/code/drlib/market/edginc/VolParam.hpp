//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolParam.hpp
//
//   Description : Interface for Parametrized Vol Surface
//
//   Author      : Regis Guichard
//
//   Date        : 02 Mai 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDR_VOL_PARAM_HPP
#define EDR_VOL_PARAM_HPP

#include "edginc/Lattice.hpp"
#include "edginc/FXAsset.hpp"
#include "edginc/Correlation.hpp"
#include "edginc/VolSurface.hpp"
#include "edginc/VolRequestDVF.hpp"
#include "edginc/Asset.hpp"

DRLIB_BEGIN_NAMESPACE

struct MARKET_DLL SImpV{
    double vol;
    double d_dT;
    double d_dStrike;
    double d2_dStrike2;
};
typedef vector<SImpV> SImpVArray;
typedef refCountPtr<SImpVArray> SImpVArraySP;


typedef CSlice<SImpV, SImpVArray, SImpVArraySP> ImpVSlice;
typedef CLattice<SImpV, SImpVArray, SImpVArraySP> ImpVSurf;


/* Originally a typical CVolParam (ie a derived instance) would
    contain the relevant vol. This has since been changed so the vol
    is passed as a parameter in the ComputeImpVol and
    spotVolSurfaceFromStrikes methods. The problem with the original
    is that this results in, for example, the CVolBaseParamSurface
    containing a CVolParam which has a reference back to the original
    CVolBaseParamSurface. This ran into memory problems (whilst
    cloning/shifting) when the top vol went out of scope. This was
    fixed by the CVolParam taking a copy of the
    CVolBaseParamSurface. However, this produced hideous problems when
    we needed to tweak the backbone for theta - two problems: (i)
    drilling down to the backbone (fixable) and (ii) how to clone the
    chain of CVolParam since we would no longer be able to just copy
    the reference in the clone method on CVolBaseParamSurface */

class MARKET_DLL CVolParam: public CObject{
public:
    /** Currently object is derived from CObject which may be useful if we
        ever want to serialise instances of this class - currently though
        doesn't really add anything. Note that clone() does work properly
        for this class at present */
    static CClassConstSP const TYPE;

    virtual ~CVolParam();

    /** Class that captures additional requirements for calculating
        implied vols for forward starting instruments. The idea was to 
        use this class to cope with any changes to the data needed for
        the calculation. Currently, though, this is only done in one
        place. */
    class MARKET_DLL FwdStart{
    public:
        ~FwdStart();

        /** Constructor for not forward starting */
        explicit FwdStart();

        /** Constructor if forward starting. Uses asset to get spot and
            fwd at start date */
        explicit FwdStart(const DateTime&          valueDate,
                          const DateTime&          startDate,
                          const TimeMetricConstSP& timeMetric,
                          const CAsset*            asset);

        /** Constructor - used info from volRequest and asset etc */
        explicit FwdStart(const DateTime&          valueDate,
                          const CVolRequestDVF*    volRequest,
                          const TimeMetricConstSP& timeMetric,
                          const CAsset*            asset);
 
        /** Constructor with everything specified explicitly */
        explicit FwdStart(const DateTime&          valueDate,
                          const DateTime&          startDate,
                          const TimeMetricConstSP& timeMetric,
                          double                   spot,
                          double                   fwdAtStart);

        /** Should a forward starting adjustment be made */
        bool       isFwdStarting() const;

        /** Only valid if isFwdStarting() is true. Returns the
            instrument's start date */
        const DateTime&  getStartDate() const;
        /** Only valid if isFwdStarting() is true. Returns value date */
        const DateTime&  getValueDate() const;
        /** Only valid if isFwdStarting() is true. Returns the time
            metric */
        const TimeMetricConstSP& getTimeMetric() const;
        /** Only valid if isFwdStarting() is true. Returns the asset's
            (whose vol this is) spot price */
        double     getSpot() const;
        /** Only valid if isFwdStarting() is true. Returns the asset's
            (whose vol this is) spot price at the start date */
        double     getFwdAtStart() const;
    private:
        bool     fwdStarting;
        DateTime          valueDate;
        DateTime          startDate;
        TimeMetricConstSP metric;
        double            spot;
        double            fwdAtStart;

        void checkFwdStarting(const string& method) const;
    };

    /** Given Lattice of Spots and maturity dates along the Lattice axis,
        returns the implied vols.
        The output CLatticeDouble must be of the right shape.
        For performance reasons, this is the one you have to do */
    virtual void ComputeImpVol(const CVolBase*          vol,
                               const CLatticeDouble&    strikes,
                               const DateTimeArray&     maturities,
                               CLatticeDouble&          impV) const = 0;

    /** Given a maturity date and a strike, returns the implied vol.
        The default implementation is just a wrapper around the
        CLattice version. */
    virtual double ComputeImpVol(const CVolBase*  vol,
                                 double           strike,
                                 const DateTime&  maturity) const;

    /** Computes forward starting implied vols (from start date to each
        supplied maturity). Note default implementation
        provided. If strikesRelToSpotAtStart is true then strikes are
        percentages relative to spot at start */
    virtual void computeFwdStartImpVol(
        const CVolBase*          vol,
        const FwdStart&          fwdStart,
        const CLatticeDouble&    strikes,
        bool                     strikesRelToSpotAtStart,
        const DateTimeArray&     maturities,
        CLatticeDouble&          impV) const;

    /** Create a VolSurface along the given strikes */
    virtual VolSurface* spotVolSurfaceFromStrikes(
        const CVolBase*       vol,
        const CDoubleArray&   strikes) const = 0;

protected:
    CVolParam(const CClassConstSP& clazz);

private:
    static void load(CClassSP& clazz);
    CVolParam(const CVolParam &rhs);
    CVolParam& operator=(const CVolParam& rhs);
};

typedef smartPtr<CVolParam> CVolParamSP;
typedef smartConstPtr<CVolParam> CVolParamConstSP;


DRLIB_END_NAMESPACE
#endif 
