//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : IndexSpecFX.hpp
//
//   Description : specification for FX underlying/payoff index
//
//----------------------------------------------------------------------------

#ifndef INDEX_SPEC_FX_HPP
#define INDEX_SPEC_FX_HPP

#include "edginc/IndexSpec.hpp"
#include "edginc/FXAsset.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL IndexSpecFX : public IndexSpec
{
public:
    static CClassConstSP const TYPE;

    /************************ exported fields ************************/
    FXAssetWrapper factor; // market factor this spec is based on
    // after setup(), factor will be NULL iff sameCurve

    /************************ transient fields ************************/
protected:
    string riskCurve, baseCurve;
public:
    bool sameCurve;

    /*************************** methods ***************************/
    /* IProdCreator:: */
    virtual void setup(const IModel* model, const MarketData* market);
    /* IProdCreator:: */
    virtual double getValue(DateTime date, CashflowInfo &cfi ) const;

    virtual IMarketFactorConstSP getFactor() const { return factor.getSP(); }
    virtual IMarketFactorSP getFactor() { return factor.getSP(); }

    IndexSpecFX(const string &name, const FXAssetWrapper &factor) 
        : IndexSpec(TYPE, name), factor(factor), sameCurve(false) {}

    IndexSpecFX(const string &riskCurve, const string &baseCurve) 
        : IndexSpec(TYPE, ""), riskCurve(riskCurve), 
        baseCurve(baseCurve), sameCurve(riskCurve==baseCurve) {}

    /** IndexSpec override to process static future estimations */
    double pastValue(
	    const DateTime&             sampleDate,
	    const ObservationType*      obsType,
	    const ObservationSource*    source,
	    const FixingType*           fixType,
	    const IObservationOverride* overrides,
	    const SamplingConvention*   sampleRule) const;

protected:
    IndexSpecFX(const CClassConstSP &type = TYPE) : IndexSpec(type), sameCurve(false) {}
    virtual AssetHistoryConstSP getAssetHistory(const string &source) const;

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor() { return new IndexSpecFX(); }
};

typedef smartPtr<IndexSpecFX> IndexSpecFXSP;
typedef smartConstPtr<IndexSpecFX> IndexSpecFXConstSP;
typedef array<IndexSpecFXSP, IndexSpecFX> IndexSpecFXArray;

DRLIB_END_NAMESPACE

#endif
