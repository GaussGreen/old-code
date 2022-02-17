//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PDFDefaultLNStrike.hpp
//
//   Description : Implementation of PDFCalculator for parameterised vols
//
//   Author      : Mark A Robson
//
//   Date        : 24 March 2002
//
//
//----------------------------------------------------------------------------

#ifndef EDR_PDFPARAMLNSTRIKE_HPP
#define EDR_PDFPARAMLNSTRIKE_HPP

#include "edginc/PDFDefaultLNStrike.hpp"
#include "edginc/VolBase.hpp"
#include "edginc/Asset.hpp"
#include "edginc/PDFRequestLNStrike.hpp"
#include "edginc/VolParam.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL PDFParamLNStrike : public PDFDefaultLNStrike {
public:
    static CClassConstSP const TYPE;

    virtual ~PDFParamLNStrike();

    /** Constructor - note fwd starting not supported currently */
    PDFParamLNStrike(
        const DateTime&                  valueDate,
        const CAssetConstSP&             asset,
        const TimeMetricConstSP&         metric, // to turn vols into vars
        const CVolBaseConstSP&           vol,
        const CVolParamConstSP&          paramVol,
        const PDFRequestLNStrikeConstSP& lnRequest);

protected:
    /** Caclulates spreads using specified lo and hi strikes. Overrides
        default implementation  in PDFDefaultLNStrike to use the 
        ComputeImpVol on the VolParam */
    virtual void spreads(const CLatticeDouble&    loStrikes,
                         const CLatticeDouble&    hiStrikes,
                         const DateTimeArray&     maturities,
                         const DoubleArray&       fwds, // at maturities
                         CLatticeDouble&          spread) const;

private:
    PDFParamLNStrike();
    class Helper;
    friend class Helper;

    // fields
    TimeMetricConstSP           metric;
    CVolBaseConstSP             vol;
    CVolParamConstSP            paramVol;
};

DRLIB_END_NAMESPACE
#endif
