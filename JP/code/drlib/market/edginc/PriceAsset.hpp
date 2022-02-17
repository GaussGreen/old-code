//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PriceAsset.hpp
//
//   Description : PriceAsset interface
//
//   Author      : Mark A Robson
//
//   Date        : 14 May 2004
//
//
//----------------------------------------------------------------------------

#ifndef EDR_PRICEASSET_HPP
#define EDR_PRICEASSET_HPP
#include "edginc/GeneralAsset.hpp"
#include "edginc/OutputName.hpp"

DRLIB_BEGIN_NAMESPACE
class CVolProcessedBS;
class CVolRequestLN;
class PDFCalculator;
class PDFRequest;

/** This class keeps additional information about the sensitive strike
    request. For example a protected component in a basket with vol
    surface should only calculate sensitive strikes for the protected
    forwards. This currently looks a bit OTT, but there is some
    potential for further complications, for example Futures as under-
    lyings, funds, etc. */
class MARKET_DLL SensitiveStrikeDescriptor {
public:
    /** tell sensitiveStrikes method whether to ignore the volRequest
        that's being passed around */
    bool    forwardOnly;

    /** constructor */
    SensitiveStrikeDescriptor(); // in MarketFactor.cpp
};


/** A PriceAsset covers equity/fx type assets (essentially an interface
    version of CAsset) */
class MARKET_DLL IPriceAsset: public virtual IGeneralAsset{
public:
    static CClassConstSP const TYPE; // in MarketFactor.cpp

    IPriceAsset(); // in MarketFactor.cpp

    virtual ~IPriceAsset();

    /** returns the asset name without the effect of any current treatment 
        ie the name that would have been obtained if currency treatment had 
        been CCY_TREATMENT_VANILLA. Default implementation returns getName() */
    virtual string getTrueName() const = 0;

    //// Returns the ccy treatment for the asset
    virtual string getCcyTreatment() const = 0;

    /** Returns an LN processed vol - which combines the vol market
        data with the instrument data in the volRequest */
    virtual CVolProcessedBS* getProcessedVol(
        const CVolRequestLN* volRequest) const = 0;
 
    /** Calculates the expected number of shares at given fwdDate,
        based on $1 notional. By default, related convexity adjustment
        is TRUE */
    virtual double expNumberShares( const DateTime& valueDate, 
                                    const DateTime& fwdDate,
                                    const bool& convAdjust) const = 0;

    /** Returns the sensitive strikes for this asset using the data in
        the vol request */
    virtual void getSensitiveStrikes(
        const CVolRequest*               volRequest,
        OutputNameConstSP                outputName,
        const SensitiveStrikeDescriptor& sensStrikeDesc,
        DoubleArraySP                    sensitiveStrikes) const = 0;

    /** return a pdf calculator provided our vol supports it */
    virtual PDFCalculator* pdfCalculator(const PDFRequest* request) const = 0;

private:
    static void load(CClassSP& clazz); // in MarketFactor.cpp
};

DRLIB_END_NAMESPACE
#endif
