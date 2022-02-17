//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FXVol.hpp
//
//   Description : Class to take ordinary equity vol and make it appear as
//                 an fx vol
//
//   Author      : Mark A Robson
//
//   Date        : 5 Dec 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDR_FXVOL_HPP
#define EDR_FXVOL_HPP

#include "edginc/FXVolBase.hpp"
#include "edginc/FXVega.hpp"
#include "edginc/FXVegaPointwise.hpp"
#include "edginc/PDFCalculator.hpp"

DRLIB_BEGIN_NAMESPACE

/** Class to take ordinary equity vol and make it appear as
    an fx vol. The key aspect of which is turning of the equity based vega
    tweaks and turning on the fx vega based tweaks */
class MARKET_DLL FXVol: public FXVolBase,
             virtual public FXVega::IShift,
             virtual public FXVegaPointwise::IShift,
             virtual public IPDFCalculator {
public:
    static CClassConstSP const TYPE;

    /** Returns name of vol */
    virtual string getName() const;

    /** Combines market and instrument data together to give a
        Processed Vol */
    virtual IVolProcessed* getProcessedVol(const CVolRequest* volRequest,
                                           const CAsset*      asset) const;

    /** Combines market and instrument data together to give a
        Processed Vol. Here the processed volatility is a processed
        struck volatility ie it reflects the combination of this
        CVolBase together with the supplied FX asset and the
        correlation between this CVolBase and the vol of the
        FX. */
    virtual IVolProcessed* getProcessedVol(
        const CVolRequest* volRequest,
        const CAsset*      eqAsset,
        const FXAsset*     fxAsset,
        const Correlation* eqFXCorr) const;

    /** populate from market cache - default implementation provided */
    virtual void getMarket(const IModel* model, const MarketData* market);

    virtual PDFCalculator* getPDFCalculator(
        const PDFRequest* request,
        const CAsset*     asset) const;

    virtual ~FXVol();

    //// constructor
    FXVol(const CVolBase* vol);

   /** Returns name identifying vol for fx vega  */
    virtual string sensName(FXVega* shift) const;
    /** Shifts the object using given shift */
    virtual bool sensShift(FXVega* shift);

    /** Returns name identifying vol for fx vega pointwise */
    virtual string sensName(FXVegaPointwise* shift) const;
    /** Returns the array of expiries (ie maturities/benchmark dates) that
        need to be tweaked for this vol */
    virtual ExpiryArrayConstSP sensExpiries(FXVegaPointwise* shift) const;
    /** Shifts the object using given shift */
    virtual bool sensShift(FXVegaPointwise* shift);

private:
    FXVol(const FXVol& rhs);
    FXVol& operator=(const FXVol& rhs);

    FXVol();
    static IObject* defaultFXVol();
    static void load(CClassSP& clazz);

    /// fields
    string          name;    // name of the vol
    CVolBaseWrapper vol;
};


DRLIB_END_NAMESPACE

#endif
