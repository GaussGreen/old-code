//-----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : CDSIndexParSpreads.hpp
//
//   Description : CDS par spreads for an index
//                 No fields or methods are present at the moment (ie, it is
//                 identical to the parent class CDSParSpreadsBase) but this
//                 class has been created to avoid interface changes if any
//                 such fields or methods are required in the future
//
//   Date        : 2 Sept 2005
//
//-----------------------------------------------------------------------------

#ifndef QLIB_CDSINDEXPARSPREADS_HPP
#define QLIB_CDSINDEXPARSPREADS_HPP

#include "edginc/Object.hpp"
#include "edginc/CDSParSpreadsBase.hpp"

DRLIB_BEGIN_NAMESPACE

/** Holds the current par spreads for index CDSs */
class MARKET_DLL CDSIndexParSpreads: public CDSParSpreadsBase {
public:
    static CClassConstSP const TYPE;

    //possible values for the representation field
    static const string BOOTSTRAPPABLE;
    static const string INDEX_QUOTES;

    virtual ~CDSIndexParSpreads();

    //// checks parameters immediately after object is constructed
    void validatePop2Object();

    //overload parent to validate quoted fees
    virtual void getMarket(const IModel* model, const MarketData *market);

    //overloads CDSParSpreadsBase to validate representation
    virtual DefaultRatesSP defaultRates() const;

    /** Allow setting the par spreads. Used, e.g., to allow interpolating
     * the original par spreads */
    void setParSpreadCurve(ParSpreadCurveWrapper curveWrapper);

    /** Used in the calculation of index basis
     *  calculate the upfront so that par CDS (with quoted fee as upfront) = 0 */
    DoubleArraySP getUpfrontFees() const;

    /** Used in the calculation of index basis
     *  return the contract spreads */
    DoubleArraySP getContractSpreads() const;

private:
    CDSIndexParSpreads();

    static void load (CClassSP& clazz);
    static IObject* defaultCDSIndexParSpreads();
    //allow member functions to override
    void setRepresentation(const string& rep);

    // fields
    DoubleArraySP contractSpreads; //fixed fees specified at index roll
    string        representation;  //what the par spreads actually mean
                                   //BOOTSTRAPPABLE : a term structure, can be used for pricing
                                   //INDEX_QUOTES   : a series of individual quotes, can only be used
                                   //                 for calculating upfronts
};

typedef smartConstPtr<CDSIndexParSpreads> CDSIndexParSpreadsConstSP;
typedef smartPtr<CDSIndexParSpreads>      CDSIndexParSpreadsSP;
typedef MarketWrapper<CDSIndexParSpreads> CDSIndexParSpreadsWrapper;

DRLIB_END_NAMESPACE

#endif
