//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : GenericNAssetAndImnt.hpp
//
//   Description : Base class for N underlying generic instruments that also
//                 have another instrument as an underlying
//                 If you want to book something like a trail fee
//                 in Pyramid 'Generics' this is the mandatory starting point
//
//   Author      : Andrew J Swain
//
//   Date        : 10 April 2002
//
//
//----------------------------------------------------------------------------

#ifndef EDR_GENERICNASSETANDIMNT_HPP
#define EDR_GENERICNASSETANDIMNT_HPP

#include "edginc/Instrument.hpp"
#include "edginc/InstrumentSettlement.hpp"
#include "edginc/Asset.hpp"
#include "edginc/Model.hpp"
#include "edginc/Theta.hpp"

DRLIB_BEGIN_NAMESPACE


/** Base class for N underlying generic instruments that also
    have another instrument as an underlying
    If you want to book something like a trail fee
    in Pyramid 'Generics' this is the mandatory starting point  */
class PRODUCTS_DLL GenericNAssetAndImnt: public CInstrument,
                            virtual public Theta::Shift {
public:
    static CClassConstSP const TYPE;

    /** Implementation of abstract method in Instrument */
    virtual DateTime getValueDate()const;

    /** Validate instrument having aquired market data */
    void Validate();

    /** Get the asset and discount market data */
    void GetMarket(const IModel*          model, 
                   const CMarketDataSP    market);

    /** roll through time (setting historic values) */
    virtual bool sensShift(Theta* theta);

    /** Returns the name of the instrument's discount currency. */
    virtual string discountYieldCurveName() const;


protected:
    GenericNAssetAndImnt(CClassConstSP clazz);
    
    // fields
    DateTime                valueDate;
    InstrumentSettlementSP  instSettle;       // instrument settlement details 
    CAssetWrapperArray      assets;           // the underlyings
    YieldCurveWrapper       discount;         // ccy to discount payoff 
    CInstrumentSP           imnt;             // underlying instrument
    IModelSP                imntModel;        // and how to price it

private:
    friend class GenericNAssetAndImntHelper;

    GenericNAssetAndImnt(); // not implemented
    GenericNAssetAndImnt(const GenericNAssetAndImnt& rhs);
    GenericNAssetAndImnt& operator=(const GenericNAssetAndImnt& rhs);
};

typedef smartConstPtr<GenericNAssetAndImnt> GenericNAssetAndImntConstSP;
typedef smartPtr<GenericNAssetAndImnt> GenericNAssetAndImntSP;

DRLIB_END_NAMESPACE
#endif

