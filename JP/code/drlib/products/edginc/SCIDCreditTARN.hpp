//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : SCIDCreditTARN.hpp
//
//   Description : A pilot pricer for SCID model 
//
//   Author      : Adrian Bozdog
//
//----------------------------------------------------------------------------

#ifndef QR_SCIDCreditTARN_HPP
#define QR_SCIDCreditTARN_HPP

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/SCID.hpp"
#include "edginc/Model.hpp"
#include "edginc/MaturityPeriod.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(SCIDCreditTARN)
FORWARD_DECLARE_WRAPPER(YieldCurve)

class PRODUCTS_DLL SCIDCreditTARN :    public CInstrument, 
                     public virtual SCID::IIntoProduct {
public:
    static CClassConstSP const TYPE;

    virtual ~SCIDCreditTARN();
    /*=========================================================================
     * I/CInstrument Interface
     *=======================================================================*/
    virtual void GetMarket(const IModel* model, const CMarketDataSP);
    virtual void Validate();
	virtual DateTime getValueDate() const { return valueDate; } ;
	virtual string discountYieldCurveName() const { return discount.getName(); };

    /*=========================================================================
     * Pricing model interface implementations
     *=======================================================================*/
    virtual SCID::IProduct* createProduct(SCID* model) const;

protected:
    SCIDCreditTARN();

private:
    SCIDCreditTARN(const SCIDCreditTARN& rhs);
    SCIDCreditTARN& operator=(const SCIDCreditTARN& rhs);

    /*=========================================================================
     * DATA FIELDS
     *=======================================================================*/
    /**Trade value date*/
    DateTime					    valueDate;
    YieldCurveWrapper               discount;        
	int								nbDefaultedNames;
	IntArray						Triggers;
	DateTimeArray					TriggerDates;

    /*=========================================================================
     * FOR REFLECTION
     *=======================================================================*/
    static void load(CClassSP& clazz);
    static IObject* defaultSCIDCreditTARN();

    /*=========================================================================
     * OTHER METHODS
     *=======================================================================*/
    /**Returns the closed form price of the instrument*/
    void priceClosedForm(
        CResults* results, 
        Control* control, 
        SCID* model);

    /*=========================================================================
     * FRIENDS
     *=======================================================================*/
    friend class SCIDCreditTARNHelper;
    friend class SCIDCreditTARNSCID;
};

DRLIB_END_NAMESPACE
#endif




