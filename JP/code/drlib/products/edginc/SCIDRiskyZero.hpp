//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : SCIDRiskyZero.hpp
//
//   Description : A pilot pricer for SCID model 
//
//   Author      : Adrian Bozdog
//
//----------------------------------------------------------------------------

#ifndef QR_SCIDRiskyZero_HPP
#define QR_SCIDRiskyZero_HPP

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/SCID.hpp"
#include "edginc/Model.hpp"
#include "edginc/MaturityPeriod.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(SCIDRiskyZero)
FORWARD_DECLARE_WRAPPER(YieldCurve)

class PRODUCTS_DLL SCIDRiskyZero :    public CInstrument, 
                     public virtual SCID::IIntoProduct {
public:
    static CClassConstSP const TYPE;

    virtual ~SCIDRiskyZero();
    /*=========================================================================
     * I/CInstrument Interface
     *=======================================================================*/
    virtual void GetMarket(const IModel* model, const CMarketDataSP);
    virtual void Validate();
    virtual DateTime getValueDate() const;
    virtual string discountYieldCurveName() const;

    /*=========================================================================
     * Pricing model interface implementations
     *=======================================================================*/
    virtual SCID::IProduct* createProduct(SCID* model) const;

protected:
    SCIDRiskyZero();

private:
    SCIDRiskyZero(const SCIDRiskyZero& rhs);
    SCIDRiskyZero& operator=(const SCIDRiskyZero& rhs);

	void GetParSpread(MaturityPeriodSP& freq, DoubleMatrix &modelPS, DoubleMatrix &marketPS, SCIDparametersSP &sCIDparamSP );

    /*=========================================================================
     * DATA FIELDS
     *=======================================================================*/
    /**Trade value date*/
    DateTime        valueDate;
   
    YieldCurveWrapper               discount;        
	DateTimeArray					maturities;
	DayCountConventionSP			dcc;
	DateTime						forwardDate;

    /*=========================================================================
     * FOR REFLECTION
     *=======================================================================*/
    static void load(CClassSP& clazz);
    static IObject* defaultSCIDRiskyZero();

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
    friend class SCIDRiskyZeroHelper;
    friend class SCIDRiskyZeroSCID;
};

DRLIB_END_NAMESPACE
#endif




