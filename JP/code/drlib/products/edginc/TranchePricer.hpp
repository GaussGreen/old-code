//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : TranchePricer.hpp
//
//   Description : A pilot pricer for SCID model 
//
//   Author      : Adrian Bozdog
//
//----------------------------------------------------------------------------

#ifndef QR_TRANCHEPRICER_HPP
#define QR_TRANCHEPRICER_HPP

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/SCID.hpp"
#include "edginc/SCIDtree.hpp"
#include "edginc/Model.hpp"
#include "edginc/MaturityPeriod.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(TranchePricer)
FORWARD_DECLARE(ICDS)
FORWARD_DECLARE_WRAPPER(ICDSParSpreads)
FORWARD_DECLARE_WRAPPER(YieldCurve)

class PRODUCTS_DLL TranchePricer :    public CInstrument, 
					 public virtual SCIDtree::IIntoProduct,
//					 public virtual SCID2::IIntoProduct,
                     public virtual SCID::IIntoProduct {
public:
    static CClassConstSP const TYPE;

    virtual ~TranchePricer();
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
//	virtual SCID2::IProduct* createProduct(SCID2* model) const;
	virtual SCIDtree::IProduct* createProduct(SCIDtree* model) const;

protected:
    TranchePricer();

private:
    TranchePricer(const TranchePricer& rhs);
    TranchePricer& operator=(const TranchePricer& rhs);

    /*=========================================================================
     * DATA FIELDS
     *=======================================================================*/
    /**Trade value date*/
    DateTime						valueDate;
    
    YieldCurveWrapper               discount;        
	int								coupon;
	DayCountConventionSP			dcc;

	DoubleArray						kmin, kmax;
	DateTimeArray					maturities;


	bool							forwardPricer;
	DateTime						forwardDate;
	

    /*=========================================================================
     * FOR REFLECTION
     *=======================================================================*/
    //static void load(CClassSP& clazz);
    //static IObject* defaultTranchePricer();

    /*=========================================================================
     * OTHER METHODS
     *=======================================================================*/
    /**Returns the closed form price of the instrument*/
    void priceClosedForm(CResults* results, Control* control, SCID* model);
//	void priceClosedForm(CResults* results, Control* control, SCID2* model);
	void priceClosedForm(CResults* results, Control* control, SCIDtree* model);

    /*=========================================================================
     * FRIENDS
     *=======================================================================*/
    /** for reflection **/
    friend class TranchePricerHelper;
    friend class TranchePricerSCID;
	friend class TranchePricerSCIDtree;
//	friend class TranchePricerSCID2;

};

DRLIB_END_NAMESPACE
#endif




