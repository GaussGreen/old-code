#ifndef QR_TranchePricerCalibrator_HPP
#define QR_TranchePricerCalibrator_HPP

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/SCIDtree.hpp"
#include "edginc/Model.hpp"
#include "edginc/MaturityPeriod.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(TranchePricerCalibrator)
FORWARD_DECLARE(ICDS)
FORWARD_DECLARE_WRAPPER(ICDSParSpreads)
FORWARD_DECLARE_WRAPPER(YieldCurve)

class PRODUCTS_DLL TranchePricerCalibrator :    public CInstrument, 
												public virtual SCIDtree::IIntoProduct
{
public:
    static CClassConstSP const TYPE;

    virtual ~TranchePricerCalibrator();
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
	virtual SCIDtree::IProduct* createProduct(SCIDtree* model) const;

protected:
    TranchePricerCalibrator();
 
private:
    TranchePricerCalibrator(const TranchePricerCalibrator& rhs);
    TranchePricerCalibrator& operator=(const TranchePricerCalibrator& rhs);

	void priceClosedForm(CResults* results, Control* control, SCIDtree* model);

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
	DoubleMatrix					calibrateParSpread;  
	double							worldDx;
	DoubleArray						smoothness;


    /*=========================================================================
     * FRIENDS
     *=======================================================================*/
    /** for reflection **/
    friend class TranchePricerCalibratorHelper;
	friend class TranchePricerCalibratorSCIDtree;
};

DRLIB_END_NAMESPACE
#endif




