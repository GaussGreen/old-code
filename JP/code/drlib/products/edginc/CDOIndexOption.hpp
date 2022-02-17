//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : Index Option.hpp
//
//   Description : A pilot pricer for SCID model 
//
//   Author      : Adrian Bozdog
//
//----------------------------------------------------------------------------

#ifndef QR_CDOINDEXOPTION_HPP
#define QR_CDOINDEXOPTION_HPP

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/SCID.hpp"
#include "edginc/Model.hpp"
#include "edginc/MaturityPeriod.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(CDOIndexOption)
FORWARD_DECLARE_WRAPPER(YieldCurve)

class PRODUCTS_DLL CDOIndexOption :    public CInstrument, 
                     public virtual SCID::IIntoProduct {
public:
    static CClassConstSP const TYPE;

    virtual ~CDOIndexOption();
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
    
    /*=========================================================================
     * Instrument specific functions
     *=======================================================================*/
    //specifically designed for TrancheOption which inherits from this class
    virtual void setFastMCParameters(SCIDparametersSP &sCIDparamSP){}
    //TODO:we need to add SCIDparameters into the model. Thus, the model parameter will be not needed here for TrancheOptions as it should be 
    virtual double getLossAndLegs(long expiryIndexInSim, DateTimeArray & simDates, DateTimeArray& maturities, DateTime& expiry, DoubleArray& RA, DoubleArray& DL, SCID * model); 
    double payoff(double dl, double rd, double loss, double strike);

protected:
	CDOIndexOption(CClassConstSP clazz = TYPE);
    /*=========================================================================
     * OTHER METHODS
     *=======================================================================*/
    /**Returns the closed form price of the instrument*/
    void priceClosedForm(
        CResults* results, 
        Control* control, 
        SCID* model);

    /*=========================================================================
     * DATA FIELDS
     *=======================================================================*/
    /**Trade value date*/
    DateTime						valueDate;
    
    YieldCurveWrapper               discount;        
	int								coupon;
	DayCountConventionSP			dcc;

    /** Expiries for tranche maturities**/
    ExpiryArrayConstSP  trancheMaturities;
    
	DateTimeArray					expiries;
	DoubleArray				    	strikes;
    bool                            flagLoss;

private:
    CDOIndexOption(const CDOIndexOption& rhs);
    CDOIndexOption& operator=(const CDOIndexOption& rhs);

    /*=========================================================================
     * FRIENDS
     *=======================================================================*/
    /** for reflection **/
    friend class CDOIndexOptionHelper;
    friend class CDOIndexOptionSCID;
};

DRLIB_END_NAMESPACE
#endif




