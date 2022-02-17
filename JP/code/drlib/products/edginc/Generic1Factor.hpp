//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Generic1Factor.hpp
//
//   Description : Base class for 1 factor generic instruments
//
//   Author      : Stephen Hope
//
//   Date        : 4 Sep 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDR_GENERIC_1_FACTOR_HPP
#define EDR_GENERIC_1_FACTOR_HPP

#include "edginc/Class.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/Asset.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/InstrumentSettlement.hpp"
#include "edginc/SensitiveStrikes.hpp"
#include "edginc/Theta.hpp"
#include "edginc/CreditSupport.hpp"

DRLIB_BEGIN_NAMESPACE

class Generic1FactorCreditSupport;

class PRODUCTS_DLL Generic1Factor: public CInstrument,
                      virtual public Theta::Shift,
                      virtual public ISensitiveStrikes{
public:
    static CClassConstSP const TYPE;

    virtual ~Generic1Factor();

    virtual DateTime getValueDate()const;

    /** Get the asset and discount market data */
    virtual void GetMarket(const IModel*          model, 
                           const CMarketDataSP    market);
    
    /** Get the sensitive strike for this volRequest.
        Will be called from the product that inherits 
        from this class not from the infrastructure **/
    void getSensStrikes(const OutputNameConstSP&         outputName,
                        const CVolRequest*               volRequest,
                        const SensitiveStrikeDescriptor& sensStrikeDesc,
                        const DoubleArraySP&             sensitiveStrikes);
    
    /** indicates whether VEGA_MATRIX is sensible for this instrument. The
        implementation here will only return false if the model is 
        a MonteCarlo one and this instrument implements it together with the
        MC vol normal log interface. Derived classes may well need to
        override this */
    bool avoidVegaMatrix(const IModel* model);


    /** See comment for avoidVegaMatrix - returns all strikes on the
        vol surface to which this instrument is sensitive provided the
        model is a ModelCarlo and the instrument supports
        this. Derived classes may well need to override this */
    DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                      const IModel*      model);

    //// roll through time (adjusts valueDate, initialSpot, and fwdStarting
    //// as appropriate
    bool sensShift(Theta* theta);

    /** Returns the name of the instrument's discount currency. */
    virtual string discountYieldCurveName() const;

protected:
    Generic1Factor(CClassConstSP clazz);

    /** Do some asset specific validation */
    void validate();

    /** Adds the DELAY_PRICE and FWD_AT_MAT if requested */
    void addRequests(Control* control,
                     Results* results,
                     double fairValue,
                     const DateTime& maturityDate)const;
    
    DateTime                valueDate;
    DateTime                startDate;        /* When it starts */
    bool                    fwdStarting;      /* Indicates if fwd starting */
    bool                    oneContract;      /* FALSE - fixed notional
                                                 TRUE - one contract */
    double                  notional;         /* Size of deal */
    double                  initialSpot;      
    InstrumentSettlementSP  instSettle;       /* Instrument settlement details */
    InstrumentSettlementSP  premiumSettle;    /* When premium is paid */

    CAssetWrapper           asset;            /* The underlying */
    YieldCurveWrapper       discount;         /* Ccy to discount payoff */
    string                  ccyTreatment;     /* None(vanilla), protected or struck */

private:
    friend class Generic1FactorHelper;
    friend class Generic1FactorCreditSupport;

    Generic1Factor(); // not implemented
    Generic1Factor(const Generic1Factor& rhs);
    Generic1Factor& operator=(const Generic1Factor& rhs);
};

typedef smartPtr<Generic1Factor> Generic1FactorSP;


// virtual class to access private member of the Generic1Factor for credit support
class PRODUCTS_DLL Generic1FactorCreditSupport : public CreditSupport
{
public:
    /** return asset */
    virtual CAssetSP getAsset() const;

    /** return instrument ccy ISO code */
    virtual string getInstCcyCode() const;

    virtual Generic1Factor *getInst() const = 0;

protected:
    Generic1FactorCreditSupport();
};


DRLIB_END_NAMESPACE
#endif

