//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Generic1FactorCredit.hpp
//
//   Description : Base class for 1 factor generic credit instruments
//
//   Author      : Stephen Hope
//
//   Date        : 27 Sep 2002
//
//
//----------------------------------------------------------------------------

#ifndef EDR_GENERIC_1_FACTOR_CREDIT_HPP
#define EDR_GENERIC_1_FACTOR_CREDIT_HPP

#include "edginc/Class.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/CDSParSpreads.hpp"
#include "edginc/FirmAsset.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/InstrumentSettlement.hpp"
#include "edginc/SensitiveStrikes.hpp"

DRLIB_BEGIN_NAMESPACE


class PRODUCTS_DLL Generic1FactorCredit: public CInstrument,
                            virtual public Theta::Shift,
                            virtual public ISensitiveStrikes{
public:
    static CClassConstSP const TYPE;

    virtual DateTime getValueDate()const;

    /** Get the asset and discount, and par curve market data */
    virtual void GetMarket(const IModel*          model, 
                           const CMarketDataSP    market);
    
    /** Get the sensitive strike for this volRequest.
        Will be called from the product that inherits 
        from this class not from the infrastructure **/
    void getSensStrikes(const OutputNameConstSP&         outputName,
                        const CVolRequest*               volRequest,
                        const SensitiveStrikeDescriptor& sensStrikeDesc,
                        const DoubleArraySP&             sensitiveStrikes);
    
    /** indicates whether VEGA_MATRIX is sensible for this instrument. */
    bool avoidVegaMatrix(const IModel* model);


    /** See comment for avoidVegaMatrix - returns all strikes on the
        vol surface to which this instrument is sensitive */
    DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                      const IModel*      model);

    //// roll through time 
    bool sensShift(Theta* theta);

    /** Returns the name of the instrument's discount currency. */
    virtual string discountYieldCurveName() const;

protected:
    Generic1FactorCredit(CClassConstSP clazz);
    
    /* Returns the cdsParSpreads */
    const ICDSParSpreads* getICDSParSpreads() const;

    /** Do some asset specific validation */
    void validate()const;

    DateTime                valueDate;
    bool                    oneContract;      /* FALSE - fixed notional
                                                 TRUE - one contract */
    double                  notional;         /* Size of deal */
    string                  ccyTreatment;     /* None(vanilla), protected or struck */
    InstrumentSettlementSP  instSettle;       /* Instrument settlement details */
    InstrumentSettlementSP  premiumSettle;    /* When premium is paid */
    CAssetWrapper           asset;            /* The credit underlying */
    ICDSParSpreadsWrapper   cdsParSpreads;    /* The CDS par curve */
    YieldCurveWrapper       discount;         /* Ccy to discount payoff */

private:
    friend class Generic1FactorCreditHelper;

    Generic1FactorCredit(); // not implemented
    Generic1FactorCredit(const Generic1FactorCredit& rhs);
    Generic1FactorCredit& operator=(const Generic1FactorCredit& rhs);
};

typedef smartPtr<Generic1FactorCredit> Generic1FactorCreditSP;

DRLIB_END_NAMESPACE
#endif
