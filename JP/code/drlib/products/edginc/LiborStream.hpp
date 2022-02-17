//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : LiborStream.hpp
//
//   Description : LiborStream instrument
//
//   Author      : Stephen Hope
//
//   Date        : 24 May 2001
//
//
//----------------------------------------------------------------------------

#ifndef LIBORSTREAM_HPP
#define LIBORSTREAM_HPP
#include "edginc/Instrument.hpp"
#include "edginc/ClosedForm.hpp"
#include "edginc/FloatRate.hpp"
#include "edginc/PayStream.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/Theta.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/LiborStreamCreditSupport.hpp"
#include "edginc/ITaxableInst.hpp"
#include "edginc/FXAsset.hpp"

DRLIB_BEGIN_NAMESPACE

/** LiborStream instrument - present value of future cash flows */

class PRODUCTS_DLL LiborStream: public CInstrument, 
                   public ClosedForm::IIntoProduct,
                   public LastSensDate,                  
                   public Theta::Shift,
                   public CreditSupport::Interface,
                   public ITaxableInst::Basic,
                   public ITaxableInst::WithCoupons {
public:
    static CClassConstSP const TYPE;
    friend class LiborStreamHelper;
    friend class LiborStreamClosedForm;
    friend class LiborStreamCreditSupport;

    /** copy market data relevant to the instrument */
    virtual void GetMarket(const IModel*, const CMarketDataSP);

    /** instrument validation */
    virtual void Validate();

    /** Implementation of ClosedForm::IntoProduct interface */
    virtual ClosedForm::IProduct* createProduct(ClosedForm* model) const;
    
    /** Shifts the object using given shift. */
    virtual bool sensShift(Theta* shift);

    /** what's today ? */
    virtual DateTime getValueDate() const;
    
    /** when to stop tweaking */
    virtual DateTime endDate(const Sensitivity* sensControl) const;
    
    virtual CreditSupportSP createCreditSupport(CMarketDataSP market);

    /** Support for ITaxableInstBasic */
    virtual const DateTime getFinalPaymentDate() const;

    /** Support for ITaxableInstWithCoupons */
    virtual CashFlowArrayConstSP getCoupons() const;

    /** Returns the name of the instrument's discount currency. */
    virtual string discountYieldCurveName() const;

private:
    LiborStream();
    LiborStream(const LiborStream& rhs);
    LiborStream& operator=(const LiborStream& rhs);

    void price(Control* control, CResults* results) const;
    void requests(Control*        control, 
                  CResults*       results,
                  bool            isCallable,
                  const DateTime& cutoff,
                  double          callFwdRate) const;

    DateTime            valueDate;
    bool                isCallable;
    DateTime            callDate;
    int                 callDateOffset;
    YieldCurveWrapper   discount;
    PayStreamSP         payStream;
    CashFlowArraySP     fees;
	string				ccyTreatment;		// optional

	YieldCurveWrapper	underlyingCurve;	// optional
	FXAssetSP			fxAsset;			// transient and optional
};

typedef smartPtr<LiborStream> LiborStreamSP;

DRLIB_END_NAMESPACE
#endif
