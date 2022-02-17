//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : AUD3YBondFuture.hpp
//
//   Description : AUD3YBondFuture instrument
//
//   Author      : Keiji Kitazawa
//
//   Date        : 27 Jun 2002
//
//
//----------------------------------------------------------------------------

#ifndef AUD3YBondFuture_HPP
#define AUD3YBondFuture_HPP
#include "edginc/ClosedForm.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/BondParams.hpp"
#include "edginc/GenericSimpleIR.hpp"

DRLIB_BEGIN_NAMESPACE

/** AUD3YBondFuture instrument - present value of future cash flows */

class PRODUCTS_DLL AUD3YBondFuture: public GenericSimpleIR, 
                   public ClosedForm::IIntoProduct,
                   public LastSensDate{
public:
    static CClassConstSP const TYPE;
    friend class AUD3YBondFutureHelper;
    friend class AUD3YBondFutureClosedForm;

    /** copy market data relevant to the instrument */
    virtual void GetMarket(const IModel*, const CMarketDataSP);

    /** instrument validation */
    virtual void Validate();

    /** Implementation of ClosedForm::IntoProduct interface */
    virtual ClosedForm::IProduct* createProduct(ClosedForm* model) const;
        
    /** when to stop tweaking */
    virtual DateTime endDate(const Sensitivity* sensControl) const;
    
private:
    AUD3YBondFuture();
    AUD3YBondFuture(const AUD3YBondFuture& rhs);
    AUD3YBondFuture& operator=(const AUD3YBondFuture& rhs);

    void price(Control* control, CResults* results) const;
    void requests(Control*        control, 
                  CResults*       results,
                  double          value) const;

    DateTime            startDate;
    BondParamsSP        PseudoBond ;

};

DRLIB_END_NAMESPACE
#endif
