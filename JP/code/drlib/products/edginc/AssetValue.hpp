//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : AssetValue.hpp
//
//   Description : Asset instrument
//
//   Author      : Andrew J Swain
//
//   Date        : 7 September 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDR_ASSETVALUE_HPP
#define EDR_ASSETVALUE_HPP
#include "edginc/Class.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/ClosedFormLN.hpp"
#include "edginc/Theta.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/Asset.hpp"
#include "edginc/AssetValueCreditSupport.hpp"

DRLIB_BEGIN_NAMESPACE

/** Asset instrument */
class PRODUCTS_DLL AssetValue: public CInstrument, 
                  public CClosedFormLN::IIntoProduct,
                  public LastSensDate, 
                  public Theta::IShift,
                  public CreditSupport::Interface
{
public:
    static CClassConstSP const TYPE;

    /** instrument validation */
    virtual void Validate();

    /** input data validation */
    virtual void validatePop2Object();

    /** retrieve market data needed by Vanilla - just valueDate, asset and
        discount yield curve */
    void GetMarket(const IModel*          model, 
                   const CMarketDataSP    market);

    /** Implementation of ClosedForm::IntoProduct interface */
    virtual CClosedFormLN::IProduct* createProduct(CClosedFormLN* model) const;

   /** Returns rolls value date and sets initial spot for Theta,
       return true if sub objects need to be tweaked */
    virtual bool sensShift(Theta* shift);

    /** when to stop tweaking */
    virtual DateTime endDate(const Sensitivity* sensControl) const;

    virtual DateTime getValueDate() const;

    virtual CreditSupportSP createCreditSupport(CMarketDataSP market){
                return CreditSupportSP(new AssetValueCreditSupport(this, market));}

    /** Returns the name of the instrument's discount currency */
    virtual string discountYieldCurveName() const;

private:
    friend class AssetValueHelper;
    friend class AssetValueClosedForm;
    friend class AssetValueCreditSupport;

    AssetValue();
    AssetValue(const AssetValue& rhs);
    AssetValue& operator=(const AssetValue& rhs);

protected:
    AssetValue(CClassConstSP clazz);

    CAssetWrapper asset;
    DateTime      valueDate;
};

typedef smartPtr<AssetValue> AssetValueSP;

DRLIB_END_NAMESPACE
#endif
