//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : AssetValueCreditSupport.hpp
//
//   Description : Credit support object for AssetValue
//
//   Author      : Ning Shen
//
//   Date        : 19 Sept 2002
//
//
//----------------------------------------------------------------------------

#ifndef ASSETVALUE_CREDIT_HPP
#define ASSETVALUE_CREDIT_HPP

#include "edginc/CreditSupport.hpp"

DRLIB_BEGIN_NAMESPACE

class AssetValue;

/** AssetValue credit support object  */
class PRODUCTS_DLL AssetValueCreditSupport : virtual public CreditSupport
{
public:

    /** preprocess instrument for a given set of path dates */
    virtual void preProcess(const DateTimeArray& dates, 
							const DoubleArray& atmFwd, 
							const DoubleArray& atmVar);

    /** return asset */
    virtual CreditUndSP getUnderlier() const;

    /** calculate values for a given path */
    virtual void calcPathValues(DoubleArray& results, const DateTimeArray& dates, 
                    const double* spots, double spotRef);

    /** return model for this instrument */
    virtual IModelSP getModel();

    /** return instrument ccy ISO code */
    virtual string getInstCcyCode() const;

	/** return instrument's last exposure date */
    virtual DateTime getInstLastExposureDate() const;

    AssetValueCreditSupport(CInstrument*, CMarketDataSP market);

private:
    IModelSP model;
    
    AssetValue* instrAVOriginal;

    smartPtr<AssetValue> instrAV;

    DoubleArray FwdCache;
    DoubleArray FwdDeltaCache;
};

DRLIB_END_NAMESPACE

#endif
