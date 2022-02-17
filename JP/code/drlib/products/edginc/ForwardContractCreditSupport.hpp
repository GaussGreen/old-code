//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ForwardContractCreditSupport.hpp
//
//   Description : Credit support object for ForwardContract
//
//   Author      : Ning Shen
//
//   Date        : 10 June 2002
//
//
//----------------------------------------------------------------------------

#ifndef FORWARDCONTRACT_CREDIT_HPP
#define FORWARDCONTRACT_CREDIT_HPP

#include "edginc/CreditSupport.hpp"

DRLIB_BEGIN_NAMESPACE

class CForwardContract;

/** ForwardContract credit support object  */
class PRODUCTS_DLL ForwardContractCreditSupport : virtual public CreditSupport
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

    ForwardContractCreditSupport(CInstrument*, CMarketDataSP market);

private:
    IModelSP model;
    
    CForwardContract* fwdOriginal;

    smartPtr<CForwardContract> fwdContract;

    DoubleArray FwdCache;
    DoubleArray FwdDeltaCache;
};

DRLIB_END_NAMESPACE

#endif
