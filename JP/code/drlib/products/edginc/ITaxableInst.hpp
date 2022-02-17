//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ITaxableInst.hpp
//
//   Description : Allow an instrument to specify cash flows and hence be
//                 fully supported for tax treatment
//
//   Date        : Oct 2003
//
//
//----------------------------------------------------------------------------

#ifndef EDR_ITAXABLE_INST_HPP
#define EDR_ITAXABLE_INST_HPP

#include "edginc/Object.hpp"
#include "edginc/CashFlow.hpp"

DRLIB_BEGIN_NAMESPACE

class PRODUCTS_DLL ITaxableInst {
public:
    /** Basic info needed for tax treatment 
        - ability to price (implicit in being an instrument), and
        - final pay date (indicating end of tax liability) */
    class PRODUCTS_DLL Basic {
    public:
        virtual const DateTime getFinalPaymentDate() const = 0;
        
        virtual ~Basic() {};
    };

    /** An instrument may (i.e. this is optional)
        provide more info : what cash flows it generates */
    class PRODUCTS_DLL WithCoupons {
    public:
        virtual CashFlowArrayConstSP getCoupons() const = 0;
        
        virtual ~WithCoupons() {};
    };
};

DRLIB_END_NAMESPACE

#endif

