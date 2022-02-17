//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ClosedFormLN.hpp
//
//   Description : Closed Form Log Normal Algorithm (asks instrument to do it)
//
//   Author      : Mark A Robson
//
//   Date        : 15 Jan 2001
//
//
//----------------------------------------------------------------------------

#ifndef CLOSED_FORM_LN_HPP
#define CLOSED_FORM_LN_HPP
#include "edginc/ModelLN.hpp"
#include "edginc/IRVegaPointwise.hpp"

DRLIB_BEGIN_NAMESPACE

/** Implementation of Model algorithm where the algorithm is specific to
    the instrument and yet is based on a log normal methodology */
class PRODUCTS_DLL CClosedFormLN: public CModelLN,
    /* should stop deriving from ISensitivePoints once all the ir products
       switch to using ClosedFormIRLN */
                     public virtual IRVegaPointwise::ISensitivePoints{
public:
    static CClassConstSP const TYPE;
    friend class CClosedFormLNHelper;

    /** the class that the product must be able to create */
    class PRODUCTS_DLL IProduct{
    public:
        IProduct();
        virtual void price(CClosedFormLN*    model,
                           Control*          control, 
                           CResults*         results) const = 0;
        virtual ~IProduct();
    };

    /** interface that the instrument must implement */
    class PRODUCTS_DLL IIntoProduct: virtual public CModel::IModelIntoProduct {
    public:
        friend class CClosedFormLNHelper;
        static CClassConstSP const TYPE;
        IIntoProduct();
        virtual ~IIntoProduct();
        virtual IProduct* createProduct(CClosedFormLN* model) const = 0;
    };

    /** Constructor takes type of vol to use */
    CClosedFormLN(const string& volType,
                  bool          allowNegativeFwdVar 
                        = CModelLN::allowNegativeFwdVar_default);

    /** defaults type of vol to use to VolatilityBS */
    CClosedFormLN();

    /** calculate single price and store result in CResult */
    virtual void Price(CInstrument*  instrument, 
                       CControl*     control, 
                       CResults*     results);
    
    /** Essentially relies on instrument implementing ISensitiveIRVolPoints.
        If not returns null. */
    virtual IRGridPointAbsArraySP getSensitiveIRVolPoints(
        OutputNameConstSP  outputName,
        const CInstrument* inst) const;
private:
    CClosedFormLN(const CClosedFormLN &rhs);
    CClosedFormLN& operator=(const CClosedFormLN& rhs);
};

DRLIB_END_NAMESPACE
#endif
