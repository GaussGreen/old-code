//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ClosedForm.hpp
//
//   Description : Closed Form Algorithm (asks instrument to do it)
//
//   Author      : Andrew J Swain
//
//   Date        : 16 February 2001
//
//
//----------------------------------------------------------------------------

#ifndef CLOSEDFORM_HPP
#define CLOSEDFORM_HPP
#include "edginc/Model.hpp"

DRLIB_BEGIN_NAMESPACE

/** Implementation of Model algorithm where the algorithm is specific to
    the instrument and doesn't depend on any vo assumptions.
    What we're looking at here are things like cashflows, equities and 
    forwards
*/
class RISKMGR_DLL ClosedForm: public CModel {
public:
    static CClassConstSP const TYPE;
    friend class ClosedFormHelper;

    /** the class that the product must be able to create */
    class RISKMGR_DLL IProduct {
    public:
        virtual void price(ClosedForm* model,
                           Control*    control, 
                           CResults*   results) const = 0;

        virtual ~IProduct();
    };

    /** interface that the instrument must implement */
    class RISKMGR_DLL IIntoProduct: virtual public CModel::IModelIntoProduct {
    public:
        friend class ClosedFormHelper;
        static CClassConstSP const TYPE;
        virtual IProduct* createProduct(ClosedForm* model) const = 0;
    };

    /** Simple constructor */
    ClosedForm();

    /** calculate single price and store result in CResult */
    virtual void Price(CInstrument*  instrument, 
                       CControl*     control, 
                       CResults*     results);
    
    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this model
     *
     * Returns riskMappingIrrelevant, because ClosedForm models by definition
     * don't depend on any vol, hence a fortiori don't involve parametric vol.
     * See IModel::wantsRiskMapping().
     */

    virtual IModel::WantsRiskMapping wantsRiskMapping() const;

private:
    ClosedForm(const ClosedForm &rhs);
    ClosedForm& operator=(const ClosedForm& rhs);
};

DRLIB_END_NAMESPACE
#endif
