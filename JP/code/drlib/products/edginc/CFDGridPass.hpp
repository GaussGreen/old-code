//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CFDGridPass.hpp
//
//   Description : FD associate a product to a model
//
//
//----------------------------------------------------------------------------

#ifndef CFDGridPass_HPP
#define CFDGridPass_HPP
#include "edginc/Model.hpp"

DRLIB_BEGIN_NAMESPACE


class PRODUCTS_DLL CFDGridPass: public CModel{
public:
    static CClassConstSP const TYPE;
    friend class CFDGridPassHelper;

    /** associate a product to a model */
  //  virtual void GetProduct(CInstrument *inst);

    
    class PRODUCTS_DLL IProduct{
    public:
        virtual void price(CFDGridPass*    model,
                           Control*          control, 
                           CResults*         results) const = 0;
    };

    /** interface that the instrument must implement */
    class PRODUCTS_DLL IIntoProduct{
    public:
        friend class CFDGridPassHelper;
        static CClassConstSP const TYPE;
        virtual IProduct* createProduct(CFDGridPass* model) const = 0;
    };

    CFDGridPass(const string& volType);
    /** Simple constructor */
    CFDGridPass();

  
    /** calculate single price and store result in CResult */
    virtual void Price(CInstrument*  instrument, 
                       CControl*     control, 
                       CResults*     results);

    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this model
     *
     * Returns riskMappingIrrelevant.  See IModel::wantsRiskMapping().
     */

    virtual IModel::WantsRiskMapping wantsRiskMapping() const;

	string volType;
	int    iMax;
	int    nMax; 
	int    divType;
    
private:

    CFDGridPass(const CFDGridPass &rhs);
    CFDGridPass& operator=(const CFDGridPass& rhs);

};


DRLIB_END_NAMESPACE
#endif
