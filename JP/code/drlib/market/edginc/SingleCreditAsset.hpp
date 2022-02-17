//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Description : Empty class that just derives from CreditAsset. See 
//                 CreditAsset.hpp for details
//
//   Date        : March 2005
//
//----------------------------------------------------------------------------

#ifndef EDR_SINGLE_CREDIT_ASSET_HPP
#define EDR_SINGLE_CREDIT_ASSET_HPP

#include "edginc/CreditAsset.hpp"

DRLIB_BEGIN_NAMESPACE

//------------------
// SingleCreditAsset
//------------------
class MARKET_DLL SingleCreditAsset: public CreditAsset {
public:
    static CClassConstSP const TYPE;

    virtual ~SingleCreditAsset();
    
private:
    friend class StopCompilerWarning;
    SingleCreditAsset();
    SingleCreditAsset(const SingleCreditAsset& rhs); // don't use
    SingleCreditAsset& operator=(const SingleCreditAsset& rhs); // don't use
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();
};

typedef smartPtr<SingleCreditAsset> SingleCreditAssetSP;
typedef smartConstPtr<SingleCreditAsset> SingleCreditAssetConstSP;

DRLIB_END_NAMESPACE

#endif
