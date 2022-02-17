//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ModelLN.hpp
//
//   Description : Abstract Log Normal Algorithm (selects type of vol)
//
//   Author      : Mark A Robson
//
//   Date        : 22 Mar 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_MODEL_LN_HPP
#define EDG_MODEL_LN_HPP
#include "edginc/Model.hpp"
#include "edginc/MarketDataFetcherLN.hpp"

DRLIB_BEGIN_NAMESPACE

/** Abstract Model algorithm which chooses a Log Normal volatility */
class PRODUCTS_DLL CModelLN: public CModel{
public:
    static CClassConstSP const TYPE;

    static const bool allowNegativeFwdVar_default;

    /** Override default createMDF in order to set the right MDF */
    virtual MarketDataFetcherSP createMDF() const;

    bool negativeFwdVarAllowed() const;

    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this model
     *
     * Returns riskMappingIrrelevant, because we price using a plain old vol
     * surface.  See IModel::wantsRiskMapping().
     */

    virtual IModel::WantsRiskMapping wantsRiskMapping() const;

protected:
    /** constructor takes type of vol to use */
    CModelLN(const CClassConstSP& clazz, 
             const string&        volType,
             bool                 allowNegativeFwdVar = allowNegativeFwdVar_default);

    /** defaults type of vol to use to VolatilityBS */
    CModelLN(const CClassConstSP& clazz);
    
    /** return the volType in use */
    string getVolType() const;

private:
    CModelLN(const CModelLN &rhs);
    CModelLN& operator=(const CModelLN& rhs);
    static void load(CClassSP& clazz);

    string volType;
    bool   allowNegativeFwdVar;

    mutable MarketDataFetcherLNSP mdf;  // not in registration $unregistered
};

DRLIB_END_NAMESPACE
#endif
