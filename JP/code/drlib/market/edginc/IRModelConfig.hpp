//----------------------------------------------------------------------------
//
//   Group       : QR Rates
//
//   Filename    : IRModelConfig.hpp
//
//   Description : ir model configuration base class
//
//   Author      : Anwar E Sidat
//
//   Date        : 18-Aug-2006
//
//----------------------------------------------------------------------------
#ifndef QLIB_IR_MODEL_CONFIG_HPP
#define QLIB_IR_MODEL_CONFIG_HPP

#include "edginc/MarketObject.hpp"
#include "edginc/MarketTable.hpp"

DRLIB_BEGIN_NAMESPACE

/** inherit from market object purely for convenience for now. Easily retrieve data from market env.
    However, there should be proper model configuration storage place, not in market */
class MARKET_DLL IRModelConfig : public MarketObject
{
public:
    static CClassConstSP const TYPE;

    virtual ~IRModelConfig();

    /** overrides clone */
    IObject* clone() const;

protected:
    IRModelConfig(const CClassConstSP& clazz);
    IRModelConfig(const IRModelConfig& irv);
    IRModelConfig& operator=(const IRModelConfig& irv);

private:
    static void load(CClassSP& clazz);
};


// Support for smart pointer, wrapper and market table
typedef smartConstPtr<IRModelConfig> IRModelConfigConstSP;
typedef smartPtr<IRModelConfig> IRModelConfigSP;
typedef array<IRModelConfigSP, IRModelConfig> IRModelConfigArray;
typedef MarketWrapper<IRModelConfig> IRModelConfigWrapper;
typedef smartPtr<IRModelConfigWrapper> IRModelConfigWrapperSP;
typedef array<IRModelConfigWrapperSP,IRModelConfigWrapper> IRModelConfigWrapperArray;
typedef MarketTable<IRModelConfig> IRModelConfigTable;
typedef MarketWrapper<IRModelConfigTable> IRModelConfigTableWrapper;

#ifndef QLIB_IR_MODEL_CONFIG_CPP
    EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<IRModelConfig>);
    EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<IRModelConfig>);
    EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<IRModelConfig>);
    EXTERN_TEMPLATE(class MARKET_DLL_SP MarketTable<IRModelConfig>);
    EXTERN_TEMPLATE(class MARKET_DLL_SP MarketWrapper<IRModelConfigTable>);
#else
    INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<IRModelConfig>);
    INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<IRModelConfig>);
    INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<IRModelConfig>);
    INSTANTIATE_TEMPLATE(class MARKET_DLL MarketTable<IRModelConfig>);
    INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<IRModelConfigTable>);
#endif

DRLIB_END_NAMESPACE

#endif
