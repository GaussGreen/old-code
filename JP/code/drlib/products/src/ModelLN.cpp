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

#include "edginc/config.hpp"
#include "edginc/ModelLN.hpp"
#include "edginc/VolatilityBS.hpp"

DRLIB_BEGIN_NAMESPACE

const bool CModelLN::allowNegativeFwdVar_default = false;

/** constructor takes type of vol to use */
CModelLN::CModelLN(const CClassConstSP& clazz, 
                   const string&        volType,
                   bool                 allowNegativeFwdVar):
CModel(clazz), 
volType(volType), 
allowNegativeFwdVar(allowNegativeFwdVar), 
mdf(0){}

/** defaults type of vol to use to VolatilityBS */
CModelLN::CModelLN(const CClassConstSP& clazz):
CModel(clazz), 
volType(IVolatilityBS::TYPE->getName()), 
allowNegativeFwdVar(allowNegativeFwdVar_default), 
mdf(0){}

/** Override default createMDF in order to set the right MDF */
MarketDataFetcherSP CModelLN::createMDF() const{
    return MarketDataFetcherSP(new MarketDataFetcherLN(volType));
}

IModel::WantsRiskMapping CModelLN::wantsRiskMapping() const {
    return riskMappingIrrelevant;
}

/** return the volType in use */
string CModelLN::getVolType() const {
    return volType;
}


bool CModelLN::negativeFwdVarAllowed() const{
    return allowNegativeFwdVar;
}

/** Invoked when Class is 'loaded' */
void CModelLN::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(CModelLN, clazz);
    SUPERCLASS(Model);
    FIELD(volType, "Type of vol to use");
    FIELD(allowNegativeFwdVar, "If true, allow negative fwd var; don't, otherwise");
    FIELD_MAKE_OPTIONAL(allowNegativeFwdVar);
}


CClassConstSP const CModelLN::TYPE = CClass::registerClassLoadMethod(
    "ModelLN", typeid(CModelLN), load);

DRLIB_END_NAMESPACE
