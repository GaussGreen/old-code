#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#define QLIB_CREDITINDEXSPREADPOINTWISE_CPP
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/CreditIndexSpreadPointwise.hpp"
#include "edginc/CreditIndex.hpp"
#include "edginc/TRACE.hpp"
#include "edginc/SimpleTweakNameResolver.hpp"

DRLIB_BEGIN_NAMESPACE

CreditIndexSpreadPointwise::CreditIndexSpreadPointwise(): CObject(TYPE) {}
CreditIndexSpreadPointwise::~CreditIndexSpreadPointwise() {}

static void CreditIndexSpreadPointwise_load(CClassSP& clazz) {
    REGISTER(CreditIndexSpreadPointwise, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(DefaultConstructor<CreditIndexSpreadPointwise>::iObject);
}

CClassConstSP const CreditIndexSpreadPointwise::TYPE = CClass::registerClassLoadMethod("CreditIndexSpreadPointwise", typeid(CreditIndexSpreadPointwise), CreditIndexSpreadPointwise_load);

RiskProperty_TYPES(CreditIndexSpreadPointwise)

//implementation of riskMappingMatrix for the CreditIndexSpreadPointwise tweak
//overrides default method in IRiskProperty
//this is exactly the same as the parallel method
//the risk mapping matrix caters for both cases
template <>
RiskMappingMatrixConstSP RiskProperty<CreditIndexSpreadPointwise>::riskMappingMatrix(
    IObjectConstSP world,
    OutputNameConstSP name) const
{
    static const string method = "RiskProperty<CreditIndexSpreadPointwise>::riskMappingMatrix";

    TRACE_METHOD;
    TRACE("Attempting to find CreditIndex with name " + name->toString());
    // find the credit index in the world
    SimpleTweakNameResolver stnr(name);
    SensMgrConst smc(world);
    IObjectConstSP theIndexObject;
    try {
        theIndexObject = smc.theFirst(CreditIndex::TYPE, &stnr);
    } catch (exception&) {
        TRACE("...failed");
        throw ModelException(method, "Failed to find index " + name->toString());
    }

    //validate
    if (!(theIndexObject.get()))
    {
        TRACE("...failed");
        throw ModelException(method, "Failed to find index " + name->toString());
    }

    TRACE("...succeeded");
    // and now ask it for a risk mapping matrix
    // since it is the credit index that has knowledge
    // of the names & expiries to map to

    CreditIndexConstSP theIndex = CreditIndexConstSP::dynamicCast(theIndexObject);
    return theIndex->getRiskMappingMatrix();
}

DRLIB_END_NAMESPACE
