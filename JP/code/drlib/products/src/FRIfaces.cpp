//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FRIfaces.cpp
//
//   Description : For registration
//
//   Date        : 3 April 2002
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/FRIfaces.hpp"
DRLIB_BEGIN_NAMESPACE

FRIfaces::IRValue::~IRValue(){}
FRIfaces::IHoldsState::~IHoldsState(){}
FRIfaces::ILValue::~ILValue(){}
FRIfaces::ILConstValue::~ILConstValue(){}
FRIfaces::IRValueDouble::~IRValueDouble(){}
FRIfaces::IRValueInt::~IRValueInt(){}
FRIfaces::IRValueBool::~IRValueBool(){}
FRIfaces::IRValueDate::~IRValueDate(){}
FRIfaces::IRValueDoubleArray::~IRValueDoubleArray(){}
FRIfaces::IRValueIntArray::~IRValueIntArray(){}
FRIfaces::IRValueBoolArray::~IRValueBoolArray(){}
FRIfaces::IRValueSchedule::~IRValueSchedule(){}
FRIfaces::IRValueTabulatedFunc::~IRValueTabulatedFunc(){}
FRIfaces::ILValueDouble::~ILValueDouble(){}
FRIfaces::ILValueDoubleArray::~ILValueDoubleArray(){}
FRIfaces::ILValueIntArray::~ILValueIntArray(){}
FRIfaces::ILValueBoolArray::~ILValueBoolArray(){}
FRIfaces::ILValueInt::~ILValueInt(){}
FRIfaces::ILValueBool::~ILValueBool(){}
FRIfaces::ILValueDate::~ILValueDate(){}
FRIfaces::IRValueExpression::~IRValueExpression(){}
FRIfaces::ILValueExpression::~ILValueExpression(){}
FRIfaces::IProductView::~IProductView(){}
FRIfaces::IAssignment::~IAssignment(){}
FRIfaces::ITimePointRules::~ITimePointRules(){}
FRIfaces::IAlgorithm::~IAlgorithm(){}
FRIfaces::ITimePointNoOp::~ITimePointNoOp(){}

/** Invoked when Class is 'loaded' - needed to use Addin feature in load */
#define MY_LOAD(funcName,className) static void funcName(CClassSP& clazz){\
    REGISTER_INTERFACE(className, clazz);\
    EXTENDS(IObject);\
    clazz->setPublic(); /* make visible to EAS/spreadsheet */\
}

MY_LOAD(load1, FRIfaces::IRValueExpression);
CClassConstSP const FRIfaces::IRValueExpression::TYPE = 
CClass::registerInterfaceLoadMethod(
    "FRIfaces::IRValueExpression", typeid(FRIfaces::IRValueExpression), 
    load1);

// work around for msvc 7 bug
typedef FRIfaces::IRValueExpressionArray FRIfacesIRValueExpressionArray;
DEFINE_TEMPLATE_TYPE_WITH_NAME("FRIfaces::IRValueExpressionArray", FRIfacesIRValueExpressionArray);


/** Invoked when Class is 'loaded' - needed to use Addin feature in load */
static void loadILValueExpression(CClassSP& clazz) {
    REGISTER_INTERFACE(FRIfaces::ILValueExpression, clazz);
    EXTENDS(FRIfaces::IRValueExpression);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

CClassConstSP const FRIfaces::ILValueExpression::TYPE = 
CClass::registerInterfaceLoadMethod(
    "FRIfaces::ILValueExpression", typeid(FRIfaces::ILValueExpression), 
    loadILValueExpression);

// work around for msvc 7 bug
typedef FRIfaces::ILValueExpressionArray FRIfacesILValueExpressionArray;
DEFINE_TEMPLATE_TYPE_WITH_NAME("FRIfaces::ILValueExpressionArray", FRIfacesILValueExpressionArray);

MY_LOAD(load2, FRIfaces::IAlgorithm);
CClassConstSP const FRIfaces::IAlgorithm::TYPE = 
CClass::registerInterfaceLoadMethod(
    "FRIfaces::IAlgorithm", typeid(FRIfaces::IAlgorithm), load2);

MY_LOAD(load3, FRIfaces::ITimePointRules);
CClassConstSP const FRIfaces::ITimePointRules::TYPE = 
CClass::registerInterfaceLoadMethod(
    "FRIfaces::ITimePointRules", typeid(FRIfaces::ITimePointRules), load3);

// work around for msvc 7 bug
typedef FRIfaces::ITimePointRulesArray FRIfacesITimePointRulesArray;
DEFINE_TEMPLATE_TYPE_WITH_NAME("FRIfaces::ITimePointRulesArray", FRIfacesITimePointRulesArray);

MY_LOAD(load4, FRIfaces::ITimePointNoOp);
CClassConstSP const FRIfaces::ITimePointNoOp::TYPE = 
CClass::registerInterfaceLoadMethod(
    "FRIfaces::ITimePointNoOp", typeid(FRIfaces::ITimePointNoOp), load4);

// work around for msvc 7 bug
typedef FRIfaces::ITimePointNoOpArray FRIfacesITimePointNoOpArray;
DEFINE_TEMPLATE_TYPE_WITH_NAME("FRIfaces::ITimePointNoOpArray", FRIfacesITimePointNoOpArray);

DRLIB_END_NAMESPACE

