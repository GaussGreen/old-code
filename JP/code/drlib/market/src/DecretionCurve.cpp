//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : DecretionCurve.cpp
//
//   Description : A decretion curve base for ABS
//
//   Author      : Keith Jia
//
//   Date        : 03 January 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Class.hpp"
#include "edginc/DecretionCurve.hpp"
#include "edginc/Actual365F.hpp"
#include "edginc/Maths.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

//------------------------------------------
// DecretionCurve methods
//------------------------------------------

DecretionCurve::DecretionCurve(CClassConstSP clazz) : MarketObject(clazz)
{}

DecretionCurve::~DecretionCurve() 
{}

//---------------------------------------------------------
// load
//---------------------------------------------------------

void DecretionCurve::load(CClassSP& clazz){
    clazz->setPublic();
    REGISTER(DecretionCurve, clazz);
    SUPERCLASS(MarketObject);
    IMPLEMENTS(IDecretionCurve);
    FIELD(name, "name");
}

//---------------------------------------------------------
// MarketObject methods
//---------------------------------------------------------
string DecretionCurve::getName() const
{
    return name;
}

//---------------------------------------------------------
// static variables
//---------------------------------------------------------
CClassConstSP const DecretionCurve::TYPE = 
    CClass::registerClassLoadMethod("DecretionCurve", typeid(DecretionCurve), load);

DEFINE_TEMPLATE_TYPE(DecretionCurveWrapper);


DRLIB_END_NAMESPACE
