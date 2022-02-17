//----------------------------------------------------------------------------
//
//   Group       : QR Rates
//
//   Filename    : VolRequestCalib.cpp
//
//   Description : Vol request for calibration index.
//
//   Author      : Anwar E Sidat
//
//   Date        : 06-Nov-2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_VolRequestCalib_CPP
#include "edginc/VolRequestCalib.hpp"

DRLIB_BEGIN_NAMESPACE


VolRequestCalib::~VolRequestCalib(){}

VolRequestCalib::VolRequestCalib(const CClassConstSP& clazz)
    : CVolRequest(clazz){}

VolRequestCalib::VolRequestCalib(IRVol::CalibType calibType, ExpirySP calibTenor, CDoubleSP calibVolOverride)
    : CVolRequest(TYPE), m_calibType(calibType), m_calibTenor(calibTenor), m_calibVolOverride(calibVolOverride)
{
}

void VolRequestCalib::load(CClassSP& clazz)
{
    clazz->setPublic();
    clazz->setDescription("Vol request for calibration index.");
    REGISTER(VolRequestCalib, clazz);
    SUPERCLASS(CVolRequest);
}

CClassConstSP const VolRequestCalib::TYPE = CClass::registerClassLoadMethod(
    "VolRequestCalib", typeid(VolRequestCalib), VolRequestCalib::load);

// initialise type for array of CVolRequestCalib
DEFINE_TEMPLATE_TYPE_WITH_NAME("VolRequestCalibArray", VolRequestCalibArray);

DRLIB_END_NAMESPACE
