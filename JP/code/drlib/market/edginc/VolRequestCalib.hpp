//----------------------------------------------------------------------------
//
//   Group       : QR Rates
//
//   Filename    : VolRequestCalib.hpp
//
//   Description : Vol request for calibration index.
//
//   Author      : Anwar E Sidat
//
//   Date        : 06-Nov-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_VolRequestCalib_HPP
#define QLIB_VolRequestCalib_HPP

#include "edginc/VolRequest.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/IRVol.hpp"

DRLIB_BEGIN_NAMESPACE

class IRVol;

/** VolRequestCalib class.
 */
class MARKET_DLL VolRequestCalib : public CVolRequest
{
public:
    static CClassConstSP const TYPE;

    virtual ~VolRequestCalib();
    VolRequestCalib(IRVol::CalibType calibType, ExpirySP calibTenor, CDoubleSP calibVolOverride);
    
    const IRVol::CalibType calibType() const { return m_calibType; }
    ExpirySP calibTenor() const { return m_calibTenor; }
    CDoubleSP calibVolOverride() const { return m_calibVolOverride; }

protected:
    VolRequestCalib(const CClassConstSP& clazz);
    VolRequestCalib(const VolRequestCalib& irv);
    VolRequestCalib& operator=(const VolRequestCalib& irv);

    //Fields
    string           name;         /* optional - but recommended.  */
    StringArray      style;        // for each 'style' there is a DoubleArray
    StringArray      paramLabel;   // name for param column there is an optional label
    DoubleArrayArray params;

private:
    static void load(CClassSP& clazz);
    
    IRVol::CalibType m_calibType;
    ExpirySP         m_calibTenor;        // tenor or fixed maturity date
    CDoubleSP        m_calibVolOverride;
};


// Support for smart pointers and array types
typedef smartConstPtr<VolRequestCalib> VolRequestCalibConstSP;
typedef smartPtr<VolRequestCalib> VolRequestCalibSP;
typedef array<VolRequestCalibSP, VolRequestCalib> VolRequestCalibArray;

#ifndef QLIB_VolRequestCalib_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<VolRequestCalib>);
EXTERN_TEMPLATE(class MARKET_DLL array<VolRequestCalibSP _COMMA_ VolRequestCalib>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<VolRequestCalib>);
INSTANTIATE_TEMPLATE(class MARKET_DLL array<VolRequestCalibSP _COMMA_ VolRequestCalib>);
#endif

typedef smartConstPtr<VolRequestCalibArray> VolRequestCalibArrayConstSP;
typedef smartPtr<VolRequestCalibArray> VolRequestCalibArraySP;

DRLIB_END_NAMESPACE

#endif
