//----------------------------------------------------------------------------
//
//   Group       : QR Rates
//
//   Filename    : VolCalibInterp.hpp
//
//   Description : Simulates a VolCalibInterp for testing purposes!
//
//   Author      : Anwar E Sidat
//
//   Date        : 07-Oct-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_VolCalibInterp_HPP
#define QLIB_VolCalibInterp_HPP

#include "edginc/Addin.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/BenchmarkDate.hpp"
#include "edginc/Surface.hpp"
#include "edginc/IRVolBase.hpp"
#include "edginc/IRVol.hpp"
#include "edginc/VolProcessedCalib.hpp"


DRLIB_BEGIN_NAMESPACE

/** VolCalibInterp Class.
 */
class MARKET_DLL VolCalibInterp : public CObject
{
public:
    static CClassConstSP const TYPE;

    VolCalibInterp();
    virtual ~VolCalibInterp();

    /** Returns name of model. */
    virtual string getName() const;

    /** overrides default */
    virtual void validatePop2Object();

    /** goVolCalibInterp */
    CDoubleMatrixSP goVolCalibInterp();

    /** static interpVolsForCalibration */
    static VolProcessedCalibSP interpVolsForCalibration(
        const ExpiryArray&    xArray,       // Expiries
        const ExpiryArray&    yArray,       // Maturity Tenors
        const CDoubleMatrix&  zMatrix,      // 2D Matrix
        DateTime              baseDate,     // BaseDate
        IRVol::CalibType      calibType,
        ExpirySP              calibTenor,
        CDoubleSP             calibVolOverride = CDoubleSP());
  
protected:

    VolCalibInterp(const CClassConstSP& clazz);
    VolCalibInterp(const VolCalibInterp& irv);
    VolCalibInterp& operator=(const VolCalibInterp& irv);

    //Fields
    IRVolBaseWrapper vol;
    IRVol::CalibType calibType;
    ExpirySP         calibTenor;
    CDoubleSP        calibVolOverride;

private:

    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new VolCalibInterp(); }
};

DRLIB_END_NAMESPACE

#endif
