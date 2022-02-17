//----------------------------------------------------------------------------
//
//   Group       : QR - Credit Hybrids
//
//   Filename    : WeightedInstrumentShift.hpp
//
//   Description : Sensitivity to a shift in instrument weight 
//                 (implemented by WeightedInstrumentCollection)
//
//   Author      : Linus Thand
//
//   Date        : 13 July 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_INSTRUMENTWEIGHTSHIFT_H
#define QLIB_INSTRUMENTWEIGHTSHIFT_H

/* A sensitivity to a shift in instrument weight (as in WeightedInstrumentCollection),
 * used by managed trades through ImpliedScalarShiftMulti.
 * Inheriting from AllNamesRiskPropertySensitivity, so it's a declarative greek.
 * Implementing IScalarPerNameShift, so it can be used for solving.
 */

#include "edginc/WeightedInstrumentTweak.hpp"
#include "edginc/AllNamesRiskPropertySensitivity.hpp"
#include "edginc/ScalarRiskPropertySensitivity.hpp"

DRLIB_BEGIN_NAMESPACE

class RISKMGR_DLL WeightedInstrumentShift : public AllNamesRiskPropertySensitivity, 
                                            virtual public IScalarPerNameShift {
 
 public:
    static CClassConstSP const TYPE;
    static const double DEFAULT_SHIFT;
  
    //// To implement IScalarPerNameShift
    IHypothesis::AlternateWorldSP appliedTo(OutputNameConstSP name,
                                            double shiftSize,
                                            IObjectSP world);
   
    OutputNameArrayConstSP allNames(const IObject* object) const;
    
    WeightedInstrumentShift(double shiftSize = DEFAULT_SHIFT);
 private:
    ~WeightedInstrumentShift();
    WeightedInstrumentShift(const WeightedInstrumentShift& rhs); 
    WeightedInstrumentShift& operator=(const WeightedInstrumentShift& rhs);
    static const double SENSITIVITY_UNIT;
    static const string NAME;
    static void load(CClassSP& clazz);
    Deriv deriv() const;
    double getShiftSize() const;
    
    // Field
    IntArray instrumentList; 
};

FORWARD_DECLARE(WeightedInstrumentShift)

DRLIB_END_NAMESPACE

#endif
