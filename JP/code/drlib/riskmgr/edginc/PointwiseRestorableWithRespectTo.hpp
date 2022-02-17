/**
 * \file PointwiseRestorableWithRespectTo.hpp
 *
 *
 * $History$
 *
 */

#ifndef EDG_RISKMGR_POINTWISERESTORABLEWITHRESPECTTO_H
#define EDG_RISKMGR_POINTWISERESTORABLEWITHRESPECTTO_H

#include "edginc/PointwiseTweakableWithRespectTo.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Interface for (market data) objects which can be subjected to an in-place
 * tweak/restore of a certain type
 *
 * See GenericVectorShfit and GenericScalarShift for info about the purpose of
 * this interface.  It is an explicit representation of the implicit signature
 * common to classes like Delta::IRestore, CorrelationSqueeze::IRestore etc.
 *
 * TWEAK might for example be ParSpreadRhoPointwiseTweak.
 */

template <class TWEAK>
struct PointwiseRestorableWithRespectTo: public PointwiseTweakableWithRespectTo<TWEAK> {
  static CClassConstSP const TYPE;

  virtual void sensRestore(TWEAK *shift) = 0;

  static void load(CClassSP &clazz) {
    REGISTER_INTERFACE(PointwiseRestorableWithRespectTo, clazz);
    EXTENDS(PointwiseTweakableWithRespectTo<TWEAK>);
  }
};

DRLIB_END_NAMESPACE

#endif
