/**
 * \file RestorableWithRespectTo.hpp
 *
 *
 * $History$
 *
 */

#ifndef EDG_RISKMGR_RESTORABLEWITHRESPECTTO_H
#define EDG_RISKMGR_RESTORABLEWITHRESPECTTO_H

#include "edginc/Class.hpp"
#include "edginc/TweakableWithRespectTo.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Interface for (market data) objects which can be subjected to an in-place
 * tweak/restore of a certain type
 *
 * See GenericScalarShift for info about the purpose of this interface.  It
 * is an explicit representation of the implicit signature common to classes
 * like Delta::IRestore, CorrelationSqueeze::IRestore etc.
 *
 * TWEAK might for example be ParSpreadRhoParallelTweak.
 */

template <class TWEAK>
struct RestorableWithRespectTo: public TweakableWithRespectTo<TWEAK> {
  static CClassConstSP const TYPE;

  virtual void sensRestore(TWEAK *shift) = 0;

  static void load(CClassSP &clazz) {
    REGISTER_INTERFACE(RestorableWithRespectTo, clazz);
    EXTENDS(TweakableWithRespectTo<TWEAK>);
  }
};

DRLIB_END_NAMESPACE

#endif
