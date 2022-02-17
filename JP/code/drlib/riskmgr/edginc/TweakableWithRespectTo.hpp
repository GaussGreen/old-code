/**
 * \file TweakableWithRespectTo.hpp
 *
 *
 * $History$
 *
 */

#ifndef EDG_RISKMGR_TWEAKABLEWITHRESPECTTO_H
#define EDG_RISKMGR_TWEAKABLEWITHRESPECTTO_H

#include "edginc/Class.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Interface for (market data) objects which can be subjected to a tweak of a
 * certain type
 *
 * See GenericScalarShift for info about the purpose of this interface.  It
 * is an explicit representation of the implicit signature common to classes
 * like Delta::IShift, CorrelationSqueeze::IShift etc.
 *
 * For example: ParSpreadCurve implements
 * TweakableWithRespectTo<ParSpreadRhoParallelTweak>.
 */

template <class TWEAK>
struct TweakableWithRespectTo: public virtual IObject {
  static CClassConstSP const TYPE;

  virtual bool sensShift(TWEAK *shift) = 0;
  virtual string sensName(TWEAK *shift) const = 0;

  virtual ~TweakableWithRespectTo() {}

  static void load(CClassSP &clazz) {
      REGISTER_INTERFACE(TweakableWithRespectTo, clazz);
      EXTENDS(IObject);
  }
};

DRLIB_END_NAMESPACE

#endif
