/**
 * \file TweakableWith.hpp
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
 * NOTE NOTE NOTE: this interface may steadily be replaced with
 * ITweakableWithRespectTo, as sensitivities are ported over to the
 * "declarative" framework.  For an overview see IRiskQuantityFactory.
 *
 * For example: CurrencyBasis implements
 * TweakableWith<CurrencyBasisSpreadLevelTweak>.
 */

template <class TWEAK>
struct TweakableWith: public virtual IObject {
  static CClassConstSP const TYPE;

  virtual bool sensShift(TWEAK *shift) = 0;
  virtual string sensName(TWEAK *shift) const = 0;

  virtual ~TweakableWith() {}

  static void load(CClassSP &clazz) {
      REGISTER_INTERFACE(TweakableWith, clazz);
      EXTENDS(IObject);
  }
};

DRLIB_END_NAMESPACE

#endif
