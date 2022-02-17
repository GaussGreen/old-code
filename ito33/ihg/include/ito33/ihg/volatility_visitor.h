/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/volatility_visitor.h
// Purpose:     Visitor for Volatility-derived classes
// Author:      Yann d'halluin
// Created:     2004-18-07
// RCS-ID:      $Id: volatility_visitor.h,v 1.4 2005/09/05 14:00:43 zhang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/ihg/volatility_visitor.h
    @brief Visitor class accepting different Volatility kinds.

    This is the central part of the implementation of the visitor pattern for
    the volatility. To do something different depending on the exact type of
    a volatility object you have to define a new class inheriting from this one
    and do whatever is required in its methods.
 */

#ifndef _ITO33_IHG_VOLATILITY_VISITOR_H_
#define _ITO33_IHG_VOLATILITY_VISITOR_H_

namespace ito33
{

namespace ihg
{

  class VolatilityFlat;
  class VolatilityCombo;
  class VolatilityTimeOnly;
  class VolatilityPower;
  class VolatilityTanh;

/**
    Volatility visitor.

    An object of a class derived from this one must be passed by the user code
    to any functions which may access heterogeneous collections of Volatility,
    e.g. XML::Reader in IHG currently.

    Using a visitor might be less usual than querying the data directly for
    but it allows us to not lose the type information about the volatility
    and keep the code maintainable and extensible.
 */
class VolatilityVisitor
{
public:
   /// Virtual dtor for any base class
  virtual ~VolatilityVisitor() { }

  /// Flat volatility
  virtual void OnVolatilityFlat(const VolatilityFlat& volFlat) = 0; 

  
  /// Time only volatility
  virtual void OnVolatilityTimeOnly(const VolatilityTimeOnly& volTimeOnly) = 0;

  /// power volatility
  virtual void OnVolatilityPower(const VolatilityPower& volPower) = 0;

  /// Tanh volatility
  virtual void OnVolatilityTanh(const VolatilityTanh& volTanh) = 0;

  // -----  (semi) internal volatility classes  ------------------------
  /// Volatility combo
  virtual void OnVolatilityCombo(const VolatilityCombo&)
  {
    FAIL ( " Not supported yet. Please implement" );
  }
};

} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_VOLATILITY_VISITOR_H_
