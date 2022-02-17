/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/constconstraint.h
// Purpose:     a generic constant constraint class
// Author:      (z)
// Created:     2003/08/14
// RCS-ID:      $Id: constconstraint.h,v 1.5 2006/06/16 18:00:51 dave Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/constconstraint.h
    @brief The declaration of the ConstConstraint class.

    The base class for a single (typically min or max) constant constraint.
*/

#ifndef _ITO33_PRICING_CONSTCONSTRAINT_H_
#define _ITO33_PRICING_CONSTCONSTRAINT_H_

#include "ito33/debug.h"

#include "ito33/pricing/constraints.h"

namespace ito33
{

namespace pricing
{
/// The base class for a single (typically min or max) constant constraint.
class ConstConstraint : public Constraints
{
public:

  /// default constructor
  ConstConstraint() : m_dValue(0), m_dTol(1.e-7), m_bOn(false) {}

  /**
    Constructor.

    @param dValue the constraint value
    @param dTolerance the relative tolerance
    */
  ConstConstraint(double dValue, double dTolerance = 1.e-7);

  // See base class
  virtual void Get(double* pdPrices, int* piFlags, 
    size_t nNbValues) const;

  /// turn the constraint off
  void TurnOff() { m_bOn = false; }

  /**
    update the constraint value

    @param dValue the new constraint value
    */
  void Update(double dValue)
  {
    m_dValue = dValue;
    m_bOn = true;
  }

  /**
    set the relative tolerance

    @param dTolerance the new relative tolerance
    */
  void SetTolerance(double dTolerance)
  {
    ASSERT_MSG(dTolerance >= 0,
               "The tolerance must be positive or zero!");

    m_dTol = dTolerance;
  }

  /// Check if the constraint is on or off
  bool IsOn() const { return m_bOn; }

protected:
  
  /// Constraint value
  double m_dValue;

  /// Relative tolerance
  double m_dTol;

  /// Whether the constraint is active
  bool m_bOn;
};

} // namespace pricing

} // namespace ito33

#endif
