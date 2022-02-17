/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/minconstraint.h
// Purpose:     a min constraint class
// Author:      David Pooley
// Created:     2003/08/14
// RCS-ID:      $Id: minconstraint.h,v 1.9 2006/06/16 18:00:51 dave Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/minconstraint.h
    @brief The declaration of the minimum constraints class.

    Classes to define minimum constraints (such as an American option).
*/

#ifndef _ITO33_PRICING_MINCONSTRAINT_H_
#define _ITO33_PRICING_MINCONSTRAINT_H_

#include "ito33/pricing/constraint.h"

namespace ito33 
{

namespace pricing
{

/** 
    The declaration of a minimum constraint class.

    MinConstraint defines a minimum constraint (such as an American
    option).
*/
class MinConstraint : public Constraint
{

public:

  /**
    Default constructor.
    
    Default tolerance of 1e-12.

    @param iFlagValue the value to set/check in the flag array (defaults to 1)
    @param bClearBeforeApply to clear or not the flag array before applying
           the constraint (defaults to true)
   */
  MinConstraint(int iFlagValue = 1, bool bClearBeforeApply = true) 
    : Constraint(iFlagValue, bClearBeforeApply) 
  { }

  /**
    Constructor by an array of constraint values.

    @param pdValues pointer holds the constraint values
    @param nNbValues the number of the constraint values
    @param dTolerance the relative tolerance used when a price value is
                    compared to the constraint value
    @param nMaxSize the number to define the real size of the containted
                    constraint pointer. We tolerate a input value which is 
                    less than nNbValues.
    @param iFlagValue the value to set/check in the flag array (defaults to 1)
    @param bClearBeforeApply to clear or not the flag array before applying
           the constraint (defaults to true)
    */
  MinConstraint(double *pdValues, size_t nNbValues, 
      double dTolerance = 1.e-12, size_t nMaxSize = 0,
      int iFlagValue = 1, bool bClearBeforeApply = true) 
    : Constraint(pdValues, nNbValues, dTolerance, nMaxSize,
                 iFlagValue, bClearBeforeApply) 
  { }

  // See base class 
  void Apply(double* pdPrices, int* piFlags, size_t nNbValues) const;

  // See base class
  virtual void ApplyWithoutTolerance(double* pdPrices, int* piFlags, 
    size_t nNbValues) const;

  // See base class
  virtual void ApplyPenalty(const double* pdPrices, int* piFlags, 
                            double* pdConstraints, size_t nNbValues) const;

};

} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_MINCONSTRAINT_H_

