/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/maxconstraint.h
// Purpose:     a max constraint class
// Author:      Nabil
// Created:     2003/12/12
// RCS-ID:      $Id: maxconstraint.h,v 1.7 2006/06/16 18:00:51 dave Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/maxconstraint.h
    @brief The declaration of the maximum constraints class.

    Class to define maximum constraint
*/

#ifndef _ITO33_PRICING_MAXCONSTRAINT_H_
#define _ITO33_PRICING_MAXCONSTRAINT_H_

#include "ito33/pricing/constraint.h"

namespace ito33 
{

namespace pricing
{

/** 
    The declaration of a maximum constraint class.

    MaxConstraint defines a maximum constraint
*/
class MaxConstraint : public Constraint
{

public:

  /**
    Default constructor.
    
    Default tolerance of 1e-12.

    @param iFlagValue the value to set/check in the flag array (defaults to 1)
    @param bClearBeforeApply to clear or not the flag array before applying
           the constraint (defaults to true)
   */
  MaxConstraint(int iFlagValue = 1, bool bClearBeforeApply = true) 
    : Constraint(iFlagValue, bClearBeforeApply) { }

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
  MaxConstraint(double *pdValues, size_t nNbValues, 
      double dTolerance = 1.e-12, size_t nMaxSize = 0,
      int iFlagValue = 1, bool bClearBeforeApply = true) 
    : Constraint(pdValues, nNbValues, dTolerance, nMaxSize, 
                 iFlagValue, bClearBeforeApply) {}

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

#endif // #ifndef _ITO33_PRICING_MAXCONSTRAINT_H_

