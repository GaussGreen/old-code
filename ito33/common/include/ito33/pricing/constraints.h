/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/constraints.h
// Purpose:     a generic constraints class
// Author:      David Pooley
// Created:     2003/08/14
// RCS-ID:      $Id: constraints.h,v 1.12 2006/06/16 18:00:51 dave Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/constraints.h
    @brief The declaration of the constraints class.

    Definition of the constraint interface to the numerical solving
    classes.  Single constraints (e.g. American) should be derived from
    this class.  Contracts with multiple constraints (e.g. convertible bonds)
    should derive a constraint handler from this class (Constraints). The
    constraint handler can then use single constraints derived from this
    class.  See the file "pricer/docs/constraints.tex" for more information.   
 */

#ifndef _ITO33_PRICING_CONSTRAINTS_H_
#define _ITO33_PRICING_CONSTRAINTS_H_

#include "ito33/common.h"
#include <cstddef>

namespace ito33 
{

namespace pricing
{

/** 
   Base constraints class.

   Definition of the interface between the nonlinear numerical
   routines and the contract constraints.  Classes derived from 
   Constraints can represent a single constraint (e.g. American) 
   or can represent a constraint handle that manages other
   single constraints.
 */
class Constraints
{
public:
 
  /**
    Default ctor.
   */
  Constraints() { }

  /**
    Virtual dtor for base class.
   */
  virtual ~Constraints() { }

  /**
    Apply the constraint with a tolerance. 
    
    Updates both pdPrices and piFlags.

    @param pdPrices (input&output) the price values to be adjusted by the
                    constraints. Their values are modified when the constraints
                    are applied.
    @param piFlags  (output) indicates whether the constraint value
                     has been applied to the price value and the constraint
                     type.
    @param nNbValues (input) the size of pdPrices and piFlags. That is,
                     these two arrays must contain nNbValues values.

    PROMISE:
    - Set piFlag[i] to 0 if the constraint at
      node i is NOT active.
    - Set piFlag[i] to a non-zero integer if the constraint at
      node i is active.
    - Set pdPrice[i] to the constrained value at node i 
      if the constraint at node i is active.
    - Leave pdPrice[i] unchanged if the constraint at node i
      is not active.
    */
  virtual void Apply(double* pdPrices, int* piFlags, 
                     size_t nNbValues) const = 0;

  /**
    Apply the constraint without a tolerance. 
    
    Updates both pdPrices and piFlags using a tolerance of zero.

    @param pdPrices (input&output) the price values to be adjusted by the
                    constraints. Their values are modified when the constraints
                    are applied.
    @param piFlags  (output) indicates whether the constraint value
                     has been applied to the price value and the constraint
                     type.
    @param nNbValues (input) the size of pdPrices and piFlags. That is,
                     these two arrays must contain nNbValues values.

    PROMISE:
    - Set piFlag[i] to 0 if the constraint at
      node i is NOT active.
    - Set piFlag[i] to a non-zero integer if the constraint at
      node i is active.
    - Set pdPrice[i] to the constrained value at node i 
      if the constraint at node i is active.
    - Leave pdPrice[i] unchanged if the constraint at node i
      is not active.
    */
  virtual void ApplyWithoutTolerance(double* pdPrices, int* piFlags, 
                     size_t nNbValues) const = 0;

  /**
    Apply the constraint for the penalty method. 
    
    Updates both pdConstraints and piFlags without using any tolerances.
    This function is used by the penalty solver.

    @param pdPrices (input) the price values to be compared with the 
                     constraints.
    @param piFlags  (output) indicates whether the constraint value
                     should be applied and also the constraint type.
    @param pdConstraints  (output) the constraint values at nodes where the 
                          constraint is active.
    @param nNbValues (input) the size of pdPrices, piFlags and pdConstraints.
                     These three arrays must contain nNbValues values.

    PROMISE:
    - Set piFlag[i] to 0 if the constraint at node i is NOT active.
    - Set piFlag[i] to a non-zero integer if the constraint at
      node i is active.
    - Set pdConstraints[i] to the constrained value at node i 
      if the constraint at node i is active.
    - Leave pdPrices[i] unchanged for all i.
    */
  virtual void 
    ApplyPenalty(const double* pdPrices, int* piFlags, 
                 double* pdConstraints, size_t nNbValues) const = 0;

  /**
    Gets the constraint. 
    
    Updates pdPrices. piFlags may be updated when the input value is invalid.     

    @param pdPrices (output) gets out the constraint value according to the
                    piFlags value.
    @param piFlags  (input&output) Its value at each node indicates whether 
                    the constraint value should be gotten out at this node.
    @param nNbValues (input) the size of pdPrices and piFlags. That is,
                     these two arrays must contain nNbValues values.

    PROMISE:
    - modify piFlag[i] to zero only when its original value doesn't correspond
      to the state of the constraint at node i (For example, the constraint is
      inactive while the input flag value is non-zero).
    - Set pdPrice[i] to the constrained value at node i
      if piFlag[i] is non-zero.
    - Leave pdPrice[i] unchanged if piFlag[i] is zero.
    */
  virtual void Get(double* pdPrices, int* piFlags, 
                   size_t nNbValues) const = 0;
};


/**
  Apply the constraints to a greek. 
  
  Updates pdGreeks.
  
  @param pdGreeks (input&output) the greek values to be adjusted by the
                  constraints. Their values are modified when the constraints
                  are applied. For the cb, when the constraints are applied,
                  all the greeks become zero.
  @param piFlagConstraints  (input) indicates whether the constraint value
                   has been applied to the price value and the constraint
                   type.
  @param nNbValues (input) the size of pdGreeks and piFlagConstraints. That is,
                   these two arrays must contain nNbValues values.
*/
void ApplyConstraintsToGreek(double* pdGreeks, const int* piFlagConstraints, 
                             size_t nNbValues);

} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_CONSTRAINTS_H_

