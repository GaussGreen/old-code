/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/constraint.h
// Purpose:     a generic constraint class
// Author:      David Pooley
// Created:     2003/08/14
// RCS-ID:      $Id: constraint.h,v 1.13 2006/06/16 18:00:51 dave Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/constraint.h
    @brief The declaration of the constraint class.

    The base class for a single (typically min or max) constraint.  It is
    derived from Constraints so that it can be called from the solver
    code, but is still pure virtual.
 */

#ifndef _ITO33_PRICING_CONSTRAINT_H_
#define _ITO33_PRICING_CONSTRAINT_H_

#include "ito33/array.h"

#include "ito33/pricing/constraints.h"

namespace ito33
{

namespace pricing
{

/** 
  The declaration of the constraint class.

  The base class for a single (typically min or max) constraint which contains
  an array of constraint values.  It is derived from Constraints so that it can
  be called from the solver code, but is still pure virtual.
 */
class Constraint : public Constraints
{
public:

  /**
    Default constructor.
    
    Default values for all parameters.
   */
  Constraint() 
    : m_nNbValues(0), m_nMaxSize(0),
      m_dTol(1.e-12), m_bOn(false),
      m_iFlagValue(1), m_bClearBeforeApply(true)
  {
  }

  /**
    Constructor accepting flag value and clear status.
    
    Primarily used for 'Constraints' derived classes that manage multiple
    'Constraint' classes.  Default tolerance of 1e-12.

    @param iFlagValue the value to set/check in the flag array 
    @param bClearBeforeApply to clear or not the flag array before applying
           the constraint
   */
  Constraint(int iFlagValue, bool bClearBeforeApply) 
    : m_nNbValues(0), m_nMaxSize(0),
      m_dTol(1.e-12), m_bOn(false),
      m_iFlagValue(iFlagValue), m_bClearBeforeApply(bClearBeforeApply)
  {
  }

  /**
    Constructor by an array of constraint values.

    @param pdValues pointer holds the constraint values
    @param nNbValues the number of the constraint values
    @param dTolerance the relative tolerance used when a price value is
                      compared to the constraint value
    @param nMaxSize the number to define the real size of the containted
                    constraint pointer. We tolerate an input value which is 
                    less than nNbValues.
    @param iFlagValue the value to set/check in the flag array (defaults to 1)
    @param bClearBeforeApply to clear or not the flag array before applying
           the constraint (defaults to true)
   */
  Constraint(const double *pdValues, size_t nNbValues, 
             double dTolerance, size_t nMaxSize,
             int iFlagValue = 1, bool bClearBeforeApply = true);

  /**
     Virtual dtor for base class.
   */
  virtual ~Constraint()
  {
  }

  // See base class
  virtual void Apply(double* pdPrices, int* piFlags, 
                     size_t nNbValues) const = 0;

  // See base class
  virtual void ApplyWithoutTolerance(double* pdPrices, int* piFlags, 
                                     size_t nNbValues) const = 0;

  // See base class
  virtual void ApplyPenalty(const double* pdPrices, int* piFlags, 
                            double* pdConstraints, size_t nNbValues) const = 0;

  // See base class
  virtual void Get(double* pdPrices, int* piFlags, size_t nNbValues) const;

  /**
     Turns the constraint off.
   */
  void TurnOff() { m_bOn = false; }

  /**
    Turns the constraint on.
    
    The user must pay attention with function. When the function is called,
    the constraint array is assumed to have been initialized at least once.
   */
  void TurnOn();

  /**
     Checks if the constraint is on or off.
   */
  bool IsOn() const { return m_bOn; }

  /**
    Updates contained constraint values.
    
    @param pdValues pointer holds the new constraint values
    @param nNbValues the number of constraint values
   */
  void Update(double *pdValues, size_t nNbValues);

  /**
    Sets the relative tolerance.

    @param dTolerance (input) the relative tolerance

    PROMISE: update the contained absolute tolerance array 
   */
  void SetTolerance(double dTolerance);

  /**
    Gets the pointer to constraint values.

    @return pointer to constraint values
    */
  const double* GetValues() const;

  /**
    Gets the number of nodes.

    @return number of nodes
    */
  size_t GetNbValues() const;


protected:

  /// The number of stored constraint values
  size_t m_nNbValues;
  
  /**
    The maximum size of the m_pdValues array

    To avoid the re-allocation of m_pdValues when Update() function is called.
   */
  size_t m_nMaxSize;

  /// The stored constraint values
  Array<double> m_pdValues;

  /// Relative tolerance for comparison
  double m_dTol;

  /// Absolute tolerance for comparison
  Array<double> m_pdTols;

  /// Boolean to indicate if the constraint is on
  bool m_bOn;

  /// The value to set in the flag array (defaults to 1)
  const int m_iFlagValue;

  /// Whether or not to clear flags before applying (defaults to true)
  const bool m_bClearBeforeApply;

private:

  NO_COPY_CLASS(Constraint);

};


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_CONSTRAINT_H_

