/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/constraint_cbcall.h
// Purpose:     cb call constraint class, is pseudo-maxconstraint
// Author:      ZHANG Yunzhi
// Created:     2005/01/14
// RCS-ID:      $Id: constraint_cbcall.h,v 1.5 2006/06/16 18:00:51 dave Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/constraint_cbcall.h
    @brief The declaration of the cb call constraint class.

    Class to define cb call constraint
*/

#ifndef _ITO33_PRICING_CONSTRAINT_CBCALL_H_
#define _ITO33_PRICING_CONSTRAINT_CBCALL_H_

#include "ito33/pricing/constraint.h"

namespace ito33 
{

namespace pricing
{

/// ConstraintCBCall defines cb call constraint
class ConstraintCBCall : public Constraint
{

public:

  /**
    Default constructor.
    
    Default tolerance of 1e-12.

    @param iFlagValue the value to set/check in the flag array (defaults to 1)
    @param bClearBeforeApply to clear or not the flag array before applying
           the constraint (defaults to true)
   */
  ConstraintCBCall(int iFlagValue = 1, bool bClearBeforeApply = true) 
    : Constraint(iFlagValue, bClearBeforeApply) 
  { }

  /**
    Constructor by an array of constraint values.

    @param pdValues pointer holds the constraint values
    @param nNbValues the number of the constraint values
    @param nIdxStartConversion index of the point where the 
                               instrument is called to be converted.
    @param dTolerance the relative tolerance used when a price value is
                    compared to the constraint value
    @param nMaxSize the number to define the real size of the containted
                    constraint pointer. We tolerate a input value which is 
                    less than nNbValues.
    @param iFlagValue the value to set/check in the flag array (defaults to 1)
    @param bClearBeforeApply to clear or not the flag array before applying
           the constraint (defaults to true)
   */  
  ConstraintCBCall(double *pdValues,
                   size_t nNbValues,
                   size_t nIdxStartConversion,
                   double dTolerance = 1.e-12,
                   size_t nMaxSize = 0,
                   int iFlagValue = 1, 
                   bool bClearBeforeApply = true) 
    : Constraint(pdValues, nNbValues, dTolerance, nMaxSize, iFlagValue, 
                 bClearBeforeApply),
      m_nIdxStartConversion(nIdxStartConversion)
  {
    if (m_nIdxStartConversion > m_nNbValues)
      m_nIdxStartConversion = m_nNbValues;
  }

  // See base class
  void Apply(double* pdPrices, int* piFlags, size_t nNbValues) const;
  
  // See base class
  virtual void ApplyWithoutTolerance(double* pdPrices, int* piFlags, 
    size_t nNbValues) const;

  // See base class
  virtual void ApplyPenalty(const double* pdPrices, int* piFlags, 
                            double* pdConstraints, size_t nNbValues) const;

  /**
    Updates constraint data.

    @param pdValues pointer holds the constraint values
    @param nNbValues the number of the constraint values
    @param nIdxStartConversion index of the point where the 
                               instrument is called to be converted.
   */
  void Update(double *pdValues,
              size_t nNbValues,
              size_t nIdxStartConversion);

private:

  /// Index of the point where the called instrument is to be converted.
  size_t m_nIdxStartConversion;
};

} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_CONSTRAINT_CBCALL_H_

