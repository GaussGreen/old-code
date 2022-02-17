/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/cbconstraints.h
// Purpose:     a cb constraints class
// Author:      Nabil
// Created:     2003/12/12
// RCS-ID:      $Id: cbconstraints.h,v 1.18 2006/07/03 16:08:15 nabil Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/cbconstraints.h
    @brief The declaration of the CBConstraints class.
  
 */

#ifndef _ITO33_PRICING_CBCONSTRAINTS_H_
#define _ITO33_PRICING_CBCONSTRAINTS_H_

#include "ito33/constants.h"
#include "ito33/pricing/constraints.h"
#include "ito33/pricing/minconstraint.h"
#include "ito33/pricing/constraint_cbcall.h"

namespace ito33 
{

namespace pricing
{

/** 
   CB constraints class.

   Definition of the interface between the nonlinear numerical
   routines and the cb constraints. The CBConstraints class is derived
   from the generic class Constraints.
 */
class CBConstraints : public Constraints
{
public:

  /**
    The possible constraint states.
   */
  enum State
  {
    State_none,
    State_put,
    State_conv,
    State_call,
    
    State_Max
  };
 
  /**
    Ctor.

    Construct the individual constraints with the appropriate flag values.
    Also, indicate that the constraints do not clear the flags so the same 
    flag array can be passed to all the constraints in sequence.
   */
  CBConstraints()
    : m_put(State_put, false),
      m_call(State_call, false),
      m_conv(State_conv, false)
  {
  } 

  /**
    Dtor.
   */
  virtual ~CBConstraints() 
  { 
  }

  // See base class
  virtual void Apply(double* pdPrices, int* piFlags, size_t nNbValues) const;
  
  // See base class
  virtual void ApplyWithoutTolerance(double* pdPrices, int* piFlags, 
                                     size_t nNbValues) const;

  // See base class
  virtual void ApplyPenalty(const double* pdPrices, int* piFlags, 
                            double* pdConstraints, size_t nNbValues) const;

  // See base class
  virtual void Get(double* pdPrices, int* piFlags,
                   size_t nNbValues) const;  

  /**
     Updates the put constraint values.

     @param pdValues the new put consraint values
     @param nNbValues the number of values
   */
  void UpdatePut(double *pdValues, size_t nNbValues);

  /**
     Updates the conversion constraint values.

     @param pdValues the new conversion consraint values
     @param nNbValues the number of values
   */
  void UpdateConv(double *pdValues, size_t nNbValues);

  /**
     Updates the call constraint values.

     @param pdValues the new call consraint values
     @param nNbValues the number of values
     @param nIdxStartConversion index of the point where the 
                                instrument is called to be converted.
   */
  void UpdateCall(double *pdValues,
                  size_t nNbValues,
                  size_t nIdxStartConversion);

  /**
    Turns off the call constraint.
   */
  void DisableTheCall() { m_call.TurnOff(); }
  
  /**
    Turns off the conversion constraint.
   */
  void DisableTheConv() { m_conv.TurnOff(); }
  
  /**
    Turns off the put constraint.
   */
  void DisableThePut() { m_put.TurnOff(); }

  
protected:

  /// The put constraint
  MinConstraint m_put;

  /// The conversion constraint
  MinConstraint m_conv;

  /// The call constraint
  ConstraintCBCall m_call;

private:

  NO_COPY_CLASS(CBConstraints);

};

} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_CBCONSTRAINTS_H_

