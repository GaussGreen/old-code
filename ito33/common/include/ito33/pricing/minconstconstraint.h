/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/minconstconstraint.h
// Purpose:     a min constant constraint class
// Author:      (z)
// Created:     2003/10/16
// RCS-ID:      $Id: minconstconstraint.h,v 1.6 2006/06/16 18:00:51 dave Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/minconstconstraint.h
    @brief The declaration of the minimum constant constraint class.

    Class to define minimum constant constraints.
*/

#ifndef _ITO33_PRICING_MINCONSTCONSTRAINT_H_
#define _ITO33_PRICING_MINCONSTCONSTRAINT_H_

#include "ito33/pricing/constconstraint.h"

namespace ito33
{

namespace pricing
{

/// defines a minimum constant constraint
class MinConstConstraint : public ConstConstraint
{
public:

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

#endif
