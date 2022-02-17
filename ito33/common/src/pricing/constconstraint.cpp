/////////////////////////////////////////////////////////////////////////////
// Name:        pricing/constconstraint.cpp
// Purpose:     a constant constraint class
// Author:      (z)
// Created:     2003/10/16
// RCS-ID:      $Id: constconstraint.cpp,v 1.2 2004/10/05 09:13:47 pedro Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file pricing/constconstraint.cpp
    @brief The implementation of the constant constraint class.
*/

#include "ito33/pricing/constconstraint.h"

using ito33::pricing::ConstConstraint;

void ConstConstraint::Get(double* pdPrices, int* piFlags, 
    size_t nNbValues) const
{
  size_t 
    nIdx;

  if (m_bOn)
  {
    for (nIdx = 0; nIdx < nNbValues; nIdx++)
      if (piFlags[nIdx])
        pdPrices[nIdx] = m_dValue;
  }
  else
    for (nIdx = 0; nIdx < nNbValues; nIdx++)
      piFlags[nIdx] = 0;

}



