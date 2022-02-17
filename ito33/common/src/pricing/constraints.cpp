/////////////////////////////////////////////////////////////////////////////
// Name:        pricing/constraints.cpp
// Purpose:     a generic constraints class
// Author:      Nabil
// Created:     2004/01/02
// RCS-ID:      $Id: constraints.cpp,v 1.4 2005/09/15 16:15:46 nabil Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file pricing/constraints.cpp
*/

#include "ito33/pricing/constraints.h"
#include "ito33/autoptr.h"

namespace ito33
{

  // implement the autoptrdeleter of the Constraints class
  ITO33_IMPLEMENT_AUTOPTR(pricing::Constraints);

} // namespace ito33

void ito33::pricing::ApplyConstraintsToGreek(double* pdGreeks, 
                             const int* piFlagConstraints, 
                             size_t nNbValues)
{
  size_t
    nIdxSpot;

  for(nIdxSpot = 0; nIdxSpot < nNbValues; nIdxSpot++)
  {
    if(piFlagConstraints[nIdxSpot])
      pdGreeks[nIdxSpot] = 0;
  }
}
