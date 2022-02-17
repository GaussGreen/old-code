/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/problemtype.h
// Purpose:     Enumeration for the problem type (backward, forward)
// Author:      ICARE
// Created:     2003/10/27
// RCS-ID:      $Id: problemtype.h,v 1.3 2004/10/05 09:13:39 pedro Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/problemtype.h
    @brief Enumeration for the problem type (backward, forward)

    Define the type of problem. e.g forward or backward equations
 */

#ifndef _ITO33_PRICING_PROBLEMTYPE_H_
#define _ITO33_PRICING_PROBLEMTYPE_H_

namespace ito33 
{

namespace pricing
{

/// Enumeration for the type of problem (forward, backward)
enum ProblemType
{
  ProblemType_Forward,
  ProblemType_Backward,

  ProblemType_Max
};

} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_PROBLEMTYPE_H_

