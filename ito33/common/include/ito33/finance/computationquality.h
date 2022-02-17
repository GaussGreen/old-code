/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/computationquality.h
// Purpose:     ComputationQuality enum values
// Author:      Vadim Zeitlin (extracted from qualitycontrol.h)
// Created:     2004-09-03
// RCS-ID:      $Id: computationquality.h,v 1.3 2006/01/09 10:26:00 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/computationquality.h
    @brief Constants defining the quality of computation.
 */
#ifndef _ITO33_FINANCE_COMPUTATIONQUALITY_H_
#define _ITO33_FINANCE_COMPUTATIONQUALITY_H_

namespace ito33
{

namespace finance
{

/// Quality of numerical computation.
enum ComputationQuality
{
  ComputationQuality_Standard,

  ComputationQuality_High,

  ComputationQuality_VeryHigh_But_MuchSlower

  #ifndef __CPP2ANY__
  ,
  /// noexport
  ComputationQuality_Max
  #endif
};

} // namespace finance

} // namespace ito33

#endif // _ITO33_FINANCE_COMPUTATIONQUALITY_H_
