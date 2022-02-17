/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/qualitycontrol.h
// Purpose:     Class defines user end parameters for quality control
// Author:      ZHANG Yunzhi
// Created:     2004-05-15
// RCS-ID:      $Id: qualitycontrol.h,v 1.17 2005/03/31 10:24:34 pedro Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/qualitycontrol.h
    @brief Class defines user end parameters for quality control
 */
#ifndef _ITO33_FINANCE_QUALITYCONTROL_H_
#define _ITO33_FINANCE_QUALITYCONTROL_H_

#include "ito33/finance/computationquality.h"

#include "ito33/dlldecl.h"

namespace ito33
{

namespace finance
{

/**
    Determines the speed/precision tradeoff of the computations.
 */
class ITO33_DLLDECL QualityControl
{
public:
  
  /**
      Default constructor sets the computation quality to standard.
   */
  QualityControl() : m_computationQuality(ComputationQuality_Standard) 
  {
  }

  // default copy ctor, operator=() and dtor are ok

  /**
     The requested computation quality

     @param computationQuality user requested computation quality
   */
  void SetComputationQuality(ComputationQuality computationQuality)
  {
    m_computationQuality = computationQuality;
  }

  /**
     The requested computation quality

     @return user requested computation quality
   */
  ComputationQuality GetComputationQuality() const
  {
    return m_computationQuality;
  }

private:

  /// user requested quality of numerical computation
  ComputationQuality m_computationQuality;

}; // class QualityControl


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_QUALITYCONTROL_H_

