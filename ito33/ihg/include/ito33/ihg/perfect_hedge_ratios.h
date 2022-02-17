/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/perfect_hedge_ratios.h
// Purpose:     PerfectHedgeRatios class
// Author:      ZHANG
// Created:     2005/02/03
// RCS-ID:      $Id: perfect_hedge_ratios.h,v 1.8 2006/07/11 14:34:05 nabil Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/ihg/common.h"

namespace ito33
{

namespace ihg
{


/**
   PerfectHedgeRatios class
    
   @nocreate
 */
class ITO33_IHG_DLLDECL PerfectHedgeRatios
{
public:

  /// constructor
  PerfectHedgeRatios() : m_dDefaultHedgeRatio(0),
                         m_dUnderlyingHedgeRatio(0),
                         m_dFXHedgeRatio(0)
  {
  }

  /**
     Default hedge ratio.

     @return Default hedge ratio
   */
  double GetDefaultHedgeRatio() const { return m_dDefaultHedgeRatio; }

  /**
     Underlying hedge ratio.

     @return Underlying hedge ratio
   */
  double GetUnderlyingHedgeRatio() const { return m_dUnderlyingHedgeRatio; }

  /**
     FX hedge ratio.

     @return FX hedge ratio
   */
  double GetFXHedgeRatio() const { return m_dFXHedgeRatio; }

  /**
     @internal  
     Underlying hedge ratio.

     @param dHedgeRatio Underlying hedge ratio

     @noexport
   */
  void SetUnderlyingHedgeRatio(double dHedgeRatio)
  {
    m_dUnderlyingHedgeRatio = dHedgeRatio;
  }

  /**
     @internal
     FX hedge ratio. Should be used only for cross-currency cases.

     @param dHedgeRatio FX hedge ratio

     @noexport
   */
  void SetFXHedgeRatio(double dHedgeRatio)
  {
    m_dFXHedgeRatio = dHedgeRatio;
  }

  /**
     @internal  
     Default hedge ratio.

     @param dHedgeRatio Default hedge ratio

     @noexport
   */
  void SetDefaultHedgeRatio(double dHedgeRatio)
  {
    m_dDefaultHedgeRatio = dHedgeRatio;
  }


private:

  double m_dDefaultHedgeRatio;
  
  double m_dUnderlyingHedgeRatio;

  double m_dFXHedgeRatio;
};


} // namespace ihg

} // namespace ito33

