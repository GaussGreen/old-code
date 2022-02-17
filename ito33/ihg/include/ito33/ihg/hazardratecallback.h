/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/hazardratecallback.h
// Purpose:     Callback hazard rate class.
// Author:      Nabil
// Created:     2003/11/06
// RCS-ID:      $Id: hazardratecallback.h,v 1.20 2005/12/02 15:52:57 zhang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////


#ifndef _ITO33_IHG_HAZARDRATECALLBACK_H_
#define _ITO33_IHG_HAZARDRATECALLBACK_H_

#include "ito33/ihg/callbackfunction.h"
#include "ito33/ihg/hazardrate.h"

namespace ito33
{

namespace ihg
{

 class HazardRateVisitor;

/**
    Callback hazard rate class, using a user defined function of type 
    CallBackFunction

    @nocreate
 */
class HazardRateCallBack : public HazardRate
{
public:
  
  /**
    usuel constructor

    @param Func function of type CallBackFunction
    @param iUserData user data
    */
  HazardRateCallBack(CallBackFunction Func, int iUserData)
  {
    m_FuncHazardRate = Func;
    m_iUserData = iUserData;
  }

  /// dtor
  ~HazardRateCallBack() { }
  
  void GetHazardRates(double dTime,
                      const double *pdS,
                      double *pdValues, 
                      size_t nNumber) const
  {
    m_FuncHazardRate(dTime, pdS, pdValues, nNumber, m_iUserData);
  }
    
  virtual void Dump(ito33::XML::Tag& tagParent) const;

  virtual void Visit(HazardRateVisitor& visitor) const;
  
  virtual 
    void GetModelParameters(finance::ModelParametersConsumer&) const
  {
    // can't do anything
  }

  
private:

  /// user defined function
  CallBackFunction m_FuncHazardRate;
  
  /// user data
  int m_iUserData;
  
}; // class HazardRateCallBack


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_HAZARDRATECALLBACK_H_
