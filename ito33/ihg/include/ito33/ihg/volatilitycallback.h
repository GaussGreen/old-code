/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/volatiliycallback.h
// Purpose:     implementation for callback volatility class
// Author:      Wang
// Created:     2004/03/17
// RCS-ID:      $Id: volatilitycallback.h,v 1.14 2006/08/20 09:47:41 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/ihg/volatilitycallback.h
 */

#ifndef _ITO33_IHG_VOLATILITYCALLBACK_H_
#define _ITO33_IHG_VOLATILITYCALLBACK_H_

#include "ito33/ihg/callbackfunction.h"
#include "ito33/ihg/volatility.h"

namespace ito33
{

namespace ihg
{

  class VolatilityVisitor;

/**
    Base class from which callback implementations of Volatility should be
    derived.

    We currently have VolatilityCallBack which uses a C callback but a better
    solution in C++ is to derive a new class from this one and implement its
    DoGetValues() method. This is more flexible and typesafe than using C
    callback and doesn't require using "user data" to do anything nontrivial.
 */
class VolatilityCallBackBase : public Volatility
{
public:
  /**
      Default ctor.

      The shift parameter is used by Perturb().

      @param dShift the amount by which to offset the values returned by
                    GetVols()
   */
  VolatilityCallBackBase(double dShift = 0.) : m_dShift(dShift)
  {
  }

  // default dtor is ok

  /**
     Get the values of the volatility at a given time for an array of spots.

     This method forwards to DoGetValues() after taking into account a possible
     shift.
   */
  virtual void GetVols(double dTime,
                       const double *pdS,
                       double *pdVols,
                       size_t nNbS) const;

  /**
      Return a copy of this object shifted by the specified amount.

      The implementation of this method uses Clone() pure virtual which must be
      implemented in the derived class.
   */
  virtual shared_ptr<Volatility> Perturb(double dShift) const
  {
    return shared_ptr<Volatility>(Clone(dShift));
  }

  virtual void Dump(ito33::XML::Tag& tagParent) const;

  virtual void Visit(VolatilityVisitor& visitor) const;


protected:

  // this must implemented for Perturb() to work: simply return a new object
  virtual VolatilityCallBackBase *Clone(double dShift) const = 0;

  // and this is the real work function which must be implemented in the
  // derived class
  virtual void DoGetValues(double dTime,
                           const double *pdS,
                           double *pdVols,
                           size_t nNbS) const = 0;


  // the shift from Perturb()
  const double m_dShift;

  NO_COPY_CLASS(VolatilityCallBackBase);
};

/// call back volatility function
class VolatilityCallBack : public VolatilityCallBackBase
{

public:
  /**
     Ctor constructs a volatility by using a callback function

     The end user should always use the default value 0 for dShift.

     @param Func the user supplied call back function
     @param iUserData the data that the user needs for the call back function
     @param dShift the amount that the volatilty got shifted
   */
  VolatilityCallBack(CallBackFunction Func, int iUserData, double dShift = 0.)
     : VolatilityCallBackBase(dShift),
       m_FuncVolatility(Func),
       m_iUserData(iUserData)
  {
  }

  // Default dtor is ok

  // GetVolsSquared: the same as in base class.

protected:
  // implement base class pure virtuals
  virtual VolatilityCallBackBase *Clone(double dShift) const
  {
    return new VolatilityCallBack(m_FuncVolatility, m_iUserData, dShift);
  }

  virtual void DoGetValues(double dTime,
                           const double *pdS,
                           double *pdVols,
                           size_t nNbS) const
  {
    m_FuncVolatility(dTime, pdS, pdVols, nNbS, m_iUserData);
  }

private:
  CallBackFunction m_FuncVolatility;

  int m_iUserData;

  NO_COPY_CLASS(VolatilityCallBack);
}; // class VolatilityCallBack


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_VOLATILITYCALLBACK_H_

