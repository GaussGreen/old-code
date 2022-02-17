/**********************************************************************
 * Name:        ito33/license.h
 * Purpose:     class for license issue to hide CLic class
 * Author:      (z)
 * Created:     2004/6/30
 * RCS-ID:      $Id: license.h,v 1.7 2006/06/12 17:03:23 zhang Exp $
 * Copyright:   (c) 2003 Trilemma LLP
 **********************************************************************/

#ifndef _ITO33_LICENSE_H_
#define _ITO33_LICENSE_H_

#include "ito33/useexception.h"
#include "ito33/debug.h"

// define HAVE_FLEXLM for using FlexLM in any case
#define HAVE_FLEXLM  
#include "ito33/protection.h"

namespace ito33
{

/**
  License class. 

  We'd like to hide CLic class. The difference of these two classes is
  that this one can throw exception when the license requirement is not
  met. That is why the constructor of this class takes error as parameter.

  The difficulty is that an Error object is global variable. So the user
  can't use directly a global object of this class (however, we can do it
  with CLic). ihg/src/common/license.cpp is an example for using ito33::License
  class.
  */
class License
{
public:
  /**
    creates a license object by its feature name, version and error object
    showing no-license for this feature.

    @param pcFeature name of the feature to be protected
    @param pcVersion name of the version
    @param errCode no-license error code
    @param errMsg readable error message for no-license
    */
  License(const char* pcFeature, 
          const char* pcVersion,
          int errCode,
          const std::string& errMsg)
    : m_lic(pcFeature, pcVersion), 
      m_errCode(errCode),
      m_errMsg(errMsg)
  {
    #ifndef NDEBUG
      m_bCanCallIsOK = false;
    #endif
  }

  /**
  The Go function tries to get a license for the feature. If it returns a
  value Factor which is 1 when Go() succeeds. The user of this class can
  multiply his results by Factor to better protect his calculus.
  */
  void Go(void* p)
  {
    double *pd = static_cast<double*>(p);
    *pd = m_lic.Go();

    #ifndef NDEBUG
      m_bCanCallIsOK = true;
    #endif
  }

  /**
  This function throws exception when it is about a unregistered copy.

  REQUIRE: The function can ONLY be used after Go() has been called!
  */
  void Check()
  {
    ASSERT_MSG(m_bCanCallIsOK, "Calling IsOk() before calling Go().");

    if(!m_lic.IsOk())
    {
      m_lic.ErrorMsg();
      throw ito33::Exception(m_errCode, m_errMsg, "", 0, "");
    }
  }

  /**
  This function throws exception when it is about a unregistered copy.

  The function does exactly Go() + Check(). So it is slower than
  Check() function.
  */
  void Verify()
  {
    if(!m_lic.Verify())
    {
      m_lic.ErrorMsg();
      throw ito33::Exception(m_errCode, m_errMsg, "", 0, "");
    }

    #ifndef NDEBUG
      m_bCanCallIsOK = true;
    #endif
  }

private:
  CLic m_lic;
  int m_errCode;
  std::string m_errMsg;

  #ifndef NDEBUG
    bool m_bCanCallIsOK;
  #endif 

  NO_COPY_CLASS(License);
};


} // namespace ito33


#endif // _ITO33_LICENSE_H_
