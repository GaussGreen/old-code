/**********************************************************************
 * Name:        ito33/protection.h
 * Purpose:     class for license issue
 * Author:      (z)
 * Created:     10 march.03
 * RCS-ID:      $Id: protection.h,v 1.32 2006/05/27 20:13:16 zhang Exp $
 * Copyright:   (c) 2003 Trilemma LLP
 **********************************************************************/

#ifndef _ITO33_PROTECTION_H_
#define _ITO33_PROTECTION_H_

/**
    \file   ito33/protection.h
    \brief  Licence verification classes.

    This header always defines the CLic class but the class may be empty if the
    protection is disabled which happens if ITO33_NO_PROTECTION is predefined
    which, in turn, is done by default in debug builds (NDEBUG not defined).
    However it is also possible to predefine _ITOSELFPROTECT macro which
    overrides this and which does activate licence checking no matter what.
 */

#include "ito33/common.h"
#include "ito33/getenv.h"

/**
    @def ITO33_NO_PROTECTION

    Define this macro to disable licence verification.

    This is done by default in debug builds.
 */
#ifndef NDEBUG
  #define ITO33_NO_PROTECTION
#endif

/**
    @def _ITOSELFPROTECT

    Define this macro to force the licence checking code to be enabled.
 */
#ifndef _ITOSELFPROTECT
  #ifndef ITO33_NO_PROTECTION
      #define _ITOSELFPROTECT
  #endif
#endif 

#ifdef _ITOSELFPROTECT

#include "ito33/common.h"

#include "ito33/beforestd.h"
#include <string>
#include "ito33/afterstd.h"
#include "ito33/thread.h"

// we can use either FlexLM or our own old protection scheme, HAVE_FLEXLM
// defines whether to use FlexLM
//
// under Unix it is set by configure so it is enough to include the config.h
// file where it is defined but under Windows we have no way to autodetect
// this, so assume we do have FlexLM by default
#if !defined(HAVE_ITO33_CONFIG_H) && defined(_WIN32) && !defined(HAVE_FLEXLM)
    #define HAVE_FLEXLM
#endif

#ifdef HAVE_FLEXLM
    #include <lmpolicy.h>

    // include the required FlexLM libs implicitly
    #ifdef _MSC_VER
        #pragma comment(lib, "lmgr_md.lib")
        #pragma comment(lib, "flock_md.lib")
        #pragma comment(lib, "lm_new_md.lib")
        #pragma comment(lib, "netapi32.lib")
        #pragma comment(lib, "comctl32.lib")
        #pragma comment(lib, "oldnames.lib")
        #pragma comment(lib, "wsock32.lib")
    #endif // _MSC_VER
#endif // HAVE_FLEXLM

namespace ito33
{

/**
    CLic class manages the license for a specified feature. The usage of this 
    class is multi-thread safe.
 */
class CLic
{
public:

#  ifdef HAVE_FLEXLM
  /// TODO
  CLic(const char* _pcFeature, const char* _pcVersion):
      mdFact(1.3), mbOk(false), mbFirstTime(true), mbFirstTimeMsg(true),
      mstrFeature(_pcFeature), mstrVersion(_pcVersion)
        {
        }
  /// TODO
  ~CLic()
    {
    if(mbOk)
      lp_checkin(mp_lphandle);
    }
#  endif

  /**
      Tries to get a license for the feature. If it returns a value Factor
      which is 1 when Go() succeeds. The user of this class can multiply his
      results by Factor to better protect his calculus.
   */
  double Go(void);

  /**
      Returns false when it is about a unregistered copy.

      The function can ONLY be used after Go() has been called!
   */
  bool IsOk(void);

  /**
      Returns false when it is about a unregistered copy.

      The function does exactly Go() + IsOk(). So it is slower than IsOk()
      function.
   */
  bool Verify(void);

  /**
      This function displays the error message. The function is usually called
      when IsOk() return false
   */
  int ErrorMsg(void);

private:
  CriticalSection mobjSC;
  std::string mstrFeature;
  bool mbOk;
  bool mbFirstTime;
  bool mbFirstTimeMsg;
  double mdFact;

#  ifdef HAVE_FLEXLM
    LP_HANDLE *mp_lphandle;
    std::string mstrVersion;
#  endif

  NO_COPY_CLASS(CLic);
};

inline double CLic::Go()
{
  Lock<CriticalSection> lock(mobjSC);

  if(mbFirstTime)
  {
    int iFact = 2;

    // the directory containing the registration information
    std::string strRegDir(GetEnv("ITO33InstDir"));

    if ( strRegDir.empty() )
    {
#ifdef _WIN32
      strRegDir = "c:\\ito33";
#else // !_WIN32
      strRegDir = "/etc/ito33";
#endif
    }

    // waiting for the FileName class...
#ifdef _WIN32
    strRegDir += '\\';
#else // !_WIN32
    strRegDir += '/';
#endif

    strRegDir += "reg";

#ifdef HAVE_FLEXLM
    /*
        We use LM_MANUAL_HEARTBEAT here to prevent FlexLM from using a
        background thread for automatic heartbeat notifications because they
        result in mysterious crashes when the DLL is being unloaded. This is
        presumably a bug in FlexLM because it happens after our code had
        been unloaded but the exact details are unknown. Disabling the heart
        beat thread is by far the simplest thing to do and we don't need it
        anyhow [for now] so we don't lose anything.
     */
    iFact = lp_checkout
            (
              LPCODE,
              LM_RESTRICTIVE | LM_MANUAL_HEARTBEAT,
              const_cast<char *>(mstrFeature.c_str()),
              const_cast<char *>(mstrVersion.c_str()),
              1,
              const_cast<char *>(strRegDir.c_str()),
              &mp_lphandle
            );
#endif // HAVE_FLEXLM/!HAVE_FLEXLM

    mdFact = 0.1 * iFact + 1.;

    mbFirstTime = false;
    if(iFact == 0)
      mbOk = true;
  }

  return mdFact;
}


inline bool CLic::IsOk()
{
  return mbOk;
}


inline bool CLic::Verify()
{
  Go();
  return mbOk;
}


inline int CLic::ErrorMsg()
{
  Lock<CriticalSection> lock(mobjSC);

  if(mbFirstTimeMsg)
  {
#   ifdef HAVE_FLEXLM
    {
      lp_perror(mp_lphandle, "Failed to acquire the license");
    }
#   endif
    
    mbFirstTimeMsg = false;

    return 0;
  }

  return 1;
}


} // namespace ito33

#else // !_ITOSELFPROTECT

// Declare a trivial CLic class just to avoid peppering the user code with
// #ifdef _ITOSELFPROTECT
namespace ito33
{

class CLic
{
public:

  CLic(const char *, const char *) { }
  static bool Verify() { return true; }  
  static bool IsOk() { return true; }  
  static int ErrorMsg() {return 0;}
  static double Go() {return 1;}
};

} // namespace ito33

#endif // _ITOSELFPROTECT/!_ITOSELFPROTECT

#endif // _ITO33_PROTECTION_H_
