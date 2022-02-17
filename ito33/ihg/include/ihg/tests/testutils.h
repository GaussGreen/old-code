/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/testutils.h
// Purpose:     Helper functions for coreinterface tests
// Created:     2006/03/23
// RCS-ID:      $Id: testutils.h,v 1.3 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_IHG_TESTS_TESTUTILS_H_
#define _ITO33_IHG_TESTS_TESTUTILS_H_

#include "ito33/dlldecl.h"
#include "ito33/sharedptr.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL SessionData;
}

namespace ihg
{

shared_ptr<finance::SessionData> MakeSessionData(double dSpot);



} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_TESTS_TESTUTILS_H_
