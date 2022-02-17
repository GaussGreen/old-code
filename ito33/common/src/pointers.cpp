/////////////////////////////////////////////////////////////////////////////
// Name:        pointers.cpp
// Purpose:     non-inline parts of AutoPtr, SharedPtr &c
// Author:      Vadim Zeitlin
// Created:     2004-07-21
// RCS-ID:      $Id: pointers.cpp,v 1.4 2005/05/25 14:48:28 zeitlin Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"

using namespace ito33;
using namespace ito33::Policy;
using namespace ito33::Private;

// the ctor and dtor for SharedPtrCount must be out of line to avoid problems
// with mixing DLL heaps

/* static */
SharedPtrCount<Policy::MTUnsafe> *MTUnsafe::CreateSharedCount()
{
  return new SharedPtrCount<Policy::MTUnsafe>();
}

/* static */
void MTUnsafe::Destroy(Private::SharedPtrCount<MTUnsafe> *p)
{
  delete p;
}

/* static */
SharedPtrCount<Policy::MTSafe> *MTSafe::CreateSharedCount()
{
  return new SharedPtrCount<Policy::MTSafe>();
}

/* static */
void MTSafe::Destroy(Private::SharedPtrCount<MTSafe> *p)
{
  delete p;
}

