/////////////////////////////////////////////////////////////////////////////
// Name:        com/coclass.cpp
// Purpose:     COM coclasses support stuff
// Author:      Vadim Zeitlin
// Created:     29.01.03
// RCS-ID:      $Id: coclass.cpp,v 1.5 2004/10/05 09:13:42 pedro Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/com/coclass.h"

using namespace ito33::COM;

CoClassInfo *CoClassInfo::ms_first = NULL;

/* static */
const CoClassInfo *CoClassInfo::Find(const CLSID& clsid)
{
  // do linear search in the linked list (ok here as the list is small)
  //
  // NB: we don't use std::find() here because this would force us to
  //     define iterators for our linked list (and we can't use std::list,
  //     of course, because it allocates memory which we don't want to do in
  //     these static objects) and it just seems as too much trouble for too
  //     little gain
  const CoClassInfo *coclass = CoClassInfo::GetFirst();
  while ( coclass )
  {
    if ( coclass->GetClsid() == clsid )
    {
      // found the right coclass
      break;
    }

    coclass = coclass->GetNext();
  }

  return coclass;
}

