/////////////////////////////////////////////////////////////////////////////
// Name:        pricing/forwardmeshmanager_fix.cpp
// Purpose:     base mesh manager for forward PDE problems
// Author:      Wang
// Created:     2004/03/10
// RCS-ID:      $Id: forwardmeshmanager_fix.cpp,v 1.2 2004/10/05 09:13:47 pedro Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////


#include "ito33/pricing/forwardmeshmanager_fix.h"

using ito33::pricing::ForwardMeshManagerFix;

void ForwardMeshManagerFix::SetupMe()
{
  ForwardMeshManager::SetupMe();

  ConstructSpaceMesh();
}
