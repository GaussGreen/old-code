/////////////////////////////////////////////////////////////////////////////
// Name:        pricing/backwardmeshmanager_fix.cpp
// Purpose:     base mesh manager for backward PDE problems
// Author:      Wang
// Created:     2004/02/11
// RCS-ID:      $Id: backwardmeshmanager_fix.cpp,v 1.7 2004/10/05 09:13:46 pedro Exp $
// Copyright:   (c) 2003-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////


#include "ito33/pricing/backwardmeshmanager_fix.h"

using ito33::pricing::BackwardMeshManagerFix;

void BackwardMeshManagerFix::SetupMe()
{
  BackwardMeshManager::SetupMe();

  ConstructSpaceMesh();
}
