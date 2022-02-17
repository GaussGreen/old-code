/////////////////////////////////////////////////////////////////////////////
// Name:        pricing/backwardmeshmanager_withspacemesh.cpp
// Purpose:     base mesh manager for backward PDE problems with space mesh
// Created:     2004/09/02
// RCS-ID:      $Id: backwardmeshmanager_withspacemesh.cpp,v 1.4 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/pricing/backwardmeshmanager_withspacemesh.h"

using ito33::pricing::BackwardMeshManagerWithSpaceMesh;

void BackwardMeshManagerWithSpaceMesh::SetupMe()
{
  BackwardMeshManager::SetupMe();

  ConstructSpaceMesh();
}
