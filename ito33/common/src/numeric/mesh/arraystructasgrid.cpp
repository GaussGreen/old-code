/////////////////////////////////////////////////////////////////////////////
// Name:        numeric/mesh/arraystructasgrid.cpp
// Purpose:     class for array structured as grid with change of variable.
// Author:      Nabil
// Created:     2003/12/04
// RCS-ID:      $Id: arraystructasgrid.cpp,v 1.3 2004/10/05 09:13:46 pedro Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////
#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/autoptr.h"

#include "ito33/numeric/mesh/arraystructasgrid.h"

using ito33::AutoPtr;

using ito33::numeric::mesh::ArrayStructAsGrid;

// implement the AutoPtrDeleter for ArrayStructAsGrid class
namespace ito33
{

  ITO33_IMPLEMENT_AUTOPTR(ArrayStructAsGrid<double>);

}
