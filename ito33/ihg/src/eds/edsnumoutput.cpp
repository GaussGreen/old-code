/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/eds/edsnumoutput.cpp
// Purpose:     implementation of  numeric output class for EDS
// Created:     2005/01/26
// RCS-ID:      $Id: edsnumoutput.cpp,v 1.3 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"

#include "ito33/numeric/domain_fixedspacemesh.h"
#include "ito33/numeric/surfacegeneral.h"

#include "ito33/pricing/edsparams.h"
#include "ito33/pricing/edsmeshmanager.h"

#include "ihg/edsnumoutput.h"
#include "ihg/edsinstdata.h"

#include "ito33/ihg/modeloutput.h"

using namespace ito33::numeric;
using namespace ito33;
using namespace ito33::ihg;

// implement the AutoPtrDeleter for EDSNumOutput
namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(ihg::EDSNumOutput);
}

void EDSNumOutput::Init(EDSInstData &instdata)
{
  BackwardNumOutput::InitWithFixMesh(instdata, m_params, false);
}

void EDSNumOutput::SaveSurface(BackwardInstData& instdata, double dTime)
{
  ASSERT_MSG
  (
    dynamic_cast<DomainFixedSpaceMesh*>(m_pDomain.get()), 
    "The type of Domain in CDSNumOutput should be DomainFixedSpaceMesh."
  );

  DomainFixedSpaceMesh& 
    domainFixed = static_cast<DomainFixedSpaceMesh&>(*m_pDomain);

  BackwardNumOutput::SaveSurfaceWithFixMesh(instdata, domainFixed, dTime);
}


// Return the requested data to the user
shared_ptr<ihg::ModelOutput> EDSNumOutput::GetOutput()
{
  // Construct the output class
  shared_ptr<ihg::ModelOutput> pOutput( new ihg::ModelOutput );

  BackwardNumOutput::GetOutput( pOutput.get() );

  return pOutput;
}
