/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/cds/cdsnumoutput.cpp
// Purpose:     implementation of cds numeric output class 
// Author:      Based on ICARE version
// RCS-ID:      $Id: cdsnumoutput.cpp,v 1.29 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/numeric/interpolation.h"
#include "ito33/numeric/deltagamma.h"
#include "ito33/numeric/domain_fixedspacemesh.h"
#include "ito33/numeric/surfacegeneral.h"

#include "ito33/pricing/cdsparams.h"
#include "ito33/pricing/cdsmeshmanager.h"

#include "ito33/ihg/modeloutput.h"

#include "ihg/cdsnumoutput.h"
#include "ihg/cdsinstdata.h"
#include "ihg/backwardinstdata.h"

using ito33::shared_ptr;
using ito33::AutoPtr;

using namespace ito33::numeric;

using namespace ito33::ihg;

// implement the AutoPtrDeleter for CDSNumOutput
namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(ihg::CDSNumOutput);
}

void CDSNumOutput::Init(CDSInstData &instdata)
{
  BackwardNumOutput::InitWithFixMesh(instdata, m_params, false);
}

void CDSNumOutput::SaveSurface(BackwardInstData& instdata, double dTime)
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
shared_ptr<ModelOutput> CDSNumOutput::GetOutput()
{
  // Construct the output class
  shared_ptr<ModelOutput> pOutput (new ihg::ModelOutput());

  BackwardNumOutput::GetOutput(pOutput.get());

  return pOutput;
}
