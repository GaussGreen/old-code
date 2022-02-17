/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/parbond/parbondnumoutput.cpp
// Purpose:     implementation of parbond numeric output class 
// Author:      Based on ICARE version
// RCS-ID:      $Id: parbondnumoutput.cpp,v 1.3 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/modeloutput.h"

#include "ito33/numeric/interpolation.h"
#include "ito33/numeric/deltagamma.h"
#include "ito33/numeric/domain_fixedspacemesh.h"
#include "ito33/numeric/surfacegeneral.h"

#include "ito33/pricing/parbondparams.h"
#include "ito33/pricing/parbondmeshmanager.h"

#include "ihg/parbondnumoutput.h"
#include "ihg/parbondinstdata.h"
#include "ihg/backwardinstdata.h"

#include "ito33/ihg/modeloutput.h"

using ito33::shared_ptr;
using ito33::AutoPtr;

using namespace ito33::numeric;
using namespace ito33;
using namespace ito33::ihg;

// implement the AutoPtrDeleter for ParBondNumOutput
namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(ihg::ParBondNumOutput);
}

void ParBondNumOutput::Init(ParBondInstData &instdata)
{
  BackwardNumOutput::InitWithFixMesh(instdata, m_params, false);
}

void ParBondNumOutput::SaveSurface(BackwardInstData& instdata, double dTime)
{
  ASSERT_MSG
  (
    dynamic_cast<DomainFixedSpaceMesh*>(m_pDomain.get()), 
    "The type of Domain in ParBondNumOutput should be DomainFixedSpaceMesh."
  );

  DomainFixedSpaceMesh& 
    domainFixed = static_cast<DomainFixedSpaceMesh&>(*m_pDomain);

  BackwardNumOutput::SaveSurfaceWithFixMesh(instdata, domainFixed, dTime);
}


// Return the requested data to the user
shared_ptr<ModelOutput> ParBondNumOutput::GetOutput()
{
  // Construct the output class
  shared_ptr<ModelOutput> pOutput (new ModelOutput);

  BackwardNumOutput::GetOutput(pOutput.get());

  return pOutput;
}
