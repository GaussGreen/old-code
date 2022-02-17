/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/varianceswap/varianceswapnumoutput.cpp
// Purpose:     implementation of VarianceSwapNumOutput class 
// Created:     2006/02/21
// RCS-ID:      $Id: varianceswapnumoutput.cpp,v 1.3 2006/08/20 09:31:05 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"

#include "ito33/numeric/domain_fixedspacemesh.h"

#include "ito33/pricing/varianceswapparams.h"

#include "ito33/ihg/modeloutput.h"

#include "ihg/backwardnumoutput.h"
#include "ihg/varianceswapnumoutput.h"
#include "ihg/varianceswapinstdata.h"

using ito33::shared_ptr;
using ito33::AutoPtr;

using namespace ito33::numeric;

using namespace ito33::ihg;

// implement the AutoPtrDeleter for VarianceSwapNumOutput
namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(ihg::VarianceSwapNumOutput);
}

void VarianceSwapNumOutput::Init(VarianceSwapInstData& instdata)
{  
  BackwardNumOutput::InitWithFixMesh(instdata, m_params, false);
}

void VarianceSwapNumOutput::SaveSurface(BackwardInstData& instdata, 
                                        double dTime)
{
  ASSERT_MSG
  (
    dynamic_cast<DomainFixedSpaceMesh*>(m_pDomain.get()), 
    "The type of Domain in VarianceSwapNumOutput should be DomainFixedSpaceMesh."
  );

  DomainFixedSpaceMesh& 
    domainFixed = static_cast<DomainFixedSpaceMesh&>(*m_pDomain);

  BackwardNumOutput::SaveSurfaceWithFixMesh(instdata, domainFixed, dTime);
}


// Return the requested data to the user
shared_ptr<ModelOutput> 
VarianceSwapNumOutput::GetModelOutput()
{
  // Construct the output class
  shared_ptr<ModelOutput> pOutput(new ModelOutput);

  BackwardNumOutput::GetOutput(pOutput.get());

  return pOutput;
}
