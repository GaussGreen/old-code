/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/onetouch/onetouchnumoutput.cpp
// Purpose:     implementation of  numeric output class for OneTouch
// Created:     2006/08/11
// RCS-ID:      $Id: onetouchnumoutput.cpp,v 1.2 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"

#include "ito33/numeric/domain_fixedspacemesh.h"
#include "ito33/numeric/surfacegeneral.h"

#include "ito33/pricing/onetouchparams.h"
#include "ito33/pricing/onetouchmeshmanager.h"

#include "ihg/onetouchnumoutput.h"
#include "ihg/onetouchinstdata.h"

#include "ito33/ihg/modeloutput.h"

// implement the AutoPtrDeleter for OneTouchNumOutput
namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(ihg::OneTouchNumOutput);
}

namespace ito33
{

namespace ihg
{
  using namespace numeric;

void OneTouchNumOutput::Init(OneTouchInstData& instdata)
{
  BackwardNumOutput::InitWithFixMesh(instdata, m_params, false);
}

void OneTouchNumOutput::SaveSurface(BackwardInstData& instdata, double dTime)
{
  ASSERT_MSG
  (
    dynamic_cast<DomainFixedSpaceMesh*>(m_pDomain.get()), 
    "The type of Domain in OneTouchNumOutput should be DomainFixedSpaceMesh."
  );

  DomainFixedSpaceMesh& 
    domainFixed = static_cast<DomainFixedSpaceMesh&>(*m_pDomain);

  BackwardNumOutput::SaveSurfaceWithFixMesh(instdata, domainFixed, dTime);
}


// Return the requested data to the user
shared_ptr<ihg::ModelOutput> OneTouchNumOutput::GetOutput()
{
  // Construct the output class
  shared_ptr<ihg::ModelOutput> pOutput( new ihg::ModelOutput );

  BackwardNumOutput::GetOutput( pOutput.get() );

  return pOutput;
}

} // namespace ihg

} // namespace ito33
