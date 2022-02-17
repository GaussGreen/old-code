/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/option/optionnumoutput.cpp
// Purpose:     implementation of OptionNumOutput class 
// Author:      Based on ICARE version
// RCS-ID:      $Id: optionnumoutput.cpp,v 1.46 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"

#include "ito33/numeric/interpolation.h"
#include "ito33/numeric/deltagamma.h"
#include "ito33/numeric/domain_fixedspacemesh.h"
#include "ito33/numeric/surfacegeneral.h"
#include "ito33/numeric/surfaceflag.h"

#include "ito33/finance/exercisetype.h"

#include "ito33/pricing/optionparams.h"
#include "ito33/pricing/optionmeshmanager.h"

#include "ito33/ihg/modeloutput.h"

#include "ihg/backwardnumoutput.h"
#include "ihg/optionnumoutput.h"
#include "ihg/optioninstdata.h"

using ito33::shared_ptr;
using ito33::AutoPtr;

using namespace ito33::numeric;

using namespace ito33::ihg;

// implement the AutoPtrDeleter for OptionNumOutput
namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(ihg::OptionNumOutput);
}

void OptionNumOutput::Init(OptionInstData& instdata)
{  
  // Check if instdata has constraint flags
  bool bHasConstraintFlags = false;
  if (instdata.m_pConstraints)
    bHasConstraintFlags = true;

  BackwardNumOutput::InitWithFixMesh
                     (instdata, m_params, bHasConstraintFlags);
}

void OptionNumOutput::SaveSurface(BackwardInstData& instdata, double dTime)
{
  ASSERT_MSG
  (
    dynamic_cast<DomainFixedSpaceMesh*>(m_pDomain.get()), 
    "The type of Domain in OptionNumOutput should be DomainFixedSpaceMesh."
  );

  DomainFixedSpaceMesh& 
    domainFixed = static_cast<DomainFixedSpaceMesh&>(*m_pDomain);

  BackwardNumOutput::SaveSurfaceWithFixMesh(instdata, domainFixed, dTime);
}


// Return the requested data to the user
shared_ptr<ito33::ihg::ModelOutput> OptionNumOutput::GetModelOutput()
{
  // Construct the output class
  shared_ptr<ihg::ModelOutput> pOutput (new ihg::ModelOutput);

  BackwardNumOutput::GetOutput(pOutput.get());

  return pOutput;
}
