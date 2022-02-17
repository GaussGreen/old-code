/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/mandatoryparams.cpp
// Purpose:     Mandatory params class
// Created:     2004/08/20
// RCS-ID:      $Id: mandatoryparams.cpp,v 1.15 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2004 - 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"

#include "ito33/pricing/mandatory.h"
#include "ito33/pricing/mandatoryparams.h"
#include "ito33/pricing/cbconstraints.h"

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::pricing;
using namespace ito33::numeric;

namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(pricing::MandatoryParams);
}

MandatoryParams* MandatoryParams::Clone() const
{
  // Copy the underlying contract. The new params class will manage the
  // memory
  AutoPtr<pricing::Mandatory> clonedMandatory( new Mandatory(m_mandatory) );

  // Construct and setup the cloned cb params
  MandatoryParams* pClonedParams = new MandatoryParams(clonedMandatory);

  pClonedParams->SetNumParams(m_pNumParams);
  pClonedParams->SetMeshParams(m_pMeshParams);
  pClonedParams->SetYieldCurve( GetYieldCurve() );
  pClonedParams->SetYieldCurveForMesh( GetYieldCurveForMesh() );
  pClonedParams->SetForeignCurve( GetForeignCurve() );
  pClonedParams->SetDividends( GetDividends() );
  pClonedParams->SetValuationTime( GetValuationTime() );
  pClonedParams->SetSpotSharePrice( GetSpotSharePrice() );  

  return pClonedParams;
}
