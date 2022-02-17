/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/theoreticalmodel.cpp
// Purpose:     common TheoreticalModel methods implementation
// Author:      Vadim Zeitlin
// Created:     Dec 6, 2003
// RCS-ID:      $Id: theoreticalmodel.cpp,v 1.25 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2003 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/beforestd.h"
#include <map>
#include <utility>
#include "ito33/afterstd.h"

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"
#include "ito33/dateutils.h"
#include "ito33/debugparameters.h"

#include "ito33/finance/qualitycontrol.h"
#include "ito33/finance/theoreticalmodel.h"
#include "ito33/finance/derivative.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/computationalflags.h"

#include "ito33/numeric/numparams.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/common.h"

using namespace ito33;
using namespace ito33::finance;

// ============================================================================
// implementation
// ============================================================================

/// Default ctor
TheoreticalModel::TheoreticalModel() : m_pQualityControl(new QualityControl),
                                       m_pDebug(new DebugParameters),
                                       m_bHasExternalFlags(false)
{
}

TheoreticalModel::~TheoreticalModel()
{
}

void TheoreticalModel::SetExternalFlagsToDefaults()
{
  SetExternalFlags( make_ptr(new ComputationalFlags) );
}

void TheoreticalModel::SetFlagsFrom
     ( const Derivative& derivative ) const
{
  m_pComputFlags = derivative.GetComputationalFlags();

  if ( !m_pComputFlags )
    m_pComputFlags = make_ptr(new ComputationalFlags);
}

std::string TheoreticalModel::GetDebugOutputFile() const
{
  return m_pDebug->GetDebugOutputFile("ito.xml");
}

numeric::NumParams* 
TheoreticalModel::GetNumParams(const Derivative& derivative) const
{
  return new numeric::NumParams
             (
               *m_pQualityControl,
                GetDoubleFrom( derivative.GetMaturityDate() ) 
                - GetDoubleFrom( derivative.GetSessionData()->GetValuationDate() )
             );
}

void 
TheoreticalModel::GetModelParameters(ModelParametersConsumer&) const
{
  FAIL("hg is not supported yet");
}

TheoreticalModel* TheoreticalModel::DeepCopy() const
{
  FAIL("hg is not supported yet");

  return NULL;
}

void TheoreticalModel::EnableDebugOutput(bool bEnable)
{
  m_pDebug->isEnabled = bEnable;
}

bool TheoreticalModel::IsDebugOutputEnabled() const
{
  return m_pDebug->isEnabled;
}

void TheoreticalModel::SetDebugOutputFile(const std::string& filename)
{
  m_pDebug->filename = filename;
  EnableDebugOutput();
}
