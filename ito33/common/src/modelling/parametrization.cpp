/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/parametrization.cpp
// Purpose:     Implementation of shared ptr for Parametrization class
// Created:     2005/07/29
// RCS-ID:      $Id: parametrization.cpp,v 1.9 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/debugparameters.h"

#include "ito33/finance/parametrization.h"
#include "ito33/finance/calibrationprocess.h"

namespace ito33
{

namespace finance
{

Parametrization::Parametrization() 
  : m_pProcess(new CalibrationProcess),
    m_pDebug(new DebugParameters)
{
}

Parametrization::~Parametrization() 
{
}

void Parametrization::SetDebugOutputFile(const std::string& filename)
{
  m_pDebug->filename = filename;
  EnableDebugOutput();
}

bool Parametrization::IsDebugOutputEnabled() const 
{
  return m_pDebug->isEnabled;
}

void Parametrization::EnableDebugOutput(bool bEnable) 
{
  m_pDebug->isEnabled = bEnable;
}

void Parametrization::SetCalibrationControl(CalibrationControl* pControl)
{
  m_pProcess->SetControl(pControl);
}

} // namespace finance

} // namespace ito33
