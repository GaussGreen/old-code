/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/calibratorgeneral.cpp
// Purpose:     implementation of general calibration for all models
// Created:     2005/07/04
// RCS-ID:      $Id: calibratorgeneral.cpp,v 1.17 2006/08/15 17:10:51 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "ito33/sharedptr.h"
#include "ito33/error.h"

#include "ito33/finance/error.h"
#include "ito33/finance/derivative.h"
#include "ito33/finance/derivatives.h"
#include "ito33/finance/option.h"
#include "ito33/finance/forwardoption.h"
#include "ito33/finance/basket_visitor.h"

#include "ito33/pricing/translator.h"
#include "ito33/pricing/calibratorgeneral.h"

extern const ito33::finance::Error ITO33_CALIBRATION_CANCELLED;

using namespace ito33::finance;

namespace ito33
{

namespace pricing
{


void CalibratorGeneral::InitializeArrays(const Translator& translator)
{
  // Save the number of unknowns.  Needed by the objective function
  m_nNbUnknowns = translator.GetNbParameters();

  // Start setting up the data needed by the optimization routine.
  m_pdX.resize(m_nNbUnknowns);
  m_pdFinalX.resize(m_nNbUnknowns);

  m_pdLowerBounds = translator.GetLowerBounds();
  m_pdUpperBounds = translator.GetUpperBounds();

  // Initial guess for the calibration.
  // The order of this constuction must be 'undone' by the objective function.
  // This is made easier by the translator class
  translator.GetParameters(&m_pdX[0]);
}


void 
CalibratorGeneral::ConstructDerivativeLists(const Derivatives& derivatives,
                                            bool bUseForward)
{
  BasketVisitor visitor(derivatives);

  visitor.EnableForwardForOption(bUseForward);

  visitor.Run();

  m_pForwardOption = visitor.GetForwardOption();

  m_pDerivatives = visitor.GetGenericDerivatives();

} // ConstructDerivativeLists

void CalibratorGeneral::ThrowCalibrationCancelledException()
{
  throw EXCEPTION(ITO33_CALIBRATION_CANCELLED);
}

} // namespace pricing

} // namespace ito33
