/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/common/computesensitivity.cpp
// Purpose:     implement sensitivity computation by finite difference
// Created:     2005/05/17
// RCS-ID:      $Id: computesensitivity.cpp,v 1.10 2006/08/22 18:08:31 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/computationalflags.h"
#include "ito33/finance/qualitycontrol.h"

#include "ito33/hg/theoreticalmodel.h"
#include "ito33/hg/underlyingprocess.h"

#include "hg/translator.h"
#include "hg/computesensitivity.h"


namespace ito33
{

namespace hg
{

std::vector<double> 
ComputeSensitivity(const TheoreticalModel& model, 
                   const finance::Derivative& derivative,
                   double dShift,
                   const boost::optional<double>& unshiftedPrice)
{
  // Compute all sensitivities
  shared_ptr<UnderlyingProcess> pProcess = model.GetUnderlyingProcess();
  size_t nNbParams = pProcess->GetNbParameters();
  std::vector<bool> pbSensitivityFlags( nNbParams, true );

  return ComputeSpecifiedSensitivities(pProcess, model.GetQualityControl(),
                                       derivative, dShift, 
                                       unshiftedPrice, pbSensitivityFlags);
}


std::vector<double> 
ComputeSpecifiedSensitivities(
                   const shared_ptr<UnderlyingProcess>& pUnderlyingProcess,
                   const shared_ptr<finance::QualityControl>& pQualityControl,
                   const finance::Derivative& derivative,
                   double dShift, 
                   const boost::optional<double>& unshiftedPrice, 
                   const std::vector<bool>& pbSensitivityFlags)
{
  // Create a model for pricing
  shared_ptr<TheoreticalModel> 
    myModel( new TheoreticalModel(pUnderlyingProcess) );
  
  myModel->SetExternalFlagsToDefaults();

  // respect the internal flags
  shared_ptr<finance::ComputationalFlags>
    pFlags( derivative.GetComputationalFlags() );

  if ( pFlags )
  {
    myModel->GetComputationalFlags()->SetSolverType( pFlags->GetSolverType() );
  
    myModel->GetComputationalFlags()->SetDiscretizationMethod
       ( pFlags->GetDiscretizationMethod() );

    // no need to respect SensitivityMethod flag since we are not going to
    // compute sensitivity by PDE anyway
  }

  // Set the quality control (assume same as used for unshifted price)
  myModel->SetQualityControl( pQualityControl );

  // Use the same mesh for the shifted and unshifted prices
  myModel->SetUnderlyingProcessForMesh( pUnderlyingProcess );

  // Compute unshifted price if not already set
  double dPrice = unshiftedPrice 
                ? *unshiftedPrice : myModel->Compute(derivative)->GetPrice();

  // Shift the parameters, using translator to make things easy
  Translator translator( *( myModel->GetUnderlyingProcess() ) );
  translator.SetFlags(pbSensitivityFlags);
  std::vector<double> parameters( translator.GetParameters() );
  
  std::vector<double> sensitivities( parameters.size() );
  for (size_t nIdx = 0; nIdx < parameters.size(); nIdx++)
  {
    // Setup the perturbed parameters
    std::vector<double> newParameters(parameters);
    newParameters[nIdx] += dShift;

    // Get the perturbed underlying process
    shared_ptr<UnderlyingProcess> pUP = translator(newParameters);

    myModel->SetUnderlyingProcess(pUP);

    // Do the finite difference
    double dPerturbedPrice = myModel->Compute(derivative)->GetPrice();

    sensitivities[nIdx] = (dPerturbedPrice - dPrice) / dShift;  
  }

  return sensitivities;
}


} // namespace hg

} // namespace ito33
