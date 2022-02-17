/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/common/theoreticalmodel.cpp
// Purpose:     HG TheoreticalModel implementation
// Created:     2005/01/13
// RCS-ID:      $Id: theoreticalmodel.cpp,v 1.31 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/constants.h"
#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/error.h"
#include "ito33/finance/derivative.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/qualitycontrol.h"
#include "ito33/finance/modeloutput.h"
#include "ito33/finance/yieldcurve.h"

#include "ito33/hg/error.h"
#include "ito33/hg/theoreticalmodel.h"

#include "ito33/hg/version.h"

#include "ito33/beforestd.h"
#include <map>
#include <memory>
#include <utility>
#include <fstream>
#include "ito33/afterstd.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/common.h"
#include "hg/xml/common.h"

#include "hg/numoutput.h"
#include "hg/computesensitivity.h"

extern const ito33::finance::Error ITO33_INVALID_UNDERLYINGPROCESS;

extern const ito33::hg::Error ITO33_HG_COMPUTEFUNCTION, 
                              ITO33_HG_SHARPERATIO;

namespace ito33
{

  using finance::Derivative;

namespace hg
{
  
#if defined(_MSC_VER) && !defined(_CPPRTTI)
  #error "This project must be compiled with RTTI (/GR)!"
#endif


TheoreticalModel::TheoreticalModel
                  (const shared_ptr<UnderlyingProcess>& pUnderlyingProcess)
                  : m_pUnderlyingProcess(pUnderlyingProcess),
                   m_bHasExternalProcessForMesh(false),
                   m_dSharpeRatio(0)
{
  SetDebugOutputFile("hg.xml");
  EnableDebugOutput(false);

  CHECK_COND(pUnderlyingProcess, ITO33_INVALID_UNDERLYINGPROCESS);
}


void TheoreticalModel::SetSharpeRatio(double dSharpeRatio)
{
  CHECK_COND(dSharpeRatio >= 0, ITO33_HG_SHARPERATIO);

  m_dSharpeRatio = dSharpeRatio;
}


// ----------------------------------------------------------------------------
// types
// ----------------------------------------------------------------------------

// this map allows us to find the pricing function for the given model and
// derivative
typedef std::map< TypeInfo, TheoreticalModel::ModelPriceFunction > 
        PriceFunctions;

// ----------------------------------------------------------------------------
// functions
// ----------------------------------------------------------------------------

// return the price map (this is a degenerate singleton object)
static PriceFunctions& GetPriceFunctionsMap()
{
  static PriceFunctions s_priceFunctions;

  return s_priceFunctions;
}

// ============================================================================
// implementation
// ============================================================================
void TheoreticalModel::Dump(XML::Tag& tagParent) const
{
  XML::Tag tagParam(XML_TAG_MODEL, tagParent);

  m_pUnderlyingProcess->Dump(tagParam);

  tagParam.Element(XML_TAG_HG_SR)(m_dSharpeRatio);
}

// ----------------------------------------------------------------------------
// dispatching over derivatives
// ----------------------------------------------------------------------------

/* static */ void
TheoreticalModel::DoRegisterPriceFunction(const TypeInfo& tiDerivative,
                                          ModelPriceFunction pf)
{
  // TODO: check if it already exists in the map!
  GetPriceFunctionsMap().insert(std::make_pair(tiDerivative, pf));
}

TheoreticalModel::ModelPriceFunction 
TheoreticalModel::FindPriceFunction(const finance::Derivative& derivative) const
{
  PriceFunctions::key_type key(TypeInfo(typeid(derivative)));
  PriceFunctions::const_iterator p = GetPriceFunctionsMap().find(key);

  return p == GetPriceFunctionsMap().end() ? 0 : p->second;
}

shared_ptr<finance::ModelOutput>
TheoreticalModel::ComputeAll(const Derivative& derivative, 
                             ModelPriceFunction pf) const
{  
  if ( !m_bHasExternalFlags )
    SetFlagsFrom(derivative);

  // Compute the price and greeks which we can compute by PDE
  shared_ptr<finance::ModelOutput> pModelOutput = (this->*pf)(derivative);

  shared_ptr<finance::ModelOutput> pModelOutputNew;

  // Create a clone of TM to price with for the greeks computation.
  // This avoid the recompute of the greeks already computed by PDE.
  shared_ptr<TheoreticalModel> pTM ( Clone() );
  
  // Explicitly defines the underlying process for mesh
  pTM->m_pUnderlyingProcessForMesh = m_pUnderlyingProcessForMesh;

  // Sets the same quality control
  pTM->m_pQualityControl = m_pQualityControl; 
  
  pTM->SetExternalFlagsToDefaults();

  // respect the internal flags
  pTM->m_pComputFlags->SetSolverType( m_pComputFlags->GetSolverType() );
  
  pTM->m_pComputFlags->SetDiscretizationMethod
                       ( m_pComputFlags->GetDiscretizationMethod() );

  pTM->m_pComputFlags->SetUseAnalyticSolvers
                       ( m_pComputFlags->GetUseAnalyticSolvers() );

  pTM->m_pComputFlags->SetUseSimilarityReductions
                       ( m_pComputFlags->GetUseSimilarityReductions() );

  // no need to respect SensitivityMethod flag since we are not going to compute
  // sensitivity by PDE anyway

  // internal flags
  // If price surface was computed, we have to compute it also for the greeks
  if ( m_pComputFlags->GetComputeSurface() )
    pTM->m_pComputFlags->SetComputeSurface(true);
  
  // If price at analysis date was computed, we have to compute also
  // the greeks at this date.
  if ( pModelOutput->HasPriceAtAnalysisDate() )
    pTM->m_pComputFlags->SetAnalysisDate( m_pComputFlags->GetAnalysisDate() );

  if ( m_pComputFlags->GetComputeRho() )
  { 
    double dYCShift = SHIFT;
    double dInverseYCShift = 1. / dYCShift;

    // Compute the underlying rho ***************************

    //shared_ptr<finance::Issuer> 
    //  pIssuer = derivative.GetSessionData()->GetEquity()->GetIssuer();
    
    // Save the origin YC
    shared_ptr<finance::YieldCurve> 
      pYC = derivative.GetSessionData()->GetYieldCurve();

    // Perturb the YC
    shared_ptr<finance::YieldCurve> pYCNew ( pYC->Perturb(dYCShift) );

    derivative.GetSessionData()->SetYieldCurve(pYCNew);

    // Compute the price with the perturbed YC
    pModelOutputNew = (pTM.get()->*pf)(derivative);

    // Compute the underlying rho results (scalar and surfaces)
    pModelOutput->SetUnderlyingRhoResults(pModelOutputNew, dInverseYCShift);

    // Restore the origin YC
    derivative.GetSessionData()->SetYieldCurve(pYC);

    //Compute the derivative rho *******************************

    if ( !derivative.IsCrossCurrency() )
      pModelOutput->SetRhoResults(pModelOutputNew, dInverseYCShift);
    else
    {
      shared_ptr<finance::Numeraire> pNumeraire = derivative.GetNumeraire();
      
      ASSERT_MSG( pNumeraire, "The numeraire (currency) for "
                              "the derivative must be defined in order to "
                              "compute the derivative rho if "
                              "the instrument is cross-currency.");
      // Save the original YC
      shared_ptr<finance::RateData> 
        pRateData = derivative.GetSessionData()->GetRateData();

      shared_ptr<finance::YieldCurve> 
        pYC = pRateData->GetYieldCurve(pNumeraire);

      // Perturb the YC
      shared_ptr<finance::YieldCurve> pYCNew ( pYC->Perturb(dYCShift) );

      pRateData->SetYieldCurve(pNumeraire, pYCNew);

      // Compute the price with the perturbed YC
      pModelOutputNew = (pTM.get()->*pf)(derivative);

      // Compute the derivative rho results (scalar and surfaces)
      pModelOutput->SetRhoResults(pModelOutputNew, dInverseYCShift);

      // Restore the origin YC
      pRateData->SetYieldCurve(pNumeraire, pYC);
    }
  }

  
  // If sensitivities are needed but were not computed by PDE or adjoint.
  shared_ptr<NumOutput>
    pNumOutput( static_pointer_cast<NumOutput>(pModelOutput->GetNumOutput()) );

  if ( ( m_pComputFlags->AreAllSensitivitiesActivated()
         || m_pComputFlags->HasSensitivityFlags() )
      && ! pNumOutput->HasSensitivities() )
  {
    double dSensShift = SHIFT;

    std::vector<double> pdSensitivityByFD;
    if ( m_pComputFlags->AreAllSensitivitiesActivated() )
    {
      // Compute all sensitivities
      pdSensitivityByFD = ComputeSensitivity(*this, 
                                             derivative, 
                                             dSensShift, 
                                             pModelOutput->GetPrice() );
    }
    else
    {
      // Compute partial sensitivities
      pdSensitivityByFD = ComputeSpecifiedSensitivities(
                                      m_pUnderlyingProcess, 
                                      m_pQualityControl,
                                      derivative, 
                                      dSensShift, 
                                      pModelOutput->GetPrice(),
                                      m_pComputFlags->GetSensitivityFlags() );
    }


    pNumOutput->SetSensitivities( pdSensitivityByFD );  
  }

  return pModelOutput;
}

shared_ptr<finance::ModelOutput>
TheoreticalModel::Compute
(const Derivative& derivative, TheoreticalModel::ModelPriceFunction pf) const
{
  // Validate at first the derivative
  derivative.ValidateAll();

  // Sets the yc for mesh equal to the yc before the compute process.
  derivative.GetSessionData()->SetYieldCurveForMesh( 
    derivative.GetSessionData()->GetYieldCurve() );

  // Set the underlying process for mesh before the compute process.
  // Could also call GetUnderlyingProcessForMesh().
  if ( !m_bHasExternalProcessForMesh )
    m_pUnderlyingProcessForMesh = m_pUnderlyingProcess;

  CHECK_COND(pf, ITO33_HG_COMPUTEFUNCTION);

  if ( IsDebugOutputEnabled() )
  {
    // dump everything we can before calling the pricing function in case it
    // fails
    std::ofstream ofs(GetDebugOutputFile().c_str());

    XML::RootTag tagRoot(XML_TAG_HG_ROOT, ofs);
    tagRoot.precision(10);
    tagRoot.Attr(XML_ATTR_ROOT_VERSION, ITO33_HG_SERIES_DOT_STRING);

    derivative.GetSessionData()->Dump(tagRoot);

    // dump the Model itself
    Dump(tagRoot);

    tagRoot.Element(XML_TAG_DERIVATIVES).Element(derivative);

    // now do call it
    shared_ptr<finance::ModelOutput> pMO( ComputeAll(derivative, pf) );

    // and dump its output
    if ( pMO )
    {
      XML::Tag tagOutput(XML_TAG_OUTPUT, tagRoot);
      pMO->Dump(tagOutput);
    }

    return pMO;
  }
  else // simply return the result
  {
    return ComputeAll(derivative, pf);
  }
}

shared_ptr<finance::ModelOutput>
TheoreticalModel::Compute(const Derivative& derivative) const
{
  ModelPriceFunction pf = FindPriceFunction(derivative);
  
  return Compute(derivative, pf);
}

TheoreticalModel*
TheoreticalModel::Clone() const
{
  TheoreticalModel *pTM = new TheoreticalModel(m_pUnderlyingProcess);

  return pTM;
}


} // namespace hg

} // namespace ito33
