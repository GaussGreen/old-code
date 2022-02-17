/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/theoreticalmodel.cpp
// Purpose:     IHG TheoreticalModel implementation
// Created:     2004/03/26
// RCS-ID:      $Id: theoreticalmodel.cpp,v 1.67 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/constants.h"
#include "ito33/sharedptr.h"
#include "ito33/autoptr.h"
#include "ito33/useexception.h"
#include "ito33/debugparameters.h"

#include "ito33/finance/error.h"
#include "ito33/finance/derivative.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/qualitycontrol.h"
#include "ito33/finance/modelparametersconsumer.h"
#include "ito33/finance/yieldcurve.h"
#include "ito33/finance/modeloutput.h"

#include "ito33/numeric/exception.h"

#include "ito33/ihg/error.h"
#include "ito33/ihg/volatility.h"
#include "ito33/ihg/hazardrate.h"
#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/version.h"

#include "ihg/impliedvol.h"
#include "ihg/numoutput.h"

#include "ito33/beforestd.h"
#include <map>
#include <memory>
#include <utility>
#include <fstream>
#include "ito33/afterstd.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/underlyingprocess.h"

#include "ihg/xml/common.h"

extern const ito33::finance::Error ITO33_INVALID_UNDERLYINGPROCESS;

extern const ito33::ihg::Error ITO33_IHG_CANT_COMPUTE,
                               ITO33_IHG_IMPLIEDVOL;

namespace ito33
{

  using finance::Derivative;

namespace ihg
{
  
#if defined(_MSC_VER) && !defined(_CPPRTTI)
  #error "This project must be compiled with RTTI (/GR)!"
#endif

bool g_bComputeWithSameMesh = false;

// ----------------------------------------------------------------------------
// types
// ----------------------------------------------------------------------------

// this map allows us to find the pricing function for the given model and
// derivative
typedef std::map< TypeInfo, TheoreticalModel::PriceFunction > 
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

TheoreticalModel::TheoreticalModel(): finance::TheoreticalModel()
{
  m_pUnderlyingProcess = shared_ptr<UnderlyingProcess>(new UnderlyingProcess());

  SetDebugOutputFile("ihg.xml");
  EnableDebugOutput(false);
}

TheoreticalModel::TheoreticalModel(const shared_ptr<UnderlyingProcess>& 
                                   pUnderlyingProcess)
                 :finance::TheoreticalModel(), 
                  m_pUnderlyingProcess(pUnderlyingProcess)
{ 
  CheckAll();

  SetDebugOutputFile("ihg.xml");
  EnableDebugOutput(false);
}
  
TheoreticalModel::TheoreticalModel(const shared_ptr<Volatility>& pVolatility, 
                    const shared_ptr<HazardRate>& pHazardRate)
                 :finance::TheoreticalModel()
{
  m_pUnderlyingProcess = shared_ptr<UnderlyingProcess>
                          (  
                            new UnderlyingProcess(pVolatility, pHazardRate) 
                          );

  CheckAll();

  SetDebugOutputFile("ihg.xml");
  EnableDebugOutput(false);
}

void TheoreticalModel::Dump(XML::Tag& tagParent) const
{
  XML::Tag tagParam(XML_TAG_IHG_MODEL, tagParent);

  // dump the underlying process
  GetUnderlyingProcess()->Dump(tagParam);
}

// ----------------------------------------------------------------------------
// dispatching over derivatives
// ----------------------------------------------------------------------------

/* static */ void
TheoreticalModel::DoRegisterPriceFunction(const TypeInfo& tiDerivative,
                                          PriceFunction pf)
{
  // TODO: check if it already exists in the map!
  GetPriceFunctionsMap().insert
    (std::make_pair(tiDerivative, pf));
}

TheoreticalModel::PriceFunction TheoreticalModel::FindPriceFunction
            (const finance::Derivative& derivative) const
{
  PriceFunctions::key_type key(TypeInfo(typeid(derivative)));
  PriceFunctions::const_iterator p = GetPriceFunctionsMap().find(key);

  return p == GetPriceFunctionsMap().end() ? 0 : p->second;
}

void TheoreticalModel::CheckAll() const
{
  CheckUnderlyingProcess();

  m_pUnderlyingProcess->CheckAll();
}

void TheoreticalModel::CheckUnderlyingProcess() const
{
  CHECK_COND(m_pUnderlyingProcess, ITO33_INVALID_UNDERLYINGPROCESS);
}

void TheoreticalModel::SetUnderlyingProcess(
                        const shared_ptr<UnderlyingProcess>& pUnderlyingProcess)
{    
  m_pUnderlyingProcess = pUnderlyingProcess;

  CheckUnderlyingProcess();
}

void TheoreticalModel::SetVolatility(const shared_ptr<Volatility>& pVolatility)
{
  m_pUnderlyingProcess->SetVolatility(pVolatility);
}

void TheoreticalModel::SetHazardRate(const shared_ptr<HazardRate>& pHazardRate)
{
  m_pUnderlyingProcess->SetHazardRate(pHazardRate);
}

const shared_ptr<UnderlyingProcess>& 
TheoreticalModel::GetUnderlyingProcess() const
{
  return m_pUnderlyingProcess;
}

const shared_ptr<UnderlyingProcess>& 
TheoreticalModel::GetUnderlyingProcessForMesh() const
{
  ASSERT_MSG( m_pUnderlyingProcessForMesh, 
              "Invalid underlying process for mesh in TheoreticalModel "
              "class" );
  
  m_pUnderlyingProcessForMesh->CheckAll();

  return m_pUnderlyingProcessForMesh;
}

const shared_ptr<Volatility>& TheoreticalModel::GetVolatility() const
{
  CheckUnderlyingProcess();

  return m_pUnderlyingProcess->GetVolatility();
}

const shared_ptr<HazardRate>& TheoreticalModel::GetHazardRate() const
{
  CheckUnderlyingProcess();

  return m_pUnderlyingProcess->GetHazardRate();
}

shared_ptr<finance::ModelOutput>
TheoreticalModel::ComputeAll(const Derivative& derivative, 
                             PriceFunction pf) const
{  
  if ( !m_bHasExternalFlags )
    SetFlagsFrom(derivative);

  // Compute the price and greeks which we can compute by PDE
  shared_ptr<finance::ModelOutput> pModelOutput( (this->*pf)(derivative) );

  shared_ptr<finance::ModelOutput> pModelOutputNew;

  // Create a clone of TM to price with for the greeks computation.
  // This avoid the recompute of the greeks already computed by PDE.
  shared_ptr<TheoreticalModel> pTM ( Clone() );
  pTM->SetExternalFlagsToDefaults();

  // respect the internal flags
  pTM->m_pComputFlags->SetSolverType( m_pComputFlags->GetSolverType() );
  
  pTM->m_pComputFlags->SetDiscretizationMethod
                       ( m_pComputFlags->GetDiscretizationMethod() );

  // no need to respect SensitivityMethod flag since we are not going to compute
  // sensitivity by PDE anyway

  // Sets the same quality control
  pTM->m_pQualityControl = m_pQualityControl;

  // use the same underlying process for mesh
  pTM->m_pUnderlyingProcessForMesh = m_pUnderlyingProcessForMesh;

  // If price surface was computed, we have to compute it also for the greeks
  if ( m_pComputFlags->GetComputeSurface() )
    pTM->m_pComputFlags->SetComputeSurface(true);
  
  // If price at analysis date was computed, we have to compute also
  // the greeks at this date.
  if ( pModelOutput->HasPriceAtAnalysisDate() )
    pTM->m_pComputFlags->SetAnalysisDate( m_pComputFlags->GetAnalysisDate() );

  // FXDelta computation
  if ( derivative.IsCrossCurrency() )
  {
    double dFXRateShift = SHIFT;
    double dInverseFXRateShift = 1. / dFXRateShift;

    // Perturb the FX rate
    derivative.PerturbFXRate( dFXRateShift );

    // Compute the price with the perturbed FX rate
    pModelOutputNew = (pTM.get()->*pf)(derivative);

    pModelOutput->SetFXDeltaResults(pModelOutputNew, dInverseFXRateShift);

    // Restore the FX rate
    derivative.PerturbFXRate( -dFXRateShift );
  }

  // Rho computations
  if ( m_pComputFlags->GetComputeRho() )
  { 
    double dYCShift = SHIFT;
    double dInverseYCShift = 1. / dYCShift;

    // Compute the underlying rho ***************************
  
    // Save the original YC
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
      shared_ptr<finance::YieldCurve> pYCNew( pYC->Perturb(dYCShift) );

      pRateData->SetYieldCurve(pNumeraire, pYCNew);

      // Compute the price with the perturbed YC
      pModelOutputNew = (pTM.get()->*pf)(derivative);

      // Compute the derivative rho results (scalar and surfaces)
      pModelOutput->SetRhoResults(pModelOutputNew, dInverseYCShift);

      // Restore the origin YC
      pRateData->SetYieldCurve(pNumeraire, pYC);
    }
  }

  //If vega needed and not computed by PDE.
  if ( m_pComputFlags->GetComputeVega() && !pModelOutput->HasVega() )
  {
    double dShiftVol = SHIFT;
    double dInverseShiftVol = 1. / dShiftVol;
    
    // Save the origin volatility
    shared_ptr<Volatility> 
      pVolatility = pTM->GetVolatility();

    // Perturb the volatility
    pTM->SetVolatility( pTM->GetVolatility()->Perturb(dShiftVol) );

    // Compute the price with the perturbed volatility
    pModelOutputNew = (pTM.get()->*pf)(derivative);

    // Compute the vega results (scalar and surfaces)
    pModelOutput->SetVegaResults(pModelOutputNew, dInverseShiftVol);

    // Restore the origin volatility
    pTM->SetVolatility( pVolatility );
  }

  return pModelOutput;
}

shared_ptr<finance::ModelOutput>
TheoreticalModel::Compute(const Derivative& derivative) const
{
  // Validate at first the derivative
  derivative.ValidateAll();

  if ( !g_bComputeWithSameMesh )
  {
    // Sets the yc for mesh to the yc before the compute process.    
    shared_ptr<finance::SessionData> pSessionData = derivative.GetSessionData();
    pSessionData->SetYieldCurveForMesh( pSessionData->GetYieldCurve() );

    // Sets the volatility for mesh at the volatility before the compute process.
    m_pUnderlyingProcessForMesh = GetUnderlyingProcess();
  }

  ASSERT_MSG( m_pUnderlyingProcessForMesh, "Underlying process for mesh not "
              "defined.");

  // Compute

  PriceFunction pf = FindPriceFunction(derivative);

  CHECK_COND(pf, ITO33_IHG_CANT_COMPUTE);

  if ( IsDebugOutputEnabled() )
  {
    // dump everything we can before calling the pricing function in case it
    // fails
    std::ofstream ofs(GetDebugOutputFile().c_str());
    XML::RootTag tagRoot(XML_TAG_IHG_ROOT, ofs);
    tagRoot.precision(10);
    tagRoot.Attr(XML_ATTR_IHG_ROOT_VERSION, ITO33_IHG_VERSION_DOT_STRING);

    // Dump the session data
    derivative.GetSessionData()->Dump(tagRoot);

    // dump the Model itself
    Dump(tagRoot);

    tagRoot.Element(XML_TAG_IHG_DERIVATIVES).Element(derivative);

    // now do call it
    shared_ptr<finance::ModelOutput> pMO ( ComputeAll(derivative, pf) );

    // and dump its output
    if ( pMO )
    {
      XML::Tag tagOutput(XML_TAG_IHG_OUTPUT, tagRoot);
      pMO->Dump(tagOutput);
    }

    return pMO;
  }
  else // simply return the result
  {
    return ComputeAll(derivative, pf);
  }
}

double
TheoreticalModel::ComputeImpliedBrownianVol(const Derivative& derivative) const
{
  if ( IsDebugOutputEnabled() )
  {
    // dump everything we can before calling the pricing function in case it
    // fails
    std::ofstream ofs(GetDebugOutputFile().c_str());
    XML::RootTag tagRoot(XML_TAG_IHG_ROOT, ofs);
    tagRoot.precision(10);
    tagRoot.Attr(XML_ATTR_IHG_ROOT_VERSION, ITO33_IHG_VERSION_DOT_STRING);

    {
      derivative.GetSessionData()->Dump(tagRoot);
    }

    // dump only the hazard rate as the
    // volatility is not defined 
    {
      XML::Tag tagParam(XML_TAG_UNDERLYING_PROCESS, tagRoot);  
      GetUnderlyingProcess()->GetHazardRate()->Dump(tagParam);
    }
    
    tagRoot.Element(XML_TAG_IHG_DERIVATIVES).Element(derivative);
 
  }

  ImpliedVol impliedVol(derivative, m_pUnderlyingProcess->GetHazardRate(), 
                        m_pQualityControl);
 
  double dImpliedVol;

  try
  {
    dImpliedVol = impliedVol.Compute();
  }
  catch(const numeric::Exception&)
  {
    throw EXCEPTION(ITO33_IHG_IMPLIEDVOL);
  }

  return dImpliedVol;
}

TheoreticalModel* 
TheoreticalModel::Clone() const
{
  TheoreticalModel* pTmp = new TheoreticalModel(m_pUnderlyingProcess);

  return pTmp;
}

void 
TheoreticalModel::GetModelParameters(
                    finance::ModelParametersConsumer& visitor) const
{
  GetVolatility()->GetModelParameters( visitor );
  GetHazardRate()->GetModelParameters( visitor );
}

} // namespace ihg

} // namespace ito33
