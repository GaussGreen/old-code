/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/modeloutput.cpp
// Purpose:     implementation of non inline methods of finance::ModelOutput
// Created:     April 28, 2003
// RCS-ID:      $Id: modeloutput.cpp,v 1.19 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/error.h"
#include "ito33/finance/modeloutput.h"

#include "ito33/numeric/surfacedouble.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/modeloutput.h"

extern const ito33::finance::Error ITO33_OUTPUT_NOT_AVAILABLE;

namespace ito33
{
  
namespace finance
{

/* static */
void ModelOutput::ThrowValueNotAvailable()
{
  throw EXCEPTION_MSG
        (
          ITO33_OUTPUT_NOT_AVAILABLE,
          TRANS("Required scalar value not available!")
        );
}
            
/* static */
void ModelOutput::ThrowDataAtAnalysisDateNotAvailable()
{
  throw EXCEPTION_MSG
        (
          ITO33_OUTPUT_NOT_AVAILABLE,
          TRANS("Required value at analysis date not available!")
        );
}

/* static */
void ModelOutput::ThrowSurfaceNotAvailable()
{
  throw EXCEPTION_MSG
        (
          ITO33_OUTPUT_NOT_AVAILABLE,
          TRANS("Required surface not available!")
        );
}

void ModelOutput::SetRhoResults( const shared_ptr<ModelOutput>& pModelOutputNew, 
                                 double dInverseYCShift)
{  
  double 
    dRho = ( pModelOutputNew->GetPrice() - m_dPrice ) * dInverseYCShift;

  SetRho(dRho);

  if ( m_bHasSurface )
  {
    shared_ptr<numeric::SurfaceDouble>
      pNumRhoSurface = m_pPriceSurface->GetImpl()
                                      ->ComputeFiniteDifference
                                        ( 
                                          *( pModelOutputNew->GetPriceSurface()
                                                          ->GetImpl() ),
                                          dInverseYCShift
                                        );
    
    SharedSurface pRhoSurface (new SurfaceDouble(pNumRhoSurface));

    SetRhoSurface(pRhoSurface);
  }

  if ( m_bHasDataAtAnalysisDate )
  {    
    const Values& shiftedValues = pModelOutputNew->GetPricesAtAnalysisDate();

    size_t nNbValues = m_pdPricesAtAnalysisDate.size();
    ASSERT_MSG(nNbValues == shiftedValues.size(),
               "Array size mismatch when computing rho");

    Values pdRhoValues(nNbValues);
    for (size_t nIdx = 0; nIdx < nNbValues; nIdx++)
      pdRhoValues[nIdx] = 
                        (shiftedValues[nIdx] - m_pdPricesAtAnalysisDate[nIdx]) 
                        * dInverseYCShift;

    SetRhosAtAnalysisDate(pdRhoValues);
  }
}

void ModelOutput::SetUnderlyingRhoResults
     ( const shared_ptr<ModelOutput>& pModelOutputNew, double dInverseYCShift )
{  
  double dUnderlyingRho = ( pModelOutputNew->GetPrice() - m_dPrice )
                        * dInverseYCShift;

  SetUnderlyingRho(dUnderlyingRho);

  if ( m_bHasSurface )
  {
    shared_ptr<numeric::SurfaceDouble>
      pNumUnderlyingRhoSurface = m_pPriceSurface->GetImpl()
                                      ->ComputeFiniteDifference
                                        ( 
                                          *( pModelOutputNew->GetPriceSurface()
                                                          ->GetImpl() ),
                                          dInverseYCShift
                                        );
    
    SharedSurface 
      pUnderlyingRhoSurface( new SurfaceDouble(pNumUnderlyingRhoSurface) );

    SetUnderlyingRhoSurface(pUnderlyingRhoSurface);
  }

  if ( m_bHasDataAtAnalysisDate )
  {    
    const finance::Values&
      shiftedValues = pModelOutputNew->GetPricesAtAnalysisDate();

    size_t nNbValues = m_pdPricesAtAnalysisDate.size();
    ASSERT_MSG(nNbValues == shiftedValues.size(),
               "Array size mismatch when computing rho");

    Values pdUnderlyingRhoValues(nNbValues);
    for (size_t nIdx = 0; nIdx < nNbValues; nIdx++)
      pdUnderlyingRhoValues[nIdx] = 
                        (shiftedValues[nIdx] - m_pdPricesAtAnalysisDate[nIdx]) 
                        * dInverseYCShift;

    SetUnderlyingRhosAtAnalysisDate(pdUnderlyingRhoValues);
  }
}

void 
ModelOutput::SetVegaResults( const shared_ptr<ModelOutput>& pModelOutputNew, 
                             double dInverseShiftVol )
{  
  double 
    dVega = ( pModelOutputNew->GetPrice() - m_dPrice ) * dInverseShiftVol;

  SetVega(dVega);

  if ( m_bHasSurface )
  {
    shared_ptr<numeric::SurfaceDouble>
      pNumVegaSurface = m_pPriceSurface->GetImpl()
                                      ->ComputeFiniteDifference
                                        ( 
                                          *( pModelOutputNew->GetPriceSurface()
                                                          ->GetImpl() ),
                                          dInverseShiftVol
                                        );
    
    SharedSurface pVegaSurface (new finance::SurfaceDouble(pNumVegaSurface));

    SetVegaSurface(pVegaSurface);
  }

  if ( m_bHasDataAtAnalysisDate )
  {    
    const finance::Values&
      shiftedValues = pModelOutputNew->GetPricesAtAnalysisDate();

    size_t nNbValues = m_pdPricesAtAnalysisDate.size();
    ASSERT_MSG(nNbValues == shiftedValues.size(),
               "Array size mismatch when computing vega");

    Values pdVegaValues(nNbValues);
    for (size_t nIdx = 0; nIdx < nNbValues; nIdx++)
      pdVegaValues[nIdx] = 
                        (shiftedValues[nIdx] - m_pdPricesAtAnalysisDate[nIdx]) 
                        * dInverseShiftVol;

    SetVegasAtAnalysisDate(pdVegaValues);
  }
}

void 
ModelOutput::SetFXDeltaResults( const shared_ptr<ModelOutput>& pModelOutputNew, 
                                double dInverseShift )
{  
  double 
    dFXDelta = ( pModelOutputNew->GetPrice() - m_dPrice ) * dInverseShift;

  SetFXDelta(dFXDelta);
}

void ModelOutput::Dump(ito33::XML::Tag& tagParent) const
{
  tagParent.Element(XML_TAG_OUTPUT_PRICE)(m_dPrice);
  tagParent.Element(XML_TAG_OUTPUT_DELTA)(m_dDelta);
  tagParent.Element(XML_TAG_OUTPUT_GAMMA)(m_dGamma);

  if ( HasTheta() )
    tagParent.Element(XML_TAG_OUTPUT_THETA)(m_dTheta);

  if ( HasFXDelta() )
    tagParent.Element(XML_TAG_OUTPUT_FXDELTA)(m_dFXDelta);

  if ( HasVega() )
    tagParent.Element(XML_TAG_OUTPUT_VEGA)(m_dVega);
  if ( HasRho() )
    tagParent.Element(XML_TAG_OUTPUT_RHO)(m_dRho);
  if ( HasUnderlyingRho() )
    tagParent.Element(XML_TAG_OUTPUT_UNDERLYING_RHO)(m_dUnderlyingRho);
  if ( HasFugit() )
    tagParent.Element(XML_TAG_OUTPUT_FUGIT)(m_dFugit);
}

} // namespace finance

} // namespace ito33
