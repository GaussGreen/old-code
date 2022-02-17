/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/copytheoreticalmodel.cpp
// Purpose:     implement copy function of TheoreticalModel object
// Author:      ZHANG Yunzhi
// Created:     2005/09/05
// RCS-ID:      $Id: copytheoreticalmodel.cpp,v 1.8 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/ihg/theoreticalmodel.h"

#include "ito33/ihg/volatility.h"
#include "ito33/ihg/hazardrate.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/volatilitypower.h"
#include "ito33/ihg/volatilitytanh.h"
#include "ito33/ihg/volatilitytimeonly.h"
#include "ito33/ihg/volatilitytimeonly.h"
#include "ito33/ihg/hazardratecombo.h"
#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/hazardratepower.h"
#include "ito33/ihg/hrspotcomponentpower.h"

#include "ito33/ihg/volatility_visitor.h"
#include "ito33/ihg/hazardrate_visitor.h"

namespace ito33
{

namespace ihg
{

namespace
{

// implement visitor of theoretical model to copy a new TM object
class TMCreator : public VolatilityVisitor,
                  public HazardRateVisitor
{
public:

  TheoreticalModel* Copy(const TheoreticalModel* source)
  {
    source->GetVolatility()->Visit(*this);
    source->GetHazardRate()->Visit(*this);
    
    shared_ptr<UnderlyingProcess>
      pUP( new UnderlyingProcess(m_pVolatility, m_pHazardRate) );    
    
    pUP->SetPostDefaultVolatility( 
      source->GetUnderlyingProcess()->GetPostDefaultVolatility() );

    TheoreticalModel 
      *pTM = new TheoreticalModel( pUP ); 

    return pTM;
  }

  virtual void OnVolatilityFlat(const VolatilityFlat& volFlat)
  { 
    m_pVolatility = shared_ptr<Volatility>
                    ( new VolatilityFlat(volFlat.GetValue() ) );
  }

  virtual void OnVolatilityPower(const VolatilityPower& volPower)
  {
    m_pVolatility = shared_ptr<Volatility>
                    ( new VolatilityPower( volPower.GetAlpha(),
                                            volPower.GetBeta(),
                                            volPower.GetS0() ) );
  }
  
  virtual void OnVolatilityTanh(const ihg::VolatilityTanh& volTanh)
  {
    m_pVolatility = shared_ptr<Volatility>
                    ( new VolatilityTanh( volTanh.GetLeft(),
                                          volTanh.GetRight(),
                                          volTanh.GetScale(),
                                          volTanh.GetS0() ) );
  }

  virtual void OnVolatilityTimeOnly
                (const VolatilityTimeOnly& volTimeOnly)
  {
    m_pVolatility 
        = shared_ptr<Volatility>
          ( new VolatilityTimeOnly( volTimeOnly.GetDates(),
                                    volTimeOnly.GetTimeComponentValues()));
  }

  /// Hazard rate combo
  virtual void OnHazardRateCombo(const HazardRateCombo& hrCombo)
  {
    ASSERT ( dynamic_cast<HRSpotComponentPower*>
                              ( hrCombo.GetSpotComponent().get() ) );
    const HRSpotComponentPower*
            pSC = static_cast<HRSpotComponentPower*>
                              ( hrCombo.GetSpotComponent().get() );
    shared_ptr<SpotComponent>
      pSpotComponent( new HRSpotComponentPower(pSC->GetBeta(), pSC->GetS0()) );
    m_pHazardRate = shared_ptr<HazardRate>
                    ( new HazardRateCombo
                          ( pSpotComponent,
                            hrCombo.GetDates(),
                            hrCombo.GetTimeComponentValues() ) );
  }

  /// Time only hazard rate
  virtual void OnHazardRateTimeOnly(const HazardRateTimeOnly &hrTimeOnly)
  {
    m_pHazardRate = shared_ptr<HazardRate>
                    ( new HazardRateTimeOnly
                          ( hrTimeOnly.GetDates(),
                            hrTimeOnly.GetTimeComponentValues() ) );
  }

  //power hazard rate
  virtual void OnHazardRatePower(const HazardRatePower &hrPower)
  {
    m_pHazardRate = shared_ptr<HazardRate>
                    ( new HazardRatePower
                          ( hrPower.GetAlpha(),
                            hrPower.GetBeta(),
                            hrPower.GetS0() ) );
  }
  
private:
  shared_ptr<Volatility> m_pVolatility;
  shared_ptr<HazardRate> m_pHazardRate;
};

} // end of anonymous namespace


TheoreticalModel*
TheoreticalModel::DeepCopy() const
{
  // we use Visit function of Volatility and HR (where comes TMCreator)
  // instead of adding Clone() which leads to interface change.
  // Moreover, Clone() doesn't know how to handle CallBack classes
  TMCreator c;
  return c.Copy(this);
}

} // namespace ihg

} // namespace ito33
