/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/mandatory_payoff_structure.h
// Purpose:     payoff structure of general mandatory
// Author:      ZHANG Yunzhi
// Created:     2005/03/24
// RCS-ID:      $Id: mandatory_payoff_structure.h,v 1.10 2006/03/22 13:10:23 yann Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @brief The declaration of the mandatory payoff structure. 
   @todo The interaction between new share and averaging probably needs 
          more work.
 */
         
#ifndef _ITO33_PRICING_MANDATORY_PAYOFF_STRUCTURE_H_
#define _ITO33_PRICING_MANDATORY_PAYOFF_STRUCTURE_H_

namespace ito33
{

namespace pricing
{

 
/// The declaration of the mandatory payoff structure.
class MandatoryPayoffStructure
{
public:
  
  // default constructor required by vector<MandatoryPayoffStructure>
  MandatoryPayoffStructure() : m_dDownsideConversionRatio(-1),
                               m_dConversionRatio(-1.),
                               m_dStockAverage(-1.)
  {}

  MandatoryPayoffStructure
      (
        double dDownsideConversionRatio,
        double dUpsideBaseConversionRatio,
        double dLowerStrike,
        double dHigherStrike
      )
      : m_dDownsideConversionRatio(dDownsideConversionRatio),
        m_dUpsideBaseConversionRatio(dUpsideBaseConversionRatio),
        m_dLowerStrike(dLowerStrike),
        m_dHigherStrike(dHigherStrike),
        m_dConversionRatio(-1.),
        m_dStockAverage(-1.)
  {
    m_dResetLevel = dLowerStrike * dDownsideConversionRatio;
    m_dUpsideCorrection
      = dHigherStrike * dUpsideBaseConversionRatio - m_dResetLevel;
  }

  MandatoryPayoffStructure( double dCapPrice,
                double dMinRatio, double dMaxRatio)
              : m_dDownsideConversionRatio(dMaxRatio),
                m_dUpsideBaseConversionRatio(dMinRatio),
                m_dLowerStrike(dCapPrice / dMaxRatio),
                m_dHigherStrike(dCapPrice / dMinRatio),
                m_dResetLevel(dCapPrice),
                m_dUpsideCorrection(0),
                m_dConversionRatio(-1.),
                m_dStockAverage(-1.)
  {
  }

  MandatoryPayoffStructure(double dCapPrice, double dRatio)
    : m_dDownsideConversionRatio(dRatio),
      m_dUpsideBaseConversionRatio(0),
      m_dLowerStrike(dCapPrice / dRatio),
      m_dHigherStrike(1.e30),
      m_dResetLevel(dCapPrice),
      m_dUpsideCorrection(0),
      m_dConversionRatio(-1.),
      m_dStockAverage(-1.)
  {
  }

  MandatoryPayoffStructure(double dRatio)
    : m_dDownsideConversionRatio(dRatio),
      m_dUpsideBaseConversionRatio(dRatio),
      m_dLowerStrike(0),
      m_dHigherStrike(0),
      m_dResetLevel(0),
      m_dUpsideCorrection(0),
      m_dConversionRatio(-1.),
      m_dStockAverage(-1.)
  {
  }

  // Compute the conversion ratio based on the current stock value
  double GetConversionRatio(double dS) const
  {
    if ( dS <= m_dLowerStrike )
      return m_dDownsideConversionRatio;
   
    if ( dS >= m_dHigherStrike )
      return m_dUpsideBaseConversionRatio - m_dUpsideCorrection / dS;
 
    return m_dResetLevel / dS;
  }

  // calculates the conversion values on given spots
  void GetValues(const double *pdS, double *pdValues, size_t nNbS, 
                 const double *pdNewSharePrices) const
  {
    GetValues(pdS, pdValues, nNbS, pdNewSharePrices, 1.0);
  }

  // calculates the conversion values on given spots
  void GetValues(const double *pdS, double *pdValues, size_t nNbS, 
                 const double *pdNewSharePrices, double dFXRate) const
  {
    ASSERT_MSG(m_dDownsideConversionRatio > 0, "bad data");

    double dDownsideConversionRatio   = m_dDownsideConversionRatio * dFXRate;
    double dUpsideBaseConversionRatio = m_dUpsideBaseConversionRatio * dFXRate;
    double dResetLevel                = m_dResetLevel * dFXRate;
    double dUpsideCorrection          = m_dUpsideCorrection * dFXRate;
 
    //We are considering the case where the stock
    //average is used instead of the stock.
    if ( m_dStockAverage >= 0 )
    {  
      
      if ( m_dStockAverage <= m_dLowerStrike )    
      {
        for (size_t n = 0; n < nNbS; n++)
          pdValues[n] = pdNewSharePrices[n] * dDownsideConversionRatio;       
      }
      else if ( m_dStockAverage >= m_dHigherStrike)       
      {
        for (size_t n = 0; n < nNbS; n++)       
          pdValues[n] = pdNewSharePrices[n] * (dUpsideBaseConversionRatio        
           - dUpsideCorrection / m_dStockAverage);  
      }
      else       
      {  
        for (size_t n = 0; n < nNbS; n++)       
          pdValues[n] = pdNewSharePrices[n] * dResetLevel / m_dStockAverage;  
      }
     
      return;
    } // StockAveraging
  
    // We are considering the case where the average conversion ratio is used
    if ( m_dConversionRatio >= 0.0 )
    {
      for (size_t n = 0; n < nNbS; n++)
        pdValues[n]  = pdNewSharePrices[n] * m_dConversionRatio * dFXRate;

      return;
    }

    // classic way
    size_t n = 0;
    for (n = 0; n < nNbS && pdS[n] <= m_dLowerStrike; n++)
      pdValues[n] = pdNewSharePrices[n] * dDownsideConversionRatio;

    for (; n < nNbS && pdS[n] < m_dHigherStrike; n++)
      pdValues[n] = pdNewSharePrices[n] * dResetLevel / pdS[n];

    for (; n < nNbS; n++)
      pdValues[n] = pdNewSharePrices[n] * (dUpsideBaseConversionRatio
                  - dUpsideCorrection / pdS[n]);
  }

  /// The minimum conversion ratio
  double GetMinConversionRatio() const
  {
    return m_dResetLevel / m_dHigherStrike;
  }

  /// the higher strike
  double GetHigherStrike() const
  {
    return m_dHigherStrike;
  }
   
  // indicate that the Average of the stock should be used
  // instead of the stock value for conversion  
  void SetStockAverage(double dStockAverage)
  {
    m_dStockAverage = dStockAverage;
  }
   
  // indicate that the average conversion ratio must be used
  // instead of the stock value for conversion 
  void SetConversionRatioAverage(double dConversionRatio)
  {  
    m_dConversionRatio = dConversionRatio; 
  }


protected:

  double
    m_dDownsideConversionRatio,
    m_dUpsideBaseConversionRatio,
    m_dLowerStrike,
    m_dHigherStrike,
    m_dResetLevel,
    m_dUpsideCorrection;

  double m_dStockAverage;
  double m_dConversionRatio;

};


} // namespace pricing

} // namespace ito33

#endif  //  #ifndef _ITO33_PRICING_MANDATORY_PAYOFF_STRUCTURE_H_
