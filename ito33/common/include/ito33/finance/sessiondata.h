/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/sessiondata.h
// Purpose:     Class containing current market information for pricing
// Created:     2006/03/16
// RCS-ID:      $Id: sessiondata.h,v 1.5 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/sessiondata.h
    @brief SessionData class.

    Class containing current market information required for pricing.  
    Current information includes things such as interest rates,  
    foreign exchange rates and the spot price.
 */

#ifndef _ITO33_FINANCE_SESSIONDATA_H_
#define _ITO33_FINANCE_SESSIONDATA_H_

#include "ito33/sharedptr.h"
#include "ito33/dlldecl.h"
#include "ito33/date.h"

#include "ito33/finance/ratedata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/spotfxrates.h"

namespace ito33
{

namespace XML
{
  class Tag;
}

namespace finance
{

  class ITO33_DLLDECL Dividends;
  class ITO33_DLLDECL YieldCurve;
  class ITO33_DLLDECL Numeraire;
  

/**
    The SessionData class.
    
    Class containing current market information required for pricing.  
    Current information includes things such as interest rates,  
    foreign exchange rates and the spot price.
 */
class ITO33_DLLDECL SessionData
{
public:

  /**
      Ctor accepting rate data, an equity, and the valuation date.

      @param pRateData valid (ie non @c NULL) rate data to use
      @param pEquity a valid (ie non @c NULL) equity to use
      @param valuationDate valid valuation date
   */
  SessionData(const shared_ptr<RateData>& pRateData,
              const shared_ptr<Equity>& pEquity,
              Date valuationDate);

  /**
      Dtor is not virtual: this class should not be used polymorphically.
   */
  ~SessionData() {}

  /**
      @name Modifiers.
   */
  //@{

  /**
      The valuation date to be used for price computations.

      This date is used for pricing, calibration and hedging.

      @param dValuationDate The valuation date used for price computations
   */
  void SetValuationDate(Date dValuationDate)
  {
    m_dValuationDate = dValuationDate;
  }

  /**
      The data associated with rates (exchange rates, yield curves).

      In general, each rate is associated with a currency.

      @param pRateData The rate data
   */
  void SetRateData(const shared_ptr<RateData>& pRateData)
  {
    m_pRateData = pRateData;
  }

  /**
      The equity to use for price computations.

      @param pEquity The equity to use for price computations
   */
  void SetEquity(const shared_ptr<Equity>& pEquity)
  {
    m_pEquity = pEquity;
  }

  //@}


  /**
      @name Accessors.
   */
  //@{

  /**
      The valuation date to be used for price computations.
      
      @return valuation date for price computations
   */
  Date GetValuationDate() const 
  { 
    return m_dValuationDate; 
  }

  /**
      The data associated with rates (exchange rates, yield curves).

      In general, each rate is associated with a currency.

      @return rate data
   */
  const shared_ptr<RateData>& GetRateData() const 
  { 
    return m_pRateData; 
  }

  /**
      The equity to use for price computations.      

      @return equity used for pricing
   */
  const shared_ptr<Equity>& GetEquity() const 
  { 
    return m_pEquity; 
  }

  //@}

  /// @name noexport helper functions for numerical use.
  //@{

  /**
      @internal
      @brief The spot of the underlying equity.

      @return the current spot of the underlying equity

      @noexport
   */
  double GetSpotSharePrice() const;

  
  /**
      @internal
      @brief The previous share price of the underlying equity.

      @return the previous share price of the underlying equity

      @noexport
   */
  double GetPreviousSharePrice() const;


  /**
      @internal
      @brief The collection of the discrete dividends of the underlying.

      @return the discrete dividends of the underlying equity

      @noexport
   */
  const shared_ptr<Dividends>& GetDividends() const;

  /**
      @internal
      @brief The currency of the underlying equity.

      @return the currency of the underlying equity

      @noexport
   */
  const shared_ptr<Numeraire>& GetNumeraire() const;

  /**
      @internal
      @brief The yield curve for the currency of the underlying equity.

      @return the yield curve for the currency of the underlying equity

      @noexport
   */
  const shared_ptr<YieldCurve>& GetYieldCurve() const;

  /**
      @internal
      @brief The yield curve for the currency of the underlying equity.

      @param pYC the yield curve for the currency of the underlying equity

      @noexport
   */
  void SetYieldCurve(const shared_ptr<YieldCurve>& pYC) const;

  /**
      @internal
      @brief The yield curve for mesh construction for the underlying equity.

      @return the yield curve for mesh construction for the underlying equity

      @noexport
   */
  const shared_ptr<YieldCurve>& GetYieldCurveForMesh() const;

  /**
      @internal
      @brief The yield curve for mesh construction for the underlying equity.

      @param pYC The yield curve for mesh construction for the underlying

      @noexport
   */
  void SetYieldCurveForMesh(const shared_ptr<YieldCurve>& pYC) const;

  /**
      @internal
      @brief The continuous dividend curve of the underlying equity.
      
      Sometimes referred to as the borrow curve or foreign curve, although
      this function has nothing to do with foreign exchanges.
    
      @return the continuous dividend curve of the underlying equity

      @noexport
   */
  const shared_ptr<YieldCurve>& GetForeignCurve() const;

  /**
      @internal
      @brief The spot FX rate between two currencies.
      
      This is the number of units of base currency needed to purchase 
      one unit of the foreign currency (indirect quotation). 

      @param pNumeraire1 the foreign currency
      @param pNumeraire2 the base currency
      @return the spot FX rate between two currencies.

      @noexport
   */
  double GetSpotFXRate(const shared_ptr<Numeraire>& pNumeraire1, 
                       const shared_ptr<Numeraire>& pNumeraire2) const;

  //@}   

  /**
      @internal
      @brief Dumps all data stored in this object in XML format.

      This method is usually called by the function doing the pricing,
      calibration &c but can also be called "manually" if needed.

      @param tagParent the parent tag under which our tag(s) should be created
      
      @noexport
   */
  void Dump(XML::Tag& tagParent) const;


private:

  /// the rate data
  shared_ptr<RateData> m_pRateData;

  /// the underlying equity of the instruments using this session data
  shared_ptr<Equity> m_pEquity;

  /// the pricing date 
  Date m_dValuationDate;

};

// ----------------------------------------------------------------------------
// inline functions
// ----------------------------------------------------------------------------

inline double SessionData::GetSpotSharePrice() const
{
  return GetEquity()->GetSpotSharePrice();
}

inline double SessionData::GetPreviousSharePrice() const
{
  return GetEquity()->GetPreviousSharePrice();
}

inline const shared_ptr<Dividends>& SessionData::GetDividends()  const
{
  return GetEquity()->GetDividends();
}

inline const shared_ptr<Numeraire>& SessionData::GetNumeraire() const
{
  return GetEquity()->GetNumeraire();
}

inline const shared_ptr<YieldCurve>& SessionData::GetYieldCurve() const
{
  return m_pRateData->GetYieldCurve( GetEquity()->GetNumeraire() );
}

inline void 
SessionData::SetYieldCurve(const shared_ptr<YieldCurve>& pYC) const
{
  m_pRateData->SetYieldCurve( GetEquity()->GetNumeraire(), pYC );
}

inline const shared_ptr<YieldCurve>& SessionData::GetYieldCurveForMesh() const
{
  return m_pRateData->GetYieldCurveForMesh( GetEquity()->GetNumeraire() );
}

inline void 
SessionData::SetYieldCurveForMesh(const shared_ptr<YieldCurve>& pYC) const
{
  m_pRateData->SetYieldCurveForMesh( GetEquity()->GetNumeraire(), pYC );
}

inline const shared_ptr<YieldCurve>& SessionData::GetForeignCurve() const
{
  return GetEquity()->GetBorrowCurve();
}

inline double 
SessionData::GetSpotFXRate(const shared_ptr<Numeraire>& pNumeraire1, 
                           const shared_ptr<Numeraire>& pNumeraire2) const
{
  return m_pRateData->GetSpotFXRates()->GetFXRate(pNumeraire1, pNumeraire2);
}


} // namespace finance

} // namespace ito33

#endif // _ITO33_FINANCE_SESSIONDATA_H_
