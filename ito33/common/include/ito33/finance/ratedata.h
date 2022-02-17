/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/ratedata.h
// Purpose:     Class containing all rate related information for pricing
// Created:     2006/03/17
// RCS-ID:      $Id: ratedata.h,v 1.6 2006/08/21 14:27:26 zhang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/ratedata.h
    @brief Class containing all rate related information

    Rate information includes yield curves and exchange rates.
 */

#ifndef _ITO33_FINANCE_RATEDATA_H_
#define _ITO33_FINANCE_RATEDATA_H_

#include "ito33/beforestd.h"
#include <map>
#include <string>
#include "ito33/afterstd.h"

#include "ito33/sharedptr.h"
#include "ito33/dlldecl.h"

namespace ito33
{

namespace XML
{
  class Tag;
}


namespace finance
{

class ITO33_DLLDECL YieldCurve;
class ITO33_DLLDECL MoneyMarket;
class ITO33_DLLDECL Numeraire;
class ITO33_DLLDECL SpotFXRates;

/**
    Class containing rate related information needed for pricing.
    
    Includes yieldcurves and exchange rates.
 */
class ITO33_DLLDECL RateData
{

public:
  
  /// typedef to simplify the code
  typedef std::string NumeraireCode; 

  /// The container type used for yield curve storage
  typedef std::map< NumeraireCode, shared_ptr<MoneyMarket> > YCElements;  

  /**
      Constructs a RateData object.

      Initially empty.  Rate data is currency (numeraire) specific, so 
      must be initialized via set functions.
   */
  RateData();


  /**
      @name Modifiers.
   */
  //@{

  /**
      The YieldCurve of the specified currency (numeraire).

      @param pNumeraire the currency
      @param pYieldCurve the yield curve of the currency
   */
  void SetYieldCurve(const shared_ptr<Numeraire>& pNumeraire, 
                     const shared_ptr<YieldCurve>& pYieldCurve);
  
  /**
      The YieldCurve of the specified currency (numeraire).

      @param currency the currency code
      @param pYieldCurve the yield curve of the currency
      
      @noexport
   */
  void SetYieldCurve(const NumeraireCode& currency, 
                     const shared_ptr<YieldCurve>& pYieldCurve);

  /**
      The spot foreign exchange rates.
     
      By using a SpotFXRate container, all rates can be updated in one call.

      @param pSpotFXRates the foreign exchange rates
   */
  void SetSpotFXRates(const shared_ptr<SpotFXRates>& pSpotFXRates);

  //@}


  /**
      @name Accessors.
   */
  //@{

  /**
      The YieldCurve of specified currency (numeraire).

      Throws an exception if the yield curve is not set.

      @param pNumeraire the currency of the desired yieldcurve
      @return the yieldcurve of the specified currency
   */
  const shared_ptr<YieldCurve>& 
  GetYieldCurve(const shared_ptr<Numeraire>& pNumeraire);

  /**
      The spot foreign exchange rates.     

      @return the foreign exchange rates
   */
  const shared_ptr<SpotFXRates>& GetSpotFXRates() const
  {
    return m_pSpotFXRates;
  }
 
  //@}

  /**
      @internal
      @brief Yield curve for mesh construction

      @param pNumeraire the currency
      @param pYieldCurve the yield curve for mesh construction of the currency 

      @noexport
   */
  void SetYieldCurveForMesh(const shared_ptr<Numeraire>& pNumeraire, 
                            const shared_ptr<YieldCurve>& pYieldCurve);
 
  /**
      @internal
      @brief Yield curve for mesh construction

      @param pNumeraire the currency
      @return the yieldcurve for mesh construction of the currency 

      @noexport
   */
  const shared_ptr<YieldCurve>& 
  GetYieldCurveForMesh(const shared_ptr<Numeraire>& pNumeraire);

  /**
      @internal

      @noexport
   */
  void Dump(XML::Tag& tagParent) const;

  /**
      @iternal
      @brief Gets collection of yield curves.

      @return collection of yield curves.

      @noexport
   */
  const YCElements& GetAll() const
  {
    return m_MoneyMarkets;
  }

private:

  /// the foreign exchange rates
  shared_ptr<SpotFXRates> m_pSpotFXRates;
 
  /// storage for the MoneyMarkets
  YCElements m_MoneyMarkets;

}; // class RateData


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_RATEDATA_H_
