/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/spotfxrates.h
// Purpose:     Class for current foreign exchange rate matrix
// Author:      Wang
// Created:     2004/09/02
// RCS-ID:      $Id: spotfxrates.h,v 1.18 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/finance/spotfxrates.h
   @brief Class for current foreign exchange rate matrix.

   @sa Numeraire
 */

#ifndef _ITO33_FINANCE_SPOTFXRATES_H_
#define _ITO33_FINANCE_SPOTFXRATES_H_

#include "ito33/beforestd.h"
#include <map>
#include "ito33/afterstd.h"

#include "ito33/dlldecl.h"
#include "ito33/sharedptr.h"
#include "ito33/string.h"

namespace ito33
{

namespace XML
{
  class Tag;
}

namespace finance
{

class ITO33_DLLDECL Numeraire;

/**
   Class for foreign exchange rate matrix.  
 */
class ITO33_DLLDECL SpotFXRates 
{

public:
  
  /// typedef to simplify the code
  typedef std::string NumeraireCode; 

  /// typedef to simplify the code
  typedef std::pair<NumeraireCode, NumeraireCode> NumerairePair;

  /// The container type used for spot FX rate storage
  typedef std::map<NumerairePair, double> Elements;

  /**
     Constructs a SpotFXRates object.

     Initially empty.  Exchange rate data is currency (numeraire) 
     specific, so initialization must be done via set functions.
   */
  SpotFXRates() { }

  /**
      @name Modifiers.
   */
  //@{

  /**
     The spot FX rate between two currencies.
     
     This is the number of units of base currency needed to purchase 
     one unit of the foreign currency (indirect quotation).

     @param pNumeraire1 the foreign currency
     @param pNumeraire2 the base currency
     @param dRate the spot FX rate between two currencies.
   */
  void SetFXRate(const shared_ptr<Numeraire>& pNumeraire1, 
                 const shared_ptr<Numeraire>& pNumeraire2,
                 double dRate);

  //@}


  /**
      @name Accessors.
   */
  //@{

  /**
     The spot FX rate between two currencies.
     
     This is the number of units of base currency needed to purchase 
     one unit of the foreign currency (indirect quotation).
     Throws an exception if the FX rate is not set.

     @param pNumeraire1 the foreign currency
     @param pNumeraire2 the base currency
     @return the spot FX rate between two currencies.
   */
  double GetFXRate(const shared_ptr<Numeraire>& pNumeraire1, 
                   const shared_ptr<Numeraire>& pNumeraire2);
  
  //@}

  /**
     @internal
     @brief Dumps all data stored in this object in XML format.

     @param tagParent the parent tag under which our tag(s) should be created
      
     @noexport
   */
  XML::Tag Dump(XML::Tag& tagParent) const;


private:

  /// storage for the spot exchange rates
  Elements m_FXRates;

}; // class SpotFXRates


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_SPOTFXRATES_H_

