/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/conversionpricereset.h
// Purpose:     Characteristics of a single reset date 
//              that is applied to the conversion price
// Author:      Yann and David
// Created:     2004/10/20
// RCS-ID:      $Id: conversionpricereset.h,v 1.9 2005/04/21 17:19:48 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/bondlike/conversionpricereset.h
    @brief characteristics of a single conversion price reset date
 */

#ifndef _ITO33_FINANCE_BONDLIKE_CONVERSIONPRICERESET_H_
#define _ITO33_FINANCE_BONDLIKE_CONVERSIONPRICERESET_H_

#include "ito33/sharedptr.h"
#include "ito33/date.h"

#include "ito33/dlldecl.h"


namespace ito33
{

namespace XML
{
  class Tag;
}

namespace finance
{

/**
   Characteristics of a single reset date.
 */
class ITO33_DLLDECL ConversionPriceReset
{
public:
  
  /**
     Constructor.

     @param resetDate the date of the reset
     @param dFloor The floor conversion price. It is expressed as a 
            percentage of the initial or prevailing conversion price
   */
  ConversionPriceReset(Date resetDate, double dFloor);

  /// virtual dtor for base class
  virtual ~ConversionPriceReset() { }


  /**
     Gets the date at which the conversion price is reset

     @return The date at which the conversion price is reset
   */
  Date GetDate() const { return m_resetDate; }

  /**
     The cap conversion price. It is expressed as a percentage of the 
     prevailing conversion price

     @return The cap 
   */
  double GetCap() const { return m_dCap; }

  /**
     The floor conversion price. It is expressed as a percentage of the 
     initial or prevailing conversion price

     @return The floor rate
   */
  double GetFloor() const { return m_dFloor; }

  /**
     The multiplier used to compute new conversion prices.
     @param dMultiplier The multiplier when computing the new price
   */
  void SetMultiplier(double dMultiplier);

  /**
     The multiplier used to compute new conversion prices.

     @return The multiplier used to compute new conversion prices
   */
  double GetMultiplier() const { return m_dMultiplier; }

  /**
     The cap conversion price. It is expressed as a percentage of the 
     prevailing conversion price

     @param dCap the cap conversion price
   */
  void SetCap(double dCap);

  /**
     @internal
     @brief Writes myself to parent tag

     @param tagParent parent tag
     @return the tag created for the reset terms

     @noexport
   */
  XML::Tag Dump(ito33::XML::Tag& tagParent) const;
 

protected:

  /// reset dates
  Date m_resetDate;
  
  /// cap 
  double m_dCap;

  /// floor 
  double m_dFloor;

  /// multiplier for computing new conversion ratios
  double m_dMultiplier;

}; // class ConversionPriceReset


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_BONDLIKE_CONVERSIONPRICERESET_H_
