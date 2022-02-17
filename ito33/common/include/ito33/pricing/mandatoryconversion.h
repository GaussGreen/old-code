/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/mandatoryconversion.h
// Purpose:     Mandatory conversion class
// Author:      Wang
// Created:     2004/08/16
// RCS-ID:      $Id: mandatoryconversion.h,v 1.22 2006/03/22 13:10:23 yann Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_PRICING_MANDATORYCONVERSION_H_
#define _ITO33_PRICING_MANDATORYCONVERSION_H_

#include "ito33/date.h"

#include "ito33/finance/bondlike/triggeraspercentageof.h"

#include "ito33/numeric/mesh/roots.h"

#include "ito33/pricing/conversionprovisions.h"
#include "ito33/pricing/mandatory_payoff_structure.h"


namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL PEPSLike;
  class ITO33_DLLDECL GeneralizedPEPSLike;
  class ITO33_DLLDECL PERCSLike;
}

namespace pricing
{

class CBCalls;

class MandatoryConversion : public ConversionProvisions
{
public:

  MandatoryConversion(): ConversionProvisions() {  }

  MandatoryConversion(const finance::PEPSLike& peps);

  MandatoryConversion(const finance::GeneralizedPEPSLike& peps);

  MandatoryConversion(const finance::PERCSLike& percs);

  // Default dtor is ok
  

  /// @name implement virtual functions
  //@{
  virtual void ComputeRoots(double dTime, 
                            size_t &nNbRoots, numeric::mesh::Root *pRoots, 
                            bool bPlus = true);

  virtual bool GetGrossParities
       (const double* pdS, size_t nNbS, 
        const double* pdNewSharePrices, double* pdValues) const;
  
  virtual double 
  GetConversionPrice
  (finance::TriggerAsPercentageOf triggerAsPercentageOf) const;

  virtual double 
  GetConversionPrice
  (double dTime, finance::TriggerAsPercentageOf triggerAsPercentageOf, 
   bool bPlus = true) const;
   
  //@}

  /// @name helper function to access the mandatory payoff structure
  //@{
  double GetConversionRatio(double dS) const;

  //@}

  /// @name Setters
  //@{

  /**
     Indicates that the Average of the stock should be used
     instead of the stock value for conversion.
   */
  void SetStockAverage(double dStockAverage);
   
  /**
     Indicates that the average conversion ratio must be used
     instead of the stock value for conversion
   */
  void SetConversionRatioAverage(double dConversionRatio);

  virtual void SetRatios(double /* dRatio */)
  {
    FAIL("Mandatory Conversion: SetRatios function should not be called.");
  }


  //@}

private:

  void ComputeRootsWithCallCash(const CBCalls& calls, double dTime, 
                            size_t &nNbRoots, numeric::mesh::Root *pRoots, 
                            bool bPlus = true);

  void 
  ComputeRootsWithCallTrigger
  (double dTriggerRate, finance::TriggerAsPercentageOf asPercentagOf, 
   size_t &nNbRoots, numeric::mesh::Root *pRoots, 
   double dTime, bool bPlus = true);

  std::vector<MandatoryPayoffStructure> m_pPeriods;

  double m_dConversionPrice;

  double m_dRatioForRootComputation;

}; // MandatoryConversion


} // namespace pricing

} // namespace ito33

#endif  //  #ifndef _ITO33_PRICING_MANDATORYCONVERSION_H_
