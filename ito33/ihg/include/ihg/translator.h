/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/translator.h
// Purpose:     translate between array of parameter values and an ihg model
// Created:     2005/07/04
// RCS-ID:      $Id: translator.h,v 1.8 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/translator.h
    @brief translate between an array of parameter values and an ihg model
 */

#ifndef _IHG_TRANSLATOR_H_
#define _IHG_TRANSLATOR_H_

#include "ito33/vector.h"
#include "ito33/sharedptr.h"

#include "ito33/pricing/translator.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL Derivatives;
}

namespace ihg
{

  enum VolType
  {
    VolType_flat,
    VolType_power,
    VolType_tanh
  };

  enum HRType
  {
    HRType_flat,
    HRType_power,
  };

  enum SpotType
  {
    SpotType_power
  };

  class TheoreticalModel;
  class Volatility;
  class HazardRate;

/**
   Class helps to translate between a parameter array and an ihg model
 */
class Translator : public pricing::Translator
{
  
public:

  /**
     ctor takes the volatility type, hazard rate type, and session (for spot).
     
     @param volType The type of volatility to use
     @param hrType The type of hazard rate to use
     @param dSpot The current spot share price
   */
  Translator(VolType volType, HRType hrType, double dSpot);

  /**
     ctor takes the volatility and spot component for hazard rate.
     
     @param volType The type of volatility to use
     @param spotType The type of hazard rate spot component to use
     @param derivatives The derivative list (for access to maturity dates)
     @param dSpot The current spot share price
   */
  Translator(VolType volType, 
             SpotType spotType, 
             finance::Derivatives derivatives,
             double dSpot);


  /**
     Translate the parameters into an IHG model.

     @param pdX the known parameters with a well defined order
     @return a new IHG model corresponding to the parameters
   */
  shared_ptr<ihg::TheoreticalModel> operator()(const std::vector<double>& pdX);

  /**
     Translate the parameters into an IHG model.

     @param pdX the known parameters with a well defined order
     @return a new IHG model corresponding to the parameters
   */
  shared_ptr<ihg::TheoreticalModel> operator()(const double* pdX);

  /**
     Get calibration parameters from the underlying model structure.

     @param pdX The array to be filled by the calibration parameters.
   */
  void GetParameters(double* pdX) const;

  /**
     Get lower bounds for the parameters.

     @return Lower bounds for the model parameters
  */
  std::vector<double> GetLowerBounds() const;

  /**
     Get upper bounds for the parameters.

     @return Upper bounds for the model parameters
  */
  std::vector<double> GetUpperBounds() const;


protected:

  /**
   Initialize the volatility params (including bounds).

   @param nCounter Position in arrays to start updating values
  */
  void InitVolParams(size_t nCounter);

  /**
   Initialize the hazard rate params (including bounds).

   @param nCounter Position in arrays to start updating values
  */
  void InitHRParams(size_t nCounter);

  /**
   Initialize the spot hazard rate params (including bounds).

   @param nCounter Position in arrays to start updating values
  */
  void InitSpotHRParams(size_t nCounter);

  /**
   Initialize the time component params (including bounds).

   @param nCounter Position in arrays to start updating values
  */
  void InitTimeComponentParams(size_t nCounter, size_t nNbParams);

  /**
   Extract unique maturity dates from the derivative list

   @param pDerivatives the list of derivatives
   */
  std::vector<Date> GetMaturityDates(finance::Derivatives derivatives);

  /** 
     Restore the volatility.

     @param pdX the parameter values to construct the model
     @param nFlagCounter current index into flag array
     @param nXCounter current index into pdX array
     @return The appropriate volatility class
   */
  shared_ptr<Volatility> RestoreVolatility(const double* pdX,
                                          size_t& nFlagCounter, 
                                          size_t& nXCounter) const;

  /** 
     Restore the hazard rate.

     @param pdX the parameter values to construct the model
     @param nFlagCounter current index into flag array
     @param nXCounter current index into pdX array
     @return The appropriate hazard rate class
   */
  shared_ptr<HazardRate> RestoreHazardRate(const double* pdX,
                                          size_t& nFlagCounter, 
                                          size_t& nXCounter) const;

  /** 
     Restore hazard rate combo.

     @param pdX the parameter values to construct the model
     @param nFlagCounter current index into flag array
     @param nXCounter current index into pdX array
     @return The appropriate hazard rate class
   */
  shared_ptr<HazardRate> RestoreComboHR(const double* pdX,
                                       size_t& nFlagCounter, 
                                       size_t& nXCounter) const;



  /// The type of volatility used in the model
  VolType m_VolType;

  /// The type of hazard rate used in the model
  HRType m_HRType;

  /// The spot type for a hazard rate combo
  SpotType m_spotType;

  /// The parameters used to construct the volatility and hazard rate
  std::vector<double> m_pdParams;

  /// The lower bounds of the parameters
  std::vector<double> m_pdLowerBounds;

  /// The upper bounds of the parameters
  std::vector<double> m_pdUpperBounds;

  /// The dates/times used for time component of a combo hazard rate
  std::vector<Date> m_TimeComponentDates;

  /// The current spot (used by several volatility and hazard rate types)
  double m_dSpot;
  
private:

  NO_COPY_CLASS(Translator);

}; // class Translator


} // namespace ihg

} // namespace ito33

#endif // _IHG_TRANSLATOR_H_
