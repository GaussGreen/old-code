/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/voltimeonlycalibrator_useoptions.h
// Purpose:     class for time only volatility calibration using options   
// Author:      Wang, David
// Created:     2004/05/19
// RCS-ID:      $Id: voltimeonlycalibrator_useoptions.h,v 1.8 2006/08/20 09:47:41 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/ihg/voltimeonlycalibrator_useoptions.h
   @brief class for time only volatility calibration using options

   @todo some comment need to be corrected.
 */


#ifndef _ITO33_IHG_VOLTIMEONLYCALIBRATOR_USEOPTIONS_H_
#define _ITO33_IHG_VOLTIMEONLYCALIBRATOR_USEOPTIONS_H_

#include "ito33/vector.h"
#include "ito33/sharedptr.h"
#include "ito33/array.h"

#include "ito33/finance/option.h"
#include "ito33/finance/modeloutput.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/volatilitytimeonly.h"
#include "ito33/ihg/volatilitycombo.h"
#include "ito33/ihg/modeloutput.h"


namespace ito33
{

namespace ihg
{

/// class for time only volatility calibration using options
class VolTimeOnlyCalibratorUseOptions 
{

public:

  typedef std::vector< shared_ptr<finance::Option> > OptionElements;

  /**
     Ctor takes a list of options to be calibrated

     @param pOptions A list of shared pointer of options, required to be 
                 ordered by maturity
   */
  VolTimeOnlyCalibratorUseOptions(const OptionElements &pOptions,
                                  shared_ptr<ihg::HazardRate>& pHazardRate);

  /**
     Run the calibrator

     @return the calibrated time only volatility pointer
   */
  virtual shared_ptr<VolatilityTimeOnly> Calibrate();

  /**
     Set an external volatility to be combined with the time only values.

     @param pVol The external volatility to be combined with internal time only
  */
  void SetExternalVolatility(shared_ptr<Volatility>& pVol)
  {
    m_pExternalVol = pVol;
  }


  /**
     Helper operator to be used for numerical calibrator

     @return the difference of the theoritical price and the market price
             for the current option
   */
  double operator()(double dVolatility)
  {
    m_pdVolatilities[m_nIdxOption] = dVolatility;
 
    shared_ptr<VolatilityTimeOnly>
      pVolatility = new VolatilityTimeOnly
                        (
                          m_pdMaturities.Get(), m_pdVolatilities.Get(), 
                          m_nIdxOption + 1
                        );

    // check if we have an external volatility to use during calibration
    if ( !m_pExternalVol )
      m_theoreticalModel.SetVolatility(pVolatility);
    else
      m_theoreticalModel.SetVolatility
      ( new VolatilityCombo(pVolatility, m_pExternalVol) );

    shared_ptr<ihg::ModelOutput> 
      output = m_theoreticalModel.PriceOption(*m_pOptions[m_nIdxOption]);

    double dScale = m_dMarketPrice;
    if ( m_dMarketPrice < 1.0 )
      dScale = 1.0;

    return (output->GetPrice() - m_dMarketPrice)/dScale; 
  }


private:

  /// Model used during calibration. This class updates the hazard rate
  ihg::TheoreticalModel m_theoreticalModel;

  /// a list of options to be calibrated
  OptionElements m_pOptions;

  /// the market price of the current option to be calibrated  
  double m_dMarketPrice;

  /// helper for index of the current option
  size_t m_nIdxOption;

  /// The number of the maturities
  size_t m_nNbMaturities;

  /// the calibrated hazard rates
  Array<double> m_pdVolatilities;

  /// the maturities of the option
  Array<double> m_pdMaturities;

  /// External volatility to be combined with the time only vol
  shared_ptr<Volatility> m_pExternalVol;

  NO_COPY_CLASS(VolTimeOnlyCalibratorUseOptions);

}; // class VolTimeOnlyCalibratorUseOptions


} // namespace ihg

} // namespace ito33


#endif // #ifndef _ITO33_IHG_VOLTIMEONLYCALIBRATOR_USEOPTIONS_H_

