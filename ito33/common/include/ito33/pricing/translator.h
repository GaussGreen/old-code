/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/translator.h
// Purpose:     translate between an array of parameter into model params
// Created:     2005/07/04
// RCS-ID:      $Id: translator.h,v 1.2 2006/06/01 21:13:51 dave Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/translator.h
    @brief translate between an array of parameter into model params
 */

#ifndef _ITO33_PRICING_TRANSLATOR_H_
#define _ITO33_PRICING_TRANSLATOR_H_

#include "ito33/vector.h"
#include "ito33/sharedptr.h"


namespace ito33
{

namespace pricing
{


/**
   Class helps to translate between a parameter array and models
 */
class Translator
{
  
public:

  /**
     Empty constructor.
   */
  Translator() {}

  virtual ~Translator() {}

  
  /**
     Get calibration parameters from the underlying model structure.

     @param pdX The array to be filled by the calibration parameters.
   */
  virtual void GetParameters(double* pdX) const = 0;

  /**
     Get lower bounds for the parameters.

     @return Lower bounds for the model parameters
  */
  virtual std::vector<double> GetLowerBounds() const = 0;

  /**
     Get upper bounds for the parameters.

     @return Upper bounds for the model parameters
  */
  virtual std::vector<double> GetUpperBounds() const = 0;

  /**
     Sets the translatation flag for each parameter. 

     @param flags The translation flags for each parameter
   */
  void SetFlags(const std::vector<bool>& flags);

  /**
     Gets the translation flags for each parameter. 

     The flags can then be used by the pricing process to determine which
     sensitivity is required.
   */
  const std::vector<bool>& GetFlags() const
  {
    return m_flags;
  }

  /**
     Get the number of active parameters needed by the model

     @return The number of active parameters
   */
  size_t GetNbParameters() const;

  /**
     Get calibration parameters from the underlying model structure.

     @return The calibration parameters as a vector
   */
  std::vector<double> GetParameters() const;

  /**
     Apply the lower and upper bounds to the specified parameters.

     @param pdX the original parameters

     @return the parameters with the bounds applied (if neccessary)
  */
  std::vector<double> ApplyBounds(const double* pdX);

protected:

  /// The number of flags, equivalent to total number of model parameters
  size_t m_nNbFlags;

  /// Flags indicating which params are frozen
  std::vector<bool> m_flags;

 
}; // class Translator


} // namespace pricing

} // namespace ito33

#endif // _ITO33_PRICING_TRANSLATOR_H_
