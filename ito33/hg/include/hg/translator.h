/////////////////////////////////////////////////////////////////////////////
// Name:        hg/translator.h
// Purpose:     translate to/from parameter values and an underlying process
// Created:     2005/07/04
// RCS-ID:      $Id: translator.h,v 1.8 2006/08/19 23:46:52 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/translator.h
    @brief translate to/from parameter values and an underlying process
 */

#ifndef _HG_TRANSLATOR_H_
#define _HG_TRANSLATOR_H_

#include "ito33/vector.h"
#include "ito33/sharedptr.h"

#include "ito33/pricing/translator.h"
#include "ito33/hg/underlyingprocess.h"

namespace ito33
{

namespace hg
{


/**
   Class helps to translate between a parameter array and an underlying process

   Only one order is implemented so far.
 */
class Translator : public pricing::Translator
{
  
public:

  /**
     ctor takes the structure of an well defined underlying process.
     
     @param pUP The underlying process defining the structure of jumps
   */
  Translator(const UnderlyingProcess& underlyingProcess);

  /**
     Translate the parameters into an underlying process.

     @param pdX the known parameters with a well defined order
     @return a new underlying process corresponding to the parameters
   */
  shared_ptr<UnderlyingProcess> operator()(const std::vector<double>& pdX);

  /**
     Translate the parameters into an underlying process.

     @param pdX the known parameters with a well defined order
     @return a new underlying process corresponding to the parameters
   */
  shared_ptr<UnderlyingProcess> operator()(const double* pdX);
 
  /**
     Get calibration parameters from the underlying process structure.

     @param pdX The array to be filled by the calibration parameters.
   */
  void GetParameters(double* pdX) const;

  /**
     Get calibration parameters from the underlying model structure.

     @return The calibration parameters as a vector
   */
  std::vector<double> GetParameters() const;

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

  const UnderlyingProcess& m_underlyingProcess;
  
private:

  NO_COPY_CLASS(Translator);

}; // class Translator


} // namespace hg

} // namespace ito33

#endif // _HG_TRANSLATOR_H_
