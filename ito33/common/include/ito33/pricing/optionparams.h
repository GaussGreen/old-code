/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/optionparams.h
// Purpose:     option params class
// Created:     2004/02/11
// RCS-ID:      $Id: optionparams.h,v 1.20 2006/08/19 22:18:03 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/optionparams.h
    @brief option params class

    Implementation of the params class for options.
 */

#ifndef _ITO33_PRICING_OPTIONPARAMS_H_
#define _ITO33_PRICING_OPTIONPARAMS_H_

#include "ito33/autoptr.h"
#include "ito33/list.h"

#include "ito33/pricing/option.h"
#include "ito33/pricing/params.h"

namespace ito33 
{
 
namespace finance
{
  class Payoff;

  class ITO33_DLLDECL SessionData;
}

namespace pricing
{


/// Pricing parameters for European and American options
class OptionParams : public Params
{

public: 

  /**
      Ctor for option contract objet. 

      @param option reference to option contract objet
   */
  OptionParams(Option& option) : Params(option),
                                 m_option(option),
                                 m_clonedOption(0)
  { 
  }

  /**
      Ctor by option contract objet. This Ctor is essentially the same
      as above, except that the memory for the contract class will
      be managed internally through the autoptr.

      @param Option autoptr to Option contract objet
   */
  OptionParams(AutoPtr<Option> option) 
    : Params(*option),
      m_option(*option), 
      m_clonedOption(option)
  {  
  }

  /**
      Ctor by option contract and common financial and numerical datas. 

      @param option reference to Option contract objet
      @param sessionData reference to financial session data
      @param pNumParams numeric parameters
      @param pMeshParams parameters for mesh builder
   */
  OptionParams(Option& option,
               const finance::SessionData& sessionData,
               const shared_ptr<numeric::NumParams>& pNumParams,
               const shared_ptr<numeric::MeshParams>& pMeshParams)
    : Params(option, sessionData, pNumParams, pMeshParams),
      m_option(option), m_clonedOption(0)
  {
  }

  virtual ~OptionParams() { }

  virtual void Init();

  /// Get a reference to the option contract
  Option& GetOption() const 
  { 
    return m_option; 
  }

  /// Get the payoff pointer
  const finance::Payoff* GetPayoff() const 
  { 
    return m_pPayoff.get(); 
  }
  
  /// Set the payoff pointer for non-standard payoff conditions
  void SetPayoff(const shared_ptr<finance::Payoff>& pPayoff)
  {
    m_pPayoff = pPayoff;
  }

  /// If the object is cloned, then the clone needs to manage memory
  AutoPtr<OptionParams> Clone();


  /// returns the strike of the option
  virtual double GetStrike() const
  {
    return m_option.GetStrike();
  }

protected:

  virtual void ConstructPayoff();

  /// Helper for initial values and constraints
  shared_ptr<finance::Payoff> m_pPayoff;

  Option& m_option;

  /// If the object is cloned, then the clone needs to manage memory
  AutoPtr<Option> m_clonedOption;

private:

  NO_COPY_CLASS(OptionParams);

}; // class OptionParams


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_OPTIONPARAMS_H_
