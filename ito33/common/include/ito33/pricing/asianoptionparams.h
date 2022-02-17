/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/asianoptionparams.h
// Purpose:     Asian option params class
// Author:      ITO 33 Canada
// Created:     April 6, 2005
// RCS-ID:      $Id: asianoptionparams.h,v 1.9 2006/08/19 22:18:03 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/asianoptionparams.h
    @brief asian option params class

    Implementation of the params class for asian options.
 */

#ifndef _ITO33_PRICING_ASIANOPTIONPARAMS_H_
#define _ITO33_PRICING_ASIANOPTIONPARAMS_H_

#include "ito33/vector.h"
#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"
#include "ito33/list.h"

#include "ito33/pricing/asianoption.h"
#include "ito33/pricing/optionparams.h"

namespace ito33 
{
 
namespace finance
{
  class Payoff;
  class ITO33_DLLDECL SessionData;
}

namespace pricing
{

  class Model;

/// Pricing parameters for Asian Option
class AsianOptionParams : public OptionParams
{

public: 

  /**
      Ctor for asian option contract objet. 

      @param asian option reference to option contract objet
   */
  AsianOptionParams(AsianOption& asianOption) 
    : OptionParams(asianOption), 
      m_asianOption(asianOption),
      m_clonedAsianOption(0), 
      m_dStrike(asianOption.GetStrike()) ,
      m_dAverage(0.)
  { 
  }

  /**
      Ctor by asian option contract objet. This Ctor is essentially the same
      as above, except that the memory for the contract class will
      be managed internally through the autoptr.

      @param asianOption autoptr to AsianOption contract objet
   */
  AsianOptionParams(AutoPtr<AsianOption> asianOption) 
    : OptionParams(*asianOption),
      m_asianOption(*asianOption), 
      m_clonedAsianOption(asianOption),
      m_dStrike( asianOption->GetStrike() ),
      m_dAverage(0.)
  {  
  }

  /**
      Ctor by Asianoption contract and common financial and numerical datas. 

      @param asianOption reference to asian Option contract objet
      @param sessionData reference to financial SessionData
      @param pNumParams numeric parameters
      @param pMeshParams parameters for mesh builder
   */
  AsianOptionParams(AsianOption& asianOption,
                    const finance::SessionData& sessionData,
                    const shared_ptr<numeric::NumParams>& pNumParams,
                    const shared_ptr<numeric::MeshParams>& pMeshParams)
    : OptionParams(asianOption, sessionData, pNumParams, pMeshParams),
      m_asianOption(asianOption), m_clonedAsianOption(0),
      m_dStrike(asianOption.GetStrike()),m_dAverage(0.)
  {
  }

  /// destructor
  ~AsianOptionParams() { }


  /// Gets a reference to the option contract
  AsianOption& GetAsianOption() const 
  { 
    return m_asianOption; 
  }

  /// Clone
  AutoPtr<AsianOptionParams> Clone();

  /**
      Sets the strike.

      @param strike value
   */
  void SetStrike(double dStrike)
  {
    m_dStrike = dStrike;
  }

  /** 
      Gets the strike.

      @return the strike value
   */
  double GetStrike() const
  {
    return m_dStrike;
  }

  /**
      Sets the Average, this is needed for fixed Asian Strike
      option type.
 
      @param dAverage average value
   */
  void SetAverage(double dAverage)
  {
    m_dAverage = dAverage;
  }

  /**
      Gets the average of the path.

      @return average of the path
   */
  double GetAverage() const
  {
    return m_dAverage;
  }
  
  /**
      Constructs the averaging grid.

      @param pdGridY constructed grid
      @param model pricing model
      @param bHasSimilarityReduction indicate if a 
             similarity reduction is possible
   */
  void ConstructPathDepGrid(std::vector<double>& pdGridY, 
                            Model& model, bool bHasSimilarityReduction);

  /**
      Constructs the events.

      @param bHasSimilarityReduction 
             if a similarity reduction is feasible

      @return list of path dependent events
   */
  std::list< shared_ptr<PathDepEvent> > 
  ConstructEvents(bool bHasSimilarityReduction);

  /**
      Used to determined whether 
      or not we can do a similarity transform.

      @return true/false
   */
  bool IsPathDependent() const;

  /** 
      Constructs parameters vector for Asian option
      based on the Averaging grid.

      @param pdGridY averaging grid

      @return vector of auto pointers of options param
   */
  std::vector< AutoPtr<OptionParams> > 
  ConstructParam(const std::vector<double>& pdGridY);

  /**
      Gets the path to save.

      @param pdGridY averaging grid

      @return path to save
   */
  size_t GetPathToSave(const std::vector<double>& pdGridY) const;

protected:

  /// Strike of the option 
  ///  if fixed strike then strike is strike :)
  ///  if floating then strike is set to current average
  double m_dStrike;

  /// Average of the Path
  double m_dAverage;

  /// Construct payoff for fixed strike or floating strike
  void ConstructPayoff();

  /// AsianOption object
  AsianOption& m_asianOption;

  /// If the object is cloned, then the clone needs to manage memory
  AutoPtr<AsianOption> m_clonedAsianOption;

private:

  NO_COPY_CLASS(AsianOptionParams);

}; // class AsianOptionParams


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_ASIANOPTIONPARAMS_H_
