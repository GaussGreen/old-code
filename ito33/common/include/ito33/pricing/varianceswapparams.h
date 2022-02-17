/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/varianceswapparams.h
// Purpose:     variance swap params class
// Created:     2006/02/21
// RCS-ID:      $Id: varianceswapparams.h,v 1.13 2006/08/19 22:18:03 wang Exp $
// Copyright:   (c) 2006 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/varianceswapparams.h
    @brief variance swap params class

    Implementation of the params class for variance swaps.
 */

#ifndef _ITO33_PRICING_VARIANCESWAPPARAMS_H_
#define _ITO33_PRICING_VARIANCESWAPPARAMS_H_

#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"

#include "ito33/pricing/varianceswap.h"
#include "ito33/pricing/params.h"

namespace ito33 
{

namespace pricing
{


/// Pricing parameters for variance swaps.
class VarianceSwapParams : public Params
{

public: 

  /**
      Ctor by VarianceSwap contract and common financial and numerical datas.

      @param varianceSwap reference to variance swap contract object
      @param sessionData reference to financial session data
      @param pNumParams numeric parameters
      @param pMeshParams parameters for mesh builder
   */
  VarianceSwapParams(VarianceSwap& varianceSwap,
                     const finance::SessionData& sessionData,
                     const shared_ptr<numeric::NumParams>& pNumParams,
                     const shared_ptr<numeric::MeshParams>& pMeshParams)
    : Params(varianceSwap, sessionData, pNumParams, pMeshParams),
      m_varianceSwap(varianceSwap), 
      m_clonedVarianceSwap(0),
      m_dConditionalPayoff(-1)
  {
  }

  /**
      Ctor by variance swap contract object. This Ctor is essentially the 
      same as above, except that the memory for the contract class will
      be managed internally through the autoptr.

      @param varianceSwap autoptr to VarianceSwap contract object
   */
  VarianceSwapParams(AutoPtr<VarianceSwap> varianceSwap) 
    : Params(*varianceSwap),
      m_varianceSwap(*varianceSwap), 
      m_clonedVarianceSwap(varianceSwap),
      m_dConditionalPayoff(-1)
  {  
  }

  virtual ~VarianceSwapParams() { }

  /**
      Sets the average of the squared returns to date.

      @param dAvgSqrReturn the average of the squared returns
   */
  void SetAvgSqrReturn(double dAvgSqrReturn)
  {
    m_dAvgSqrReturns = dAvgSqrReturn;
  }

  /**
      Sets the spot price at the previous observation.

      @param dPreviousSpot the spot price at the previous observation
   */
  void SetPreviousSpot(double dPreviousSpot)
  {
    m_dPreviousSpot = dPreviousSpot;
  }

  /**
      Sets the conditional variance swap payoff value.

      @param dPayoff the conditional swap payoff value
   */
  void SetConditionalPayoff(double dPayoff)
  {
    m_dConditionalPayoff = dPayoff;
  }

  /**
      Gets the average of the squared returns to date.

      @return the average of the squared returns
   */
  double GetAvgSqrReturn() const
  {
    return m_dAvgSqrReturns;
  }

  /**
      Gets the spot price at the previous observation.

      @return the spot price at the previous observation
   */
  double GetPreviousSpot() const
  {
    return m_dPreviousSpot;
  }

  /**
      Gets the conditional variance swap payoff.

      @return the conditional variance swap payoff
   */
  double GetConditionalPayoff() const
  {
    return m_dConditionalPayoff;
  }

  /**
      Gets a reference to the variance swap contract.

      @return reference to the underlying variance swap contract
   */
  VarianceSwap& GetVarianceSwap() const 
  { 
    return m_varianceSwap; 
  }

  /**
      Gets the payoff.

      @return the payoff
   */
  const shared_ptr<finance::Payoff>& GetPayoff() const 
  { 
    return m_pPayoff; 
  }

  /**
      Checks if a similarity reduction can be used.

      A similarity reduction can be used to reduce the problem by
      one dimension if:
      - there are no cash dividends between valuation and maturity
      - the variance swap does not have any exotic features

      @return true if a similarity reduction is possible, false otherwise
   */
  bool HasSimilarityReduction();

  /**
      The master space grid.

      This grid is used to create the previous spot grid, and is also
      used to create the space grids for each path. 

      @param model the pricing model

      @return the grid
   */
  std::vector<double> ConstructMasterGrid(const Model& model);

  /**
      Grid for the average of the squared returns.

      @return the average squared return grid
   */
  std::vector<double> ConstructAvgSqrReturnGrid();

  /**
      Grid for the previous spots.

      @param model the pricing model
      @param bIsSimilarityReduction indicate whether or not
                 a similarity reduction can be used
      @return the grid
   */
  std::vector<double>
  ConstructPreviousSpotGrid(const Model& model, bool bIsSimilarityReduction);

  /**
      Constructs the params for each path.

      @param pdAvgSqrReturnGrid grid for the average of the squared returns
      @param pdPreviousSpotGrid grid for the previous spot prices

      @return vector of variance swap param auto pointers
   */
  std::vector< AutoPtr<VarianceSwapParams> > 
  ConstructParams(const std::vector<double>& pdAvgSqrReturnGrid, 
                  const std::vector<double>& pdPreviousSpotGrid);

  /**
      Constructs the params for each path for the conditional var swap.

      @param pdYGrid grid for count within the corridor

      @return vector of variance swap param auto pointers
   */
  std::vector< AutoPtr<VarianceSwapParams> > 
  ConstructConditionalParams(const std::vector<double>& pdYGrid);

  /**
      Gets the path to save.

      @param pParams the vector of params
      @param pdAvgSqrReturnGrid grid for the average of the squared returns
      @param pdPreviousSpotGrid grid for the previous spot prices

      @return the path to save
   */
  size_t 
  GetPathToSave(const std::vector<double>& pdAvgSqrReturnGrid, 
                const std::vector<double>& pdPreviousSpotGrid);

  /**
      Gets the path to save for the conditional variance swap.

      @param pdYGrid grid for the count within the corridor

      @return the path to save
   */
  size_t GetConditionalPathToSave(const std::vector<double>& pdYGrid);

  /**
      Constructs the events.

      @param bIsSimilarityReduction indicate whether or not
                 a similarity reduction can be used

      @return list of path dependent events
   */
  std::list< shared_ptr<PathDepEvent> > 
  ConstructEvents(bool bIsSimilarityReduction);

  /**
      Constructs the events for the fixed leg of conditional swap.

      @return list of path dependent events
   */
  std::list< shared_ptr<PathDepEvent> > ConstructConditionalEvents();

  /**
      Grid for the number of times wihtin the corridor.

      @return the grid
   */
  std::vector<double> ConstructConditionalGrid();

  /**
      Compute payoff if fixed leg of conditional swap.

      @return payoff of conditional leg
   */
  double ComputeConditionalPayoff();

  /**
      Clones this object.

      @return autoptr to cloned params
   */
  AutoPtr<VarianceSwapParams> Clone();  
  
  /**
      Used to indicate if the variance swap to be priced is a variance swap
      and has no cap.
   */
  bool IsSpecialVarianceSwap() const;

  /**
      Checks if the variance swap can be priced analytically.

      @return true if the vs can be priced analytically, false otherwise
   */
  bool IsAnalytical() const;

  /**
      Whether or not the variance swap is forward starting.

      @return true if forward starting, false otherwise
   */
  bool IsForwardStarting() const;

  // base class functions
  virtual void Init();

protected:

  /// Constructs the payoff
  void ConstructPayoff();

  /// Reference to underlying variance swap contract
  VarianceSwap& m_varianceSwap;

  /// If the object is cloned, then the clone needs to manage memory
  AutoPtr<VarianceSwap> m_clonedVarianceSwap;

  /// The average of the squared returns to date
  double m_dAvgSqrReturns;

  /// The spot price at the previous observation time
  double m_dPreviousSpot;

  /// The conditional variance swap payoff
  double m_dConditionalPayoff;

  /// Store the payoff
  shared_ptr<finance::Payoff> m_pPayoff;

private:

  NO_COPY_CLASS(VarianceSwapParams);

}; // class VarianceSwapParams


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_VARIANCESWAPPARAMS_H_
