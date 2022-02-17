/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/params.h
// Purpose:     base parameter class 
// Author:      ITO team
// Created:     2003/10/28
// RCS-ID:      $Id: params.h,v 1.83 2006/08/19 22:18:03 wang Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/params.h
    @brief base parameter class, including common parts for pricers with 
           different models.
 */

#ifndef _ITO33_PRICING_PARAMS_H_
#define _ITO33_PRICING_PARAMS_H_

#include "ito33/sharedptr.h"
#include "ito33/vector.h"

#include "ito33/pricing/contracts.h"
#include "ito33/pricing/eventmanager.h"

#include "ito33/numeric/interpolation.h"
#include "ito33/numeric/extrapolationmode.h"

namespace ito33 
{

namespace finance
{
  class Payoff;

  class ITO33_DLLDECL YieldCurve;

  class ITO33_DLLDECL Dividends;

  class ITO33_DLLDECL SessionData;
}

namespace numeric
{
  class NumParams;
  class MeshParams;

  namespace mesh
  {
    class SpecialTimes;
  }
}

namespace pricing
{

  class PathDepEvent;
  class Model;


/**
    Base parameter class, it doesn't contain information 
    related to space variable.

    The analysis time is by default set to negative, so invalid.
 */
class Params
{

public:
 
  /**
      Ctor by contracts. 
  
      Other members must be initialized by the SetXXX() functions

      @param contracts reference to Contracts
   */
  Params(Contracts& contracts); 

  /**
      Creates valid Params object.

      @param contracts reference to Contracts
      @param sessionData reference to financial Session
      @param pNumParams numeric parameters
      @param pMeshParams parameters for mesh builder
   */
  Params(Contracts& contracts,
         const finance::SessionData& sessionData,
         const shared_ptr<numeric::NumParams>& pNumParams,
         const shared_ptr<numeric::MeshParams>& pMeshParams);

  /// virtual dtor
  virtual ~Params() { } 

  /**
      Gets the current time.

      @return the current time
   */
  double GetCurrentTime() const
  {
    return m_dCurrentTime;
  }

  /**
      Sets yield curve.

      @param pYieldCurve yield curve
   */
  void SetYieldCurve(const shared_ptr<finance::YieldCurve>& pYieldCurve)
  {
    m_pYieldCurve = pYieldCurve;
  }

  /**
      Gets shared ptr to yield curve object.

      @return yield curve
   */
  const shared_ptr<finance::YieldCurve>& GetYieldCurve() const 
  { 
    return m_pYieldCurve; 
  }
   
  /**
      Sets yield curve for mesh.

      @param pYieldCurveForMesh yield curve for mesh
   */
  void SetYieldCurveForMesh(
    const shared_ptr<finance::YieldCurve>& pYieldCurveForMesh )
  {
    m_pYieldCurveForMesh = pYieldCurveForMesh;
  }

  /**
      Gets the yield curve for mesh.

      @return yield curve for mesh
   */
  const shared_ptr<finance::YieldCurve>& GetYieldCurveForMesh() const 
  { 
    return m_pYieldCurveForMesh; 
  }
  
  /**
      Sets foreign curve.

      @param pForeignCurve foreign curve    
   */
  void SetForeignCurve(const shared_ptr<finance::YieldCurve>& pForeignCurve)
  {
    m_pForeignCurve = pForeignCurve;
  }

  /**
      Gets shared ptr to foreign curve.

      @return foreign curve    
   */
  const shared_ptr<finance::YieldCurve>& GetForeignCurve() const 
  { 
    return m_pForeignCurve; 
  }

  /**
      Sets dividends.

      @param pDividends dividends
   */
  void SetDividends(const shared_ptr<finance::Dividends>& pDividends)
  {
    m_pDividends = pDividends;
  }

  /**
      Gets shared ptr to dividends.

      @return dividends
   */
  const shared_ptr<finance::Dividends>& GetDividends() const
  { 
    return m_pDividends;
  }

  /**
      Defines parameters for mesh builder.

      @param pMeshParams shared_ptr of the MeshParams object
   */
  void SetMeshParams(const shared_ptr<numeric::MeshParams>& pMeshParams)
  {
    m_pMeshParams = pMeshParams;
  }
 
  /**
      Gets mesh params.
  
      @return pointer to MeshParams
   */
  const numeric::MeshParams* GetMeshParams() const 
  { 
    return m_pMeshParams.get(); 
  }

  /**
      Defines the numerical parameters.

      @param pNumParams shared_ptr of the NumParams object
   */
  void SetNumParams(const shared_ptr<numeric::NumParams>& pNumParams)
  {
    m_pNumParams = pNumParams;
  }

  /**
      Gets NumParams.
  
      @return pointer to NumParams
   */
  const numeric::NumParams* GetNumParams() const 
  { 
    return m_pNumParams.get(); 
  }

  /**
      Sets valuation time.

      We don't save valuation date, since some problems such as 
      small convertible bond for Call Notice feature should be
      valuated at whatever time.

      @param dValuationTime 
   */
  void SetValuationTime(double dValuationTime)
  { 
    m_dValuationTime = dValuationTime; 
  }

  /**
      Gets valuation time.

      @return valuation time
   */
  double GetValuationTime() const { return m_dValuationTime; }

  /**
      Sets the stopping time.

      Stopping time is initialized to the maturity.  For problems that are
      solved in stages, this function must be called.  For example,
      solve from maturity to reset time using 2D pricer, then solve
      from reset time to valuation time using 1D solver. For the 
      second stage, this function is called with the reset time.

      @param dStoppingTime time to stop timestepping
   */
  void SetStoppingTime(double dStoppingTime)
  {
    m_dStoppingTime = dStoppingTime;
  }

  /**
      Gets stopping time.

      Usually, but not always, the same as the maturity time.

      @return the stopping time for timestepping
   */
  double GetStoppingTime() const 
  { 
    return m_dStoppingTime; 
  }

  /**
      Sets the analysis time. 

      @param dAnalysisTime the analysis time. 
   */
  void SetAnalysisTime(double dAnalysisTime);

  /**
      Gets the analysis time, negatif if the analysis date is not a valid date
      or is not between the pricing date and the maturity.

      @return the analysis date.
   */
  double GetAnalysisTime() const { return m_dAnalysisTime; }

  /**
      Sets spot share price where the scalar theoretical values are computed.

      @param dSpot initial spot
   */
  void SetSpotSharePrice(double dSpot) { m_dSpot = dSpot; }

  /**
      Gets spot share price.

      @return spot share price
   */
  double GetSpotSharePrice() const { return m_dSpot; } 

  /**
      Tells if there has any event(not only basic event) at the current time.

      Default implementation checks only the basic event.

      @return true if there has any event, false otherwise
   */
  virtual bool HasEventsNow() const 
  { 
    return m_eventManager.HasEventsNow(); 
  }

  /**
      Gets a pointer to a basic event.

      @return a pointer to a basic event if basic event exists, 0 otherwise
   */
  const Event* GetBasicEvent() const { return m_eventManager.GetEvent(); }

  /**
      Gets a pointer to a basic event but applied after constraints

      @return a pointer to a basic event if basic event exists, 0 otherwise
   */
  const Event* GetBasicEventAppliedAfterConstraints() const
  {
    return m_eventManager.GetEventAppliedAfterConstraints();
  }  

  /// @name helper functions for meshes construction
  //@{

  /**
      Gets special time points that we should put in time mesh

      Note that this is a virtual function. The default implementation requires
      that specific Contracts class implement also GetSpecialTimes() function.

      @param specialTimes (output) all special time. that means the initial
            elements in specialTimes, if any, are overwritten.
   */
  virtual void GetSpecialTimes(numeric::mesh::SpecialTimes& specialTimes) const;

  /**
      Calculates diffusion size to evaluation of the mesh domain size

      @param dSquaredTotalVol sqrt(TotalVolatility)
      @return diffusion size
   */
  double GetDiffusionSize(double dSquaredTotalVol) const;

  /**
      Calculates Convection size but only for a backward problem. the size
      is used to evaluation of the mesh domain size

      @param dSquaredTotalVol sqrt(TotalVolatility)
      @return Convection size for barkward problem
   */
  double GetBackwardConvectionSize(double dSquaredTotalVol) const;

  /**
      Calculates Convection size but only for a forward problem. the size
      is used to evaluation of the mesh domain size

      @param dSquaredTotalVol sqrt(TotalVolatility)
      @return Convection size for forward problem
   */
  double GetForwardConvectionSize(double dSquaredTotalVol) const;

  /**
      Generates a log mesh that can be used by some trivial instruments.

      @param model the pricing model
      @return the log space mesh
   */
  std::vector<double> GenerateSpaceMesh(const Model& model) const;

  //@}

  /// Intialize the params
  virtual void Init();

  /**
      Initializes state for given start time 

      @param dTime start time
   */
  virtual void SetInitialState(double dTime)
  {
    m_dCurrentTime = dTime;

    m_eventManager.SetInitialState(dTime);
  }

  /**
      Updates state for given time

      @param dTime given time where Params should be updated.
   */
  virtual void Update(double dTime)
  {
    m_dCurrentTime = dTime;

    m_eventManager.SetCurrentTime(dTime);
  }
  
  /// @name Path dependent functions
  //@{

  /**
      Gets the list of path dependent events
   */
  std::list< shared_ptr<PathDepEvent> > GetPathDepEvents()
  {
    return m_pathDepEvents; 
  }

  /**
      Sets the list of path dependent events
   */
  void SetPathDepEvents( std::list< shared_ptr<PathDepEvent> >& pathDepEvents )
  {
    m_pathDepEvents = pathDepEvents; 
  }
  
  //@}  
  
  /**
      For adjoint method in backward.

      @param pdSpots vector of spots
      @param pdMarketPrices vector of market prices
   */
  void SetObservations
  (const std::vector<double>& pdSpots, 
   const std::vector<double>& pdMarketPrices)
  {
    m_pdSpots = pdSpots;
    m_pdMarketPrices = pdMarketPrices;
  }

  /**
      Gets spots.

      @return spot vector
   */
  const std::vector<double>& GetSpots() const { return m_pdSpots; }

  /**
      Gets market prices.

      @return prices vector
   */
  const std::vector<double>& GetMarketPrices() const
  { 
    return m_pdMarketPrices; 
  }

protected:

  /// @name Input 
  //@{

  /// reference to Contracts
  Contracts& m_contracts;

  /// shared ptr to financial Yield Curve
  shared_ptr<finance::YieldCurve> m_pYieldCurve;

  /// shared ptr to financial Yield Curve for mesh
  shared_ptr<finance::YieldCurve> m_pYieldCurveForMesh;

  /// shared ptr to foreign curve
  shared_ptr<finance::YieldCurve> m_pForeignCurve;

  /// shared ptr to dividends
  shared_ptr<finance::Dividends> m_pDividends;

  /// shared ptr to NumParams
  shared_ptr<numeric::NumParams> m_pNumParams; 
  
  /// Shared ptr to MeshParams
  shared_ptr<numeric::MeshParams> m_pMeshParams; 

  /// valuation time
  double m_dValuationTime;

  /// analysis time.
  double m_dAnalysisTime;

  /**
      Stopping time for timestepping.  Usually the maturity time.  Can be 
      different if the problem is solved in different (time) stages. For 
      example, solve from maturity to reset time using a 2D pricer, then
      solve from reset time to valuation time using a 1D pricer. Cannot
      simply modify the maturity time of the contract since this can be
      used for other calculations (such as claim).
   */
  double m_dStoppingTime;

  /// initial spot
  double m_dSpot;

  //@}

  /// The current time, indicates the advancing of the pricing
  double m_dCurrentTime;

  /// A helper for event management
  EventManager m_eventManager;

  /// Helper function for backward dividend events construction
  void ConstructDividendEvents(
    numeric::ExtrapolationMode emLeft = numeric::ExtrapolationMode_Linear,
    numeric::ExtrapolationMode emRight = numeric::ExtrapolationMode_Linear,
    numeric::InterpolationMethod interpolationMethod 
                                  = numeric::InterpolationMethod_Quadratic);

  /// Helper function for forward dividend events construction.
  void ConstructForwardDividendEvents(
    numeric::ExtrapolationMode emLeft = numeric::ExtrapolationMode_Linear,
    numeric::ExtrapolationMode emRight = numeric::ExtrapolationMode_Linear,
    numeric::InterpolationMethod interpolationMethod 
                                  = numeric::InterpolationMethod_Quadratic);

  /// Path dependent events
  std::list< shared_ptr<PathDepEvent> > m_pathDepEvents;
  
  /// The spots for which we'll compute prices/sensitivities at valuation date
  std::vector<double> m_pdSpots;

  /// The observed market values at the above spots
  std::vector<double> m_pdMarketPrices;

private:

  NO_COPY_CLASS(Params);

}; // class Params


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_PARAMS_H_
