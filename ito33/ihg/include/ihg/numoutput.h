/////////////////////////////////////////////////////////////////////////////
// Name:      ihg/numoutput.h
// Purpose:   implementation of base NumOutput class 
// Author:    ZHANG Yunzhi
// RCS-ID:    $Id: numoutput.h,v 1.22 2006/03/24 10:18:27 pedro Exp $
// Copyright: (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/numoutput.h
    @brief base ihg numericial output class

    Implementation of the ihg::NumOutput class.
 */

#ifndef _IHG_NUMOUTPUT_H_
#define _IHG_NUMOUTPUT_H_

#include "ito33/vector.h"

#include "ito33/finance/numoutput.h"

#include "ito33/finance/computationalflags.h"

namespace ito33
{

namespace numeric
{
  class SurfaceGeneral;
}

namespace ihg
{


/**
    Base numercial output class, all greeks go into derived class 
    BackwardNumOutput. 
 */
class NumOutput : public finance::NumOutput
{
public:
  
  /// default ctor, initializes the final save flag to false
  NumOutput() : m_bFinalSave(false), m_bHasImpliedBrownianVol(false)
  { 
  }

  virtual ~NumOutput() { }

  /// @name Modifier methods
  //@{

  /**
      Sets the FinalSave flag.

      @param bFlag value of FinalSave flag 
   */
  void SetFinalSave(bool bFlag = true)
  {
    m_bFinalSave = bFlag;
  }
 
  /**                                                                          
      Sets the implied brownian volatility.

      @param dImpliedBrownianVol implied brownian volatility  
   */                                                                          
  void SetImpliedBrownianVol(double dImpliedBrownianVol)                       
  {        
    m_bHasImpliedBrownianVol = true;
    m_dImpliedBrownianVol = dImpliedBrownianVol;                               
  }     

  //@}

  /*
    Note, there is no SetComputationalFlags function. Please use the getter 
    below to setup individually the flags. SetComputationalFlags is dangerous
    since it can take a given flags objet which might have actived too many 
    flags that we don't need.
  */
 
  /// @name Accessors methods 
  //@{

  /**
      Gets the reference of the internal computational flags.

      @return the refrence of the computational flags
   */
  finance::ComputationalFlags& GetComputationalFlags()
  {
    return m_computationalFlags;
  }
 
  /**
      Gets the price at initial spot.

      @return the price at initial spot
   */
  double GetPrice() const { return m_dPrice; }
 
  /**
      Gets the value after default at initial spot.

      @return the value after default at initial spot
   */
  double GetValueAfterDefault() const { return m_dValueAfterDefault; }

  /**
      Gets the FinalSave flag.
    
     @return the FinalSave flag
   */
  bool GetFinalSave() const
  {
    return m_bFinalSave;
  }

  /** 
      Gets the final prices.

      Should be called only when GetFinalSave returns true.
   */
  const std::vector<double>& GetFinalPrices() const
  {
    return m_pdFinalValues;
  }

  /** 
      Gets the final computational space mesh(the space mesh at start time).     
   */
  const std::vector<double>& GetFinalMesh() const
  {
    return m_pdFinalSpots;
  }

  /**
      Gets the implied brownian volatility.
   */
  double GetImpliedBrownianVol();

  bool HasImpliedBrownianVol()
  {
    return m_bHasImpliedBrownianVol;
  }

  //@}  // name Accessor methods 

protected:

  /**
      The computational flags, although flags on greeks only make sense for
      backward problem, the flag on surface (or even flag on analysis date)
      does make sense for forward problem. 
   */
  finance::ComputationalFlags m_computationalFlags;

  /// If we need to save the final results
  bool m_bFinalSave;

  /// indicate if the implied brownian volatility has been computed
  bool m_bHasImpliedBrownianVol;

  ///implied brownian volatility
  double m_dImpliedBrownianVol;

  /// The final space mesh
  std::vector<double> m_pdFinalSpots;

  /// The final price array (presumably at valuation date) 
  std::vector<double> m_pdFinalValues;

  /// price on spot at final time
  double m_dPrice;

  /// delta on spot at final time
  double m_dDelta;

  /// gamma on spot at final time
  double m_dGamma;

  /// theta on spot at final time
  double m_dTheta;

  /// value after default on spot at final time
  double m_dValueAfterDefault;

private:

  // Explicitly forbid copy since copying a numoutput doesn't make sense
  NO_COPY_CLASS(NumOutput);

}; // class NumOutput


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_NUMOUTPUT_H_
