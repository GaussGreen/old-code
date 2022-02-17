/////////////////////////////////////////////////////////////////////////////
// Name:      hg/numoutput.h
// Purpose:   implementation of HG base NumOutput class 
// Created:   2005/01/13 
// RCS-ID:    $Id: numoutput.h,v 1.6 2006/08/18 12:42:46 wang Exp $
// Copyright: (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file hg/numoutput.h
   @brief implementation of HG base NumOutput class 
 */

#ifndef _HG_NUMOUTPUT_H_
#define _HG_NUMOUTPUT_H_

#include "ito33/finance/numoutput.h"
#include "ito33/finance/computationalflags.h"

namespace ito33
{

namespace numeric
{
  class SurfaceGeneral;
}

namespace hg
{


/**
   HG Base numercial output class, all greeks go into derived class 
   BackwardNumOutput. 
 */  
class NumOutput : public finance::NumOutput
{
public:
  
  /// default ctor, initializes the final save flag to false
  NumOutput() : m_bFinalSave(false),
                m_dPrice(0), m_dDelta(0), m_dGamma(0), m_dTheta(0) {}

  virtual ~NumOutput() { }
 
  /*
    Note, there is no SetComputationalFlags function. Please use the accessor 
    below to setup individually the flags. SetComputationalFlags is dangerous
    since it can take a given flags objet which might have actived too many 
    flags that we don't need.
  */
 
  /// @name Accessor methods 
  //@{

  /**
     Get the reference of the internal computational flags.

     @return the refrence of the computational flags
   */
  finance::ComputationalFlags& GetComputationalFlags()
  {
    return m_computationalFlags;
  }
 
  /**
     Get the price at initial spot.

     @return the price at initial spot
   */
  double GetPrice() const { return m_dPrice; }

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

     Should be called only when GetFinalSave returns true
   */
  const std::vector<double>& GetFinalPrices() const
  {
    return m_pdFinalValues;
  }

  /// Gets the final computational space mesh(the space mesh at start time)     
  const std::vector<double>& GetFinalMesh() const
  {
    return m_pdFinalSpots;
  }
  
  /// Gets the sensitivities.
  const std::vector<double>& GetSensitivities()
  {
    return m_pdSensitivities;
  }

  /// Checks whether or not the sensitivities have been set.
  bool HasSensitivities() const 
  { 
    return !m_pdSensitivities.empty(); 
  }

  //@}  // name Accessor methods 

  /// @name Modifier methods 
  //@{

  void SetSensitivities(const std::vector<double>& pdSensitivities) 
  { 
    m_pdSensitivities = pdSensitivities; 
  }

  /**
     Sets the FinalSave flag.

     @param bFlag value of FinalSave flag 
   */
  void SetFinalSave(bool bFlag = true)
  {
    m_bFinalSave = bFlag;
  }  

  //@}  // name Modifier methods 

protected:

  /**
     The computational flags, although flags on greeks only make sense for
     backward problem, the flag on surface (or even flag on analysis date)
     does make sense for forward problem. 
   */
  finance::ComputationalFlags m_computationalFlags;

  /// If we need to save the final results
  bool m_bFinalSave;

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

  /// Scalar values
  std::vector<double> m_pdSensitivities;


private:

  // Explicitly forbid copy since copying a numoutput doesn't make sense
  NO_COPY_CLASS(NumOutput);

}; // class NumOutput


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_NUMOUTPUT_H_
