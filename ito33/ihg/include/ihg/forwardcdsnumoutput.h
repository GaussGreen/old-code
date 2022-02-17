/////////////////////////////////////////////////////////////////////////////
// Name:      ihg/forwardcdsnumoutput.h
// Purpose:   NumOutput class for CDS contracts with forward PDE
// Author:    David
// RCS-ID:    $Id: forwardcdsnumoutput.h,v 1.13 2006/08/20 09:36:16 wang Exp $
// Copyright: (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/forwardcdsnumoutput.h
    @brief NumOutput class for CDS contracts with forward PDE
 */

#ifndef _IHG_FORWARDCDSNUMOUTPUT_H_
#define _IHG_FORWARDCDSNUMOUTPUT_H_

#include "ito33/autoptr.h"
#include "ito33/vector.h"
#include "ito33/list.h"

namespace ito33
{

namespace pricing
{
  class ForwardCDSParams;
}

namespace ihg
{

  class ForwardCDSInstData;
  class ModelOutput;

/**
    This class stores all pricing information for forward CDS contracts, 
    and is responsible for constructing the model output class returned
    to the user.
 */
class ForwardCDSNumOutput
{
public:
  
  ForwardCDSNumOutput(pricing::ForwardCDSParams& params) 
    : m_params(params)
  { 
  }

  // Default dtor is ok

  /**
     Return the output structure containing the price information

     @return the model output class containing requested price information
   */
  shared_ptr<ModelOutput> GetModelOutput();

  /**
     Initialize the class variables

     instdata::SetInitialState must have been called
   */
  void Init(ForwardCDSInstData& instdata);

  /**
     If required, strore the pricing data at the current timestep
   */
  void UpdateMe(ForwardCDSInstData& instdata, double dTime);

  /**
     Finalize the pricing data.

     Called after pricing is complete.  If required, any pointers to external 
     data objects (such as grids, data at the last timestep) must be copied.
   */
  void Finalize(ForwardCDSInstData& instdata);
 
  double GetPrice()
  {
    return m_CDSPrices.back();
  }

  std::vector<double> GetRecoveryTerms()
  {
    return m_pdRecoveryTerms;
  }

  std::vector<double> GetAccruedTerms()
  {
    return m_pdAccruedTerms;
  }

  std::vector<double> GetSpreadTerms()
  {
    return m_pdSpreadTerms;
  }

  std::vector<double> GetTimes()
  {
    return m_pdTimes;
  }

  std::vector<double> GetPrices()
  {
    return m_pdCDSPrices;
  }

protected:

  /// The params of the PDE
  pricing::ForwardCDSParams& m_params;

  /// The CDS prices at each timestep
  std::list<double> m_CDSPrices;
  std::vector<double> m_pdCDSPrices;
 
  /// The recovery part of the CDS price
  std::list<double> m_RecoveryTerms;
  std::vector<double> m_pdRecoveryTerms;

  /// The spread part of the CDS price
  std::list<double> m_SpreadTerms;
  std::vector<double> m_pdSpreadTerms;

  /// The accrued part of the CDS price
  std::list<double> m_AccruedTerms;
  std::vector<double> m_pdAccruedTerms;

  /// The time mesh
  std::list<double> m_Times;
  std::vector<double> m_pdTimes;

private:

  NO_COPY_CLASS(ForwardCDSNumOutput);

}; // class ForwardCDSNumOutput


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_FORWARDCDSNUMOUTPUT_H_
