/////////////////////////////////////////////////////////////////////////////
// Name:        hg/logcontract_closedform.h
// Purpose:     HG log contract pricer(analytical)
// Created:     2006/07/18
// RCS-ID:      $Id: logcontract_closedform.h,v 1.1 2006/07/19 17:39:51 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/logcontract_closedform.h
    @brief HG log contract pricer class using closed form
 */

#ifndef _HG_LOGCONTRACT_CLOSEDFORM_H_
#define _HG_LOGCONTRACT_CLOSEDFORM_H_

#include "ito33/common.h"
#include "ito33/autoptr.h"

namespace ito33
{

namespace finance 
{ 
  class ITO33_DLLDECL ComputationalFlags; 
}

namespace pricing
{
  class LogContractParams;
}

namespace hg
{

  class Model;
  class NumOutputTimeOnly;

/**
    Log contract pricer.

    Log contract prices can be computed analytically.
 */
class LogContractClosedForm
{
public:

  /**
      Constructor.

      @param params reference to variance swap pricing object
      @param model the HG model to use
      @param flags the computational flags specifying what to compute
   */
  LogContractClosedForm(pricing::LogContractParams& params, 
                        Model& model, 
                        const finance::ComputationalFlags& flags);
  
  /** 
      Prices the log contract.
      
      @return object containing all price data
   */
  AutoPtr<NumOutputTimeOnly> Price();

protected:

  /// The variance swap params being used for pricing
  pricing::LogContractParams& m_params;

  /// The HG model parameters
  Model& m_model;

  /// The computational flags (ie. what data to compute)
  const finance::ComputationalFlags& m_flags;
  
private:
  
  NO_COPY_CLASS(LogContractClosedForm);

}; // class LogContractClosedForm


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_LOGCONTRACT_CLOSEDFORM_H_
