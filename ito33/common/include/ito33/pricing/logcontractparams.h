/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/logcontractparams.h
// Purpose:     params class for log contract
// Created:     2006/07/18
// RCS-ID:      $Id: logcontractparams.h,v 1.2 2006/08/19 22:18:03 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/logcontractparams.h
    @brief params class for log contract
 */

#ifndef _ITO33_PRICING_LOGCONTRACTPARAMS_H_
#define _ITO33_PRICING_LOGCONTRACTPARAMS_H_

#include "ito33/pricing/logcontract.h"
#include "ito33/pricing/params.h"

namespace ito33 
{
 
namespace finance
{
  class ITO33_DLLDECL SessionData;
}

namespace pricing
{

class LogContract;

/// params class for log contract
class LogContractParams : public Params
{

public: 

  /**
      Ctor takes a pricing log contract objet. 

      @param logContract reference to log contract objet
   */
  LogContractParams(LogContract& logContract) 
                  : Params(logContract), m_logContract(logContract)
  { 
  }

  /**
      Ctor takes a pricing log contract objet and common financial and
      numerical datas. 

      @param logContract reference to log contract objet
      @param sessionData reference to financial session data
      @param pNumParams numeric parameters
      @param pMeshParams parameters for mesh builder
   */
  LogContractParams(LogContract& logContract,
                    const finance::SessionData& sessionData,
                    const shared_ptr<numeric::NumParams>& pNumParams,
                    const shared_ptr<numeric::MeshParams>& pMeshParams)
                  : Params(logContract, sessionData, pNumParams, pMeshParams),
                    m_logContract(logContract)
  {
  }

  virtual void Init();

  /// Get a reference to the log contract
  LogContract& GetLogContract() const 
  { 
    return m_logContract; 
  }


private:

  /// The pricing log contract of this params
  LogContract& m_logContract;

  NO_COPY_CLASS(LogContractParams);

}; // class LogContractParams


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_LOGCONTRACTPARAMS_H_
