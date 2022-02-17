/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/cboptionparams.h
// Purpose:     cb option params class
// Author:      Nabil
// Created:     2005/10/13
// RCS-ID:      $Id: cboptionparams.h,v 1.4 2006/08/19 22:18:03 wang Exp $
// Copyright:   (c) 2005 - 2006  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/cboptionparams.h
    @brief cb option params class

    Implementation of the params class for CBOption.
 */

#ifndef _ITO33_PRICING_CBOPTIONPARAMS_H_
#define _ITO33_PRICING_CBOPTIONPARAMS_H_

#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"
#include "ito33/list.h"
#include "ito33/date.h"

#include "ito33/pricing/cboption.h"
#include "ito33/pricing/cbparams.h"


namespace ito33 
{

namespace pricing
{


/// Pricing parameters for Convertible Bond Options
class CBOptionParams : public CBParams
{

public: 

  /**
      Ctor by CBOption contract and common financial and numerical datas. 

      @param cboption reference to CBOption contract objet
      @param sessionData reference to financial SessionData
      @param pNumParams numeric parameters
      @param pMeshParams parameters for mesh builder
   */
  CBOptionParams(CBOption& cboption,
                 const finance::SessionData& sessionData,
                 const shared_ptr<numeric::NumParams>& pNumParams,
                 const shared_ptr<numeric::MeshParams>& pMeshParams)
    : CBParams(cboption.GetCB(), sessionData, pNumParams, pMeshParams),
      m_cboption(cboption)
  {
  }

  /// virtual dtor for base class
  virtual ~CBOptionParams() {}

  /// @name implementation of virtual functions
  //@{

  virtual void Init();
  
  /**
      Updates time at maturity for different members (time, eventmanager etc).

      @param dTime should be the time at maturity
   */
  virtual void SetInitialState(double dTime);
  
  /**
      Updates the index of the mono date event.

      @param pCBEvent the event to treat
      @param bConstraintsNeedUpdate true if update needed    
   */
  virtual void MonoDateEventIndexTreatment
               (const CBEvent* pCBEvent, bool& bConstraintsNeedUpdate);
  
  virtual void Update(double dTime);

  //@}

  /// @name Some functions specific to the cb option
  //@{

  /**
      Updates (if necessary) the data for the cb option.

      @param dTime given time
      @param bPlus t+ or t-?
   */
  void UpdateCBOptionData(double dTime, bool bPlus);

  /// updates cb option data at current time
  void UpdateCBOptionData();

  bool InCBOptionWindow() const { return m_bInCBOptionWindow; }

  /**
      Gets the CBOptionData in the CB contract.

      @return the CBOptionData in the CB contract
   */
  CBOptionData* GetCBOptionData() const 
  { 
    return m_cboption.GetCBOptionData(); 
  }

  /**
      Gets the strike of the cb option at the current time.

      @return the strike of the cb option at the current time.
   */
  double GetCBOptionStrike() const { return m_dCBOptionStrike; }
  
  /**
      Computes the NPV of the cash flows passed as argument.

      @param cashflows the cash flow stream
      @param dTime the time
      @param nNextCashFlow the index of the next cash flow
   */
  double ComputeNPV(const CashFlows& cashflows, double dTime, 
                    size_t nNextCashFlow) const;
  
  //@}

  /// Gets a reference to the CB Option contract
  CBOption& GetCBOption() const { return m_cboption; }  

protected:
  
  CBOption& m_cboption;
  
  ///@name cb option related data
  //@{

  /// Next floating coupon index. Used for the cb option to update the strike.
  size_t m_nIdxASWNextFloatingCoupon;

  /// slope floating coupon
  double m_dASWSlopeFloatingCoupon;

  /** 
      Difference between NPV of floating legs and that of fixed legs 
      at current time (i.e that is NPVFixed - NPVFloating).
   */
  double m_dASWNPVDiff;

  /// Strike of the cb option
  double m_dCBOptionStrike;

  /// Indicates if the current time is in a cb option window
  bool m_bInCBOptionWindow;

  //@}

private:

  NO_COPY_CLASS(CBOptionParams);

}; // class CBOptionParams


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_CBOPTIONPARAMS_H_
