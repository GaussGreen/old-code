 /////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/generalizedpepslike.h
// Purpose:     PEPS-like financial class
// Created:     2004/08/17
// RCS-ID:      $Id: generalizedpepslike.h,v 1.9 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/bondlike/generalizedpepslike.h
    @brief declaration of the financial PEPS-like class   
 */

#ifndef _ITO33_FINANCE_BONDLIKE_GENERALIZED_PEPSLIKE_H_
#define _ITO33_FINANCE_BONDLIKE_GENERALIZED_PEPSLIKE_H_

#include "ito33/dlldecl.h"

#include "ito33/finance/bondlike/convertiblelike.h"
#include "ito33/finance/bondlike/generalizedpepslikecall.h"
#include "ito33/finance/bondlike/pepsaveragingperiod.h"

namespace ito33
{

namespace finance
{

class ITO33_DLLDECL BondLikeTerms;
class ITO33_DLLDECL CallSchedule;

/**
   The class describing the generalized PEPS structure.
 */
class ITO33_DLLDECL GeneralizedPEPSLike : public ConvertibleLike
{
public:
  
  /**
     Creates a PERCS-like by BondLikeTerms and GeneralizedPEPSLike payoff.

     The payoff of the generalized PEPS is defined as the following
     - When underlying share price is less than "lower strike", the conversion
       ratio is "downside conversion ratio".
     - When underlying share price S is greater than "lower strike" and less
       than "higher strike", the conversion ratio is ResetLevel/S where "Reset
       level" equal to DownsideConversionRatio * LowerStrike
     - When underlying share price S is greater than "higher strike", the 
       "upside conversion ratio" UpsideBaseConversionRatio - Cash / S. 
       In particular, we have 
       ResetLevel = UpsideBaseConversionRatio - Cash / HihgerStrike.

     @param pBondLikeTerms the terms of the PEPS that are common 
                           to a the bond terms.
     @param dDownsideConversionRatio the downside conversion ratio
     @param dLowerStrike lower strike before which downside conversion ratio
                          is applied.
     @param dHigherStrike higher strike
     @param dUpsideBaseConversionRatio upside base conversion ratio
   */
  GeneralizedPEPSLike(const shared_ptr<BondLikeTerms>& pBondLikeTerms,
                      double dDownsideConversionRatio,
                      double dLowerStrike,
                      double dUpsideBaseConversionRatio,
                      double dHigherStrike);
  
  /**
     Specifies that the security has optional conversion.
   */
  void EnableOptionalConversion()
  {
    m_bHasOptionalConversion = true;
  }

  /**
     The fixed cash call of PEPS-like instruments.

     @param pCallSchedule shared pointer to the fixed cash call.
     
     @method
   */
  void SetCallFixedCash(const shared_ptr<CallSchedule>& pCallSchedule);

  /**
     The Averaging period.

     @param pAveragingPeriod shared pointer to the averaging Period.
   */
  void 
  SetAveragingPeriod(const shared_ptr<PEPSAveragingPeriod>& pAveragingPeriod);

  /**
     The call provision of PEPS-like instruments whose value is share dependent.

     Note that this is indeed the optional conversion at issuer's option.

     @param pGeneralizedPEPSLikeCall shared pointer to GeneralizedPEPSLikeCall

     @method
   */
  void SetGeneralizedPEPSLikeCall
    (const shared_ptr<GeneralizedPEPSLikeCall>& pGeneralizedPEPSLikeCall);

  /**
      @name Methods for accessing the PEPS like instrument.
   */
  //@{

  /**
     Gets the downside conversion ratio.

     @return the downside conversion ratio
   */
  double GetDownsideConversionRatio() const
  {
    return m_dDownsideConversionRatio;
  }

  /**
     Gets the lower strike.

     @return the lower strike
   */
  double GetLowerStrike() const
  {
    return m_dLowerStrike;
  }

  /**
     Gets the higher strike.

     @return the higher strike
   */
  double GetHigherStrike() const
  {
    return m_dHigherStrike;
  }

  /**
     Gets the upside base conversion ratio.

     @return the upside base conversion ratio
   */
  double GetUpsideBaseConversionRatio() const
  {
    return m_dUpsideBaseConversionRatio;
  }

  /**
     Gets the reset level.

     @return the reset level
   */
  double GetResetLevel() const
  {
    return m_dDownsideConversionRatio * m_dLowerStrike;
  }

  /**
     Gets the cash adjustment.
    
     @return the cash adjustment
   */
  double GetCashAdjustment() const
  {
    return m_dUpsideBaseConversionRatio * m_dHigherStrike - GetResetLevel();
  }

  /**
     Gets the minimum conversion ratio.

     @return the minimum conversion ratio
   */
  double GetMinConversionRatio() const
  {
    return GetResetLevel() / m_dHigherStrike;
  }

  /**
     Specifies whether the PEPS-like security has optional conversion (fixed
     share).

     @return true if optional conversion at holder's option exists.
   */
  bool HasOptionalConversion() const
  {
    return m_bHasOptionalConversion;
  }

  
  /**
     Indicates whether the mandatory conversion involves an averaging period.

     @return true if an averaging period has been specified
   */
 bool HasAveragingPeriod() const 
  {
    if (!m_pAveragingPeriod)
      return false;

    return true;
  }

  /**
     Gets the averaging period if any, a null pointer otherwise.

     @return the averaging period
   */
  const shared_ptr<PEPSAveragingPeriod>& GetAveragingPeriod() const
  {
    return m_pAveragingPeriod;
  }

  /**
     The fixed cash call of PEPS-like instruments.
    
     @return call provision
   */
  const shared_ptr<CallSchedule>& GetCallFixedCash() const
  {
    return m_pCallSchedule; 
  }

  /**
     The optional conversion at issuer's option.

     @return optional conversion at issuer's option.
   */
  const shared_ptr<GeneralizedPEPSLikeCall>& GetGeneralizedPEPSLikeCall() const
  {
    return m_pGeneralizedPEPSLikeCall;
  }

  //@}

  void ValidateWith(const SessionData& sessionData) const;

  void Visit(DerivativeVisitor& visitor) const;
  
  void Visit(DerivativeModifyingVisitor& visitor);

  virtual XML::Tag Dump(XML::Tag& tagParent) const;


private:
    
  /// the downside conversion ratio
  double m_dDownsideConversionRatio;

  /// lower strike before which downside conversion ratio is applied
  double m_dLowerStrike;

  /// upside base conversion ratio
  double m_dUpsideBaseConversionRatio;

  /// higher strike
  double m_dHigherStrike;

  /// Is the optional conversion enabled?
  bool m_bHasOptionalConversion;

  /// fixed cash call provision
  shared_ptr<CallSchedule> m_pCallSchedule;

  /// fixed cash call provision
  shared_ptr<PEPSAveragingPeriod> m_pAveragingPeriod;

  /// fixed share call provision
  shared_ptr<GeneralizedPEPSLikeCall> m_pGeneralizedPEPSLikeCall;

  /// whether or not the call is set
  bool m_bCallSet;
}; // class PEPSLike


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_BONDLIKE_GENERALIZED_PEPSLIKE_H_
