/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/bondliketerms.h
// Purpose:     base Characteristics class for Bond and mandatory
// Created:     2004/08/13
// RCS-ID:      $Id: bondliketerms.h,v 1.34 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/bondlike/bondliketerms.h
    @brief declaration of the BondLikeTerms class
 */

#ifndef _ITO33_FINANCE_BONDLIKE_BONDLIKETERMS_H_
#define _ITO33_FINANCE_BONDLIKE_BONDLIKETERMS_H_

#include "ito33/sharedptr.h"
#include "ito33/date.h"

#include "ito33/dlldecl.h"

namespace ito33
{

namespace XML
{
  class Tag;
}

namespace finance
{

class ITO33_DLLDECL CashFlowStream;


/**
   Terms of the instument, that are common to a bond-like.  
 */
class ITO33_DLLDECL BondLikeTerms
{
public:
  
  /**
     ctor.

     @param issueDate the issue date
     @param dIssuePrice the isssue price, as a percentage of nominal. 
      (100% if no OID and no premium redemption at maturity) 
     @param maturityDate the maturity date
     @param dNominal the nominal of the security
     @param dRecoveryRate the recovery rate of the bond-like in case of default
   */
  BondLikeTerms(Date issueDate, double dIssuePrice, 
                Date maturityDate, double dNominal,
                double dRecoveryRate);

  /// virtual dtor for base class
  virtual ~BondLikeTerms() { }

  /**
     The cash distribution of the instrument (coupons or prefered dividends)
     expressed as percentage rates of nominal, not required for a zero coupon
     bond.

     @param pCashDistribution cash distribution of the instrument.
   */
  virtual void 
  SetCashDistribution(const shared_ptr<CashFlowStream>& pCashDistribution);

  /**
     Gets the issue date of the bondlike security.

     @return issue date
   */
  Date GetIssueDate() const { return m_issueDate; }

  /**
     Gets the issue price of the bondlike security, expressed as a percentage
     of the nominal.

     @return the issue price expressed as a percentage of the nominal
   */
  double GetIssuePrice() const { return m_dIssuePrice; }

  /**
     Gets the maturity date of the bondlike security.
     @return maturity date
   */
  Date GetMaturityDate() const  { return m_maturityDate; }

  /**
     Gets the nominal value of the bondlike security.

     @return nominal
   */
  double GetNominal() const { return m_dNominal; }

  /**
     Gets the recovery rate of the bondlike security. The recovery rate ranges
     from 0 to 1.

     @return the recovery rate of the security
   */
  double GetRecoveryRate() const { return m_dRecoveryRate; }

  /**
     The cash distribution of the security (coupons or prefered dividends)
     expressed as percentage rates of nominal, not required for a zero coupon
     bond.
    
     @return Cash distribution of the security
   */
  const shared_ptr<CashFlowStream>& GetCashDistribution() const 
  {
    return m_pCashDistribution;
  }  
   
  /**
     @internal
     @brief Dump to xml file.

     @param tagParent xml node
     
     @noexport
   */
  virtual XML::Tag Dump(ito33::XML::Tag& tagParent) const;
  
protected:

  /**
     Writes myself to tag root of specific bondliketerms.

     @param tagParent BondTerms etc.
   */
  void DumpMe(ito33::XML::Tag& tagParent) const;


  /// issue date
  Date m_issueDate;
  
  /// issue price
  double m_dIssuePrice;

  /// maturity date
  Date m_maturityDate;

  /// the nominal of the security
  double m_dNominal;

  /// the recovery rate of the security
  double m_dRecoveryRate;

  /// cash distribution of the security (coupons or prefered dividends)
  shared_ptr<CashFlowStream> m_pCashDistribution;

}; // class BondLikeTerms


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_BONDLIKE_BONDLIKETERMS_H_
