/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/cds.h
// Purpose:     financial cds class 
// Author:      Nabil
// Created:     2003/09/15
// RCS-ID:      $Id: cds.h,v 1.61 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 1999-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/cds.h
    @brief declaration of the financial cds class.
 */

#ifndef _ITO33_FINANCE_CDS_H_
#define _ITO33_FINANCE_CDS_H_

#include "ito33/common.h"

#include "ito33/finance/cdslike.h"

namespace ito33
{

namespace finance
{

  class ITO33_DLLDECL CashFlowStreamUniform;
  class DerivateVisitor;
  class DerivativeModifyingVisitor;

/**
   This class describes a Credit Default Swap (CDS).
 
   A Credit Default Swap (CDS) is an OTC contract that allows the buyer
   to cover itself against the risk of default of some company (usually 
   the issuer of a bond).
 */
class ITO33_DLLDECL CDS : public CDSLike
{

public:

  /**
     Ctor constructs a cds from its recovery rate and spread stream.
     
     The market price of the cds is not set.

     @param dRecoveryRate the recovery rate of the CDS
     @param pSpreadStream the spread stream of the CDS
   */
  CDS(double dRecoveryRate,
      const shared_ptr<CashFlowStreamUniform>& pSpreadStream);

  // Default dtor is ok. CDS won't be derived

  /**
      Gets the maturity date of the CDS contract which is the last
      payment date of the spread stream.

      @return maturity date of the CDS contract.
   */
  virtual Date GetMaturityDate() const;

  // implement pure virtual functions of CDSLike
  virtual double GetSpread() const;
  virtual shared_ptr<CashFlowStreamUniform> GetSpreadStream() const
  { 
    return m_pSpreadStream; 
  }
  
  // implement base class pure virtuals
  void Visit(DerivativeVisitor& visitor) const;
  void Visit(DerivativeModifyingVisitor& visitor);
  XML::Tag Dump(XML::Tag& tagParent) const;

}; // class CDS


} //namespace finance

} //namespace ito33

#endif // #ifndef _ITO33_FINANCE_CDS_H_
