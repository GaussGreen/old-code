/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/derivativevisitors/cdsvisitor.h
// Purpose:     Visitor for cds
// Author:      Vadim Zeitlin
// Created:     2004-05-11
// RCS-ID:      $Id: cdsvisitor.h,v 1.5 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/derivativevisitors/cdsvisitor.h
    @brief Visitor class accepting different Derivative kinds.

    This is the central part of the implementation of the visitor pattern for
    the derivatives. To do something different depending on the exact type of
    a Derivative object you have to define a new class inheriting from this one
    and do whatever is required in its methods.
 */

#ifndef _ITO33_FINANCE_DERIVATIVEVISITORS_CDSVISITOR_H_
#define _ITO33_FINANCE_DERIVATIVEVISITORS_CDSVISITOR_H_

#include "ito33/finance/derivative_visitor.h"
#include "ito33/finance/cds.h"

#include "ito33/dlldecl.h"

namespace ito33
{

namespace finance
{

/**
    CDS visitor 

    An object of the class must be passed by the user code
    to any functions which may access heterogeneous collections of Derivatives,
    e.g. XML::Reader in IHG currently.
 */
class ITO33_DLLDECL CDSVisitor : public DerivativeVisitor
{
public:
  
  /// Gets visited CDS, otherwise returns 0
  shared_ptr<CDS> GetCDS() { return m_pCDS; }
  
  /// Called for a CDS
  virtual void OnCDS(const CDS& contract)  
  { 
    m_pCDS = shared_ptr<CDS>(new CDS(contract));
  } 

private:

  shared_ptr<CDS> m_pCDS;
};

} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_DERIVATIVEVISITORS_CDSVISITOR_H_
