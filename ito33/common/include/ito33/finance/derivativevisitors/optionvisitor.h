/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/derivativevisitors/optionvisitor.h
// Purpose:     Visitor for Derivative-derived classes
// Author:      Vadim Zeitlin
// Created:     2004-05-11
// RCS-ID:      $Id: optionvisitor.h,v 1.5 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/derivativevisitors/optionvisitor.h
    @brief Visitor class accepting different Derivative kinds.

    This is the central part of the implementation of the visitor pattern for
    the derivatives. To do something different depending on the exact type of
    a Derivative object you have to define a new class inheriting from this one
    and do whatever is required in its methods.
 */

#ifndef _ITO33_FINANCE_DERIVATIVEVISITORS_OPTIONVISITOR_H_
#define _ITO33_FINANCE_DERIVATIVEVISITORS_OPTIONVISITOR_H_

#include "ito33/finance/derivative_visitor.h"
#include "ito33/finance/option.h"

#include "ito33/dlldecl.h"

namespace ito33
{

namespace finance
{

/**
    Option visitor

    An object of the class must be passed by the user code
    to any functions which may access heterogeneous collections of Derivatives,
    e.g. XML::Reader in IHG currently.
 */
class ITO33_DLLDECL OptionVisitor : public DerivativeVisitor
{
public:
  
  /// Gets shared_ptr to option if an option has been visited, otherwise returns 0
  shared_ptr<Option> GetOption() { return m_pOption; }

  /// Called for an Option
  virtual void OnOption(const Option& contract)  
  { 
    m_pOption = shared_ptr<Option>(new Option(contract));
  } 

private:

  shared_ptr<Option> m_pOption;

};

} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_DERIVATIVEVISITORS_OPTIONVISITOR_H_
