/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/cboptionoutput.h
// Purpose:     Output class for cb option contracts
// Created:     2006/01/03
// RCS-ID:      $Id: cboptionoutput.h,v 1.5 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/bondlike/cboptionoutput.h
    @brief Output class for cb option contracts.
 */

#ifndef _ITO33_FINANCE_BONDLIKE_CBOPTIONOUTPUT_H_
#define _ITO33_FINANCE_BONDLIKE_CBOPTIONOUTPUT_H_

#include "ito33/sharedptr.h"

#include "ito33/finance/bondlike/bondlikeoutput.h"

namespace ito33
{

namespace finance
{

/**
    ModelOutput class for cb option contracts.

    @nocreate
 */
class ITO33_DLLDECL CBOptionOutput : public ModelOutput
{

public:
  
  /// ctor
  CBOptionOutput() {}

  /// virtual dtor
  virtual ~CBOptionOutput() {}

  /**
      Gets the cb output.

      @return the cb output
   */
  shared_ptr<BondLikeOutput> GetCBOutput() const { return m_pCBOutput; }
  
  /**
      @internal
      @brief Sets the cb output.

      @param pCBOutput the cb output
      @noexport
   */
  void SetCBOutput(const shared_ptr<BondLikeOutput>& pCBOutput) 
  { 
    m_pCBOutput = pCBOutput;
  }

  virtual void Dump(ito33::XML::Tag& tagParent) const;


protected:

  /// The convertible bond output
  shared_ptr<BondLikeOutput> m_pCBOutput;

};  // class CBOptionOutput


} // namespace finance

} // namespace ito33

#endif // _ITO33_FINANCE_BONDLIKE_CBOPTIONOUTPUT_H_
