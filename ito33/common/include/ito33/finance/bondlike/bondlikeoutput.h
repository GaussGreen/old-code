/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/bondlikeoutput.h
// Purpose:     Base output class for bondlike
// Created:     September 21, 2005
// RCS-ID:      $Id: bondlikeoutput.h,v 1.4 2006/05/22 10:13:23 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/bondlike/bondlikeoutput.h
    @brief base output class for bondlike
 */

#ifndef _ITO33_FINANCE_BONDLIKE_BONDLIKEOUTPUT_H_
#define _ITO33_FINANCE_BONDLIKE_BONDLIKEOUTPUT_H_

#include "ito33/finance/modeloutput.h"

namespace ito33
{

namespace XML 
{ 
  class Tag; 
}

namespace finance
{

/**
    Bondlike output class.

    @nocreate
 */
class ITO33_DLLDECL BondLikeOutput : public ModelOutput
{
public:
  
  BondLikeOutput() : ModelOutput(), m_dBondFloor(0.) {}

  /// Dummy virtual dtor for base class
  virtual ~BondLikeOutput() { }
  
  /**
      Gets the bond floor.

      @return The bond floor
   */
  double GetBondFloor() const { return m_dBondFloor; }

  /**
      @internal
      @brief Set the bond floor

      @param dBondFloor bond floor value to set

      @noexport
   */
  void SetBondFloor(double dBondFloor)
  {
    m_dBondFloor = dBondFloor;
  }

  virtual void Dump(XML::Tag& tagParent) const;

protected:
  
  /// The bond floor
  double m_dBondFloor;

}; // class BondLikeOutput


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_BONDLIKE_BONDLIKEOUTPUT_H_
