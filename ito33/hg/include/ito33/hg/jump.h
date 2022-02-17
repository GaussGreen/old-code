/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/hg/jump.h
// Purpose:     Definition for a jump between non default regimes
// Created:     2005/01/13
// RCS-ID:      $Id: jump.h,v 1.4 2005/06/02 20:37:33 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/hg/jump.h
   @brief Definition for a jump between non default regimes

   A jump has an intensity and an amplitude. 
 */
#ifndef _ITO33_HG_JUMP_H_
#define _ITO33_HG_JUMP_H_

#include "ito33/hg/common.h"

namespace ito33
{

namespace hg
{


/**
   A jump class for multi-regime models.

   @nocreate
 */
class ITO33_HG_DLLDECL Jump
{
public: 
  /** 
     ctor initializes the member fields to the given values.
     If the values are not valid, the constructor throws.

     @param dIntensity the intensity of the jump
     @param dAmplitude the amplitude of the jump      
   */
  Jump(double dIntensity, double dAmplitude);
  
  /// get the intensity of this jump
  double GetIntensity() const { return m_dIntensity; }
  
  /// get a reference of the amplitude of this jump
  double GetAmplitude() const { return m_dAmplitude; }


private :
 
  /// the intensity of this jump
  double m_dIntensity;

  /// the amplitude of this jump
  double  m_dAmplitude;

}; // class Jump


} // namespace hg

} // namespace ito33

#endif // _ITO33_HG_JUMP_H_
