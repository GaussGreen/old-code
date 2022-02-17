/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/extrapolationmode.h
// Purpose:     enum of extrapolation mode
// Author:      WANG Xuewen
// Created:     2003/07/25
// RCS-ID:      $Id: extrapolationmode.h,v 1.3 2004/10/05 09:13:38 pedro Exp $
// Copyright:   (c) 1999-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_NUMERIC_EXTRAPOLATIONMODE_H_
#define _ITO33_NUMERIC_EXTRAPOLATIONMODE_H_

namespace ito33
{

namespace numeric
{
/**

  This enum contains extrapolation modes : constant, linear or zero
  
*/

enum ExtrapolationMode
{
  ExtrapolationMode_Constant = 0,//constant extrapolation : the value will 
    //be set to the value at the boundary
  ExtrapolationMode_Linear = 1,//linear extrapolation
  ExtrapolationMode_Zero,//no extrapolation : the value will be set to 0

  ExtrapolationMode_Max //end marker
};

} // namespace numeric

} // namespace ito33

#endif // _ITO33_NUMERIC_EXTRAPOLATIONMODE_H_
