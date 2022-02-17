/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/bondlikeoutput.h
// Purpose:     base output class for all bond classes
// Author:      Vadim Zeitlin
// Created:     2005-03-25
// RCS-ID:      $Id: bondlikeoutput.h,v 1.4 2006/05/22 10:15:43 wang Exp $
// Copyright:   (c) 1999-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/ihg/bondlikeoutput.h
    @brief Base output class for all bond classes.

    Output class adding fields common to all bond-like instruments.
*/

#ifndef _ITO33_IHG_BONDLIKEOUTPUT_H_
#define _ITO33_IHG_BONDLIKEOUTPUT_H_

#include "ito33/ihg/dlldecl.h"

#include "ito33/finance/bondlike/bondlikeoutput.h"

namespace ito33
{

namespace ihg
{

/**
   ModelOutput-derived class for Bond and Convertible like.

   @nocreate
 */
class BondLikeOutput : public finance::BondLikeOutput
{
public:
  /// ctor
  BondLikeOutput() : finance::BondLikeOutput() { }

  /// virtual dtor
  virtual ~BondLikeOutput() { }
  
};  // class BondLikeOutput


} // namespace ihg

} // namespace ito33

#endif // _ITO33_IHG_BONDLIKEOUTPUT_H_

