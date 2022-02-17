/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/hg/bondlikeoutput.h
// Purpose:     base output class for all bond like classes
// Created:     2005/04/11
// RCS-ID:      $Id: bondlikeoutput.h,v 1.3 2006/05/22 10:16:59 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/hg/bondlikeoutput.h
    @brief Base output class for all bond classes.

    Output class adding fields common to all bond-like instruments.
 */

#ifndef _ITO33_HG_BONDLIKEOUTPUT_H_
#define _ITO33_HG_BONDLIKEOUTPUT_H_

#include "ito33/sharedptr.h"
#include "ito33/hg/dlldecl.h"

#include "ito33/finance/bondlike/bondlikeoutput.h"

namespace ito33
{

namespace hg
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

  /// virtual dtor for base class
  virtual ~BondLikeOutput() { }
  
};  // class BondLikeOutput


} // namespace hg

} // namespace ito33

#endif // #ifndef _ITO33_HG_BONDLIKEOUTPUT_H_
