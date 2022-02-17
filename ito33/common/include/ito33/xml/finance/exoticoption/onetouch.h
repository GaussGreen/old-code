/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/exoticoption/onetouch.h
// Purpose:     Names of elements and attributes used in XML for a one touch
// Created:     2005/07/04
// RCS-ID:      $Id: onetouch.h,v 1.5 2006/08/19 22:21:47 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/xml/finance/exoticoption/onetouch.h
   @brief Contains the elements of the OneTouch
 */

#ifndef _ITO33_XML_FINANCE_EXOTICOPTION_ONETOUCH_H_
#define _ITO33_XML_FINANCE_EXOTICOPTION_ONETOUCH_H_

#include "ito33/sharedptr.h"

#include "ito33/xml/finance/barrier.h"
#include "ito33/xml/finance/rebate.h"

/**
   @name Tag name macros
*/
//@{

#define XML_TAG_ONETOUCH_ROOT "one_touch"

#define XML_TAG_FXONETOUCH_ROOT "fx_one_touch"

#define XML_TAG_FXONETOUCH_BSBARRIER "bs_barrier"

#define XML_TAG_FXONETOUCH_REFVOL "reference_volatility"

#define XML_TAG_FXONETOUCH_QUOTE "quote"

//@}


namespace xml
{
  class node;
}

namespace ito33
{

namespace finance
{
  class OneTouch;
  class FXOneTouch;
}


namespace XML
{

  /// Restore one touch object
  bool Restore(const xml::node& node, 
               shared_ptr<finance::OneTouch>& pOneTouch);

  /// Restore fx one touch objet
  bool Restore(const xml::node& node, 
               shared_ptr<finance::FXOneTouch>& pFXOneTouch);

}
  

} // namespace ito33

#endif // #ifndef _ITO33_XML_FINANCE_EXOTICOPTION_ONETOUCH_H_
