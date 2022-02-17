/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/bondlike/putschedule.h
// Purpose:     Names of elements and attributes used in XML for a putschedule
// Author:      ZHANG Yunzhi
// Created:     2004-09-03
// RCS-ID:      $Id: putschedule.h,v 1.15 2006/08/19 22:21:47 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/bondlike/putschedule.h
    @brief Contains the names of the elements used in the XML description of a
           putschedule

    @sa ito33/xml/finance/sessiondata.h
 */

#ifndef _ITO33_XML_FINANCE_BONDLIKE_PUTSCHEDULE_H_
#define _ITO33_XML_FINANCE_BONDLIKE_PUTSCHEDULE_H_

#include "ito33/sharedptr.h"

namespace xml
{
  class node;
}

//=============================================================================
//=              readers                                                      =
//=============================================================================
namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL PutSchedule;
}

namespace XML
{

/// Restore function for put schedule
bool Restore(const xml::node& node,
             shared_ptr<finance::PutSchedule>& pPutSchedule);

}

}

/**
   @name Tag name macros
*/
//@{
#define XML_TAG_BONDLIKE_PUTSCHEDULE_ROOT            "put_schedule"

#define XML_TAG_BONDLIKE_PUTCHEDULE_PUTS             "puts"

#define XML_TAG_BONDLIKE_PUTSCHEDULE_PUT_ROOT        "put"
//@}

#endif // _ITO33_XML_FINANCE_BONDLIKE_PUTSCHEDULE_H_
