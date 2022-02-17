/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/barrier.h
// Purpose:     Names of elements and attributes used in XML for barrier
// Created:     2005/07/05
// RCS-ID:      $Id: barrier.h,v 1.4 2006/03/24 10:18:27 pedro Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/barrier.h
    @brief Contains the names of the elements used in the XML description of
           Barrier.
 */

#ifndef _ITO33_XML_FINANCE_BARRIER_H_
#define _ITO33_XML_FINANCE_BARRIER_H_

#include "ito33/finance/barriertype.h"
#include "ito33/enum_values_names.h"

/**
   @name Tag name macros
*/
//@{

#define XML_TAG_BARRIER                            "barrier"

#define XML_TAG_BARRIER_TYPE                       "barrier_type"

#define XML_BARRIER_TYPE_UPOUT                     "up_and_out"
#define XML_BARRIER_TYPE_DOWNOUT                   "down_and_out"

//@}

//=============================================================================
//=              helper for read/write                                         =
//=============================================================================
namespace ito33
{

/// Mapping betwen macro name and finance type
const EnumValuesNames<finance::BarrierType> 
g_barrierTypes[] =
{
  { XML_BARRIER_TYPE_UPOUT, finance::Barrier_UpAndOut },
  { XML_BARRIER_TYPE_DOWNOUT, finance::Barrier_DownAndOut }
};


} // namespace ito33

#endif // #ifndef _ITO33_XML_FINANCE_BARRIER_H_
