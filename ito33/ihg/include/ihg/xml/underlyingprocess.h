/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/xml/underlyingprocess.h
// Purpose:     IHG XML stuff for UnderlyingProcess
// Created:     2006/06/01
// RCS-ID:      $Id: underlyingprocess.h,v 1.2 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2004-2006 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/xml/underlyingprocess.h
    @brief IHG XML stuff for UnderlyingProcess
 */

#ifndef _IHG_XML_UNDERLYINGPROCESS_H_
#define _IHG_XML_UNDERLYINGPROCESS_H_

#include "ito33/sharedptr.h"

#include "ito33/ihg/common.h"

namespace xml { class node; }

namespace ito33
{

namespace ihg 
{
  class ITO33_IHG_DLLDECL UnderlyingProcess;

namespace XML
{


/**
    Restores an underlying process object from XML. Otherwise, a NULL pointer
    is returned.

    @param node the root tag in DOM tree where we should find the underlying
                process object.
    @return the new underlying process object
 */
 shared_ptr<UnderlyingProcess> ReadUnderlyingProcess(const xml::node &node);


} // namespace XML

} // namespace hg

} // namespace ito33

#endif // _IHG_XML_UNDERLYINGPROCESS_H_
