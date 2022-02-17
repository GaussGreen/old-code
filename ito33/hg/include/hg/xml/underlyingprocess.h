/////////////////////////////////////////////////////////////////////////////
// Name:        hg/xml/underlyingprocess.h
// Purpose:     HG XML stuff for UnderlyingProcess
// Created:     2005/04/15
// RCS-ID:      $Id: underlyingprocess.h,v 1.4 2006/08/19 23:46:52 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/xml/underlyingprocess.h
    @brief HG XML stuff for UnderlyingProcess
 */

#ifndef _HG_XML_UNDERLYINGPROCESS_H_
#define _HG_XML_UNDERLYINGPROCESS_H_

#include "ito33/sharedptr.h"

#include "ito33/hg/common.h"

#define XML_TAG_UNDERLYING_PROCESS "underlying_process"
#define XML_TAG_REGIME_NUMBER "number_of_regimes"
#define XML_TAG_VOLATILITIES "volatilities"
#define XML_TAG_VOLATILITY "volatility"
#define XML_TAG_DEFAULTINTENSITIES "default_intensities"
#define XML_TAG_DEFAULTINTENSITY "default_intensity"

#define XML_TAG_INTER_JUMPS "inter_jumps"
#define XML_TAG_JUMP_FROM "from"
#define XML_TAG_JUMP_TO "to"
#define XML_TAG_JUMPS "jumps"
#define XML_TAG_JUMP "jump"
#define XML_TAG_INTENSITY "intensity"
#define XML_TAG_AMPLITUDE "amplitude"

namespace xml { class node; }

namespace ito33
{

namespace hg 
{
  class ITO33_HG_DLLDECL UnderlyingProcess;

namespace XML
{


/**
    Restore an underlying process object from XML. Otherwise, a NULL pointer
    is returned.

    @param node the root tag in DOM tree where we should find the underlying
                process object.
    @return the new underlying process object
 */
 shared_ptr<UnderlyingProcess> ReadUnderlyingProcess(const xml::node &node);


} // namespace XML

} // namespace hg

} // namespace ito33

#endif // _HG_XML_UNDERLYINGPROCESS_H_
