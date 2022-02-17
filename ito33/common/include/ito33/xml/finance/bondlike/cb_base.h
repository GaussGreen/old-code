/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/bondlike/cb_base.h
// Purpose:     Names of elements and attributes used in XML for 
//              CBBase class
// Author:      ZHANG Yunzhi
// Created:     2004/12/05
// RCS-ID:      $Id: cb_base.h,v 1.3 2006/01/10 16:53:16 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_XML_FINANCE_BONDLIKE_CBBASE_H_
#define _ITO33_XML_FINANCE_BONDLIKE_CBBASE_H_

#include "ito33/common.h"
#include "ito33/xml/finance/bondlike/trigger.h"

namespace xml
{
  class node;
}

namespace ito33
{
namespace finance
{
  class ITO33_DLLDECL CBBase;
}

namespace XML
{
  /**
    Get optional data of CBBase class from xml node.

    @param node xml node
    @param cbBase object whose optional CBBase
              data are required.
    */
  void GetOptionalCBBaseDataFromNode
    (const xml::node& node, finance::CBBase& cbBase);
}

}

#endif // _ITO33_XML_FINANCE_BONDLIKE_CBBASE_H_
