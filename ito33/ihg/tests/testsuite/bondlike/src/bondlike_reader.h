//
// file   : bondlike_reader.h
// author : ZHANG Yunzhi
// RCS-ID:  $Id: bondlike_reader.h,v 1.4 2005/04/19 13:04:54 wang Exp $
//

#ifndef _ITO33_XML_IHG_BONDLIKE_READER_H_
#define _ITO33_XML_IHG_BONDLIKE_READER_H_

#include "ihg/xml/pricingreader.h"
#include "ito33/xml/read.h"

#include <xmlwrapp/node.h>
#include <xmlwrapp/document.h>
#include <xmlwrapp/tree_parser.h>

#define XML_IHGTEST_BASICDDATA "basic_data_test"

namespace ito33
{

namespace ihg
{

namespace XML
{
class BondLikeReader : public ihg::XML::PricingReader
{
public:
  BondLikeReader(const char *filename) : PricingReader(filename) {}

  xml::node GetMainNode() const
  {
    const xml::node& nodeRoot = GetRootNode();
    xml::node::const_iterator
      pNodeFound = nodeRoot.find(XML_IHGTEST_BASICDDATA);

    if(pNodeFound == nodeRoot.end())
    {
      // name is not exactly message but we can still (ab)use this macro here
      typedef ito33::XML::MissingNodeException Exception;
      throw EXCEPTION_MSG(nodeRoot, XML_IHGTEST_BASICDDATA);
    }
    return (*pNodeFound);
  }

};


}

}

}


#endif // #define _ITO33_XML_IHG_BONDLIKE_READER_H_
