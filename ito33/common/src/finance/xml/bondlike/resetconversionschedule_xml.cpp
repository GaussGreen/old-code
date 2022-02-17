/////////////////////////////////////////////////////////////////////////////
// Name:        resetconversionschedule_xml.cpp
// Purpose:     Restore resetconversionschedule object from XML document
// Author:      ZHANG Yunzhi
// Created:     2004-10-20
// RCS-ID:      $Id: resetconversionschedule_xml.cpp,v 1.6 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"

#include "ito33/finance/bondlike/resetconversionschedule.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/bondlike/common.h"
#include "ito33/xml/finance/bondlike/conversionpricereset.h"
#include "ito33/xml/finance/bondlike/conversion.h"
#include "ito33/xml/finance/bondlike/resetconversionschedule.h"

using namespace ito33::finance;

namespace ito33
{

namespace XML
{ 


shared_ptr<finance::ResetConversionSchedule> 
GetResetConversionScheduleFromNode(const xml::node& node)
{
  xml::node::const_iterator
    pNodeRoot = node.find(XML_TAG_BONDLIKE_RESETCONVERSIONSCHEDULE_ROOT);

  if(pNodeRoot == node.end())
  {
    typedef MissingNodeException Exception;
    throw EXCEPTION_MSG(node, XML_TAG_BONDLIKE_RESETCONVERSIONSCHEDULE_ROOT);
  }

  Date startDate 
    = GetDateFromName(*pNodeRoot, XML_TAG_BONDLIKE_RESETCONVERSIONSCHEDULE_START);

  Date endDate 
    = GetDateFromName(*pNodeRoot, XML_TAG_BONDLIKE_RESETCONVERSIONSCHEDULE_END);

  double dInitialConvPrice
    = GetDoubleFromName(*pNodeRoot, XML_TAG_BONDLIKE_RESETCONVERSIONSCHEDULE_INITIAL);
  
  double dCurrentConvPrice
    = GetDoubleFromName(*pNodeRoot, XML_TAG_BONDLIKE_RESETCONVERSIONSCHEDULE_CURRENT);



  finance::ResetFlooredBy floorby = 
                        GetEnumFromName
                          ( 
                            *pNodeRoot,
                            XML_TAG_BONDLIKE_RESETCONVERSIONSCHEDULE_FLOOREDBY,
                            SIZEOF(g_resetFlooredBy),
                            g_resetFlooredBy
                          );



  xml::node::const_iterator pNodeData
      = pNodeRoot->find(XML_TAG_BONDLIKE_RESETCONVERSIONSCHEDULE_CONVPRICERESET);

  if(pNodeData == node.end())
  {
    typedef MissingNodeException Exception;
    throw EXCEPTION_MSG(*pNodeData, XML_TAG_BONDLIKE_RESETCONVERSIONSCHEDULE_CONVPRICERESET);
  }

  shared_ptr<finance::ResetConversionSchedule>
    pResetConversionSchedule(new finance::ResetConversionSchedule
                                 (
                                   startDate, endDate, dInitialConvPrice,
                                   dCurrentConvPrice, floorby
                                 ));
 

  RestoreCommonConversionData(*pNodeRoot, *pResetConversionSchedule);
  
  double dCash
    = GetDoubleFromName(*pNodeRoot, XML_TAG_BONDLIKE_RESETCONVERSIONSCHEDULE_CASH);

  pResetConversionSchedule->SetCash(dCash);

  size_t nCounter = 0;

  for(xml::node::const_iterator pNode = pNodeData->begin();
      pNode != pNodeData->end();
      pNode++)
  {
    shared_ptr<finance::ConversionPriceReset> pReset;
    if(Restore(*pNode, pReset))
    {
      pResetConversionSchedule->AddConversionPriceReset(pReset);
      nCounter++;
    }
  }

  if(nCounter == 0)
    throw EXCEPTION_MSG(ITO33_BAD_DATA, "None ConversionPriceReset Data.");

  return pResetConversionSchedule;
}

} // namespace XML

} // namespace ito33
