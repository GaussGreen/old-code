/////////////////////////////////////////////////////////////////////////////
// Name:        putschedule_xml.cpp
// Purpose:     Restore bondliketerms objects from XML document
// Author:      ZHANG Yunzhi
// Created:     2004-09-03
// RCS-ID:      $Id: putschedule_xml.cpp,v 1.14 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"

#include "ito33/finance/bondlike/putschedule.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/bondlike/putschedule.h"
#include "ito33/xml/finance/bondlike/common.h"
#include "ito33/xml/finance/common.h"

using namespace ito33::finance;

namespace ito33
{

namespace XML
{ 
  

bool Restore(const xml::node& node,
             shared_ptr<finance::PutSchedule>& pPutSchedule)
{
  xml::node::const_iterator
    pnodeRoot = node.find(XML_TAG_BONDLIKE_PUTSCHEDULE_ROOT);

  if( pnodeRoot == node.end() )
    return false;

  const xml::node& nodeRoot = *pnodeRoot;

  pPutSchedule = shared_ptr<finance::PutSchedule>(new PutSchedule);
  
  pPutSchedule->SetKeepAccrued
      ( GetBoolFromName(nodeRoot, XML_TAG_BONDLIKE_KEEPACCRUED) );

  pPutSchedule->SetForfeitCoupon
      ( GetBoolFromName(nodeRoot, XML_TAG_BONDLIKE_FORFEITCOUPON) );

  xml::node::const_iterator pnodeFound;

  pnodeFound = nodeRoot.find(XML_TAG_BONDLIKE_PUTCHEDULE_PUTS);

  if(pnodeFound != node.end())
  {
    xml::node::const_iterator pnodePeriod;
    for(pnodePeriod = pnodeFound->begin();
        pnodePeriod != pnodeFound->end();
        pnodePeriod++)
    {
      // each node is a put period except the first and the last one

      if( !strcmp(pnodePeriod->get_name(),
                  XML_TAG_BONDLIKE_PUTSCHEDULE_PUT_ROOT))
      {
        // All puts have a date
        Date putDate = GetDateFromName(*pnodePeriod, XML_TAG_FINANCE_DATE);

        // Check for strike or yield.  Only one should be defined.
        xml::node::const_iterator pPutType;

        pPutType = pnodePeriod->find(XML_TAG_FINANCE_STRIKE);
        double dStrike = 0.0;
        if ( pPutType != pnodePeriod->end() )
          dStrike = GetDoubleFromName(*pnodePeriod, XML_TAG_FINANCE_STRIKE);

        pPutType = pnodePeriod->find(XML_TAG_FINANCE_YIELD);
        double dYield = 0.0;
        if ( pPutType != pnodePeriod->end() )
          dYield = GetDoubleFromName(*pnodePeriod, XML_TAG_FINANCE_YIELD);

        if (dStrike > 0.0)
          pPutSchedule->AddPutWithStrike(putDate, dStrike);
        else
          pPutSchedule->AddPutWithYield(putDate, dYield);
      }
    }
  }

  return true;
}


} // namespace XML

} // namespace ito33

