/////////////////////////////////////////////////////////////////////////////
// Name:        ConversionPriceReset_xml.cpp
// Purpose:     Restore ConversionPriceReset object from XML document
// Author:      ZHANG Yunzhi
// Created:     2004-10-20
// RCS-ID:      $Id: conversionprice_reset_xml.cpp,v 1.5 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"

#include "ito33/finance/bondlike/conversionpricereset.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/bondlike/conversionpricereset.h"
#include "ito33/xml/finance/bondlike/common.h"
#include "ito33/xml/finance/common.h"

using namespace ito33::finance;

namespace ito33
{

namespace XML
{ 
  

bool Restore(const xml::node& node,
             shared_ptr<finance::ConversionPriceReset>& pConversionPriceReset)
{
  if( strcmp(node.get_name(), XML_TAG_CONVPRICERESET_ROOT) )
    return false;

  pConversionPriceReset
    = make_ptr( new ConversionPriceReset
                (
                  GetDateFromName(node, XML_TAG_FINANCE_DATE),
                  GetDoubleFromName(node, XML_TAG_CONVPRICERESET_FLOORRATE)      
                ) );
 
  double dCapRate =  GetDoubleFromName(node, XML_TAG_CONVPRICERESET_CAPRATE);

  pConversionPriceReset->SetCap(dCapRate);

  double 
    dMultiplier = GetDoubleFromName(node, XML_TAG_CONVPRICERESET_MULTIPLIER);

  pConversionPriceReset->SetMultiplier(dMultiplier);

  return true;
}


} // namespace XML

} // namespace ito33

