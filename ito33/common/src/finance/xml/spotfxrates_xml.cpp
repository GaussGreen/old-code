/////////////////////////////////////////////////////////////////////////////
// Name:        spotfxrates_xml.cpp
// Purpose:     Restore SpotFXRates object from XML document
// Author:      Vadim Zeitlin
// Created:     2004-05-04
// RCS-ID:      $Id: spotfxrates_xml.cpp,v 1.6 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/spotfxrates.h"
#include "ito33/finance/numeraire.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/spotfxrates.h"
#include "ito33/xml/finance/common.h"


namespace ito33
{

namespace XML
{


shared_ptr<finance::SpotFXRates> GetSpotFXRatesFromNode(const xml::node& node)
{
  
  xml::node::const_iterator pNodeEle;

  // SpotFXRate is initially empty.  Create and fill in below.
  shared_ptr<finance::SpotFXRates> pSpotFXRates(new finance::SpotFXRates);

  for( pNodeEle = node.begin(); pNodeEle != node.end(); pNodeEle++)
  {
    if(strcmp(pNodeEle->get_name(), XML_TAG_SPOT_FX_ROOT) == 0)
    {
      shared_ptr<finance::Numeraire> pNumeraire1(new finance::Numeraire(
        GetNodeByName(*pNodeEle, 
                      XML_TAG_SPOT_FX_FOREIGN_CURRENCY).get_content() ) );

      shared_ptr<finance::Numeraire> pNumeraire2(new finance::Numeraire(
        GetNodeByName(*pNodeEle, 
                      XML_TAG_SPOT_FX_BASE_CURRENCY).get_content() ) );

      double dRate = GetDoubleFromName(*pNodeEle, XML_TAG_FINANCE_RATE);
      pSpotFXRates->SetFXRate( pNumeraire1, pNumeraire2, dRate);
    }
  }

  return pSpotFXRates;
}

} // namespace XML

} // namespace ito33

