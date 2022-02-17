/////////////////////////////////////////////////////////////////////////////
// Name:        bond_xml.cpp
// Purpose:     Restore Bond object from XML document
// Author:      Vadim Zeitlin
// Created:     2004-05-08
// RCS-ID:      $Id: convertiblelike_xml.cpp,v 1.13 2006/07/28 21:01:10 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/bondlike/convertiblelike.h"
#include "ito33/finance/bondlike/triggeraspercentageof.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/derivative.h"
#include "ito33/xml/finance/bondlike/convertiblelike.h"
#include "ito33/xml/finance/bondlike/trigger.h"


using namespace ito33;
using namespace ito33::XML;
using namespace ito33::finance;

void ito33::XML::GetOptionalConvertibleLikeDataFromNode
        (const xml::node& node, ConvertibleLike& convertibleLike)
{
  xml::node::const_iterator pNodeFound;

  GetOptionalDerivativeDataFromNode(node, convertibleLike);

  /*
    The reason that we do node.find() here is that old xml format is supported
    without any change. 
    
    Note that bConvertIntoNewShare value is always dumped. see ConvertilbeLike::Dump().
  */
  pNodeFound = node.find(XML_TAG_CONVERTIBLELIKE_NEWSHARE);
  if ( pNodeFound != node.end() )
    convertibleLike.SetConvertIntoNewShare
        ( GetBoolFromName(node, XML_TAG_CONVERTIBLELIKE_NEWSHARE) );

  pNodeFound = node.find(XML_TAG_BONDLIKE_TRIGGERINCURRENCYOF);

  // reading of TriggerInCurrency is not symetric to writing.
  // the reason is that sometimes we don't want to write this tag
  // when the security is not cross currency.
  if ( pNodeFound != node.end() )
  {
    TriggerInCurrencyOf
      triggerIn = GetEnumFromNode
                              (
                                *pNodeFound,
                                SIZEOF(g_triggerInCurrencyOfs),
                                g_triggerInCurrencyOfs
                              );

    if (triggerIn == TriggerInCurrencyOf_Underlying)
      convertibleLike.SetTriggerInCurrencyOf( triggerIn );
  }

  // Fixed FX rate
  pNodeFound = node.find(XML_TAG_BONDLIKE_FIXEDFXRATE);
  if ( pNodeFound != node.end() )
    convertibleLike.SetFixedFXRate( GetDoubleFromNode( *pNodeFound ) );

  // Fixed quanto feature
  if( ( pNodeFound = node.find(XML_TAG_CONVERTIBLELIKE_FIXEDQUANTO_ROOT) )
                  != node.end() )
  {
    convertibleLike.SetFixedQuanto
    ( 
      GetDoubleFromName(*pNodeFound, 
        XML_TAG_CONVERTIBLELIKE_FIXEDQUANTO_FXVOL), 
      GetDoubleFromName(*pNodeFound, 
        XML_TAG_CONVERTIBLELIKE_FIXEDQUANTO_CORRELATION) 
    );
  }

  pNodeFound = node.find(XML_TAG_BONDLIKE_TAPO);

  if ( pNodeFound != node.end() )
  {
    TriggerAsPercentageOf
      triggerAs = GetEnumFromNode
                  (
                    *pNodeFound,
                    SIZEOF(g_triggerAsPercentageOfs),
                    g_triggerAsPercentageOfs
                  );

    convertibleLike.SetConversionTriggerAsPercentageOf(triggerAs);  
  }

  // exchangeable feature
  if ( ( pNodeFound = node.find(XML_TAG_CONVERTIBLELIKE_EXCHANGEABLE_ROOT) )
                  != node.end() )
  {
    convertibleLike.SetExchangeable
        ( GetBoolFromName
          (*pNodeFound, XML_TAG_CONVERTIBLELIKE_EXCHANGEABLE_COD) );
  }

  // market price
  GetMarketPrice(node, convertibleLike); 
}
