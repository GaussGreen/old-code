/////////////////////////////////////////////////////////////////////////////
// Name:        volatility_xml.cpp
// Purpose:     Restore volatility object from XML document
// Author:      Yann d'Halluin
// Created:     2004-05-18
// RCS-ID:      $Id: volatility_xml.cpp,v 1.10 2006/08/20 09:31:05 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/link.h"

#include "ito33/ihg/volatility.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/volatilitypower.h"
#include "ito33/ihg/volatilitytanh.h"

#include "ihg/xml/volatility.h"

#include "ito33/xml/read.h"

using namespace ito33;
using namespace ito33::XML;

// define the field of the Factory declared in ihg/xml/volatility.h
ITO33_IMPLEMENT_THE_FACTORY(ito33::VolatilityFactory);

ITO33_FORCE_LINK_THIS_MODULE(volatility_xml);


shared_ptr<ito33::ihg::Volatility> 
ito33::ihg::XML::ReadVolatility(const xml::node& node)
{
  shared_ptr<ihg::Volatility> volatility;
  xml::node::const_iterator i;
  for ( i = node.begin(); i != node.end(); ++i ) 
  {
    volatility = make_ptr( VolatilityFactory::Create(i->get_name(), &(*i)) );
    if ( volatility )
      break;
  }
  return volatility;
}


/**
    Restore a flat volatility object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node the session tag in DOM tree
    @return the new volatility object to be deleted by the caller
 */
static
ito33::ihg::Volatility *ReadVolatilityFlat(const xml::node *pNode)
{
  const xml::node& node = *pNode;
  double dVolatilityValue = -99.0; //guarantee to throw an error if not set 
  const xml::node::const_iterator end = node.end();
  xml::node::const_iterator i = node.find(XML_TAG_VOLATILITY_FLAT);

  dVolatilityValue = GetDoubleFromNode(*i);

  return  ( new ihg::VolatilityFlat(dVolatilityValue) );
}

ITO33_DEFINE_VOLATILITY_READER(XML_TAG_VOLATILITYFLAT_ROOT, VolatilityFlat);

//---------------------------------------------------------------------------

static
ito33::ihg::Volatility *ReadVolatilityPower(const xml::node *pNode)
{

  double dAlpha = GetDoubleFromName(*pNode, XML_TAG_VOLATILITYPOWER_ALPHA); 
  double dBeta  = GetDoubleFromName(*pNode, XML_TAG_VOLATILITYPOWER_BETA);
  double dS0    = GetDoubleFromName(*pNode, XML_TAG_VOLATILITYPOWER_S0);

  return ( new ihg::VolatilityPower(dAlpha,dBeta,dS0) );
}

ITO33_DEFINE_VOLATILITY_READER(XML_TAG_VOLATILITYPOWER_ROOT, VolatilityPower);

static
ito33::ihg::Volatility *ReadVolatilityTanh(const xml::node *pNode)
{

  double dLeft  = GetDoubleFromName(*pNode, XML_TAG_VOLATILITYTANH_LEFT); 
  double dRight = GetDoubleFromName(*pNode, XML_TAG_VOLATILITYTANH_RIGHT);
  double dScale = GetDoubleFromName(*pNode, XML_TAG_VOLATILITYTANH_SCALE);
  double dS0    = GetDoubleFromName(*pNode, XML_TAG_VOLATILITYTANH_S0);

  return ( new ihg::VolatilityTanh(dLeft, dRight, dScale, dS0) );
}

ITO33_DEFINE_VOLATILITY_READER(XML_TAG_VOLATILITYTANH_ROOT, VolatilityTanh);
