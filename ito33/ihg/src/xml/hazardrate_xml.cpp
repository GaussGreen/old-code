  /////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/xml/hazardrate_xml.cpp
// Purpose:     Restore hazard rate object from XML document
// Author:      Yann d'Halluin
// Created:     2004-05-18
// RCS-ID:      $Id: hazardrate_xml.cpp,v 1.22 2006/08/20 09:31:05 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/useexception.h"
#include "ito33/vector.h"
#include "ito33/link.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/termstructurecds.h"

#include "ito33/ihg/hazardrate.h"
#include "ito33/ihg/hazardratedecay.h"
#include "ito33/ihg/hazardrateflat.h"
#include "ito33/ihg/hazardratelinear.h"
#include "ito33/ihg/hazardratepower.h"
#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/hazardratecombo.h"

#include "ihg/xml/hazardrate.h"
#include "ihg/xml/spotcomponent.h"

#include "ito33/xml/read.h"
#include "ito33/xml/read_vector.h"
#include "ito33/xml/finance/common.h"

extern const ito33::Error ITO33_UNEXPECTED;

using namespace ito33;
using namespace ito33::XML;

// define the field of the Factory declared in ihg/xml/hazardrate.h
ITO33_IMPLEMENT_THE_FACTORY(ito33::HazardRateFactory);

ITO33_FORCE_LINK_THIS_MODULE(hazardrate_xml);

shared_ptr<ito33::ihg::HazardRate> 
ito33::ihg::XML::ReadHazardRate(const xml::node& node)
{
  shared_ptr<ihg::HazardRate> hazardrate;
  xml::node::const_iterator i;
  for ( i = node.begin(); i != node.end(); ++i ) 
  {
    hazardrate = make_ptr( HazardRateFactory::Create(i->get_name(), &(*i)) );
    if ( hazardrate )
      break;
  }
  return hazardrate;
}

// ----------------------------------------------------------------------------

/**
    Restore a combo hazard rate object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param pNode the parent tag in DOM tree
    @return the new option object to be deleted by the caller
 */
static
ito33::ihg::HazardRate *ReadHazardRateCombo(const xml::node *pNode)
{
    const xml::node& node = *pNode;
    xml::node nodeTimeComp
                  (GetNodeByName(node, XML_TAG_HAZARDRATE_TIMECOMPONENT));

    // Construct and return
    return new ihg::HazardRateCombo
                  ( ihg::XML::GetSpotComponentFromNode(node),// spot component
                    GetVectorFromNode<Date>                      // time array
                                (nodeTimeComp,
                                  XML_TAG_FINANCE_DATES,
                                  XML_TAG_FINANCE_DATE),
                    GetVectorFromNode<double>                   // value array
                                (nodeTimeComp,
                                  XML_TAG_FINANCE_VALUES,
                                  XML_TAG_FINANCE_VALUE)
                    );
    
}

ITO33_DEFINE_HAZARDRATE_READER(XML_TAG_HAZARDRATECOMBO_ROOT, HazardRateCombo);

// ----------------------------------------------------------------------------


/**
    Restore a flat hazard rate object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node the parent tag in DOM tree
    @return the new option object to be deleted by the caller
 */
static
ito33::ihg::HazardRate *ReadHazardRateFlat(const xml::node *pNode)
{
  double dFlatRate = GetDoubleFromName(*pNode, XML_TAG_HAZARDRATE_FLAT);

  return  ( new ihg::HazardRateFlat( dFlatRate) );
}

ITO33_DEFINE_HAZARDRATE_READER(XML_TAG_HAZARDRATEFLAT_ROOT, HazardRateFlat);

// ----------------------------------------------------------------------------
/**
    Restore a decay hazard rate object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node the parent tag in DOM tree
    @return the new option object to be deleted by the caller
 */
static
ito33::ihg::HazardRate *ReadHazardRateDecay(const xml::node *pNode)
{
  const xml::node& node = *pNode;
 
  double dAlpha = GetDoubleFromName(node, XML_TAG_HAZARDRATE_ALPHA); 
  double dS0    =  GetDoubleFromName(node, XML_TAG_HAZARDRATE_S0);
 
  return ( new ihg::HazardRateDecay(dAlpha,dS0) );
}

ITO33_DEFINE_HAZARDRATE_READER(XML_TAG_HAZARDRATEDECAY_ROOT, HazardRateDecay);

// ----------------------------------------------------------------------------

/**
    Restore a linear hazard rate object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node the parent tag in DOM tree
    @return the new option object to be deleted by the caller
 */
static
ito33::ihg::HazardRate *ReadHazardRateLinear(const xml::node *pNode)
{
  const xml::node& node = *pNode;
 
  double dB     = GetDoubleFromName(node, XML_TAG_HAZARDRATE_B); 
  double dSlope =  GetDoubleFromName(node, XML_TAG_HAZARDRATE_SLOPE);
 
  return ( new ihg::HazardRateLinear(dB,dSlope) );
}

ITO33_DEFINE_HAZARDRATE_READER(XML_TAG_HAZARDRATELINEAR_ROOT, HazardRateLinear);

// ----------------------------------------------------------------------------

/**
    Restore a power hazard rate object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node the parent tag in DOM tree
    @return the new option object to be deleted by the caller
 */
static
ito33::ihg::HazardRate *ReadHazardRatePower(const xml::node *pNode)
{
  const xml::node& node = *pNode;

  double dBeta  = GetDoubleFromName(node, XML_TAG_HAZARDRATE_BETA);
  double dAlpha = GetDoubleFromName(node, XML_TAG_HAZARDRATE_ALPHA); 
  double dS0    =  GetDoubleFromName(node, XML_TAG_HAZARDRATE_S0);
 
  return ( new ihg::HazardRatePower(dAlpha,dBeta,dS0) );
}

ITO33_DEFINE_HAZARDRATE_READER(XML_TAG_HAZARDRATEPOWER_ROOT, HazardRatePower);

// ----------------------------------------------------------------------------

/**
    Restore a time only hazard rate object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node the parent tag in DOM tree
    @return the new option object to be deleted by the caller
 */
static
ito33::ihg::HazardRate *ReadHazardRateTimeOnly(const xml::node *pNode)
{
 return 
   new ihg::HazardRateTimeOnly
          (
            GetVectorFromNode<Date>  // time array
                        (*pNode,
                          XML_TAG_FINANCE_DATES,
                          XML_TAG_FINANCE_DATE),
            GetVectorFromNode<double>  // value array
                        (*pNode,
                          XML_TAG_FINANCE_VALUES,
                          XML_TAG_FINANCE_VALUE)
           );
}

ITO33_DEFINE_HAZARDRATE_READER
    (
      XML_TAG_HAZARDRATETIMEONLY_ROOT,
      HazardRateTimeOnly
    );
